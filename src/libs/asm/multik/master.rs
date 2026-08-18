//! One master's iterative state: the pass-0 skeleton at `k0`, the graph
//! mutated by every validation round, and the finalized output chains.

use super::bridge::{bridge_filter, bridge_kmer, split_by_bridge};
use super::graph::{
    bubble_merge, internal_repeat_bridge_split, merge_chains, progressive_filter, recompact_graph,
    remove_unsupported, tip_remover, weak_link_remover,
};
use super::schedule::pass0_opts;
use super::{MultikOptions, MultikUnitig};
use crate::libs::asm::assemble::{assemble_unitigs_from_table, compute_links, Link, Unitig};
use crate::libs::asm::table::{base_code, Kmer as TdKmer, RefineTable};
use anyhow::Result;
use std::collections::HashSet;

/// Read-only count lookup for validation steps. Implemented for a single
/// table and for the sum of two disjoint tables: a shared reads-only table
/// plus the master's own unitig-only table, whose summed counts equal the
/// joint count over both inputs (so per-master rounds reuse the shared
/// reads table instead of recounting the reads at every k).
pub(crate) trait CountView: Sync {
    fn count(&self, kmer: &TdKmer) -> u32;

    /// Whether `kmer` — already canonicalized by the caller (e.g. via
    /// [`RollCanon`]) — reaches `threshold`. Views over rolling windows
    /// override this to skip the redundant per-window canonicalization.
    fn is_solid(&self, kmer: &TdKmer, threshold: u32) -> bool {
        self.count(kmer) >= threshold
    }
}

impl CountView for RefineTable {
    fn count(&self, kmer: &TdKmer) -> u32 {
        self.get_count(kmer)
    }

    fn is_solid(&self, kmer: &TdKmer, threshold: u32) -> bool {
        self.get_count_canonical(kmer) >= threshold
    }
}

struct SumView<'a>(&'a RefineTable, &'a RefineTable);

impl CountView for SumView<'_> {
    fn count(&self, kmer: &TdKmer) -> u32 {
        self.0
            .get_count(kmer)
            .saturating_add(self.1.get_count(kmer))
    }

    fn is_solid(&self, kmer: &TdKmer, threshold: u32) -> bool {
        // Short-circuit: solid reads alone settle most windows, so the
        // master's own unitig table is only consulted for the rare rest.
        let a = self.0.get_count_canonical(kmer);
        a >= threshold || a.saturating_add(self.1.get_count_canonical(kmer)) >= threshold
    }
}

/// Rolling window holding the forward k-mer and its reverse complement in
/// lockstep: advancing costs two packed shifts and the canonical form is a
/// byte compare, instead of rebuilding the rc base-by-base (O(k)) at every
/// window position.
pub(crate) struct RollCanon {
    fw: TdKmer,
    rc: TdKmer,
}

/// Debug trace: dumps the current unitig graph (FASTA + link lines) to
/// `$ANCHR_MULTIK_TRACE_DIR/<tag>.fa` so an external tool can follow a
/// junction across the rounds. No-op when the env var is unset.
fn trace_graph(tag: &str, unitigs: &[Unitig], links: &[Vec<Link>]) {
    use std::io::Write;
    let Some(dir) = std::env::var_os("ANCHR_MULTIK_TRACE_DIR") else {
        return;
    };
    let path = std::path::Path::new(&dir).join(format!("{tag}.fa"));
    let mut f = std::fs::File::create(&path).unwrap();
    for (i, u) in unitigs.iter().enumerate() {
        writeln!(f, ">u{i} len={} cov={:.1}", u.bases.len(), u.coverage).unwrap();
        for chunk in u.bases.chunks(100) {
            writeln!(f, "{}", String::from_utf8_lossy(chunk)).unwrap();
        }
        let ls: Vec<String> = links[i]
            .iter()
            .map(|l| format!("{}{}", l.to, if l.from_rc { "-" } else { "+" }))
            .collect();
        writeln!(f, "#L {}", ls.join(" ")).unwrap();
    }
}

impl RollCanon {
    /// Window over the first `k` bases of `bases` (`bases.len() >= k`).
    pub(crate) fn new(k: usize, bases: &[u8]) -> Self {
        let mut fw = TdKmer::new(k);
        for &b in &bases[..k] {
            fw.push_right(base_code(b));
        }
        Self { rc: fw.rc(), fw }
    }

    /// Advance one base (2-bit `code` of the incoming base): the window's
    /// rc prepends the complemented base and drops its 3' end.
    pub(crate) fn push_code(&mut self, code: u8) {
        self.fw.push_right(code);
        self.rc.push_left(3 - code);
    }

    /// The canonical (lexicographically smaller) strand.
    pub(crate) fn canon(&self) -> &TdKmer {
        if self.fw.cmp_bases(&self.rc) != std::cmp::Ordering::Greater {
            &self.fw
        } else {
            &self.rc
        }
    }
}

/// One master's iterative state: the pass-0 skeleton at `k0`, the graph
/// mutated by every validation round, and the finalized output chains.
pub(crate) struct Master {
    k0: usize,
    unitigs: Vec<Unitig>,
    links: Vec<Vec<Link>>,
    branch: Vec<bool>,
    /// Low-abundance unitigs pruned by the last round's progressive filter;
    /// re-fed into the next round's graph (megahit bubble re-feeding).
    carried: Vec<Unitig>,
    carried_branch: Vec<bool>,
    /// Finalized output (set by [`Master::finalize`]).
    pub(crate) out: Vec<MultikUnitig>,
}

impl Master {
    /// Pass 0 from a prebuilt reads(-only) count table at `k0`.
    pub(crate) fn pass0(table: &RefineTable, k0: usize, opts: &MultikOptions) -> Self {
        let unitigs = assemble_unitigs_from_table(table, &pass0_opts(k0, opts));
        Self::from_unitigs(unitigs, k0)
    }

    /// Pass 0 from already-built unitigs (the FASTQ fallback path counts
    /// inside `assemble_unitigs_core` with quality gating).
    pub(crate) fn from_unitigs(unitigs: Vec<Unitig>, k0: usize) -> Self {
        let links = compute_links(&unitigs, k0);
        trace_graph(&format!("pass0_k{k0}"), &unitigs, &links);
        // Repeat fragments (a repeated element's shared k-mers fan into four
        // or more flanking unitigs at the master k) connect to four or more
        // partners. Snapshot that branch status at pass 0: a chain through
        // such a node picks one of several genomic contexts and joins
        // distant loci, so its links never participate in recompaction. The
        // flag propagates through every unitig reindexing below.
        let branch = links
            .iter()
            .map(|ls| ls.iter().map(|l| l.to).collect::<HashSet<_>>().len() >= 4)
            .collect();
        Self {
            k0,
            unitigs,
            links,
            branch,
            carried: Vec::new(),
            carried_branch: Vec::new(),
            out: Vec::new(),
        }
    }

    /// Single-round master (no larger k validates it): still cut chimeric
    /// junctions — internal k-mers without read support (the joined
    /// sequence does not exist in the reads) mark a chimeric unitig, which
    /// validation rounds would otherwise catch. `base` counts the reads
    /// (plus a guide, when the invocation is guided).
    pub(crate) fn cut(&mut self, base: &RefineTable, threshold: u32) -> Result<()> {
        let seqs: Vec<&[u8]> = self.unitigs.iter().map(|u| u.bases.as_slice()).collect();
        // Direct counting: unitig sequence is unique, so super-mer collapse
        // never engages (byte-identical table, half the sort volume).
        let unitig_table = RefineTable::build_direct_slices(&seqs, self.k0)?;
        let view = SumView(base, &unitig_table);
        remove_unsupported(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            &view,
            self.k0,
            threshold,
        );
        Ok(())
    }

    /// One validation round at `k` (k > k0): re-feed the carried unitigs,
    /// validate every link's bridge k-mer and every unitig's internal
    /// k-mers against `base` + this master's unitigs, filter links by
    /// reads bridges, recompact, and prune low-abundance unitigs into the
    /// carry. `probe` is the shared reads-only probe table; `repeat` is the
    /// shared reads-only short-k table that gates repeat bridges during
    /// recompaction.
    pub(crate) fn round(
        &mut self,
        k: usize,
        base: &RefineTable,
        probe: &RefineTable,
        probe_half: usize,
        repeat: &RefineTable,
        repeat_k: usize,
        opts: &MultikOptions,
    ) -> Result<()> {
        let t_round = std::time::Instant::now();
        self.unitigs.append(&mut self.carried);
        self.branch.append(&mut self.carried_branch);
        let t0 = std::time::Instant::now();
        let seqs: Vec<&[u8]> = self.unitigs.iter().map(|u| u.bases.as_slice()).collect();
        // Direct counting (unique sequence; see `Master::cut`).
        let unitig_table = RefineTable::build_direct_slices(&seqs, k)?;
        let view = SumView(base, &unitig_table);
        let t_count = t0.elapsed().as_secs_f64();
        let timing = std::env::var_os("ANCHR_MULTIK_TIMING").is_some();
        let threshold = opts.min_count_extend as u32;
        let k0 = self.k0;
        // 1. Cross-round link validation (solveEdges): the bridge k-mer
        // covering the junction must be solid at the current k.
        let t1 = std::time::Instant::now();
        for (i, ls) in self.links.iter_mut().enumerate() {
            ls.retain(|l| {
                let u = &self.unitigs[i];
                let v = &self.unitigs[l.to];
                // Short unitigs (shorter than the current k-1 window) cannot
                // provide a full bridge k-mer; skip their validation and let
                // the final compaction decide via actual extremity matching
                // (the link is guaranteed to share a (k0-1)-mer).
                if u.bases.len() < k - 1 || v.bases.len() < k - 1 {
                    return true;
                }
                bridge_kmer(u, v, l, k, k0).is_some_and(|km| view.count(&km) >= threshold)
            });
        }
        let t2 = std::time::Instant::now();
        // 2. Chimeric-unitig cleanup (removeUnsupportedUnitigs): every
        // internal current-k k-mer of a long-enough unitig must be solid.
        remove_unsupported(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            &view,
            k,
            threshold,
        );
        let t3 = std::time::Instant::now();
        // 2.5 Reads-bridge validation: every surviving link must have reads
        // fully covering a probe spanning the junction. Chimeric links (two
        // distant regions joined by a shared k-mer) have no bridging reads
        // and are dropped BEFORE recompaction, so the per-round merge cannot
        // fix them into the main path (prevents relocation misassemblies).
        bridge_filter(&self.unitigs, &mut self.links, probe, k0, probe_half, 2);
        let t4 = std::time::Instant::now();
        trace_graph(
            &format!("r{k0}_k{k}_after_bridge"),
            &self.unitigs,
            &self.links,
        );
        // 3. Recompact unique chains so the main path grows between rounds
        // (metaMDBG recompacts after every abundance-removal round). No
        // abundance pruning here — that is deferred to the final filter, so
        // single-genome coverage fluctuation never drops real content.
        recompact_graph(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            k0,
            Some(repeat),
            repeat_k,
        );
        let t5 = std::time::Instant::now();
        trace_graph(
            &format!("r{k0}_k{k}_after_compact"),
            &self.unitigs,
            &self.links,
        );
        // 4. Prune low-abundance branching/isolated unitigs and carry them
        // into the next round (the final round's carry becomes output).
        let (dropped, dropped_branch) = progressive_filter(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            k0,
            Some(repeat),
            repeat_k,
        );
        self.carried = dropped;
        self.carried_branch = dropped_branch;
        trace_graph(
            &format!("r{k0}_k{k}_after_prog"),
            &self.unitigs,
            &self.links,
        );
        // 4.5 Split unitigs at internal positions that have no reads support
        // (recompact_graph can fuse chimeric unitigs — the abundance filter's
        // recompaction also does this — so split them here before the next
        // round uses the unitigs as the skeleton).  probe is reads-only, so
        // unitig self-counts never mask a junction window.
        split_by_bridge(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            probe,
            k0,
            probe_half,
            threshold,
        );
        trace_graph(
            &format!("r{k0}_k{k}_after_split"),
            &self.unitigs,
            &self.links,
        );
        if timing {
            let n = self.unitigs.len();
            let bp: usize = self.unitigs.iter().map(|u| u.bases.len()).sum();
            let edges: usize = self.links.iter().map(|l| l.len()).sum();
            eprintln!(
                "master k0={k0} round k={k}: {n} unitigs, {bp} bp, {edges} edges, count {t_count:.3}s link {:.3}s unsup {:.3}s bridge {:.3}s compact {:.3}s prog {:.3}s total {:.3}s",
                t2.duration_since(t1).as_secs_f64(),
                t3.duration_since(t2).as_secs_f64(),
                t4.duration_since(t3).as_secs_f64(),
                t5.duration_since(t4).as_secs_f64(),
                t5.elapsed().as_secs_f64(),
                t_round.elapsed().as_secs_f64()
            );
        }
        Ok(())
    }

    /// Final steps after the last validation round: split unitigs at
    /// reads-unsupported windows, re-verify links, merge bubbles, clean
    /// tips/weak links, and compact validated chains into the output.
    pub(crate) fn finalize(
        &mut self,
        probe: &RefineTable,
        probe_half: usize,
        repeat: &RefineTable,
        repeat_k: usize,
        opts: &MultikOptions,
    ) -> Result<()> {
        let k0 = self.k0;
        // Split unitigs at internal positions that have no reads support
        // (the abundance filter's recompaction can fuse chimeric links into
        // a single unitig — the source of G37 relocations). Every probe
        // window of a unitig must occur in the reads at least
        // `min_count_extend` times; an unsupported window is a chimeric
        // junction and the unitig is cut there.
        split_by_bridge(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            probe,
            k0,
            probe_half,
            opts.min_count_extend as u32,
        );
        // Split unitigs at internal narrow high-coverage spikes in the
        // short-k repeat table: a repeat bridge inside a unitig (reads from
        // both copies share a few windows) fused two distant loci, and the
        // 130-mer split above misses it because the junction windows have
        // reads support. End repeats are already gated by is_repeat_bridge
        // during the final compaction, so only internal spikes are cut here.
        internal_repeat_bridge_split(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            repeat,
            repeat_k,
            k0,
        );
        // Re-verify the links recomputed by the split: the new extremities
        // may join distant regions, so every surviving link needs bridging
        // reads.
        bridge_filter(&self.unitigs, &mut self.links, probe, k0, probe_half, 2);

        // Megahit-style bubble merge: divergent paths that reconverge at the
        // same partner collapse to the highest-coverage path, so the main
        // path is not interrupted at variant/error sites. Removed
        // alternatives stay independent output (megahit writes them to
        // bubble_seq.fa).
        let variants = bubble_merge(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            k0,
            opts.merge_similar,
            opts.merge_len,
        );

        // Megahit-style cleaning on the final (compacted) unitigs: drop
        // short low-depth tips and disconnect weak links (depth-proportional,
        // see megahit tip_remover / weak_link_remover). Doing this per round
        // was too aggressive on the k0=21 fragments (G37 longest 52.8k ->
        // 32.6k); after compaction the tips are real and few.
        tip_remover(
            &mut self.unitigs,
            &mut self.links,
            &mut self.branch,
            k0 * 2,
            20.0,
        );
        weak_link_remover(&mut self.unitigs, &mut self.links, 0.05);
        trace_graph(&format!("f{k0}_premerge"), &self.unitigs, &self.links);

        // Final compaction: merge validated chains into long unitigs.
        let mut chains = merge_chains(
            &self.unitigs,
            &self.links,
            &self.branch,
            k0,
            Some(repeat),
            repeat_k,
        )?;
        chains.extend(
            std::mem::take(&mut self.carried)
                .into_iter()
                .map(|u| MultikUnitig {
                    bases: u.bases,
                    coverage: u.coverage,
                }),
        );
        chains.extend(variants.into_iter().map(|u| MultikUnitig {
            bases: u.bases,
            coverage: u.coverage,
        }));
        chains.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
        self.out = chains;
        Ok(())
    }
}
