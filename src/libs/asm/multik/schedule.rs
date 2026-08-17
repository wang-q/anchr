//! Scheduling for the multi-k iteration: `assemble_one` runs the
//! single-master skeleton, `assemble_all_masters` runs every k as a master
//! (chain-per-master lanes fed by a bounded-lookahead table builder), and
//! `auto_ks` derives the k ladder from the reads' N50.

use super::master::Master;
use super::{MultikOptions, MultikUnitig};
use crate::libs::asm::assemble::{assemble_unitigs_core, AssembleOptions};
use crate::libs::asm::table::RefineTable;
use anyhow::Result;
use rayon::prelude::*;

/// Pass-0 assemble options (BCALM graph3 semantics, same fast paths as the
/// `asm unitig` command: super-mer counting, adaptive minimizer length, and
/// the DFA-state walk engine).
pub(crate) fn pass0_opts(k0: usize, opts: &MultikOptions) -> AssembleOptions {
    AssembleOptions {
        k: k0,
        min_count_seed: opts.min_count_seed,
        parallel: opts.parallel,
        use_supermer: true,
        supermer_m: Some((12).min((5).max(k0 / 4))),
        use_dfa: true,
        ..AssembleOptions::default()
    }
}
/// Runs the single-skeleton iteration: pass 0 builds maximal unitigs at
/// `k0`, every k in `later_ks` validates the graph (bridge k-mers, internal
/// solidity, reads bridges), and validated chains are compacted into long
/// unitigs.
pub(crate) fn assemble_one(
    infiles: &[String],
    reads: &[(Vec<u8>, Vec<u8>)],
    k0: usize,
    later_ks: &[usize],
    opts: &MultikOptions,
) -> Result<Vec<MultikUnitig>> {
    let read_seqs: Vec<&[u8]> = reads.iter().map(|(s, _)| s.as_slice()).collect();
    let timing = std::env::var_os("ANCHR_MULTIK_TIMING").is_some();
    // FASTA input: count straight from the cached reads buffer (FastK-style
    // super-mer counting, no quality gating — byte-identical to the direct
    // path); FASTQ falls back to `assemble_unitigs_core`'s quality-gated
    // counting, which re-reads the infiles.
    let fasta_input = reads.iter().all(|(_, q)| q.is_empty());
    // The pass-0 table is kept when no rounds follow: `cut` validates
    // against the same k0 table, so the single-k path counts once instead
    // of twice.
    let mut k0_table: Option<RefineTable> = None;
    let mut master = if fasta_input {
        let t = std::time::Instant::now();
        let table = RefineTable::build_supermer_slices(&read_seqs, k0)?;
        if timing {
            eprintln!(
                "reads table k={k0} build: {:.3}s",
                t.elapsed().as_secs_f64()
            );
        }
        let t = std::time::Instant::now();
        let m = Master::pass0(&table, k0, opts);
        if timing {
            eprintln!("master k0={k0} pass0: {:.3}s", t.elapsed().as_secs_f64());
        }
        if later_ks.is_empty() {
            k0_table = Some(table);
        }
        m
    } else {
        let (unitigs, _) = assemble_unitigs_core(infiles, &pass0_opts(k0, opts))?;
        Master::from_unitigs(unitigs, k0)
    };

    // Reads-only probe table (probe_half*2-mers), shared by every round's
    // bridge_filter and the final split/bridge validation: its inputs never
    // change between rounds, so it is built once instead of once per call
    // (R+2 identical full counts per run before).
    let probe_half = 30;
    let probe_len = probe_half * 2;
    let t = std::time::Instant::now();
    let probe_table = RefineTable::build_supermer_slices(&read_seqs, probe_len)?;
    if timing {
        eprintln!(
            "probe table k={probe_len} build: {:.3}s",
            t.elapsed().as_secs_f64()
        );
    }

    let threshold = opts.min_count_extend as u32;
    if later_ks.is_empty() {
        let t = std::time::Instant::now();
        let rebuilt;
        let table = match k0_table.take() {
            Some(t) => t,
            // FASTQ path counted inside pass 0; rebuild for `cut`.
            None => {
                rebuilt = RefineTable::build_supermer_slices(&read_seqs, k0)?;
                rebuilt
            }
        };
        master.cut(&table, threshold)?;
        if timing {
            eprintln!("master k0={k0} cut: {:.3}s", t.elapsed().as_secs_f64());
        }
    }
    for &k in later_ks {
        let t = std::time::Instant::now();
        let table = RefineTable::build_supermer_slices(&read_seqs, k)?;
        if timing {
            eprintln!("reads table k={k} build: {:.3}s", t.elapsed().as_secs_f64());
        }
        master.round(k, &table, &probe_table, probe_half, opts)?;
    }
    let t = std::time::Instant::now();
    master.finalize(&probe_table, probe_half, opts)?;
    if timing {
        eprintln!("master k0={k0} finalize: {:.3}s", t.elapsed().as_secs_f64());
    }
    Ok(master.out)
}

/// Multi-master iteration in k-major order: every k in `ks` is a master
/// whose skeleton the larger ks validate, and the count table at each k is
/// built once and shared by every master's pass 0 and validation rounds at
/// that k (a master's own unitigs are counted separately and summed on
/// lookup, which equals the joint count). This replaces the template's one
/// invocation per master — M masters cost M pass-0 walks plus one reads
/// count per k, not one per (master, round). Each per-k table is dropped at
/// the end of its iteration, so peak memory stays at ~two tables (count +
/// probe) regardless of how many ks are in the ladder.
///
/// With `use_guide` the first master runs first (phase 1) and its validated
/// output guides every later master (megahit `seq2sdbg --contig`): each
/// contig feeds the later masters' counts and probe table as pseudo-reads
/// repeated to the solid threshold. Phase 2 rebuilds reads + guide jointly
/// per k (one direct count, arithmetically the merged counts), so the
/// guided ladder costs ~2 reads counts per k but keeps the same bounded
/// memory.
pub(crate) fn assemble_all_masters(
    reads: &[(Vec<u8>, Vec<u8>)],
    ks: &[usize],
    opts: &MultikOptions,
) -> Result<Vec<MultikUnitig>> {
    let read_seqs: Vec<&[u8]> = reads.iter().map(|(s, _)| s.as_slice()).collect();
    let probe_half = 30;
    let probe_len = probe_half * 2;
    let threshold = opts.min_count_extend as u32;
    let last_k = *ks.last().unwrap();
    let mut out: Vec<MultikUnitig> = Vec::new();

    if opts.use_guide {
        // Phase 1: the guide producer runs its full ladder. Each per-k
        // reads table is built for the round and dropped right after —
        // no cache, so memory stays at ~one table regardless of |ks|.
        let probe = RefineTable::build_supermer_slices(&read_seqs, probe_len)?;
        let mut guide_master = {
            let t = RefineTable::build_supermer_slices(&read_seqs, ks[0])?;
            Master::pass0(&t, ks[0], opts)
        };
        for &k in ks.iter().skip(1) {
            let t = RefineTable::build_supermer_slices(&read_seqs, k)?;
            guide_master.round(k, &t, &probe, probe_half, opts)?;
        }
        guide_master.finalize(&probe, probe_half, opts)?;

        // Guide pseudo-reads: each contig repeated to the solid threshold
        // (the same records `--guide-contigs` writes as a FASTA infile).
        let mut guide_out = std::mem::take(&mut guide_master.out);
        let reps = opts.min_count_seed.max(1);
        let guide_seqs: Vec<&[u8]> = guide_out
            .iter()
            .flat_map(|u| std::iter::repeat_n(u.bases.as_slice(), reps))
            .collect();

        // Phase 2: every later master counts reads + guide in one direct
        // build (equal to the merged reads-only + guide counts) and probes
        // reads + guide, matching a per-master invocation whose infiles
        // include the guide file. The joint table is dropped after its k
        // iteration — peak memory stays at ~two tables (joint + probe).
        let mut seqs_with_guide: Vec<&[u8]> = read_seqs.clone();
        seqs_with_guide.extend_from_slice(&guide_seqs);
        let probe2 = RefineTable::build_supermer_slices(&seqs_with_guide, probe_len)?;
        let mut masters: Vec<Master> = Vec::new();
        for &k in ks.iter().skip(1) {
            let table = RefineTable::build_supermer_slices(&seqs_with_guide, k)?;
            // Same pass-0/rounds overlap as the unguided lane (independent
            // states, shared read-only table).
            let (new_master, rounds_res) = rayon::join(
                || Master::pass0(&table, k, opts),
                || {
                    masters
                        .par_iter_mut()
                        .try_for_each(|m| m.round(k, &table, &probe2, probe_half, opts))
                },
            );
            rounds_res?;
            masters.push(new_master);
            if k == last_k {
                masters.last_mut().unwrap().cut(&table, threshold)?;
            }
        }
        // Finalize independent masters concurrently; outputs are appended
        // in ladder order to keep the byte stream deterministic.
        masters
            .par_iter_mut()
            .try_for_each(|m| m.finalize(&probe2, probe_half, opts))?;
        for m in &mut masters {
            out.append(&mut m.out);
        }
        out.append(&mut guide_out);
    } else {
        // Chain-per-master scheduling: the k-major loop serialized every
        // master's round(k) behind the slowest one (a per-k barrier), yet
        // the only true dependency is WITHIN a master (its rounds are
        // sequential) — masters are independent. Each master runs as its
        // own chain thread receiving the per-k reads tables (Arc, sent in
        // ladder order by the builder lane) and processes pass0 + its
        // rounds at its own pace. The builder pushes each table to its
        // consumer chains through capacity-1 channels IN CHAIN ORDER, so a
        // slow chain throttles the builder before it can run ahead: peak
        // memory stays at ~two tables + probe. The probe table is built on
        // demand (OnceLock) while the chains are still in pass 0.
        let probe_cell: std::sync::OnceLock<
            anyhow::Result<RefineTable, std::sync::Arc<anyhow::Error>>,
        > = std::sync::OnceLock::new();
        let masters: Vec<Master> = std::thread::scope(|scope| -> anyhow::Result<Vec<Master>> {
            type TableLane = std::sync::mpsc::SyncSender<std::sync::Arc<RefineTable>>;
            type ChainLane = std::sync::mpsc::Receiver<std::sync::Arc<RefineTable>>;
            let (txs, rxs): (Vec<TableLane>, Vec<ChainLane>) = ks
                .iter()
                .map(|_| std::sync::mpsc::sync_channel::<std::sync::Arc<RefineTable>>(1))
                .unzip();
            let read_seqs_ref = &read_seqs;
            let probe_cell_ref = &probe_cell;
            let builder = scope.spawn(move || -> anyhow::Result<()> {
                // Build reads tables on the worker pool with a bounded
                // lookahead window (building + buffered tables beyond `next`)
                // and forward them to the chains in ladder order: counting
                // overlaps the chains' round work without letting the
                // builder run unboundedly ahead (peak memory stays at
                // ~window tables + the live chain tables + probe).
                let (done_tx, done_rx): (
                    std::sync::mpsc::Sender<(usize, anyhow::Result<RefineTable>)>,
                    _,
                ) = std::sync::mpsc::channel();
                // `Receiver` is not `Sync`; the scope closure only recvs
                // from the calling thread, so a plain Mutex suffices.
                let done_rx = std::sync::Mutex::new(done_rx);
                let mut buffered: Vec<Option<std::sync::Arc<RefineTable>>> = vec![None; ks.len()];
                let mut next = 0usize;
                let mut spawned = 0usize;
                const BUILD_WINDOW: usize = 3;
                rayon::scope(|s| -> anyhow::Result<()> {
                    while next < ks.len() {
                        while spawned < ks.len() && spawned < next + BUILD_WINDOW {
                            let done_tx = done_tx.clone();
                            let (i, k) = (spawned, ks[spawned]);
                            let timing = std::env::var_os("ANCHR_MULTIK_TIMING").is_some();
                            s.spawn(move |_| {
                                let t = std::time::Instant::now();
                                let res = RefineTable::build_supermer_slices(read_seqs_ref, k);
                                if timing {
                                    eprintln!(
                                        "reads table k={k} build: {:.3}s",
                                        t.elapsed().as_secs_f64()
                                    );
                                }
                                let _ = done_tx.send((i, res));
                            });
                            spawned += 1;
                        }
                        let (i, res) = done_rx
                            .lock()
                            .unwrap()
                            .recv()
                            .map_err(|_| anyhow::anyhow!("multik reads table build failed"))?;
                        // Propagate the first build error seen (the failing
                        // build may complete out of ladder order).
                        buffered[i] = Some(std::sync::Arc::new(res?));
                        // Table k serves chains ks[0..=i] (pass0 at k0 == k,
                        // rounds for k0 < k). A failed send means the chain
                        // already exited (error): stop delivering; the
                        // remaining chains see recv errors and unwind.
                        while next < buffered.len() {
                            let Some(t) = buffered[next].take() else {
                                break;
                            };
                            for tx in &txs[..=next] {
                                if tx.send(t.clone()).is_err() {
                                    return Ok(());
                                }
                            }
                            next += 1;
                        }
                    }
                    Ok(())
                })
            });
            let handles: Vec<_> = rxs
                .into_iter()
                .enumerate()
                .map(|(i, rx)| {
                    scope.spawn(move || -> anyhow::Result<Master> {
                        let k0 = ks[i];
                        let timing = std::env::var_os("ANCHR_MULTIK_TIMING").is_some();
                        let t_chain = std::time::Instant::now();
                        // Probe table on first need (shared OnceLock; the
                        // build joins the ambient worker pool).
                        let probe = || -> anyhow::Result<&RefineTable> {
                            let t = std::time::Instant::now();
                            let cell = probe_cell_ref.get_or_init(|| {
                                let r =
                                    RefineTable::build_supermer_slices(read_seqs_ref, probe_len)
                                        .map_err(std::sync::Arc::new);
                                if timing {
                                    eprintln!(
                                        "probe table k={probe_len} build: {:.3}s",
                                        t.elapsed().as_secs_f64()
                                    );
                                }
                                r
                            });
                            cell.as_ref()
                                .map_err(|e| anyhow::anyhow!("probe table: {e}"))
                        };
                        let mut master: Option<Master> = None;
                        // Only the top-of-ladder master needs the final
                        // table (for `cut`); other chains drop each table
                        // as soon as the round releases it.
                        let mut last_table: Option<std::sync::Arc<RefineTable>> = None;
                        for &k in &ks[i..] {
                            let table = rx.recv().map_err(|_| {
                                anyhow::anyhow!("multik table builder exited unexpectedly")
                            })?;
                            if k == k0 {
                                let t = std::time::Instant::now();
                                master = Some(Master::pass0(&table, k, opts));
                                if timing {
                                    eprintln!(
                                        "master k0={k0} pass0: {:.3}s",
                                        t.elapsed().as_secs_f64()
                                    );
                                }
                            } else {
                                master.as_mut().unwrap().round(
                                    k,
                                    &table,
                                    probe()?,
                                    probe_half,
                                    opts,
                                )?;
                            }
                            if k0 == last_k {
                                last_table = Some(table);
                            }
                        }
                        let mut master = master.unwrap();
                        // The top-of-ladder master has no validating k: cut
                        // chimeric junctions with its own table instead.
                        if k0 == last_k {
                            let t = std::time::Instant::now();
                            master.cut(last_table.as_deref().unwrap(), threshold)?;
                            if timing {
                                eprintln!("master k0={k0} cut: {:.3}s", t.elapsed().as_secs_f64());
                            }
                        }
                        let t = std::time::Instant::now();
                        master.finalize(probe()?, probe_half, opts)?;
                        if timing {
                            eprintln!(
                                "master k0={k0} finalize: {:.3}s chain: {:.3}s",
                                t.elapsed().as_secs_f64(),
                                t_chain.elapsed().as_secs_f64()
                            );
                        }
                        Ok(master)
                    })
                })
                .collect();
            builder
                .join()
                .map_err(|_| anyhow::anyhow!("multik table builder panicked"))??;
            // Outputs are appended in ladder order to keep the byte stream
            // deterministic (each master's own result does not depend on
            // scheduling).
            let mut masters = Vec::with_capacity(handles.len());
            for h in handles {
                masters.push(
                    h.join()
                        .map_err(|_| anyhow::anyhow!("multik master chain panicked"))??,
                );
            }
            Ok(masters)
        })?;
        for mut m in masters {
            out.append(&mut m.out);
        }
    }
    Ok(out)
}
/// Read-length N50 of the input records (used to derive the k sequence).
pub(crate) fn read_n50(reads: &[(Vec<u8>, Vec<u8>)]) -> usize {
    let mut lens: Vec<usize> = reads.iter().map(|(s, _)| s.len()).collect();
    lens.sort_unstable();
    let total: usize = lens.iter().sum();
    let mut acc = 0usize;
    for &l in lens.iter().rev() {
        acc += l;
        if acc * 2 >= total {
            return l;
        }
    }
    0
}

/// Derives the master-k ladder from the read-length N50: the fixed
/// validated ladder `31,41,51,61,71,81,101,121,128,160,192` truncated at
/// `k_max = clamp(N50/2, 81, 192)` (MG1655 N50 339 -> 31..160, G37 N50 408
/// -> 31..192; 150 bp unmerged reads -> 31..81). The top is capped at 192:
/// higher master ks fragment on residual read errors (an error kills every
/// window overlapping it, and at k ~ 0.6×N50 nearly every window overlaps
/// one), which the guide cannot repair (MG1655 auto 51..251: N50 9.4K, 5
/// mis). Low-k masters start at 31 — they anchor low-complexity regions
/// the larger ks validate. Tune per dataset with an explicit `--kmer`.
pub(crate) fn auto_ks(n50: usize) -> Vec<usize> {
    if n50 == 0 {
        return Vec::new();
    }
    let k_max = (n50 / 2).clamp(81, 192);
    [31usize, 41, 51, 61, 71, 81, 101, 121, 128, 160, 192]
        .into_iter()
        .filter(|&k| k <= k_max)
        .collect()
}
