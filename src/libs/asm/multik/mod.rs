//! Multi-k iterative assembly: unitig graph cross-round validation.
//!
//! MetaMDBG-style iteration (see notes/design/asm-multik.md): pass 0
//! assembles maximal unitigs at the first k, then every later k validates
//! the previous unitig graph — each link's bridge k-mer (the current-k
//! window covering the shared (k_prev-1)-mer junction) must be solid in
//! reads + previous unitigs, and every internal k-mer of a long-enough
//! unitig must stay solid (chimeric-unitig cleanup). Surviving links are
//! compacted into longer unitigs, so validated connections grow
//! monotonically instead of being re-split by a larger k.

mod bridge;
mod graph;
mod master;
mod schedule;

use self::schedule::{assemble_all_masters, assemble_one, auto_ks, read_n50};
use super::assemble::read_records;
use anyhow::Result;
use pgr::libs::kmer::key::Kmer;

/// Options for [`assemble_multik`].
#[derive(Debug, Clone)]
pub struct MultikOptions {
    /// Increasing k-mer lengths; `ks[0]` is the master k (the skeleton,
    /// assembled by pass 0) and each later k is a slave k that validates the
    /// graph. Run one master per invocation and merge masters' outputs in
    /// the template (bash parallel). Empty means auto-derive from the
    /// read-length N50 (see [`auto_ks`]).
    pub ks: Vec<usize>,
    /// Solid k-mer threshold for pass 0 (same as `asm unitig`).
    pub min_count_seed: usize,
    /// Solid k-mer threshold for cross-round validation (metaMDBG `>= 2`).
    pub min_count_extend: usize,
    /// Bubble merge: minimum sequence similarity of the collapsed
    /// alternative paths (megahit `--merge-similar`, default 0.95).
    pub merge_similar: f64,
    /// Bubble merge: maximum alternative-path length as a multiple of the
    /// master k (megahit `--merge-level`, default 20).
    pub merge_len: usize,
    /// Worker threads for counting; `0` = rayon global pool.
    pub parallel: usize,
    /// Every k in `ks` is a master that builds its own skeleton (validated
    /// by every larger k) — the multi-master scheme the templates merge
    /// across masters. One invocation shares the reads-only count table at
    /// each k across all masters (k-major order), instead of the template
    /// running one invocation per master and recounting the reads each
    /// time.
    pub all_masters: bool,
    /// With `all_masters`: the first master's validated output guides the
    /// later masters (megahit `seq2sdbg --contig` guidance) — each contig
    /// feeds their counts as pseudo-reads repeated to the solid threshold,
    /// so a low-k master's structure carries into higher-k skeletons.
    pub use_guide: bool,
}

impl Default for MultikOptions {
    fn default() -> Self {
        Self {
            ks: Vec::new(),
            min_count_seed: 3,
            min_count_extend: 2,
            merge_similar: 0.95,
            merge_len: 20,
            parallel: 0,
            all_masters: false,
            use_guide: false,
        }
    }
}

/// One multi-k output unitig (a validated chain).
#[derive(Debug, Clone)]
pub struct MultikUnitig {
    pub bases: Vec<u8>,
    pub coverage: f32,
}

/// Assembles reads into long unitigs by iterating over increasing k.
pub fn assemble_multik(infiles: &[String], opts: &MultikOptions) -> Result<Vec<MultikUnitig>> {
    // Parse the reads once and share the buffer across pass 0, every round's
    // counting, and the probe validation (re-reading per round was the
    // dominant cost: 2R+4 full gz-decompress + parse + phred passes).
    let reads = read_records(infiles)?;
    let mut ks = if opts.ks.is_empty() {
        auto_ks(read_n50(&reads))
    } else {
        opts.ks.clone()
    };
    anyhow::ensure!(!ks.is_empty(), "at least one k-mer length is required");
    for &k in &ks {
        anyhow::ensure!(
            (1..=Kmer::MAX_K).contains(&k),
            "k-mer length must be in 1..={}, got {k}",
            Kmer::MAX_K
        );
    }
    anyhow::ensure!(
        opts.min_count_seed >= 1,
        "min-count-seed must be >= 1, got {} (0 makes every k-mer solid in pass 0 and empties --guide-contigs pseudo-reads)",
        opts.min_count_seed
    );
    anyhow::ensure!(
        opts.min_count_extend >= 1,
        "min-count-extend must be >= 1, got {} (0 silently disables all chimeric-junction validation)",
        opts.min_count_extend
    );
    ks.sort_unstable();
    ks.dedup();

    if opts.all_masters && ks.len() > 1 {
        let mut chains = assemble_all_masters(&reads, &ks, opts)?;
        chains.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
        return Ok(chains);
    }
    // Single-master-skeleton iteration: `ks[0]` builds the graph and every
    // later k validates it.
    let mut chains = assemble_one(infiles, &reads, ks[0], &ks[1..], opts)?;
    chains.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
    Ok(chains)
}

/// Auto-derived master-k sequence from the input reads' read-length N50.
/// Public for the CLI `--print-ks` helper: the template uses it to drive
/// per-master parallel runs with the fixed ladder
/// `31,41,51,61,71,81,101,121,128,160,192` truncated at
/// `clamp(N50/2, 81, 192)` (150 bp reads get 31..81, ~450 bp merged reads
/// get the full ladder up to 192), instead of hard-coding k values.
pub fn auto_ks_for_reads(infiles: &[String]) -> Result<Vec<usize>> {
    Ok(auto_ks(read_n50(&read_records(infiles)?)))
}

#[cfg(test)]
mod tests {
    use super::bridge::{bridge_kmer, kmer_from_bases};
    use super::graph::{
        bubble_merge, merge_chains, progressive_filter, recompact_graph, sequence_similarity,
    };
    use super::master::RollCanon;
    use super::schedule::auto_ks;
    use super::*;
    use crate::libs::asm::assemble::{compute_links, Link, Unitig};
    use crate::libs::asm::table::{base_code, Kmer as TdKmer};
    use pgr::libs::nt::rev_comp;
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static SEQ: AtomicUsize = AtomicUsize::new(0);

    fn temp_fasta(seqs: &[&str]) -> String {
        let mut p = std::env::temp_dir();
        let n = SEQ.fetch_add(1, Ordering::SeqCst);
        p.push(format!("anchr_multik_{}_{}.fa", std::process::id(), n));
        let mut f = std::fs::File::create(&p).unwrap();
        for (i, s) in seqs.iter().enumerate() {
            writeln!(f, ">r{i}").unwrap();
            writeln!(f, "{s}").unwrap();
        }
        p.to_str().unwrap().to_string()
    }

    /// Builds a reads set where two unitigs share a (k1-1)-mer junction and
    /// the bridge k-mer at k2 is present (or absent) in the reads.
    fn run(seqs: &[&str], ks: &[usize]) -> Vec<MultikUnitig> {
        let f = temp_fasta(seqs);
        let opts = MultikOptions {
            ks: ks.to_vec(),
            ..Default::default()
        };
        let out = assemble_multik(std::slice::from_ref(&f), &opts).unwrap();
        let _ = std::fs::remove_file(&f);
        out
    }

    #[test]
    fn rollcanon_matches_canonical() {
        // Rolling dual-strand window must agree with the per-window
        // `TdKmer::canonical` rebuild at every position (deterministic
        // pseudo-random sequence, several k values incl. non-multiples of 4).
        let mut x: u64 = 0x243F6A8885A308D3;
        let seq: Vec<u8> = (0..600)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[x as usize % 4]
            })
            .collect();
        for k in [31usize, 41, 61, 81, 101, 128, 160] {
            let mut rc = RollCanon::new(k, &seq);
            for j in 0..=seq.len() - k {
                if j > 0 {
                    rc.push_code(base_code(seq[j + k - 1]));
                }
                let mut km = TdKmer::new(k);
                for &b in &seq[j..j + k] {
                    km.push_right(base_code(b));
                }
                assert_eq!(
                    rc.canon().cmp_bases(&km.canonical()),
                    std::cmp::Ordering::Equal,
                    "k={k} j={j}"
                );
            }
        }
    }

    #[test]
    fn linear_chain_compacts() {
        // A 21-mer genome duplicated at high coverage with a clean linear
        // junction; k=21 then k=31 should compact the two unitigs into one.
        let g = b"ACGTACGTACGTACGTACGTTTTTACGTACGTACGTACGTACGT".to_vec();
        let mut reads: Vec<String> = Vec::new();
        for _ in 0..20 {
            reads.push(String::from_utf8(g.clone()).unwrap());
        }
        let out = run(
            &reads.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &[21, 31],
        );
        assert!(!out.is_empty());
        assert!(out[0].bases.len() >= g.len() - 20);
    }

    #[test]
    fn unsupported_bridge_is_dropped() {
        // u1 and u2 share a (20)-mer but no read contains the k=31 bridge;
        // pass 0 still emits the link, the k=31 round drops it, so the two
        // unitigs are not merged.
        let s = b"ACGTACGTACGTACGTACGT".to_vec(); // 20-mer shared junction
        let left = format!("AAAA{}", String::from_utf8(s.clone()).unwrap());
        let right = format!("{}CCCC", String::from_utf8(s.clone()).unwrap());
        let mut reads: Vec<String> = Vec::new();
        for _ in 0..20 {
            reads.push(left.clone());
            reads.push(right.clone());
        }
        let out = run(
            &reads.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &[21, 31],
        );
        // The junction k=31 bridge AAAA..S..CCCC does not appear in reads,
        // so the chain must stay split (no merged unitig spanning both).
        let merged = out
            .iter()
            .any(|u| u.bases.len() > left.len().max(right.len()) + 5);
        assert!(!merged);
    }

    #[test]
    fn auto_ks_matches_read_length() {
        // Unmerged 150 bp reads: ladder floor at 81.
        assert_eq!(auto_ks(150), vec![31, 41, 51, 61, 71, 81]);
        // MG1655 merged reads (N50 339): k_max = clamp(169, 81, 192) = 169.
        assert_eq!(
            auto_ks(339),
            vec![31, 41, 51, 61, 71, 81, 101, 121, 128, 160]
        );
        // G37 merged reads (N50 408) and long reads: capped at 192.
        assert_eq!(
            auto_ks(408),
            vec![31, 41, 51, 61, 71, 81, 101, 121, 128, 160, 192]
        );
        assert_eq!(
            auto_ks(15000),
            vec![31, 41, 51, 61, 71, 81, 101, 121, 128, 160, 192]
        );
        // Zero/empty input yields no ks.
        assert!(auto_ks(0).is_empty());
    }

    /// Round 1 prunes a low-abundance isolated unitig into the carry; round
    /// 2 must re-feed it so it survives as output (a silent drop here loses
    /// the low-abundance species).
    #[test]
    fn carried_unitigs_are_refed_into_next_round() {
        // Deterministic pseudo-random sequences (non-periodic, so they
        // assemble into clean linear unitigs).
        let mut x: u64 = 0x1234_5678_9abc_def0;
        let mut rand = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[x as usize % 4]
                })
                .collect()
        };
        let a = rand(250);
        let b = rand(250);
        let mut reads: Vec<String> = Vec::new();
        for _ in 0..30 {
            reads.push(String::from_utf8(a.clone()).unwrap());
        }
        for _ in 0..3 {
            reads.push(String::from_utf8(b.clone()).unwrap());
        }
        let out = run(
            &reads.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &[21, 31, 41],
        );
        // Region B must survive as an independent unitig (re-fed carry),
        // on either strand.
        let b_head = String::from_utf8_lossy(&b[..100]).into_owned();
        let b_head_rc =
            String::from_utf8_lossy(&rev_comp(&b[..100]).collect::<Vec<_>>()).into_owned();
        let survived = out.iter().any(|u| {
            u.bases.len() >= 200
                && (String::from_utf8_lossy(&u.bases).contains(&b_head)
                    || String::from_utf8_lossy(&u.bases).contains(&b_head_rc))
        });
        assert!(
            survived,
            "low-abundance isolated unitig was lost from output"
        );
    }

    #[test]
    fn read_n50_is_median_heavy() {
        let reads = vec![
            (b"A".repeat(100), Vec::new()),
            (b"T".repeat(300), Vec::new()),
            (b"C".repeat(200), Vec::new()),
        ];
        // Total 600 bp; N50 = length covering half = 300.
        assert_eq!(read_n50(&reads), 300);
    }

    fn mk_unitig(bases: &str, cov: f32) -> Unitig {
        Unitig {
            bases: bases.as_bytes().to_vec(),
            id: 0,
            coverage: cov,
            min_cov: 0,
            max_cov: 0,
            circular: false,
            abundances: Vec::new(),
        }
    }

    /// Builds a two-path bubble: source `s` (right end shares `x`) branches
    /// into `m1`/`m2` (left end `x`, right end `y`) and reconverges at `r`
    /// (left end `y`). `m1` is the high-coverage main path, `m2` the
    /// alternative.
    fn mk_bubble(m1: &str, m2: &str) -> (Vec<Unitig>, Vec<Vec<Link>>, usize, usize) {
        let x: String = "ACGTACGTACGTACGTACGTACGTACGTAC".chars().take(30).collect();
        let y: String = "TGCATGCATGCATGCATGCATGCATGCATGCAT"
            .chars()
            .take(30)
            .collect();
        let s = format!("{}{}", "A".repeat(60), x);
        let m1 = format!("{}{}{}", x, m1, y);
        let m2 = format!("{}{}{}", x, m2, y);
        let r = format!("{}{}", y, "T".repeat(60));
        let unitigs = vec![
            mk_unitig(&s, 100.0),
            mk_unitig(&m1, 60.0),
            mk_unitig(&m2, 30.0),
            mk_unitig(&r, 100.0),
        ];
        let links = compute_links(&unitigs, 31);
        (unitigs, links, 1, 2) // dominant index 1, alternative index 2
    }

    #[test]
    fn sequence_similarity_banded() {
        let a = b"ACGTACGTACGTACGTACGT";
        assert!((sequence_similarity(a, a, 0.95) - 1.0).abs() < 1e-9);
        let mut b = a.to_vec();
        b[5] = b'G';
        assert!(sequence_similarity(a, &b, 0.95) >= 0.95);
        // A 50% divergent sequence falls below the threshold.
        let c: Vec<u8> = a
            .iter()
            .map(|&x| if x == b'A' { b'T' } else { x })
            .collect();
        assert!(sequence_similarity(a, &c, 0.95) < 0.95);
    }

    /// A zero-width band (`--merge-similar 1.0`, or near-1.0 on short
    /// sequences) must still report identical sequences as 1.0 and any
    /// differing pair as 0.0 — not collapse everything to 0.0, which would
    /// silently disable bubble merging.
    #[test]
    fn sequence_similarity_zero_width_band() {
        let a = b"ACGTACGTACGTACGTACGT";
        assert!((sequence_similarity(a, a, 1.0) - 1.0).abs() < 1e-9);
        assert!((sequence_similarity(a, a, 0.99) - 1.0).abs() < 1e-9);
        let mut b = a.to_vec();
        b[5] = b'G';
        assert_eq!(sequence_similarity(a, &b, 1.0), 0.0);
        assert_eq!(sequence_similarity(a, &b, 0.99), 0.0);
    }

    #[test]
    fn bubble_merge_keeps_dominant_path() {
        // m1 (dominant) = C*40, m2 (variant) = C*38 + G A (2 substitutions
        // over 100 bp -> ~98% similar): the alternative must be dropped and
        // the main path fused through source -> m1 -> right.
        let (mut unitigs, mut links, dom, alt) =
            mk_bubble(&"C".repeat(40), &format!("{}GA", "C".repeat(38)));
        let alt_bases = unitigs[alt].bases.clone();
        let mut branch = vec![false; unitigs.len()];
        let variants = bubble_merge(&mut unitigs, &mut links, &mut branch, 31, 0.95, 20);
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].bases, alt_bases);
        assert!(!unitigs[dom].bases.is_empty());
        // The surviving graph is one unique chain source -> dom -> right.
        let chains = merge_chains(&unitigs, &links, &branch, 31, None, 0).unwrap();
        let expected = 60 + 30 + 40 + 30 + 60; // s + m1[30..] + r[30..]
        assert_eq!(chains[0].bases.len(), expected);
    }

    #[test]
    fn bubble_merge_rejects_divergent_paths() {
        // m2 = G*40: 40/100 substitutions -> similarity 0.6 < 0.95, so no
        // merge happens and both paths stay.
        let (mut unitigs, mut links, _dom, alt) = mk_bubble(&"C".repeat(40), &"G".repeat(40));
        let mut branch = vec![false; unitigs.len()];
        let variants = bubble_merge(&mut unitigs, &mut links, &mut branch, 31, 0.95, 20);
        assert!(variants.is_empty());
        assert_eq!(unitigs.len(), 4);
        assert_eq!(unitigs[alt].bases.len(), 100);
    }

    #[test]
    fn bubble_merge_rejects_long_middles() {
        // merge_len=2 -> max_len ~65 bp, middles are 100 bp: no merge.
        let (mut unitigs, mut links, _dom, _alt) = mk_bubble(&"C".repeat(40), &"C".repeat(40));
        let mut branch = vec![false; unitigs.len()];
        let variants = bubble_merge(&mut unitigs, &mut links, &mut branch, 31, 0.95, 2);
        assert!(variants.is_empty());
        assert_eq!(unitigs.len(), 4);
    }

    #[test]
    fn bubble_merge_rejects_different_sinks() {
        // A second right unitig with a different shared y-mer: the two
        // middles converge at different sinks, so no bubble.
        let x: String = "ACGTACGTACGTACGTACGTACGTACGTAC".chars().take(30).collect();
        let y1: String = "TGCATGCATGCATGCATGCATGCATGCATGCAT"
            .chars()
            .take(30)
            .collect();
        let y2: String = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".chars().take(30).collect();
        let s = format!("{}{}", "A".repeat(60), x);
        let m1 = format!("{}{}{}", x, "C".repeat(40), y1);
        let m2 = format!("{}{}{}", x, "C".repeat(40), y2);
        let r1 = format!("{}{}", y1, "T".repeat(60));
        let r2 = format!("{}{}", y2, "T".repeat(60));
        let mut unitigs = vec![
            mk_unitig(&s, 100.0),
            mk_unitig(&m1, 60.0),
            mk_unitig(&m2, 30.0),
            mk_unitig(&r1, 100.0),
            mk_unitig(&r2, 100.0),
        ];
        let mut links = compute_links(&unitigs, 31);
        let mut branch = vec![false; unitigs.len()];
        let variants = bubble_merge(&mut unitigs, &mut links, &mut branch, 31, 0.95, 20);
        assert!(variants.is_empty());
        assert_eq!(unitigs.len(), 5);
    }

    /// The bridge k-mer covering a u→v junction equals
    /// `upstream tail (k-1) + downstream continuation` (forward strand).
    #[test]
    fn bridge_kmer_forward_strand() {
        let s = "ACGTACGTACGTACGTACGT"; // shared 20-mer
        let u = mk_unitig(&format!("{}{}", "A".repeat(40), s), 30.0);
        let v = mk_unitig(&format!("{}{}", s, "C".repeat(40)), 30.0);
        let link = Link {
            to: 1,
            from_rc: false,
            to_rc: false,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("forward bridge");
        // Expected: u[30..] (A*10 + s) + v[20] (C).
        let exp = format!("{}C", String::from_utf8_lossy(&u.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// Reverse-complemented target: u's tail matches rc(v's tail), so the
    /// partner joins reverse-complemented and its continuation is rc(v)[20]
    /// (rc of the C-run = G).
    #[test]
    fn bridge_kmer_reverse_target() {
        let s = "GGATCACAGTCTACACTGCT"; // shared 20-mer (non-palindromic)
        let s_rc = "AGCAGTGTAGACTGTGATCC"; // rc(s)
        let u = mk_unitig(&format!("{}{}", "A".repeat(40), s), 30.0);
        // v stored reverse-complemented: tail = rc(s).
        let v = mk_unitig(&format!("{}{}", "C".repeat(40), s_rc), 30.0);
        let link = Link {
            to: 1,
            from_rc: false,
            to_rc: true,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("rc bridge");
        // Expected: u[30..] + rc(v)[20] = (A*10 + s) + G.
        let exp = format!("{}G", String::from_utf8_lossy(&u.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// from_rc=true: u's left end is the junction source, the partner sits
    /// upstream (v tail == u head), continuation is u's own base after the
    /// shared overlap.
    #[test]
    fn bridge_kmer_left_end_source() {
        let s = "ACGTACGTACGTACGTACGT"; // shared 20-mer
        let u = mk_unitig(&format!("{}{}", s, "A".repeat(40)), 30.0);
        let v = mk_unitig(&format!("{}{}", "C".repeat(40), s), 30.0);
        let link = Link {
            to: 1,
            from_rc: true,
            to_rc: false,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("left bridge");
        // Expected: v[30..] (C*10 + s) + u[20] (A).
        let exp = format!("{}A", String::from_utf8_lossy(&v.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// Progressive filter keeps unique-chain unitigs regardless of their
    /// (low) coverage and only prunes isolated nodes by abundance.
    #[test]
    fn progressive_filter_keeps_main_path() {
        // A (cov 3) -> B (cov 5): both unique-chain, both must survive.
        let a = mk_unitig(&"A".repeat(60), 3.0);
        let b = mk_unitig(&"C".repeat(60), 5.0);
        let mut unitigs = vec![a, b];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
        ];
        let mut branch = vec![false; unitigs.len()];
        let (dropped, _) = progressive_filter(&mut unitigs, &mut links, &mut branch, 21, None, 0);
        assert_eq!(unitigs.len(), 2, "unique chain must survive the filter");
        assert!(dropped.is_empty());
    }

    /// Isolated low-abundance unitigs are pruned and collected as dropped.
    #[test]
    fn progressive_filter_prunes_isolated() {
        let a = mk_unitig(&"A".repeat(60), 30.0);
        let iso = mk_unitig(&"C".repeat(60), 2.0);
        let mut unitigs = vec![a, iso];
        let mut links = vec![Vec::new(), Vec::new()];
        let mut branch = vec![false; unitigs.len()];
        let (dropped, _) = progressive_filter(&mut unitigs, &mut links, &mut branch, 21, None, 0);
        assert_eq!(
            unitigs.len(),
            1,
            "only the high-coverage isolated unitig survives"
        );
        assert_eq!(dropped.len(), 1);
        assert_eq!(dropped[0].bases.len(), 60);
    }

    /// Recompaction merges a unique chain into one longer unitig; stale
    /// links whose extremities were absorbed by the merge disappear (the
    /// recomputed links only connect real endpoint overlaps).
    #[test]
    fn recompact_merges_chain_and_drops_stale_links() {
        let s1 = "ACGTACGTACGTACGTACGT";
        let s2 = "GATCGATCGATCGATCGATC";
        let a = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let b = mk_unitig(&format!("{}{}{}", s1, "G".repeat(60), s2), 30.0);
        let c = mk_unitig(&format!("{}{}", s2, "C".repeat(80)), 30.0);
        let mut unitigs = vec![a, b, c];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            vec![Link {
                to: 2,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
        ];
        let mut branch = vec![false; unitigs.len()];
        recompact_graph(&mut unitigs, &mut links, &mut branch, 21, None, 0, false);
        // A(100) + B[20..](80) + C[20..](80) = 260.
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].bases.len(), 260);
    }

    /// A convergence node (two predecessors) is not swallowed into either
    /// chain: strict end uniqueness, not link iteration order, decides chain
    /// membership (SKESA's "predecessor == 1" invariant).
    #[test]
    fn merge_chains_keeps_convergence_split() {
        let s1 = "ACGTACGTACGTACGTACGT";
        // u's end == v's begin == w's end: v has two predecessors.
        let u = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let v = mk_unitig(&format!("{}{}", s1, "G".repeat(80)), 30.0);
        let w = mk_unitig(&format!("{}{}", "T".repeat(80), s1), 30.0);
        let unitigs = vec![u, v, w];
        let links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
        ];
        let branch = vec![false; unitigs.len()];
        let out = merge_chains(&unitigs, &links, &branch, 21, None, 0).unwrap();
        assert_eq!(out.len(), 3, "v must not merge with either predecessor");
        let total: usize = out.iter().map(|u| u.bases.len()).sum();
        assert_eq!(total, 300, "no unitig is fused at the junction");
        drop(links);
        drop(unitigs);
    }

    /// Recompaction applies the same strict end uniqueness: the convergence
    /// node survives as its own unitig instead of being merged into the
    /// first-visited chain.
    #[test]
    fn recompact_keeps_convergence_split() {
        let s1 = "ACGTACGTACGTACGTACGT";
        let u = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let v = mk_unitig(&format!("{}{}", s1, "G".repeat(80)), 30.0);
        let w = mk_unitig(&format!("{}{}", "T".repeat(80), s1), 30.0);
        let mut unitigs = vec![u, v, w];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
        ];
        let mut branch = vec![false; unitigs.len()];
        recompact_graph(&mut unitigs, &mut links, &mut branch, 21, None, 0, false);
        // v keeps both incoming links (it was not absorbed), so the merged
        // graph still has three segments and v's begin retains two edges.
        assert_eq!(unitigs.len(), 3);
        let v_in: usize = links
            .iter()
            .map(|ls| ls.iter().filter(|l| l.to == 1 && !l.from_rc).count())
            .sum();
        assert_eq!(v_in, 2, "convergence node must keep both predecessors");
    }
}
