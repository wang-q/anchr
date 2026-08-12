#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;

/// Parses a FASTA into (header, sequence) pairs.
fn parse_fa(data: &[u8]) -> Vec<(String, String)> {
    let mut recs = Vec::new();
    let mut cur: Option<(String, String)> = None;
    for line in std::str::from_utf8(data).unwrap().lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            cur = Some((rest.to_string(), String::new()));
        } else if let Some(c) = cur.as_mut() {
            c.1.push_str(line);
        }
    }
    if let Some(c) = cur {
        recs.push(c);
    }
    recs
}

/// A linear genome compresses into a single maximal unitig.
#[test]
fn command_asm_unitig_linear() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Random 60 bp (all 30 k-mers unique -> linear, not cyclic).
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 4 identical reads so every k-mer is solid (count 4 >= seed 3).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n>r3\n{seq}\n>r4\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert_eq!(recs.len(), 1, "expected one unitig, got {}", recs.len());
    let rc: String = seq
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect();
    assert!(recs[0].1 == seq || recs[0].1 == rc, "got {}", recs[0].1);
}

/// A bubble (two parallel paths) stays split: each branch is its own unitig
/// instead of being merged into one representative path.
#[test]
fn command_asm_unitig_keeps_branches() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    let prefix = "ACGT".repeat(15); // 60 bp shared prefix
    let window_a = "ACGT".repeat(8); // 32 bp path A
    let window_b = "ACGA".repeat(8); // 32 bp path B (variant)
    let suffix = "TGCA".repeat(15); // 60 bp shared suffix
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    // 10 reads per path: every k-mer has count 10 (solid at seed 3).
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    // Prefix, both variant windows, and suffix form separate unitigs; a
    // bubble-popping assembler would merge them into fewer, longer contigs.
    assert!(recs.len() >= 4, "expected >= 4 unitigs, got {}", recs.len());
    assert!(
        recs.iter().all(|(_, s)| s.len() < path_a.len()),
        "a unitig spans the whole bubble"
    );
}

/// Output is deterministic and non-empty on the Lambda dataset.
#[test]
fn command_asm_unitig_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("u1.fa");
    let out2 = out_dir.path().join("u2.fa");
    for out in [&out1, &out2] {
        AnchrCmd::new()
            .args(&[
                "asm",
                "unitig",
                "tests/bbtools/Lambda/R1.2k.fq.gz",
                "tests/bbtools/Lambda/R2.2k.fq.gz",
                "-o",
                out.to_str().unwrap(),
                "--kmer",
                "31",
            ])
            .assert()
            .success();
    }
    assert_eq!(std::fs::read(&out1).unwrap(), std::fs::read(&out2).unwrap());
    let recs = parse_fa(&std::fs::read(&out1).unwrap());
    assert!(!recs.is_empty());
}

/// Raising the solid threshold drops low-count k-mers (bcalm `-abundance-min`).
#[test]
fn command_asm_unitig_min_count_seed() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 2 identical reads: every k-mer has count 2 (not solid at the default
    // threshold of 3, solid at --min-count-seed 2).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    let default_out = out_dir.path().join("default.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            default_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    assert!(parse_fa(&std::fs::read(&default_out).unwrap()).is_empty());
    let strict_out = out_dir.path().join("strict.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--min-count-seed",
            "2",
        ])
        .assert()
        .success();
    assert_eq!(parse_fa(&std::fs::read(&strict_out).unwrap()).len(), 1);
}

/// A branching graph emits GFA segments and (k-1)-overlap links.
#[test]
fn command_asm_unitig_gfa() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.gfa");
    let prefix = "ACGT".repeat(15); // 60 bp shared prefix
    let window_a = "ACGT".repeat(8); // 32 bp path A
    let window_b = "ACGA".repeat(8); // 32 bp path B (variant)
    let suffix = "TGCA".repeat(15); // 60 bp shared suffix
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--gfa",
        ])
        .assert()
        .success();
    let gfa = std::fs::read_to_string(&out).unwrap();
    let mut headers = 0;
    let mut segments = 0;
    let mut links = 0;
    for line in gfa.lines() {
        match line.as_bytes().first() {
            Some(b'H') => headers += 1,
            Some(b'S') => segments += 1,
            Some(b'L') => {
                links += 1;
                assert!(line.ends_with("\t30M"), "overlap: {line}");
            }
            _ => {}
        }
    }
    assert_eq!(headers, 1);
    assert!(segments >= 4, "segments: {segments}");
    assert!(links >= 4, "links: {links}");
}

/// `--gfa` with the default `--min-contig-len` must not emit `L` edges to
/// segments that were dropped by the length filter (dangling references).
#[test]
fn command_asm_unitig_gfa_no_dangling_links() {
    // A bubble (P+X+S vs P+Y+S) yields a short unitig bounded by branches that
    // shares endpoints with two longer unitigs; the short one is filtered out
    // by min_len but must not leave dangling `L` targets.
    let reads = [
        "GAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATT",
        "TGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGA",
        "TAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATT",
        "AGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCG",
        "GTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTG",
        "CCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAAGACAGATAGT",
        "TTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGT",
        "CAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCA",
        "GTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGG",
        "GTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGA",
        "TACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGA",
        "CCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGT",
    ];
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("bubble.fa");
    let out = out_dir.path().join("out.gfa");
    let mut fa = String::new();
    for (i, r) in reads.iter().enumerate() {
        fa.push_str(&format!(">r{i}\n{r}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--gfa",
        ])
        .assert()
        .success();
    let gfa = std::fs::read_to_string(&out).unwrap();
    let mut segments = std::collections::HashSet::new();
    for line in gfa.lines() {
        if line.starts_with('S') {
            segments.insert(line.split('\t').nth(1).unwrap().to_string());
        }
    }
    let mut links = 0usize;
    for line in gfa.lines() {
        if !line.starts_with('L') {
            continue;
        }
        links += 1;
        // L: <from> <from_ori> <to> <to_ori> <cigar>.
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(
            segments.contains(fields[1]) && segments.contains(fields[3]),
            "dangling L reference: {line}"
        );
    }
    assert!(links > 0, "expected links between kept unitigs");
}

/// `--links` appends BCALM-style `L:` entries to FASTA headers.
#[test]
fn command_asm_unitig_links_header() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    let prefix = "ACGT".repeat(15);
    let window_a = "ACGT".repeat(8);
    let window_b = "ACGA".repeat(8);
    let suffix = "TGCA".repeat(15);
    let path_a = format!("{prefix}{window_a}{suffix}");
    let path_b = format!("{prefix}{window_b}{suffix}");
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">a{i}\n{path_a}\n>b{i}\n{path_b}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--links",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert!(recs.len() >= 4);
    let any_link = recs.iter().any(|(h, _)| {
        h.split_whitespace()
            .any(|f| f.starts_with("L:+:") || f.starts_with("L:-:"))
    });
    assert!(any_link, "no L: entries in headers");
}

/// `--links` with the default `--min-contig-len` must not emit `L:` header
/// entries referencing unitigs dropped by the length filter (dangling refs,
/// the FASTA counterpart of the GFA `L`-edge zero-dangling policy).
#[test]
fn command_asm_unitig_links_no_dangling() {
    let reads = [
        "GAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATT",
        "TGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGA",
        "TAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATT",
        "AGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCG",
        "GTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTG",
        "CCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACATGCGTTCGCTCTATTGACTACGACGCGCTCATTCCCTTGTCGGAGAGTTATGGAACAAGGACGCTGTCTGAGACTAGAAGACAGATAGT",
        "TTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGT",
        "CAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCA",
        "GTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGG",
        "GTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGA",
        "TACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGA",
        "CCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCTGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGT",
    ];
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("bubble.fa");
    let out = out_dir.path().join("out.fa");
    let mut fa = String::new();
    for (i, r) in reads.iter().enumerate() {
        fa.push_str(&format!(">r{i}\n{r}\n"));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--links",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    let ids: std::collections::HashSet<String> = recs
        .iter()
        .map(|(h, _)| {
            h.split(',')
                .next()
                .unwrap()
                .trim_start_matches("unitig_")
                .to_string()
        })
        .collect();
    let mut links = 0usize;
    for (h, _) in &recs {
        for f in h.split_whitespace() {
            // L:<from_ori>:<to_id>:<to_ori>.
            if let Some(rest) = f.strip_prefix("L:") {
                links += 1;
                let to_id = rest.split(':').nth(1).unwrap();
                assert!(
                    ids.contains(to_id),
                    "dangling L: reference to {to_id} in header {h}"
                );
            }
        }
    }
    assert!(links > 0, "expected L: entries between kept unitigs");
}

/// A k-mer above the 128-base key limit must fail cleanly instead of
/// panicking in `Kmer::new().expect()` (zero-panic policy).
#[test]
fn command_asm_unitig_rejects_kmer_above_limit() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "129",
        ])
        .assert()
        .failure();
}

/// `-o` must not overwrite an input file (the writer is opened before the
/// reads are consumed).
#[test]
fn command_asm_unitig_outfile_not_input() {
    let infile = "tests/bbtools/Lambda/R1.2k.fq.gz";
    AnchrCmd::new()
        .args(&["asm", "unitig", infile, "-o", infile])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));
}

/// A cyclic genome (periodic sequence) assembles into a circular unitig
/// flagged in the FASTA header.
#[test]
fn command_asm_unitig_circular() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Periodic 80 bp "genome": with k=31 the k-mer graph is a 4-kmer cycle.
    let genome = "ACGT".repeat(20);
    let mut fa = String::new();
    for i in 0..10 {
        let start = (i * 2) % 40;
        fa.push_str(&format!(">r{i}\n{}\n", &genome[start..start + 40]));
    }
    std::fs::write(&infile, fa).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "unitig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert_eq!(
        recs.len(),
        1,
        "expected one circular unitig, got {}",
        recs.len()
    );
    assert!(recs[0].0.contains("circular"), "header: {}", recs[0].0);
}
