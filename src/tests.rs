use super::*;

/// Helper: run OrfCaller on a single FASTA record and return formatted output lines.
fn run_orfm(
    seq_input: &str,
    min_length: usize,
    table_id: usize,
    position_limit: Option<usize>,
) -> String {
    let caller = OrfCaller::new(table_id, min_length, position_limit).unwrap();
    // Parse the FASTA input to get name, comment, seq
    let mut reader = needletail::parse_fastx_reader(seq_input.as_bytes()).unwrap();
    let mut output = Vec::new();
    while let Some(Ok(record)) = reader.next() {
        let (name, comment) = split_header(record.id());
        let seq = record.raw_seq();
        let orfs = caller.find_orfs(name, comment, &seq);
        for orf in orfs {
            output.push(orf.header());
            output.push(String::from_utf8(orf.protein).unwrap());
        }
    }
    output.join("\n") + if output.is_empty() { "" } else { "\n" }
}

/// Helper: run and get transcripts too.
fn run_orfm_with_transcripts(
    seq_input: &str,
    min_length: usize,
    table_id: usize,
    position_limit: Option<usize>,
) -> (String, String) {
    let caller = OrfCaller::new(table_id, min_length, position_limit).unwrap();
    let mut reader = needletail::parse_fastx_reader(seq_input.as_bytes()).unwrap();
    let mut protein_output = Vec::new();
    let mut transcript_output = Vec::new();
    while let Some(Ok(record)) = reader.next() {
        let (name, comment) = split_header(record.id());
        let seq = record.raw_seq();
        let orfs = caller.find_orfs(name, comment, &seq);
        for orf in orfs {
            transcript_output.push(orf.header());
            transcript_output.push(String::from_utf8(orf.transcript(&seq)).unwrap());
            protein_output.push(orf.header());
            protein_output.push(String::from_utf8(orf.protein).unwrap());
        }
    }
    let prot = protein_output.join("\n") + if protein_output.is_empty() { "" } else { "\n" };
    let trans = transcript_output.join("\n")
        + if transcript_output.is_empty() {
            ""
        } else {
            "\n"
        };
    (prot, trans)
}

#[test]
fn test_single_file_defaults() {
    let input = ">638202197:1-99\nATGGATGCTGAAAAAAGATTGTTCTTAAAGGCATTAAAGGAAAAGTTTGAAGAAGACCCAAGAGAAAAATACACTAAGTTCTATGTCTTTGGCGGATGG\n";
    let expected = ">638202197:1-99_1_1_1\nMDAEKRLFLKALKEKFEEDPREKYTKFYVFGGW\n";
    assert_eq!(run_orfm(input, 96, 1, None), expected);
}

#[test]
fn test_longer_min_length() {
    let input = ">638202197 NP_247840 methyl coenzyme M reductase I, subunit alpha (mcrA) [Methanocaldococcus jannaschii DSM 2661: NC_000909] (+)strand\nATGGATGCTGAAAAAAGATTGTTCTTAAAGGCATTAAAGGAAAAGTTTGAAGAAGACCCAAGAGAAAAATACACTAAGTTCTATGTCTTTGGCGGATGGAGACAGTCAGCAAGAAAAGAGAATTCGTTGAGGCAGCACAAAAAaTTAATTGAGAAGAGAGGAGGAATTCCATTTTACAACCCAGATATTGGAGTTCCATTGGGGCAGAGAAAATTAATGCCTTACAAAGTTTCAAATACAGATGCAATTGTTGAAGGGGATGACTTACACTTCATGAACAACGCTGCAATGCAGCAGTTCTGGGATGACATAAGAAGAACAGTTATCGTTGGGATGGATACAGCTCACGCTGTTCTTGAAAAGAGATTGGGGGTAGAGGTTACTCCAGAAACaATTAATGAATACATGGAAACTATTAACCACGCTCTCCCAGGAGGAGCAGTTGTTCAGGAGCACATGGTTGAGGTCCACCCAGCATTAGTCTGGGACTGTTACGCTAAGATATTCACTGGAGATGACGAATTAGCAGATGAGATTGACAAGAGGTTCTTAATTGACATTAAcAAGTTGTTCCCAGAAGAGCAGGCAGAACAaATCAAGAAGGCAATCGGTAAGAGAACATACCAAGTTTCAAGAGTTCCAACATTAGTCGGTAGAGTTTGTGATGGGGGAACAATAGCAaGATGGAGTGCTATGCAGATTGGAATGTCATTCATTACAGCTTACAAGCTCTGTGCTGGGGAGGCAGCAATTGCTGACTTCTCATACGCTGCAAAGCACGCTGATGTCATTCAGATGGCTTCATTCTTGCCAGCAAGAAGAGCAAGAGGGCCAAATGAACCAGGAGGTATCTTCTTCGGAGTCTTGGCAGATATTGTTCAAACATCAAGAGTTTCAGATGACCCAGTTGAACAGTCATTAGAGGTTGTTGCTGCTGGGGCTATGTTGTATGACCAAATCTGGTTAGGAGGATACATGTCTGGAGGAGTCGGATTTACACAGTATGCTACAGCAACCTATACAGATGACATCTTGGATGACTTCTCATACTACGGATATGACTACATAACCAAGAAATATGGAGGATGCAaCAGCGTAAAACCAACaATGGATGTTGTTGAAGATATTGCTACTGAAGTAACTTTATATGGTTTAGAGCAGTATGACACCTTCCCAGCATTGTTAGAAGACCACTTCGGAGGTTCCCAAAGAGCAGGGGTTACAGCTGCTGCAGCAGGTATTACAACTGCATTAGCTACAGGAAACTCAAACGCTGGAGTTAACGGATGGTATCTAAGCCAGATATTGCACAAaGAATACCACAGCAGATTAGGATTCTATGGTTATGACTTACAAGACCAGTGTGGAGCAGCCAACTCATTATCATTCAGAAACGATGAAGGTTCCCCATTAGAATTGAGAGGGCCTAACTATCCAAACTACGCAATGAACGTTGGTCACCAAGGAGAATATGCTGGAATTACACAGGCTGCACACTCAGCAAGAGGAGACGCATTTGCATTGAACCCATTAATTAAGGTTGCATTTGCAGACCCATCATTAGTCTTTGACTTCACACATCCAAGAAAAGAGTTTGCAAGAGGTGCTTTAAGAGAATTCGAGCCAGCTGGAGAAAGAGATCCAATCATCCCAGCTCACTAA\n";
    let output = run_orfm(input, 300, 1, None);
    assert!(output.starts_with(
        ">638202197_1_1_1 NP_247840 methyl coenzyme M reductase I, subunit alpha (mcrA)"
    ));
    assert!(output.contains("MDAEKRLFLKALKEKFEEDPREKYTKFYVFGGW"));
}

#[test]
fn test_toy_example_internal_frame2() {
    let input = ">eg\nAATGTGAA\n";
    let expected =
        ">eg_2_2_1\nM\n>eg_1_1_2\nNV\n>eg_3_3_3\nCE\n>eg_1_4_4\nHI\n>eg_2_5_5\nSH\n>eg_3_6_6\nFT\n";
    assert_eq!(run_orfm(input, 3, 1, None), expected);
}

#[test]
fn test_n_characters() {
    let input = ">eg\nTTAANA\n";
    let expected = ">eg_1_1_1\nLX\n";
    assert_eq!(run_orfm(input, 6, 1, None), expected);
}

#[test]
fn test_lower_case() {
    let input = ">eg\nTTAAaA\n";
    let expected = ">eg_1_1_1\nLK\n";
    assert_eq!(run_orfm(input, 6, 1, None), expected);
}

#[test]
fn test_position_limit_basic() {
    let input = ">eg\nTTAANAGGGGGGGGGG\n";

    let expected = ">eg_1_1_1\nLX\n";
    assert_eq!(run_orfm(input, 6, 1, Some(6)), expected);

    let expected = ">eg_1_1_1\nLX\n>eg_2_5_2\nXL\n";
    assert_eq!(run_orfm(input, 6, 1, Some(7)), expected);

    let expected = ">eg_1_1_1\nLX\n>eg_3_3_2\nXR\n>eg_2_5_3\nXL\n>eg_3_6_4\nPX\n";
    assert_eq!(run_orfm(input, 6, 1, Some(8)), expected);
}

#[test]
fn test_position_limit_full() {
    let input = ">eg\nTTAANAGGGGGGGGGG\n";
    let expected = ">eg_1_1_1\nLXGGG\n>eg_5_2_2\nXGGG\n>eg_3_3_3\nXRGG\n>eg_4_4_4\nPPPX\n>eg_2_5_5\nPPPXL\n>eg_3_6_6\nPPPX\n";
    assert_eq!(run_orfm(input, 6, 1, Some(16)), expected);
    // -l 80 should give same result (longer than sequence)
    assert_eq!(run_orfm(input, 6, 1, Some(80)), expected);
}

#[test]
fn test_bad_codons() {
    let input = ">eg\nGYATCATAGGCCAGCCGCTGTCCAGATGCACCGGTTCATCTGCGTCAGACGACGATCTTCACCCGGTAACCCCCGCCGATCACCAGATACTCGGCCTCCCCGCGCAAAGGATGGCTCGGCAGCAGATCGTTGAAGAACAGGAGCTTCACGACCGGAACTTCGGTTTCGAGGATATAGGCACCGAACTGCCCCGCCTGCGTCCGGTCAGCGGAGAAAGAAACGATGTTGTTGAGACGCACGAGGATTTCCCGTCCGTTGCCGGCCCC\n";
    let expected = ">eg_1_1_1\nXS\n>eg_2_2_2\nXH\n>eg_3_3_3\nII\n>eg_2_5_4\nMX\n>eg_3_6_5\nYD\n";
    assert_eq!(run_orfm(input, 6, 1, Some(8)), expected);
}

#[test]
fn test_transcripts() {
    let input = ">eg\nAATGTGAA\n";
    let (_, transcripts) = run_orfm_with_transcripts(input, 3, 1, None);
    let expected_lines: Vec<&str> = vec![
        ">eg_2_2_1",
        "ATG",
        ">eg_1_1_2",
        "AATGTG",
        ">eg_3_3_3",
        "TGTGAA",
        ">eg_1_4_4",
        "CACATT",
        ">eg_2_5_5",
        "TCACAT",
        ">eg_3_6_6",
        "TTCACA",
    ];
    let actual_lines: Vec<&str> = transcripts.trim().split('\n').collect();
    assert_eq!(actual_lines, expected_lines);
}

#[test]
fn test_transcript_odd_chars() {
    let input = ">eg\nGYATCATAGGCCAGCCGCTGTCCAGATGCACCGGTTCATCTGCGTCAGACGACGATCTTCACCCGGTAACCCCCGCCGATCACCAGATACTCGGCCTCCCCGCGCAAAGGATGGCTCGGCAGCAGATCGTTGAAGAACAGGAGCTTCACGACCGGAACTTCGGTTTCGAGGATATAGGCACCGAACTGCCCCGCCTGCGTCCGGTCAGCGGAGAAAGAAACGATGTTGTTGAGACGCACGAGGATTTCCCGTCCGTTGCCGGCCCC\n";
    let (_, transcripts) = run_orfm_with_transcripts(input, 6, 1, Some(8));
    let expected_lines: Vec<&str> = vec![
        ">eg_1_1_1",
        "GYATCA",
        ">eg_2_2_2",
        "YATCAT",
        ">eg_3_3_3",
        "ATCATA",
        ">eg_2_5_4",
        "ATGATN",
        ">eg_3_6_5",
        "TATGAT",
    ];
    let actual_lines: Vec<&str> = transcripts.trim().split('\n').collect();
    assert_eq!(actual_lines, expected_lines);
}

#[test]
fn test_alternative_codon_table() {
    let input = ">eg\nAATGTGAA\n";
    let output = run_orfm(input, 3, 4, None);
    let expected_lines: Vec<&str> = vec![
        ">eg_1_1_1",
        "NV",
        ">eg_2_2_2",
        "MW",
        ">eg_3_3_3",
        "CE",
        ">eg_1_4_4",
        "HI",
        ">eg_2_5_5",
        "SH",
        ">eg_3_6_6",
        "FT",
    ];
    let actual_lines: Vec<&str> = output.trim().split('\n').collect();
    assert_eq!(actual_lines, expected_lines);
}

/// Integration test: compare output against the original OrfM binary if available.
#[test]
fn test_against_original_orfm() {
    let orfm_path = std::env::var("ORFM_BIN").unwrap_or_else(|_| {
        let home = std::env::var("HOME").unwrap_or_default();
        format!("{}/git/OrfM/orfm", home)
    });
    if !std::path::Path::new(&orfm_path).exists() {
        eprintln!(
            "Skipping integration test: OrfM binary not found at {}",
            orfm_path
        );
        return;
    }

    // Test cases: (input, min_length, position_limit)
    // Only use table 1 (default) since the installed OrfM may not support -c
    let test_cases: Vec<(&str, usize, Option<usize>)> = vec![
        (">eg\nAATGTGAA\n", 3, None),
        (">eg\nTTAANA\n", 6, None),
        (">eg\nTTAAaA\n", 6, None),
        (">eg\nTTAANAGGGGGGGGGG\n", 6, Some(6)),
        (">eg\nTTAANAGGGGGGGGGG\n", 6, Some(7)),
        (">eg\nTTAANAGGGGGGGGGG\n", 6, Some(8)),
        (">eg\nTTAANAGGGGGGGGGG\n", 6, Some(16)),
        (">638202197:1-99\nATGGATGCTGAAAAAAGATTGTTCTTAAAGGCATTAAAGGAAAAGTTTGAAGAAGACCCAAGAGAAAAATACACTAAGTTCTATGTCTTTGGCGGATGG\n", 96, None),
    ];

    for (input, min_len, pos_limit) in test_cases {
        let rust_output = run_orfm(input, min_len, 1, pos_limit);

        let mut cmd = std::process::Command::new(&orfm_path);
        cmd.arg(format!("-m{}", min_len));
        if let Some(pl) = pos_limit {
            cmd.arg(format!("-l{}", pl));
        }
        cmd.stdin(std::process::Stdio::piped());
        cmd.stdout(std::process::Stdio::piped());

        let mut child = cmd.spawn().expect("Failed to run original orfm");
        {
            use std::io::Write;
            let stdin = child.stdin.as_mut().unwrap();
            stdin.write_all(input.as_bytes()).unwrap();
        }
        let output = child.wait_with_output().unwrap();
        let c_output = String::from_utf8_lossy(&output.stdout).to_string();

        assert_eq!(
            rust_output, c_output,
            "Mismatch for input={:?} min_len={} pos_limit={:?}\nRust:\n{}\nC:\n{}",
            input, min_len, pos_limit, rust_output, c_output
        );
    }
}
