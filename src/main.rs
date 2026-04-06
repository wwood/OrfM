use clap::Parser;
use fastq::Record;
use flate2::read::GzDecoder;
use std::io::{self, BufWriter, Read, Write};

#[derive(Parser)]
#[command(
    name = "orfm",
    about = "Find and translate open reading frames",
    version
)]
struct Args {
    /// Input FASTA/FASTQ file (gzipped or uncompressed). Reads from stdin if not provided.
    input: Option<String>,

    /// Minimum ORF length in nucleotides (must be a multiple of 3)
    #[arg(short = 'm', default_value = "96")]
    min_length: usize,

    /// Output nucleotide transcripts to this file
    #[arg(short = 't')]
    transcript: Option<String>,

    /// Ignore sequence beyond this position
    #[arg(short = 'l')]
    position_limit: Option<usize>,

    /// Codon table ID (1-25, NCBI standard)
    #[arg(short = 'c', default_value = "1")]
    codon_table: usize,

    /// Print stop codon (*) at the end of ORFs that are bounded by a stop codon
    #[arg(short = 'p', default_value = "false")]
    print_stop_codons: bool,

    /// Only output ORFs that are bounded by an in-frame stop codon (exclude ORFs that end at the end of the sequence without a stop codon)
    #[arg(short = 's', default_value = "false")]
    stop_codon_only: bool,

    /// Require OrfM to be at least this version (e.g. 1.4.0); exit with error if not
    #[arg(short = 'r')]
    required_version: Option<String>,
}

fn split_header(header: &[u8]) -> (&str, &str) {
    let header_str = std::str::from_utf8(header).unwrap_or("unknown");
    match header_str.find(|c: char| c.is_whitespace()) {
        Some(pos) => (&header_str[..pos], &header_str[pos + 1..]),
        None => (header_str, ""),
    }
}

/// Detect whether a file is FASTQ (vs FASTA) by peeking at the first content byte.
/// For gzipped files, decompresses just enough to read that byte.
fn is_fastq_file(path: &str) -> bool {
    let Ok(mut f) = std::fs::File::open(path) else {
        return false;
    };
    let mut magic = [0u8; 2];
    let Ok(n) = f.read(&mut magic) else {
        return false;
    };
    if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        // Gzipped — reopen and decompress to read first content byte
        let Ok(f2) = std::fs::File::open(path) else {
            return false;
        };
        let mut gz = GzDecoder::new(f2);
        let mut first = [0u8; 1];
        matches!(gz.read(&mut first), Ok(1) if first[0] == b'@')
    } else if n >= 1 {
        magic[0] == b'@'
    } else {
        false
    }
}

/// Check that a FASTQ sequence contains only valid DNA/RNA/IUPAC characters.
/// Panics with a descriptive message if an invalid character is found.
#[inline]
fn validate_fastq_seq(seq: &[u8], header: &[u8]) {
    for (i, &b) in seq.iter().enumerate() {
        if !matches!(
            b | 0x20,
            b'a' | b'c'
                | b'g'
                | b't'
                | b'u'
                | b'r'
                | b'y'
                | b's'
                | b'w'
                | b'k'
                | b'm'
                | b'b'
                | b'd'
                | b'h'
                | b'v'
                | b'n'
        ) {
            let name = std::str::from_utf8(header).unwrap_or("<non-utf8>");
            panic!(
                "Invalid character '{}' (0x{:02x}) at position {} in FASTQ record '{}'. \
                 Only DNA/RNA/IUPAC characters are allowed.",
                b as char, b, i, name
            );
        }
    }
}

/// Process a FASTQ file using the fastq crate's zero-copy callback API.
fn process_fastq_file(
    path: &str,
    caller: &orfm::OrfCaller,
    stop_codon_only: bool,
    print_stop_codons: bool,
    out: &mut BufWriter<io::StdoutLock>,
    transcript_out: &mut Option<BufWriter<std::fs::File>>,
    transcript_path: &Option<String>,
) {
    // FASTQ sequences are always single-line, so no newline normalization needed.
    // The fastq crate auto-detects gzip/lz4 compression.
    fastq::parse_path(Some(path), |parser| {
        let mut first = true;
        parser
            .each(|record| {
                let (name, comment) = split_header(record.head());
                let seq = record.seq();
                if first {
                    validate_fastq_seq(seq, record.head());
                    first = false;
                }
                let orfs = caller.find_orfs(name, comment, seq);
                for orf in orfs {
                    if stop_codon_only && !orf.has_stop_codon {
                        continue;
                    }
                    if let Some(ref mut tw) = transcript_out {
                        writeln!(tw, ">{}", orf.name()).unwrap();
                        let transcript = orf.transcript(seq);
                        tw.write_all(&transcript).unwrap();
                        writeln!(tw).unwrap();
                    }
                    writeln!(out, ">{}", orf.name()).unwrap();
                    out.write_all(&orf.protein).unwrap();
                    if print_stop_codons && orf.has_stop_codon {
                        out.write_all(b"*").unwrap();
                    }
                    writeln!(out).unwrap();
                }
                true
            })
            .unwrap_or_else(|e| {
                eprintln!("Error parsing FASTQ file '{}': {}", path, e);
                std::process::exit(1);
            });
    })
    .unwrap_or_else(|e| {
        eprintln!("Cannot open file '{}': {}", path, e);
        std::process::exit(5);
    });
    // Ensure transcript file is opened even if no records processed
    let _ = transcript_path;
}

/// Process a FASTQ stream from a reader using the fastq crate.
fn process_fastq_reader(
    reader: impl Read,
    caller: &orfm::OrfCaller,
    stop_codon_only: bool,
    print_stop_codons: bool,
    out: &mut BufWriter<io::StdoutLock>,
    transcript_out: &mut Option<BufWriter<std::fs::File>>,
) {
    let parser = fastq::Parser::new(reader);
    let mut first = true;
    parser
        .each(|record| {
            let (name, comment) = split_header(record.head());
            let seq = record.seq();
            if first {
                validate_fastq_seq(seq, record.head());
                first = false;
            }
            let orfs = caller.find_orfs(name, comment, seq);
            for orf in orfs {
                if stop_codon_only && !orf.has_stop_codon {
                    continue;
                }
                if let Some(ref mut tw) = transcript_out {
                    writeln!(tw, ">{}", orf.name()).unwrap();
                    let transcript = orf.transcript(seq);
                    tw.write_all(&transcript).unwrap();
                    writeln!(tw).unwrap();
                }
                writeln!(out, ">{}", orf.name()).unwrap();
                out.write_all(&orf.protein).unwrap();
                if print_stop_codons && orf.has_stop_codon {
                    out.write_all(b"*").unwrap();
                }
                writeln!(out).unwrap();
            }
            true
        })
        .unwrap_or_else(|e| {
            eprintln!("Error parsing FASTQ from stdin: {}", e);
            std::process::exit(1);
        });
}

/// Process FASTA/FASTQ using needletail (handles wrapped FASTA, all formats).
fn process_with_needletail(
    mut reader: Box<dyn needletail::FastxReader>,
    caller: &orfm::OrfCaller,
    stop_codon_only: bool,
    print_stop_codons: bool,
    out: &mut BufWriter<io::StdoutLock>,
    transcript_out: &mut Option<BufWriter<std::fs::File>>,
) {
    while let Some(result) = reader.next() {
        let record = result.unwrap_or_else(|e| {
            eprintln!("Error reading record: {e}");
            std::process::exit(1);
        });
        let (name, comment) = split_header(record.id());
        let seq = record.seq();
        let orfs = caller.find_orfs(name, comment, &seq);
        for orf in orfs {
            if stop_codon_only && !orf.has_stop_codon {
                continue;
            }
            if let Some(ref mut tw) = transcript_out {
                writeln!(tw, ">{}", orf.name()).unwrap();
                let transcript = orf.transcript(&seq);
                tw.write_all(&transcript).unwrap();
                writeln!(tw).unwrap();
            }
            writeln!(out, ">{}", orf.name()).unwrap();
            out.write_all(&orf.protein).unwrap();
            if print_stop_codons && orf.has_stop_codon {
                out.write_all(b"*").unwrap();
            }
            writeln!(out).unwrap();
        }
    }
}

/// Read from a non-seekable source: peek at the first byte to detect format, then process.
fn process_from_reader(
    reader: impl Read + Send + 'static,
    caller: &orfm::OrfCaller,
    stop_codon_only: bool,
    print_stop_codons: bool,
    out: &mut BufWriter<io::StdoutLock>,
    transcript_out: &mut Option<BufWriter<std::fs::File>>,
) {
    let mut first_byte = [0u8; 1];
    let mut reader = reader;
    let n = reader.read(&mut first_byte).unwrap_or(0);
    if n == 1 && first_byte[0] == b'@' {
        // FASTQ — use fastq crate
        let chained = io::Cursor::new(first_byte.to_vec()).chain(reader);
        process_fastq_reader(
            chained,
            caller,
            stop_codon_only,
            print_stop_codons,
            out,
            transcript_out,
        );
    } else {
        // FASTA (or gzipped) — use needletail
        let chained = io::Cursor::new(first_byte[..n].to_vec()).chain(reader);
        let nr = needletail::parse_fastx_reader(chained).unwrap_or_else(|e| {
            eprintln!("Cannot read input: {}", e);
            std::process::exit(5);
        });
        process_with_needletail(
            nr,
            caller,
            stop_codon_only,
            print_stop_codons,
            out,
            transcript_out,
        );
    }
}

fn main() {
    let args = Args::parse();

    if let Some(ref required) = args.required_version {
        let current = env!("CARGO_PKG_VERSION");
        let parse = |s: &str| -> (u64, u64, u64) {
            let mut p = s.split('.').filter_map(|x| x.parse::<u64>().ok());
            (
                p.next().unwrap_or(0),
                p.next().unwrap_or(0),
                p.next().unwrap_or(0),
            )
        };
        if parse(current) < parse(required) {
            eprintln!(
                "ERROR: OrfM version {} is less than required version {}",
                current, required
            );
            std::process::exit(1);
        }
    }

    let caller = match orfm::OrfCaller::new(args.codon_table, args.min_length, args.position_limit)
    {
        Ok(c) => c,
        Err(e) => {
            eprintln!("ERROR: {}", e);
            std::process::exit(1);
        }
    };

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let mut transcript_out: Option<BufWriter<std::fs::File>> =
        args.transcript.as_ref().map(|path| {
            let f = std::fs::File::create(path).unwrap_or_else(|e| {
                eprintln!("Cannot open output transcript file '{}': {}", path, e);
                std::process::exit(5);
            });
            BufWriter::new(f)
        });

    match &args.input {
        Some(path) => {
            // Check if the path is a regular file (not a FIFO/pipe from process substitution)
            let metadata = std::fs::metadata(path);
            let is_regular_file = metadata.as_ref().map_or(false, |m| m.is_file());

            if is_regular_file && is_fastq_file(path) {
                // Use the fastq crate's zero-copy parser for FASTQ input
                process_fastq_file(
                    path,
                    &caller,
                    args.stop_codon_only,
                    args.print_stop_codons,
                    &mut out,
                    &mut transcript_out,
                    &args.transcript,
                );
            } else if is_regular_file {
                // Use needletail for FASTA input (handles wrapped sequences)
                let reader = needletail::parse_fastx_file(path).unwrap_or_else(|e| {
                    eprintln!("Cannot open file '{}': {}", path, e);
                    std::process::exit(5);
                });
                process_with_needletail(
                    reader,
                    &caller,
                    args.stop_codon_only,
                    args.print_stop_codons,
                    &mut out,
                    &mut transcript_out,
                );
            } else {
                // Non-seekable input (FIFO, /dev/fd/*, etc.) — open and detect format from first byte
                let file = std::fs::File::open(path).unwrap_or_else(|e| {
                    eprintln!("Cannot open file '{}': {}", path, e);
                    std::process::exit(5);
                });
                process_from_reader(
                    file,
                    &caller,
                    args.stop_codon_only,
                    args.print_stop_codons,
                    &mut out,
                    &mut transcript_out,
                );
            }
        }
        None => {
            // stdin: detect format from first byte
            let stdin = io::stdin();
            process_from_reader(
                stdin,
                &caller,
                args.stop_codon_only,
                args.print_stop_codons,
                &mut out,
                &mut transcript_out,
            );
        }
    }
}
