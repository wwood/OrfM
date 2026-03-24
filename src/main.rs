use clap::Parser;
use std::io::{self, BufWriter, Write};

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
}

fn split_header(header: &[u8]) -> (&str, &str) {
    let header_str = std::str::from_utf8(header).unwrap_or("unknown");
    match header_str.find(|c: char| c.is_whitespace()) {
        Some(pos) => (&header_str[..pos], &header_str[pos + 1..]),
        None => (header_str, ""),
    }
}

fn main() {
    let args = Args::parse();

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

    let mut reader: Box<dyn needletail::FastxReader> = match &args.input {
        Some(path) => needletail::parse_fastx_file(path).unwrap_or_else(|e| {
            eprintln!("Cannot open file '{}': {}", path, e);
            std::process::exit(5);
        }),
        None => needletail::parse_fastx_reader(io::stdin()).unwrap_or_else(|e| {
            eprintln!("Cannot read stdin: {}", e);
            std::process::exit(5);
        }),
    };

    while let Some(Ok(record)) = reader.next() {
        let (name, comment) = split_header(record.id());
        let seq = record.raw_seq();
        let orfs = caller.find_orfs(name, comment, seq);
        for orf in orfs {
            if let Some(ref mut tw) = transcript_out {
                writeln!(tw, "{}", orf.header()).unwrap();
                let transcript = orf.transcript(seq);
                tw.write_all(&transcript).unwrap();
                writeln!(tw).unwrap();
            }
            writeln!(out, "{}", orf.header()).unwrap();
            out.write_all(&orf.protein).unwrap();
            writeln!(out).unwrap();
        }
    }
}
