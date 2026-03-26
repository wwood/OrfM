# OrfM

A simple and not slow open reading frame (ORF) caller. No bells or whistles like frameshift detection, just a straightforward goal 
of returning a FASTA file of open reading frames over a certain length from a FASTA/Q file of nucleotide sequences.

As of version 2.0, it is a pure-Rust reimplementation of the original C OrfM. The algorithm is the same as the original, but the codebase has been rewritten from scratch in Rust, with a library API added.

## Install

OrfM can be installed in different ways:

### 1) Install via shell pipe
Follow the instructions on the [releases page](https://github.com/wwood/OrfM/releases) e.g. for Linux:
```
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/wwood/OrfM/releases/download/v<version>/orfm-installer.sh | sh
```

### 2) Install from pre-compiled binaries

OrfM can be installed by downloading pre-compiled binaries available at https://github.com/wwood/OrfM/releases. Once you have downloaded the package, extract and run it e.g. for GNU/Linux:
```sh
tar xzf orfm-<version>.tar.gz
cd orfm-<version>
./orfm -h
```

### 3) Install from source (requires Rust toolchain)

```bash
cargo install orfm
```

Or build without installing:

```bash
cargo build --release
# binary at target/release/orfm
```

### 4) Install via Conda / Pixi

```bash
conda install -c bioconda orfm
```

### 5) Install with brew
Thanks to Torsten Seemann (@tseemann), OrfM can be installed through homebrew:
```
brew install brewsci/bio/orfm
```

## Usage

```
orfm [OPTIONS] [INPUT]
```

Reads from stdin if no input file is given. Accepts FASTA or FASTQ, gzipped or uncompressed.

### Options

| Flag | Description | Default |
|---|---|---|
| `-m <LENGTH>` | Minimum ORF length in nucleotides (must be a multiple of 3) | 96 |
| `-c <TABLE_ID>` | NCBI codon table for translation (1–25) | 1 |
| `-l <LENGTH>` | Ignore sequence beyond this position | none |
| `-t <FILE>` | Write nucleotide transcripts to this file | none |
| `-p` | Append `*` to proteins whose ORF is terminated by an in-frame stop codon | off |
| `-s` | Only output ORFs that are terminated by an in-frame stop codon (suppress terminal ORFs) | off |
| `-r <VERSION>` | Exit with an error if the running OrfM version is older than `VERSION` (e.g. `2.0.2`) | none |

### Examples

```bash
# Basic usage
orfm input.fasta > orfs.faa

# From gzipped FASTQ, shorter minimum ORF length
orfm -m 30 reads.fastq.gz > orfs.faa

# Pipe from stdin, write transcripts
cat input.fasta | orfm -m 60 -t transcripts.fna > orfs.faa

# Use mitochondrial codon table
orfm -c 2 mito.fasta > orfs.faa
```

### Output
The output ORFs fasta file contains any stretch of continuous codons which does not include a stop codon. 
There is no requirement for a start codon to be included in the ORF. One could say that OrfM is an ORF caller, not a gene caller (like say prodigal or genscan).

The output ORFs are named in a straitforward manner. The name of the sequence (i.e. anything before a space) is followed by `_startPosition_frameNumber_orfNumber` and then 
the comment of the sequence (i.e. anything after the space) is given after a space, if one exists. For example,
```
$ cat eg.fasta
>abc|123|name some comment
ATGTTA
$ orfm -m 3 eg.fasta
>abc|123|name_1_1_1 some comment
ML
```
The `startPosition` of reverse frames is the left-most position in the original sequence, not the codon where the ORF starts.

## Library usage

orfm can be used as a Rust library. Add to your `Cargo.toml`:

```toml
[dependencies]
orfm = '*'
```
Then to use it:
```rust
use orfm::OrfCaller;

let caller = OrfCaller::new(1, 96, None).unwrap(); // table_id, min_length, position_limit

// Iterate over ORFs from a file
for orf in caller.call_from_file("input.fasta") {
    println!("{}", orf.header());
    println!("{}", std::str::from_utf8(&orf.protein).unwrap());
}

// Or call on a single sequence
let orfs = caller.find_orfs("seq1", "", b"ATGGATGCTGAA...");
for orf in &orfs {
    let transcript = orf.transcript(b"ATGGATGCTGAA...");
    // ...
}
```

## Benchmarking

Walltime (seconds) results on a Linux x86-64 server, on 1 million random 150 bp sequences:

```
┌─────────────────┬───────────────┬──────────────────┬────────┬┬────────────────────────┐
│      Input      │ OrfM v1.4 (C) │ OrfM v2.1 (Rust) │ getorf ││ OrfM (Rust) / OrfM (C) │
├─────────────────┼───────────────┼──────────────────┼────────┼┼────────────────────────┤
│ fasta_unwrapped │          2.14 │             1.67 │   8.77 ││                  0.78x │
├─────────────────┼───────────────┼──────────────────┼────────┼┼────────────────────────┤
│ fasta_wrapped   │          2.14 │             1.72 │   8.86 ││                  0.81x │
├─────────────────┼───────────────┼──────────────────┼────────┼┼────────────────────────┤
│ fasta_gzipped   │          2.57 │             1.92 │    N/A ││                  0.74x │
├─────────────────┼───────────────┼──────────────────┼────────┼┼────────────────────────┤
│ fastq           │          2.28 │             1.71 │   9.05 ││                  0.75x │
├─────────────────┼───────────────┼──────────────────┼────────┼┼────────────────────────┤
│ fastq_gzipped   │          2.75 │             1.67 │    N/A ││                  0.61x │
└─────────────────┴───────────────┴──────────────────┴────────┴┴────────────────────────┘
```

- A ratio < 1 means OrfM v2.1 (Rust) is faster than OrfM v1.4 (C); Rust is 19–39% faster depending on input type.
- getorf (EMBOSS) does not support gzipped input (N/A). On plain FASTA/FASTQ it is ~4–5× slower than OrfM v2.1 (Rust).
- Peak RSS memory usage is similar across all programs (~85 MB).
- All replicates produce identical output (verified by `diff`).

A Snakefile is included for comparing performance against the original C OrfM and getorf (EMBOSS). It generates 1 million random 150 bp sequences in various formats, runs all three tools (3 replicates), checks output correctness, and collects wall-clock time and peak RSS memory.

```bash
pixi run snakemake -j1
cat benchmark/results.tsv
cat benchmark/correctness.txt
```

Requires the C OrfM binary at `~/git/OrfM/orfm` (path configurable in the Snakefile). getorf is installed automatically via the pixi environment (`emboss` package).

## Notes

- Exposes a Rust library API with an iterator over translated ORFs
- Supports all codon tables that OrfM supports (NCBI tables 1–25)
- Uses [needletail](https://github.com/onecodex/needletail) for sequence parsing and [aho-corasick](https://github.com/BurntSushi/aho-corasick) for stop codon detection

## License

OrfM is licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPL-3.0).

## Citation
Since the algorithm was devised for the original version of OrfM, best to cite that:

Ben J. Woodcroft, Joel A. Boyd, and Gene W. Tyson. [_OrfM: A fast open reading frame predictor for metagenomic data_](http://bioinformatics.oxfordjournals.org/content/32/17/2702). (2016). Bioinformatics. doi:10.1093/bioinformatics/btw241.
