# OrfM

A simple and not slow open reading frame (ORF) caller. No bells or whistles like frameshift detection, just a straightforward goal 
of returning a FASTA file of open reading frames over a certain length from a FASTA/Q file of nucleotide sequences.

As of version 2.0, it is a pure-Rust reimplementation of the original C OrfM, which is no longer maintained. The algorithm is the same as the original, but the codebase has been rewritten from scratch in Rust, with a library API.

## Install

OrfM can be installed in different ways:

### 1) Install from pre-compiled binaries

OrfM can be installed by downloading pre-compiled binaries available at https://github.com/wwood/OrfM/releases. Once you have downloaded the package, extract and run it e.g. for GNU/Linux:
```sh
tar xzf orfm-<version>.tar.gz
cd orfm-<version>
./orfm -h
```

### 2) Install from source (requires Rust toolchain)

```bash
cargo install orfm
```

Or build without installing:

```bash
cargo build --release
# binary at target/release/orfm
```

### 3) Install via Conda / Pixi

```bash
conda install -c bioconda orfm
```

### 4) Install with brew
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

A Snakefile is included for comparing performance against the original C OrfM. It generates random sequences, runs both tools, checks output correctness, and collects wall-clock time and memory usage.

```bash
snakemake -j4
cat benchmark/results.tsv
cat benchmark/correctness.txt
```
When I ran it, orfm-rs was the winner by ~5% in walltime:

| tool | replicate | wall_clock_s | max_rss_kb |
|------|-----------|--------------|------------|
| orfm_c | 1 | 2.11 | 16088 |
| orfm_rs | 1 | 1.99 | 15572 |
| orfm_c | 2 | 2.1 | 15432 |
| orfm_rs | 2 | 2.0 | 16116 |
| orfm_c | 3 | 2.12 | 16164 |
| orfm_rs | 3 | 2.02 | 16148 |

Requires the C OrfM binary at `~/git/OrfM/orfm` (this path can be changed in the Snakefile).

## Differences from OrfM

- Pure Rust, no C dependencies
- Exposes a library API with an iterator over translated ORFs
- Supports all codon tables that OrfM supports (NCBI tables 1–25)
- Uses [needletail](https://github.com/onecodex/needletail) for sequence parsing and [aho-corasick](https://github.com/BurntSushi/aho-corasick) for stop codon detection

## License

OrfM is licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPL-3.0).

## Citation
Since the algorithm was devised for the original version of OrfM, best to cite that:

Ben J. Woodcroft, Joel A. Boyd, and Gene W. Tyson. [_OrfM: A fast open reading frame predictor for metagenomic data_](http://bioinformatics.oxfordjournals.org/content/32/17/2702). (2016). Bioinformatics. doi:10.1093/bioinformatics/btw241.
