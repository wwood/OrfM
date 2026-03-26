# Changelog

## [2.1.0] - 2026-03-26

### Added
- Correctly handle wrapped (multi-line) FASTA sequences. Previously `raw_seq()` was
  used which includes embedded newlines, producing incorrect ORFs for wrapped input.
- FASTQ validation: the first record is checked character-by-character for
  valid IUPAC DNA/RNA characters; subsequent records detect invalid characters
  via the translation table (a codon containing a non-IUPAC byte produces the
  sentinel value `0`, triggering a panic with position and character details).
- `flate2` switched to the `zlib-ng` backend for all gzip decompression
  (needletail and the fastq crate), roughly halving wall-clock time on
  gzipped inputs compared to the previous `miniz_oxide` default.

### Changed
- FASTQ parsing switched from needletail to the `fastq` crate (zero-copy,
  callback-based). Format is auto-detected by inspecting the first byte of input;
  FASTA input continues to use needletail.
- Aho-Corasick automaton is now case-insensitive, removing a per-sequence uppercase
  copy allocation.
- Newline detection is integrated into the Aho-Corasick scan itself: `raw_seq()` is
  used optimistically (zero-copy), and if a newline is encountered during the stop
  codon search, the sequence is re-processed via `seq()`. Once a newline is detected,
  all subsequent records in the file use `seq()` directly.
- `find_orfs` is now the public single-sequence entry point and always returns
  `Vec<Orf>` (normalizing wrapped input automatically). See the `find_orfs` docs
  for the performance trade-off vs. the iterator API.
- Library API cleaned up: `version_at_least`, `get_codon_table`, `revcom_base`,
  and `is_valid_iupac` are no longer public.
- Benchmark data files moved to `benchmark/data/`.
- Benchmark summary now prints a box-drawing table with C vs Rust wall-clock
  times and ratio (median of 3 replicates).

## [2.0.2] - 2026-03-25

### Added
- `-p` flag: append `*` to proteins whose ORF is bounded by an in-frame stop codon
- `-s` flag: suppress ORFs that are not bounded by a stop codon (terminal ORFs)
- `-r VERSION` flag: exit with an error if the running OrfM is older than the required version
- `version_at_least(current, required)` public library function for semver comparison

### Breaking changes
- `Orf` struct gained a new public field `has_stop_codon: bool`. Any code constructing
  `Orf` values directly (e.g. in tests or downstream libraries) must set this field.

## [2.0.1] - 2026-03-24

### Added
- `--version` flag

## [2.0.0] - 2026-03-23

### Added
- Complete rewrite in Rust (algorithm unchanged from original C OrfM)
- Library API
