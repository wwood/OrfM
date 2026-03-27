# Changelog

## [2.1.0] - 2026-03-27

Now ~25% faster than the old C implementation of OrfM (v1.4). Cleaned up API, fixed multi-line FASTA bug.

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
- `Orf::orfm_id()` method returns the ORF identifier
  (`{source_seq_name}_{start}_{frame}_{orf_number}`).
- `Orf::name()` method returns the full ORF name (orfm_id plus comment if present).

### Fixed
- `revcom_base` now handles all IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V)
  and RNA uracil (U). Previously these were all mapped to N, losing information in
  reverse-strand transcripts.
- Parse errors in FASTA/FASTQ input now produce an error message and exit instead of
  silently truncating output. Previously `while let Some(Ok(...))` in the iterators
  and `process_with_needletail` would treat any `Err` as end-of-file.

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
- `Orf.seq_name` renamed to `Orf.source_seq_name` for clarity. The `header()` method
  has been replaced by `orfm_id()` and `name()`.
- `Orf.seq_name` and `Orf.seq_comment` changed from `String` to `Arc<str>`, reducing
  per-ORF heap allocations to a cheap atomic reference count increment.
- Library iterators (`OrfFileIter`, `OrfReaderIter`) no longer clone each `Orf` on
  iteration; they yield owned values directly via `IntoIter`.
- Library API cleaned up: `version_at_least`, `get_codon_table`, `revcom_base`,
  and `is_valid_iupac` are no longer public.
- Benchmark data files moved to `benchmark/data/`.
- Benchmark summary now prints a box-drawing table with C vs Rust wall-clock
  times and ratio (median of 3 replicates).

### Removed
- `-r VERSION` flag and `version_at_least` function (dead code).
- `Orf::header()` method — use `orfm_id()` or `name()` instead.

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
