# Changelog

## [2.0.3] - 2026-03-25

### Fixed
- Correctly handle wrapped (multi-line) FASTA sequences. Previously `raw_seq()` was
  used which includes embedded newlines, producing incorrect ORFs for wrapped input.

### Changed
- Aho-Corasick automaton is now case-insensitive, removing a per-sequence uppercase
  copy allocation.
- Newline detection is integrated into the Aho-Corasick scan itself: `raw_seq()` is
  used optimistically (zero-copy), and if a newline is encountered during the stop
  codon search, the sequence is re-processed via `seq()`. Once a newline is detected,
  all subsequent records in the file use `seq()` directly.
- `find_orfs()` now returns `Option<Vec<Orf>>` — `None` indicates the input contained
  embedded newlines and should be retried with a normalized sequence.

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
