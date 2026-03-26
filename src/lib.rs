use aho_corasick::AhoCorasick;
use std::io;

mod codon_tables;
use codon_tables::get_codon_table;

/// Compare two semver-style version strings (MAJOR.MINOR.PATCH).
/// Returns true if `current` >= `required`.
fn version_at_least(current: &str, required: &str) -> bool {
    let parse = |s: &str| -> (u64, u64, u64) {
        let mut parts = s.split('.').filter_map(|p| p.parse::<u64>().ok());
        (
            parts.next().unwrap_or(0),
            parts.next().unwrap_or(0),
            parts.next().unwrap_or(0),
        )
    };
    parse(current) >= parse(required)
}

/// A single open reading frame found in a sequence.
#[derive(Debug, Clone)]
pub struct Orf {
    /// Name of the source sequence (before first whitespace in header).
    pub seq_name: String,
    /// Comment from the source sequence header (after first whitespace, may be empty).
    pub seq_comment: String,
    /// 0-based start position in the original nucleotide sequence.
    pub start: usize,
    /// Frame number (1-6). 1-3 are forward, 4-6 are reverse complement.
    pub frame: u8,
    /// 1-based ORF counter within this sequence.
    pub orf_number: usize,
    /// Translated amino acid sequence.
    pub protein: Vec<u8>,
    /// Whether this ORF is on the reverse complement strand.
    pub is_reverse: bool,
    /// Length of the nucleotide ORF region.
    pub nuc_length: usize,
    /// Whether this ORF is bounded by an in-frame stop codon (vs. running to sequence end).
    pub has_stop_codon: bool,
}

impl Orf {
    /// Format the FASTA header for this ORF.
    pub fn header(&self) -> String {
        if self.seq_comment.is_empty() {
            format!(
                ">{}_{}_{}_{}",
                self.seq_name,
                self.start + 1,
                self.frame,
                self.orf_number
            )
        } else {
            format!(
                ">{}_{}_{}_{} {}",
                self.seq_name,
                self.start + 1,
                self.frame,
                self.orf_number,
                self.seq_comment
            )
        }
    }

    /// Get the nucleotide transcript for this ORF from the source sequence.
    pub fn transcript(&self, seq: &[u8]) -> Vec<u8> {
        let region = &seq[self.start..self.start + self.nuc_length];
        if self.is_reverse {
            region.iter().rev().map(|&b| revcom_base(b)).collect()
        } else {
            region.to_vec()
        }
    }
}

/// Configuration for ORF calling.
#[derive(Debug, Clone)]
pub struct OrfCaller {
    fast_table: Vec<u8>,
    ac: AhoCorasick,
    num_fwd_stop_codons: usize,
    num_stop_codons: usize,
    min_length: usize,
    position_limit: Option<usize>,
}

/// Reverse complement a single base.
#[inline]
fn revcom_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'N' => b'N',
        b'n' => b'n',
        _ => b'N',
    }
}

/// Split a needletail header (id()) into (name, comment).
/// kseq.h splits on the first whitespace character.
fn split_header(header: &[u8]) -> (&str, &str) {
    let header_str = std::str::from_utf8(header).unwrap_or("unknown");
    match header_str.find(|c: char| c.is_whitespace()) {
        Some(pos) => (&header_str[..pos], header_str[pos + 1..].trim_start()),
        None => (header_str, ""),
    }
}

/// Build a fast 3-byte-indexed translation lookup table.
/// Index: (b0 << 16) | (b1 << 8) | b2 for forward codons
/// Index: 128 + (b0 << 16) | (b1 << 8) | b2 for reverse complement lookup
/// All valid IUPAC DNA/RNA characters (both cases).
/// Codons containing only these bytes produce b'X' (unknown amino acid).
/// Any byte outside this set leaves the table entry at 0, triggering a panic on translation.
const IUPAC_CHARS: &[u8] = b"ACGTURYWSKMBDHVNacgturywskmbdhvn";

/// Return true if `b` is a valid IUPAC DNA/RNA character.
#[inline]
fn is_valid_iupac(b: u8) -> bool {
    matches!(
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
    )
}

fn build_fast_table(codon_table: &[u8; 64]) -> Vec<u8> {
    let table_size = (1usize << 24) + 128;
    // 0 is the sentinel for "invalid character" — never a valid amino acid output.
    let mut table = vec![0u8; table_size];

    let order = [b'A', b'C', b'G', b'T'];
    let low_order = [b'a', b'c', b'g', b't'];
    let complements = [b'T', b'G', b'C', b'A'];
    let low_complements = [b't', b'g', b'c', b'a'];

    let mut count = 0usize;
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                let aa = codon_table[count];

                // All 8 case combinations for forward
                for &bi in &[order[i], low_order[i]] {
                    for &bj in &[order[j], low_order[j]] {
                        for &bk in &[order[k], low_order[k]] {
                            let idx = ((bi as usize) << 16) | ((bj as usize) << 8) | (bk as usize);
                            table[idx] = aa;
                        }
                    }
                }

                // All 8 case combinations for reverse complement
                for &bi in &[complements[i], low_complements[i]] {
                    for &bj in &[complements[j], low_complements[j]] {
                        for &bk in &[complements[k], low_complements[k]] {
                            let idx = 128
                                + (((bi as usize) << 16) | ((bj as usize) << 8) | (bk as usize));
                            table[idx] = aa;
                        }
                    }
                }

                count += 1;
            }
        }
    }

    // Fill all IUPAC ambiguity codon combinations with b'X' (unknown amino acid).
    // This must run after the ACGT loop so pure-ACGT entries keep their real amino acid.
    // Any codon index not reachable by IUPAC bytes stays at 0 (invalid char sentinel).
    for &b0 in IUPAC_CHARS {
        for &b1 in IUPAC_CHARS {
            for &b2 in IUPAC_CHARS {
                let fwd = ((b0 as usize) << 16) | ((b1 as usize) << 8) | (b2 as usize);
                if table[fwd] == 0 {
                    table[fwd] = b'X';
                }
                let rev = 128 + (((b0 as usize) << 16) | ((b1 as usize) << 8) | (b2 as usize));
                if table[rev] == 0 {
                    table[rev] = b'X';
                }
            }
        }
    }

    table
}

/// Translate a nucleotide sequence to amino acids using the fast lookup table.
#[inline]
fn translate_forward(seq: &[u8], start: usize, nuc_len: usize, table: &[u8]) -> Vec<u8> {
    let num_codons = nuc_len / 3;
    let mut protein = Vec::with_capacity(num_codons);
    let mut pos = start;
    for _ in 0..num_codons {
        let idx =
            ((seq[pos] as usize) << 16) | ((seq[pos + 1] as usize) << 8) | (seq[pos + 2] as usize);
        let aa = table[idx];
        if aa == 0 {
            let bad = [seq[pos], seq[pos + 1], seq[pos + 2]]
                .into_iter()
                .enumerate()
                .find(|(_, b)| !is_valid_iupac(*b))
                .map(|(i, b)| (pos + i, b))
                .unwrap_or((pos, seq[pos]));
            panic!(
                "Invalid character '{}' (0x{:02x}) at sequence position {}",
                bad.1 as char, bad.1, bad.0
            );
        }
        protein.push(aa);
        pos += 3;
    }
    protein
}

/// Translate the reverse complement of a nucleotide sequence.
#[inline]
fn translate_reverse(seq: &[u8], start: usize, nuc_len: usize, table: &[u8]) -> Vec<u8> {
    let num_codons = nuc_len / 3;
    let mut protein = Vec::with_capacity(num_codons);
    let mut pos = start + nuc_len - 3;
    for _ in 0..num_codons {
        let idx = 128
            + (((seq[pos + 2] as usize) << 16)
                | ((seq[pos + 1] as usize) << 8)
                | (seq[pos] as usize));
        let aa = table[idx];
        if aa == 0 {
            let bad = [seq[pos], seq[pos + 1], seq[pos + 2]]
                .into_iter()
                .enumerate()
                .find(|(_, b)| !is_valid_iupac(*b))
                .map(|(i, b)| (pos + i, b))
                .unwrap_or((pos, seq[pos]));
            panic!(
                "Invalid character '{}' (0x{:02x}) at sequence position {}",
                bad.1 as char, bad.1, bad.0
            );
        }
        protein.push(aa);
        if pos < start + 3 {
            break;
        }
        pos -= 3;
    }
    protein
}

impl OrfCaller {
    /// Create a new ORF caller with the given codon table ID and minimum nucleotide length.
    pub fn new(
        table_id: usize,
        min_length: usize,
        position_limit: Option<usize>,
    ) -> Result<Self, String> {
        let codon_table = get_codon_table(table_id)
            .ok_or_else(|| format!("Invalid translation table ID: {}", table_id))?;

        if min_length < 3 {
            return Err("Minimum length must be at least 3".to_string());
        }
        if !min_length.is_multiple_of(3) {
            return Err("Minimum length must be a multiple of 3".to_string());
        }
        if let Some(pl) = position_limit {
            if pl < 3 {
                return Err("Position limit must be at least 3".to_string());
            }
            if min_length > pl {
                return Err("Position limit cannot be less than minimum length".to_string());
            }
        }

        let fast_table = build_fast_table(codon_table);

        // Build stop codon patterns for Aho-Corasick
        let order = [b'A', b'C', b'G', b'T'];
        let mut patterns: Vec<Vec<u8>> = Vec::new();

        // Forward stop codons
        let mut count = 0;
        for &b0 in &order {
            for &b1 in &order {
                for &b2 in &order {
                    if codon_table[count] == b'*' {
                        patterns.push(vec![b0, b1, b2]);
                    }
                    count += 1;
                }
            }
        }
        let num_fwd = patterns.len();

        // Reverse complement stop codons
        count = 0;
        for &b0 in &order {
            for &b1 in &order {
                for &b2 in &order {
                    if codon_table[count] == b'*' {
                        patterns.push(vec![revcom_base(b2), revcom_base(b1), revcom_base(b0)]);
                    }
                    count += 1;
                }
            }
        }

        let num_stop = patterns.len();

        // Newline sentinel patterns (detected during AC scan to trigger normalization)
        patterns.push(vec![b'\n']);
        patterns.push(vec![b'\r']);

        let ac = AhoCorasick::builder()
            .ascii_case_insensitive(true)
            .build(&patterns)
            .map_err(|e| format!("Failed to build Aho-Corasick automaton: {}", e))?;

        Ok(OrfCaller {
            fast_table,
            ac,
            num_fwd_stop_codons: num_fwd,
            num_stop_codons: num_stop,
            min_length,
            position_limit,
        })
    }

    /// Find all ORFs in a pre-normalized sequence (no embedded newlines).
    ///
    /// Returns `Some(orfs)` on success, or `None` if a newline was detected
    /// (caller should normalize and retry). Used internally by the iterators,
    /// which track normalization state across records for efficiency.
    fn find_orfs_internal(&self, name: &str, comment: &str, seq: &[u8]) -> Option<Vec<Orf>> {
        let seq_len = seq.len();
        if seq_len < self.min_length {
            return Some(Vec::new());
        }

        let read_limit = match self.position_limit {
            Some(pl) if pl < seq_len => pl,
            _ => seq_len,
        };

        let mut last_pos = [0usize; 6];
        last_pos[1] = 1;
        last_pos[2] = 2;
        last_pos[4] = 1;
        last_pos[5] = 2;

        // For reverse frames (indices 3-5): tracks whether the current left boundary of the
        // frame was set by an actual stop codon. On the reverse strand, reading 5'→3' goes
        // from high to low forward positions, so the 3' end of a reverse ORF is on the LEFT
        // (low forward position). A reverse ORF has a stop codon only when its left boundary
        // (last_pos[rev_idx]) was set by a previous stop, not by the initial sequence offset.
        let mut rev_left_is_stop = [false; 3];

        let mut orfs = Vec::new();
        let mut orf_counter = 1usize;

        // Search for stop codons using Aho-Corasick (case-insensitive)
        for mat in self.ac.find_overlapping_iter(&seq[..read_limit]) {
            let match_start = mat.start();
            let mod3 = match_start % 3;
            let pattern_id = mat.pattern().as_usize();

            // Newline sentinel: sequence contains embedded newlines, caller must normalize
            if pattern_id >= self.num_stop_codons {
                return None;
            }

            if pattern_id < self.num_fwd_stop_codons {
                // Forward stop codon: stop is at the 3' end of the ORF to its left.
                let orf_len = match_start - last_pos[mod3];
                if orf_len >= self.min_length {
                    let protein = translate_forward(seq, last_pos[mod3], orf_len, &self.fast_table);
                    orfs.push(Orf {
                        seq_name: name.to_string(),
                        seq_comment: comment.to_string(),
                        start: last_pos[mod3],
                        frame: (mod3 as u8) + 1,
                        orf_number: orf_counter,
                        protein,
                        is_reverse: false,
                        nuc_length: orf_len,
                        has_stop_codon: true,
                    });
                    orf_counter += 1;
                }
                last_pos[mod3] = match_start + 3;
            } else {
                // Reverse complement stop codon.
                // On the reverse strand (reading 5'→3' = high→low in forward coords), this
                // stop codon is at the 5' end of the upcoming ORF (to its LEFT in forward
                // coords). The ORF from last_pos[rev_idx]..match_start has its 3' end at
                // last_pos[rev_idx], which is the stop-codon end only if that boundary was
                // itself set by a previous stop.
                let rev_idx = mod3 + 3;
                let orf_len = match_start - last_pos[rev_idx];
                if orf_len >= self.min_length {
                    let protein =
                        translate_reverse(seq, last_pos[rev_idx], orf_len, &self.fast_table);
                    orfs.push(Orf {
                        seq_name: name.to_string(),
                        seq_comment: comment.to_string(),
                        start: last_pos[rev_idx],
                        frame: (mod3 as u8) + 4,
                        orf_number: orf_counter,
                        protein,
                        is_reverse: true,
                        nuc_length: orf_len,
                        has_stop_codon: rev_left_is_stop[mod3],
                    });
                    orf_counter += 1;
                }
                last_pos[rev_idx] = match_start + 3;
                rev_left_is_stop[mod3] = true;
            }
        }

        // Handle final ORFs in each frame (after last stop codon to end of sequence)
        let mod3 = read_limit % 3;
        // Frame offsets: how many trailing bases to subtract for each frame given read_limit % 3
        let offsets: [[usize; 6]; 3] = [
            [0, 2, 1, 0, 2, 1], // mod3 == 0
            [1, 0, 2, 1, 0, 2], // mod3 == 1
            [2, 1, 0, 2, 1, 0], // mod3 == 2
        ];
        let offs = &offsets[mod3];

        for frame_idx in 0..6 {
            let orf_len_raw = read_limit.saturating_sub(last_pos[frame_idx]);
            let orf_len = orf_len_raw.saturating_sub(offs[frame_idx]);
            // Make it a multiple of 3
            let orf_len = (orf_len / 3) * 3;
            if orf_len >= self.min_length {
                let is_reverse = frame_idx >= 3;
                let frame_num = (frame_idx as u8) + 1;
                let protein = if is_reverse {
                    translate_reverse(seq, last_pos[frame_idx], orf_len, &self.fast_table)
                } else {
                    translate_forward(seq, last_pos[frame_idx], orf_len, &self.fast_table)
                };
                // Forward terminal ORFs never have a stop at their 3' end.
                // Reverse terminal ORFs: their 3' end (left boundary in forward coords) has a
                // stop codon only if that boundary was set by a previous RC stop codon.
                let has_stop_codon = is_reverse && rev_left_is_stop[frame_idx - 3];
                orfs.push(Orf {
                    seq_name: name.to_string(),
                    seq_comment: comment.to_string(),
                    start: last_pos[frame_idx],
                    frame: frame_num,
                    orf_number: orf_counter,
                    protein,
                    is_reverse,
                    nuc_length: orf_len,
                    has_stop_codon,
                });
                orf_counter += 1;
            }
        }

        Some(orfs)
    }

    /// Find all ORFs in a sequence, normalizing it first if necessary.
    ///
    /// This is the convenient single-sequence entry point. It first attempts
    /// to process `seq` directly; if embedded newlines are detected (wrapped
    /// FASTA), it strips them and retries automatically.
    ///
    /// # Performance note
    ///
    /// When calling `find_orfs` in a loop over many sequences, each call
    /// independently detects whether the sequence needs newline-stripping.
    /// For unwrapped input this overhead is negligible (one fast memchr scan
    /// per sequence), but once a wrapped sequence has been seen all subsequent
    /// ones still pay that cost.
    ///
    /// For performance-critical multi-record processing, prefer
    /// [`OrfCaller::call_from_file`] or [`OrfCaller::call_from_reader`], which
    /// maintain a single `needs_normalize` flag across all records: after the
    /// first wrapped sequence is encountered, every subsequent record is
    /// normalized immediately without an extra scan.
    pub fn find_orfs(&self, name: &str, comment: &str, seq: &[u8]) -> Vec<Orf> {
        match self.find_orfs_internal(name, comment, seq) {
            Some(orfs) => orfs,
            None => {
                // Newline detected — strip and retry
                let normalized: Vec<u8> = seq
                    .iter()
                    .copied()
                    .filter(|&b| b != b'\n' && b != b'\r')
                    .collect();
                self.find_orfs_internal(name, comment, &normalized).unwrap()
            }
        }
    }

    /// Iterate over all ORFs from a FASTA/FASTQ file (or gzipped).
    pub fn call_from_file(&self, path: &str) -> OrfFileIter<'_> {
        OrfFileIter {
            caller: self,
            reader: needletail::parse_fastx_file(path).expect("Failed to open sequence file"),
            current_orfs: Vec::new(),
            current_idx: 0,
            needs_normalize: false,
        }
    }

    /// Iterate over all ORFs from a reader (e.g. stdin). Reader must be Send.
    pub fn call_from_reader<'a, R: io::Read + Send + 'a>(&'a self, reader: R) -> OrfReaderIter<'a> {
        OrfReaderIter {
            caller: self,
            reader: needletail::parse_fastx_reader(reader).expect("Failed to parse sequence input"),
            current_orfs: Vec::new(),
            current_idx: 0,
            needs_normalize: false,
        }
    }
}

/// Iterator over ORFs from a file.
pub struct OrfFileIter<'a> {
    caller: &'a OrfCaller,
    reader: Box<dyn needletail::FastxReader>,
    current_orfs: Vec<Orf>,
    current_idx: usize,
    /// Once true, all subsequent records use seq() instead of raw_seq().
    needs_normalize: bool,
}

impl Iterator for OrfFileIter<'_> {
    type Item = Orf;

    fn next(&mut self) -> Option<Orf> {
        loop {
            if self.current_idx < self.current_orfs.len() {
                let orf = self.current_orfs[self.current_idx].clone();
                self.current_idx += 1;
                return Some(orf);
            }

            match self.reader.next() {
                Some(Ok(record)) => {
                    let (name, comment) = split_header(record.id());
                    if self.needs_normalize {
                        let seq = record.seq();
                        self.current_orfs =
                            self.caller.find_orfs_internal(name, comment, &seq).unwrap();
                    } else {
                        match self
                            .caller
                            .find_orfs_internal(name, comment, record.raw_seq())
                        {
                            Some(orfs) => self.current_orfs = orfs,
                            None => {
                                self.needs_normalize = true;
                                let seq = record.seq();
                                self.current_orfs =
                                    self.caller.find_orfs_internal(name, comment, &seq).unwrap();
                            }
                        }
                    }
                    self.current_idx = 0;
                }
                _ => return None,
            }
        }
    }
}

/// Iterator over ORFs from a reader.
pub struct OrfReaderIter<'a> {
    caller: &'a OrfCaller,
    reader: Box<dyn needletail::FastxReader + 'a>,
    current_orfs: Vec<Orf>,
    current_idx: usize,
    /// Once true, all subsequent records use seq() instead of raw_seq().
    needs_normalize: bool,
}

impl Iterator for OrfReaderIter<'_> {
    type Item = Orf;

    fn next(&mut self) -> Option<Orf> {
        loop {
            if self.current_idx < self.current_orfs.len() {
                let orf = self.current_orfs[self.current_idx].clone();
                self.current_idx += 1;
                return Some(orf);
            }

            match self.reader.next() {
                Some(Ok(record)) => {
                    let (name, comment) = split_header(record.id());
                    if self.needs_normalize {
                        let seq = record.seq();
                        self.current_orfs =
                            self.caller.find_orfs_internal(name, comment, &seq).unwrap();
                    } else {
                        match self
                            .caller
                            .find_orfs_internal(name, comment, record.raw_seq())
                        {
                            Some(orfs) => self.current_orfs = orfs,
                            None => {
                                self.needs_normalize = true;
                                let seq = record.seq();
                                self.current_orfs =
                                    self.caller.find_orfs_internal(name, comment, &seq).unwrap();
                            }
                        }
                    }
                    self.current_idx = 0;
                }
                _ => return None,
            }
        }
    }
}

/// Write ORFs as FASTA to the given writer.
pub fn write_orfs_fasta<W: io::Write>(
    orfs: impl Iterator<Item = Orf>,
    writer: &mut W,
) -> io::Result<()> {
    for orf in orfs {
        writeln!(writer, "{}", orf.header())?;
        writer.write_all(&orf.protein)?;
        writeln!(writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
