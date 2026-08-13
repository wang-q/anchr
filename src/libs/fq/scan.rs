//! Zero-copy FASTQ record scanning over an in-memory buffer (public
//! facility, originally extracted from the QC hot path).
//!
//! Mirrors falco's record reader: scans the whole buffer with `memchr` and
//! returns slices into it, so a record costs 4 newline searches + pointer
//! arithmetic instead of per-line `Vec` allocation and `extend_from_slice`
//! copies. Fast for single-pass, read-mostly pipelines (QC, sampling).
//!
//! Limitations (documented so callers choose deliberately):
//! * FASTQ only — a record whose third line does not start with `+` is
//!   rejected (`next_record` returns `None`); multi-line FASTA is not
//!   supported. Use `pgr::libs::fmt::seq::SeqReader` when FASTA or owned
//!   records are needed.
//! * The caller owns the buffer (mmap or an in-memory Vec); for gzip input
//!   either decompress the whole file first or keep the streaming reader.

/// One FASTQ record as slices into the owning buffer.
#[derive(Debug, Clone, Copy)]
pub struct FastqRecord<'a> {
    /// Header without the leading `@`.
    pub name: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

#[inline]
fn strip_cr(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// Parse the next FASTQ record starting at `*pos` (updated past the record
/// including its trailing newline). Returns `None` at EOF.
pub fn next_record<'a>(data: &'a [u8], pos: &mut usize) -> Option<FastqRecord<'a>> {
    let rest = &data[*pos..];
    if rest.is_empty() {
        return None;
    }
    let hdr_end = memchr::memchr(b'\n', rest)?;
    let hdr = strip_cr(&rest[..hdr_end]);
    let name = hdr
        .strip_prefix(b"@")
        .or_else(|| hdr.strip_prefix(b">"))
        .unwrap_or(hdr);

    let mut p = *pos + hdr_end + 1;
    let seq_end_rel = memchr::memchr(b'\n', &data[p..])?;
    let seq = strip_cr(&data[p..p + seq_end_rel]);
    p += seq_end_rel + 1;

    let plus_end_rel = memchr::memchr(b'\n', &data[p..])?;
    let plus = strip_cr(&data[p..p + plus_end_rel]);
    if !plus.starts_with(b"+") {
        // FASTA record: qc needs quality scores, so reject it like other
        // quality-requiring consumers instead of silently misparsing.
        return None;
    }
    p += plus_end_rel + 1;

    let qual_end_rel = memchr::memchr(b'\n', &data[p..])?;
    let qual = strip_cr(&data[p..p + qual_end_rel]);
    *pos = p + qual_end_rel + 1;

    Some(FastqRecord { name, seq, qual })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_crlf_and_plain() {
        let data = b"@r1\nACGT\n+\n!!!!\n@r2\r\nTGCA\r\n+\r\n####\r\n";
        let mut pos = 0;
        let r1 = next_record(data, &mut pos).unwrap();
        assert_eq!(r1.name, b"r1");
        assert_eq!(r1.seq, b"ACGT");
        assert_eq!(r1.qual, b"!!!!");
        let r2 = next_record(data, &mut pos).unwrap();
        assert_eq!(r2.name, b"r2");
        assert_eq!(r2.seq, b"TGCA");
        assert_eq!(r2.qual, b"####");
        assert!(next_record(data, &mut pos).is_none());
    }

}
