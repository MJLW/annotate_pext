//! A deliberately small stand-in for the part of LOFTEE that matters for
//! grouping pext annotations: the HC/LC call.
//!
//! LOFTEE's `LC` verdict is dominated by `END_TRUNC` — a pLoF that lands so
//! close to the 3' end of the CDS that the truncated protein is likely still
//! functional and the transcript likely escapes nonsense-mediated decay. This
//! module reproduces that single filter and nothing else. It does NOT
//! reproduce splice-rescue detection, ancestral-allele checks, small-intron
//! filtering, GC-to-GT donor checks, or any LOFTEE *flag*. A call of `HC` here
//! means "pLoF, and not near the transcript end" — it does not mean LOFTEE
//! would have said HC.
//!
//! Three rules are available:
//!   * `LastExon`   - LC if the variant falls in the final coding exon.
//!   * `NmdEscape`  - LC if in the final coding exon, or within `nmd_window`
//!                    coding bases 5' of the final exon-exon junction (the
//!                    classic 50 bp NMD-escape rule).
//!   * `EndTrunc`   - LC if the variant falls in the terminal
//!                    `end_trunc_fraction` of the CDS. This is closest to
//!                    LOFTEE's own END_TRUNC (default 0.05).
//!
//! Non-pLoF consequences always return `NotApplicable` ("."), matching VEP,
//! which leaves `LoF` empty for anything that is not a pLoF.

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// 1-based, inclusive on both ends (GFF convention).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start) + 1
    }

    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }

    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos <= self.end
    }

    /// 0 if `pos` is inside the interval, otherwise the gap to the nearer edge.
    pub fn distance_to(&self, pos: u64) -> u64 {
        if pos < self.start {
            self.start - pos
        } else if pos > self.end {
            pos - self.end
        } else {
            0
        }
    }
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub seqid: String,
    pub strand: Strand,
    /// Ascending by start, genome order.
    pub exons: Vec<Interval>,
    /// Ascending by start, genome order.
    pub cds: Vec<Interval>,
}

impl Transcript {
    pub fn is_coding(&self) -> bool {
        !self.cds.is_empty()
    }

    pub fn cds_len(&self) -> u64 {
        self.cds.iter().map(|i| i.len()).sum()
    }

    /// CDS blocks in transcription order (5' -> 3').
    fn cds_in_tx_order(&self) -> Vec<Interval> {
        let mut blocks = self.cds.clone();
        blocks.sort_by_key(|i| i.start);
        if self.strand == Strand::Reverse {
            blocks.reverse();
        }
        blocks
    }

    /// 0-based offset of `pos` along the CDS, counted 5' -> 3'.
    ///
    /// Essential splice sites sit in introns and so have no CDS offset of
    /// their own. Rather than give up on them (they are half of the pLoF set)
    /// we snap to the nearest coding base, which is the exon boundary the
    /// splice site belongs to. This is an approximation and is the main place
    /// where these calls diverge from LOFTEE for splice variants.
    pub fn cds_offset(&self, pos: u64) -> Option<u64> {
        let blocks = self.cds_in_tx_order();
        if blocks.is_empty() {
            return None;
        }

        let mut cumulative: u64 = 0;
        let mut best: Option<(u64, u64)> = None; // (distance, offset)

        for block in &blocks {
            let distance = block.distance_to(pos);
            let within = if self.strand == Strand::Forward {
                pos.saturating_sub(block.start).min(block.len() - 1)
            } else {
                block.end.saturating_sub(pos).min(block.len() - 1)
            };
            let offset = cumulative + within;

            if best.map_or(true, |(best_distance, _)| distance < best_distance) {
                best = Some((distance, offset));
            }
            if distance == 0 {
                break;
            }
            cumulative += block.len();
        }

        best.map(|(_, offset)| offset)
    }

    /// Coding bases between `pos` and the 3' end of the CDS.
    pub fn cds_distance_to_end(&self, pos: u64) -> Option<u64> {
        let len = self.cds_len();
        self.cds_offset(pos).map(|o| len.saturating_sub(o + 1))
    }

    /// Does `pos` belong to the final coding exon? Ties between two equidistant
    /// blocks resolve 3'-ward, so a splice site on the boundary of the last
    /// coding exon counts as being in it.
    pub fn in_last_coding_exon(&self, pos: u64) -> bool {
        let blocks = self.cds_in_tx_order();
        if blocks.is_empty() {
            return false;
        }

        let last = blocks.len() - 1;
        let mut best_distance = u64::MAX;
        let mut best_is_last = false;

        for (i, block) in blocks.iter().enumerate() {
            let distance = block.distance_to(pos);
            if distance <= best_distance {
                best_distance = distance;
                best_is_last = i == last;
            }
        }

        best_is_last
    }

    /// Final coding exon, or within `window` coding bases 5' of the final
    /// coding exon-exon junction.
    pub fn escapes_nmd(&self, pos: u64, window: u64) -> bool {
        if self.in_last_coding_exon(pos) {
            return true;
        }
        let blocks = self.cds_in_tx_order();
        if blocks.len() < 2 {
            return false;
        }
        let last_len = blocks[blocks.len() - 1].len();
        match self.cds_distance_to_end(pos) {
            Some(distance) => distance < last_len + window,
            None => false,
        }
    }

    /// Position along the CDS as a fraction in (0, 1].
    pub fn cds_fraction(&self, pos: u64) -> Option<f64> {
        let len = self.cds_len();
        if len == 0 {
            return None;
        }
        self.cds_offset(pos)
            .map(|offset| (offset as f64 + 1.0) / len as f64)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rule {
    LastExon,
    NmdEscape,
    EndTrunc,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LofCall {
    Hc,
    Lc,
    /// Not a pLoF, or the transcript is absent / non-coding: VEP writes an
    /// empty LoF field here and the pext TSVs use ".".
    NotApplicable,
}

impl LofCall {
    pub fn as_str(&self) -> &'static str {
        match self {
            LofCall::Hc => "HC",
            LofCall::Lc => "LC",
            LofCall::NotApplicable => ".",
        }
    }
}

/// Terms treated as pLoF, in both VEP and bcftools csq spellings.
pub const PLOF_TERMS: [&str; 8] = [
    "transcript_ablation",
    "splice_acceptor",
    "splice_acceptor_variant",
    "splice_donor",
    "splice_donor_variant",
    "stop_gained",
    "frameshift",
    "frameshift_variant",
];

pub fn is_plof(consequence: &str) -> bool {
    consequence
        .split('&')
        .any(|term| PLOF_TERMS.contains(&term.trim()))
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Counts {
    pub hc: u64,
    pub lc: u64,
    pub not_applicable: u64,
    pub transcript_absent: u64,
    pub transcript_noncoding: u64,
}

pub struct Classifier {
    transcripts: HashMap<String, Transcript>,
    pub rule: Rule,
    pub end_trunc_fraction: f64,
    pub nmd_window: u64,
    pub counts: Counts,
}

impl Classifier {
    pub fn new(
        transcripts: HashMap<String, Transcript>,
        rule: Rule,
        end_trunc_fraction: f64,
        nmd_window: u64,
    ) -> Self {
        Classifier {
            transcripts,
            rule,
            end_trunc_fraction,
            nmd_window,
            counts: Counts::default(),
        }
    }

    pub fn transcript_count(&self) -> usize {
        self.transcripts.len()
    }

    pub fn call(&mut self, transcript_id: &str, pos: u64, consequence: &str) -> LofCall {
        if !is_plof(consequence) {
            self.counts.not_applicable += 1;
            return LofCall::NotApplicable;
        }

        let key = strip_version(transcript_id);
        let transcript = match self.transcripts.get(&key) {
            Some(t) => t,
            None => {
                self.counts.transcript_absent += 1;
                self.counts.not_applicable += 1;
                return LofCall::NotApplicable;
            }
        };

        if !transcript.is_coding() {
            self.counts.transcript_noncoding += 1;
            self.counts.not_applicable += 1;
            return LofCall::NotApplicable;
        }

        let truncating = match self.rule {
            Rule::LastExon => transcript.in_last_coding_exon(pos),
            Rule::NmdEscape => transcript.escapes_nmd(pos, self.nmd_window),
            Rule::EndTrunc => transcript
                .cds_fraction(pos)
                .map_or(false, |f| f > 1.0 - self.end_trunc_fraction),
        };

        if truncating {
            self.counts.lc += 1;
            LofCall::Lc
        } else {
            self.counts.hc += 1;
            LofCall::Hc
        }
    }
}

fn open_maybe_gzip<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read>, Box<dyn Error>> {
    let file = File::open(&path)?;
    let reader: Box<dyn Read> = match path.as_ref().extension().and_then(|e| e.to_str()) {
        Some("gz") | Some("bgz") => Box::new(MultiGzDecoder::new(file)),
        _ => Box::new(file),
    };
    Ok(reader)
}

pub fn strip_version(id: &str) -> String {
    id.split_once('.').map(|(i, _)| i).unwrap_or(id).to_string()
}

fn normalise_contig(seqid: &str) -> String {
    seqid.strip_prefix("chr").unwrap_or(seqid).to_string()
}

fn attribute<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    for field in attributes.split(';') {
        let field = field.trim();
        if let Some(rest) = field.strip_prefix(key) {
            if let Some(value) = rest.strip_prefix('=') {
                return Some(value);
            }
        }
    }
    None
}

/// Works for both Ensembl GFF3 (`Parent=transcript:ENST...`) and GENCODE GFF3
/// (`transcript_id=ENST....N`).
fn transcript_of(attributes: &str) -> Option<String> {
    if let Some(value) = attribute(attributes, "transcript_id") {
        return Some(strip_version(value));
    }
    let parent = attribute(attributes, "Parent")?;
    let parent = parent.split(',').next()?;
    let parent = parent.strip_prefix("transcript:").unwrap_or(parent);
    if parent.is_empty() {
        None
    } else {
        Some(strip_version(parent))
    }
}

/// Reads exon and CDS features. Transcript-level records are not needed: a
/// transcript is coding iff it has at least one CDS feature.
pub fn load_gff<P: AsRef<Path>>(path: P) -> Result<HashMap<String, Transcript>, Box<dyn Error>> {
    let reader = BufReader::with_capacity(1 << 20, open_maybe_gzip(path)?);
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature = fields[2];
        if feature != "exon" && feature != "CDS" {
            continue;
        }

        let transcript_id = match transcript_of(fields[8]) {
            Some(id) => id,
            None => continue,
        };

        let start: u64 = fields[3].parse()?;
        let end: u64 = fields[4].parse()?;
        let strand = if fields[6] == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        let entry = transcripts
            .entry(transcript_id.clone())
            .or_insert_with(|| Transcript {
                id: transcript_id.clone(),
                seqid: normalise_contig(fields[0]),
                strand,
                exons: Vec::new(),
                cds: Vec::new(),
            });

        let interval = Interval { start, end };
        if feature == "exon" {
            entry.exons.push(interval);
        } else {
            entry.cds.push(interval);
        }
    }

    for transcript in transcripts.values_mut() {
        transcript.exons.sort_by_key(|i| i.start);
        transcript.cds.sort_by_key(|i| i.start);
    }

    Ok(transcripts)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Forward-strand transcript, three coding exons of 100 bp each:
    /// 1000-1099, 2000-2099, 3000-3099. CDS length 300.
    fn forward() -> Transcript {
        Transcript {
            id: "ENST_FWD".to_string(),
            seqid: "1".to_string(),
            strand: Strand::Forward,
            exons: vec![
                Interval { start: 1000, end: 1099 },
                Interval { start: 2000, end: 2099 },
                Interval { start: 3000, end: 3099 },
            ],
            cds: vec![
                Interval { start: 1000, end: 1099 },
                Interval { start: 2000, end: 2099 },
                Interval { start: 3000, end: 3099 },
            ],
        }
    }

    /// Same coordinates, reverse strand: transcription runs 3099 -> 1000.
    fn reverse() -> Transcript {
        let mut t = forward();
        t.id = "ENST_REV".to_string();
        t.strand = Strand::Reverse;
        t
    }

    #[test]
    fn cds_offset_runs_five_to_three() {
        let t = forward();
        assert_eq!(t.cds_offset(1000), Some(0));
        assert_eq!(t.cds_offset(1099), Some(99));
        assert_eq!(t.cds_offset(2000), Some(100));
        assert_eq!(t.cds_offset(3099), Some(299));
        assert_eq!(t.cds_len(), 300);
    }

    #[test]
    fn cds_offset_respects_reverse_strand() {
        let t = reverse();
        assert_eq!(t.cds_offset(3099), Some(0));
        assert_eq!(t.cds_offset(3000), Some(99));
        assert_eq!(t.cds_offset(1000), Some(299));
    }

    #[test]
    fn last_coding_exon_is_strand_aware() {
        assert!(forward().in_last_coding_exon(3050));
        assert!(!forward().in_last_coding_exon(1050));
        assert!(reverse().in_last_coding_exon(1050));
        assert!(!reverse().in_last_coding_exon(3050));
    }

    #[test]
    fn splice_sites_snap_to_the_nearest_coding_base() {
        // Donor site immediately after exon 2 of the forward transcript.
        let t = forward();
        assert_eq!(t.cds_offset(2100), Some(199));
        assert!(!t.in_last_coding_exon(2100));
        // Acceptor site immediately before the final exon.
        assert_eq!(t.cds_offset(2999), Some(200));
        assert!(t.in_last_coding_exon(2999));
    }

    #[test]
    fn nmd_escape_covers_the_last_fifty_bases_of_the_penultimate_exon() {
        let t = forward();
        // 100 bp final exon + 50 bp window => last 150 coding bases.
        assert!(t.escapes_nmd(3000, 50));
        assert!(t.escapes_nmd(2099, 50)); // 1 base from the junction
        assert!(t.escapes_nmd(2051, 50));
        assert!(!t.escapes_nmd(2000, 50)); // 100 bases from the junction
    }

    #[test]
    fn end_trunc_uses_the_terminal_fraction_of_cds() {
        let t = forward();
        assert_eq!(t.cds_fraction(3099), Some(1.0));
        // 5% of 300 bp = the final 15 bases, i.e. offsets 285..=299.
        assert!(t.cds_fraction(3085).unwrap() > 0.95);
        assert!(t.cds_fraction(3084).unwrap() <= 0.95);
    }

    #[test]
    fn only_plof_consequences_are_called() {
        assert!(is_plof("stop_gained"));
        assert!(is_plof("splice_donor_variant&intron_variant"));
        assert!(is_plof("splice_acceptor"));
        assert!(!is_plof("missense_variant"));
        assert!(!is_plof("splice_region_variant&intron_variant"));
    }

    #[test]
    fn classifier_returns_dot_for_unknown_transcripts() {
        let mut c = Classifier::new(HashMap::new(), Rule::LastExon, 0.05, 50);
        assert_eq!(c.call("ENST_MISSING", 1000, "stop_gained"), LofCall::NotApplicable);
        assert_eq!(c.counts.transcript_absent, 1);
    }

    #[test]
    fn classifier_calls_hc_and_lc() {
        let mut map = HashMap::new();
        map.insert("ENST_FWD".to_string(), forward());
        let mut c = Classifier::new(map, Rule::LastExon, 0.05, 50);
        assert_eq!(c.call("ENST_FWD.3", 3050, "stop_gained"), LofCall::Lc);
        assert_eq!(c.call("ENST_FWD", 1050, "stop_gained"), LofCall::Hc);
        assert_eq!(c.call("ENST_FWD", 1050, "missense_variant"), LofCall::NotApplicable);
    }
}
