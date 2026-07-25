use std::str::FromStr;

use strum::{Display, EnumProperty, EnumString};

use crate::gtex_table::GTExTable;

// Serializes from VEP and BCSQ consequences, respectively. Some consequences are VEP or BCSQ specific and
// will only have a single serialization, others may be identical in serialization between both VEP and BCSQ,
// so these also only have a single serialization.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, EnumString, Display, EnumProperty, Hash)]
pub enum ConsequenceType {
    #[strum(serialize = "transcript_ablation", props(protein_affecting = true))] // VEP-specific
    TranscriptAblation,

    #[strum(
        serialize = "splice_acceptor_variant",
        serialize = "splice_acceptor",
        props(protein_affecting = true)
    )]
    SpliceAcceptor,

    #[strum(
        serialize = "splice_donor_variant",
        serialize = "splice_donor",
        props(protein_affecting = true)
    )]
    SpliceDonor,

    #[strum(serialize = "stop_gained", props(protein_affecting = true))]
    StopGained,

    #[strum(
        serialize = "frameshift_variant",
        serialize = "frameshift",
        props(protein_affecting = true)
    )]
    Frameshift,

    #[strum(serialize = "stop_lost", props(protein_affecting = true))]
    StopLost,

    #[strum(serialize = "start_lost", props(protein_affecting = true))]
    StartLost,

    #[strum(
        serialize = "transcript_amplification",
        props(protein_affecting = true)
    )] // VEP-specific
    TranscriptAmplification,

    #[strum(serialize = "feature_elongation", props(protein_affecting = true))]
    FeatureElongation,

    #[strum(serialize = "feature_truncation", props(protein_affecting = true))] // VEP-specific
    FeatureTruncation,

    #[strum(serialize = "inframe_insertion", props(protein_affecting = true))]
    InframeInsertion,

    #[strum(serialize = "inframe_deletion", props(protein_affecting = true))]
    InframeDeletion,

    #[strum(serialize = "inframe_altering", props(protein_affecting = true))] // BCSQ-specific
    InframeAltering,

    #[strum(
        serialize = "missense_variant",
        serialize = "missense",
        props(protein_affecting = true)
    )]
    Missense,

    #[strum(
        serialize = "protein_altering_variant",
        props(protein_affecting = true)
    )] // VEP-specific
    ProteinAltering,

    #[strum(
        serialize = "splice_donor_5th_base_variant",
        props(protein_affecting = true)
    )] // VEP-specific
    SpliceDonorBase5,

    #[strum(
        serialize = "splice_region_variant",
        serialize = "splice_region",
        props(protein_affecting = true)
    )]
    SpliceRegion,

    #[strum(
        serialize = "splice_donor_region_variant",
        props(protein_affecting = true)
    )] // VEP-specific
    SpliceDonorRegion,

    #[strum(
        serialize = "splice_polypyrimidine_tract_variant",
        props(protein_affecting = true)
    )]
    SplicePolypyrimidineTract,

    #[strum(
        serialize = "incomplete_terminal_codon_variant",
        props(protein_affecting = true)
    )]
    IncompleteTerminalCodon,

    #[strum(
        serialize = "start_retained_variant",
        serialize = "start_retained",
        props(protein_affecting = true)
    )]
    StartRetained,

    #[strum(
        serialize = "stop_retained_variant",
        serialize = "stop_retained",
        props(protein_affecting = true)
    )]
    StopRetained,

    #[strum(
        serialize = "synonymous_variant",
        serialize = "synonymous",
        props(protein_affecting = true)
    )]
    Synonymous,

    #[strum(
        serialize = "coding_sequence_variant",
        serialize = "coding_sequence",
        props(protein_affecting = true)
    )]
    CodingSequence,

    #[strum(serialize = "mature_miRNA_variant", props(protein_affecting = false))] // VEP-specific
    MatureMiRNA,

    #[strum(
        serialize = "5_prime_UTR_variant",
        serialize = "5_prime_utr",
        props(protein_affecting = false)
    )]
    UTR5Prime,

    #[strum(
        serialize = "3_prime_UTR_variant",
        serialize = "3_prime_utr",
        props(protein_affecting = false)
    )]
    UTR3Prime,

    #[strum(
        serialize = "non_coding_transcript_exon_variant",
        props(protein_affecting = false)
    )] // VEP-specific
    NonCodingExon,

    #[strum(
        serialize = "intron_variant",
        serialize = "intron",
        props(protein_affecting = false)
    )]
    Intron,

    #[strum(serialize = "NMD_transcript_variant", props(protein_affecting = false))]
    // VEP-specific
    NMDTranscript,

    #[strum(
        serialize = "non_coding_transcript_variant",
        serialize = "non_coding",
        props(protein_affecting = false)
    )]
    NonCoding,

    #[strum(
        serialize = "coding_transcript_variant",
        props(protein_affecting = false) // Can't be confidently said to be protein affecting
    )] // VEP-specific
    Coding,

    #[strum(serialize = "upstream_gene_variant", props(protein_affecting = false))] // VEP-specific
    UpstreamGene,

    #[strum(
        serialize = "downstream_gene_variant",
        props(protein_affecting = false)
    )] // VEP-specific
    DownstreamGene,

    #[strum(serialize = "TFBS_ablation", props(protein_affecting = false))]
    TFBSAblation,

    #[strum(serialize = "TFBS_amplification", props(protein_affecting = false))]
    TFBSAmplification,

    #[strum(
        serialize = "TF_binding_site_variant",
        props(protein_affecting = false)
    )]
    TFBindingSite,

    #[strum(
        serialize = "regulatory_region_ablation",
        props(protein_affecting = false)
    )]
    RegulatoryRegionAblation,

    #[strum(
        serialize = "regulatory_region_amplification",
        props(protein_affecting = false)
    )]
    RegulatoryRegionAmplification,

    #[strum(serialize = "regulatory_region_variant", props(protein_coding = false))]
    RegulatoryRegion,

    #[strum(
        serialize = "intergenic_variant",
        serialize = "intergenic",
        props(protein_affecting = false)
    )]
    Intergenic,

    #[strum(serialize = "sequence_variant", props(protein_affecting = false))]
    Sequence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Consequence {
    pub consequence_type: ConsequenceType,
    pub protein_affecting: bool,
    pub transcript_id: String,
    pub gene_id: String,
    pub group_columns: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    TooFewFields {
        csq: String,
        expected: usize,
        found: usize,
    },
    UnsupportedConsequence {
        consequence: String,
    },
    TranscriptNotInGTEx {
        transcript: String,
    },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::TooFewFields {
                csq,
                expected,
                found,
            } => write!(
                f,
                "Too few CSQ fields: expected at least {expected}, found {found} in \"{csq}\""
            ),
            ParseError::UnsupportedConsequence { consequence } => {
                write!(
                    f,
                    "Cannot find consequence {consequence}. See 'https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html' for valid VEP consequences and 'https://samtools.github.io/bcftools/bcftools.html' for valid BCSQ consequences (\"Supported consequence types\")."
                )
            }
            ParseError::TranscriptNotInGTEx { transcript } => write!(
                f,
                "Could not map ENST to ENSG because the transcript '{transcript}' is not in the provided GTEx table."
            ),
        }
    }
}

impl std::error::Error for ParseError {}

impl Consequence {
    fn strip_ensembl_version<S: AsRef<str>>(s: S) -> String {
        let s = s.as_ref();
        match s.split_once('.') {
            Some((id, _)) => id.to_string(),
            None => s.to_string(),
        }
    }

    pub fn parse_from_vep(csq: &str, table: &GTExTable) -> Result<Self, ParseError> {
        // Gene_id (ensembl or symbol) is added to group columns as well to assist in deduplication logic
        const CONSEQUENCE: usize = 1;
        const FEATURE: usize = 6;
        // const LOFTEE: usize = 24;
        // const GROUP_COLUMNS: [usize; 2] = [CONSEQUENCE, LOFTEE];

        let fields: Vec<&str> = csq.split('|').collect();
        const MIN_LEN: usize = 25;
        if fields.len() < MIN_LEN {
            return Err(ParseError::TooFewFields {
                csq: csq.to_string(),
                expected: MIN_LEN,
                found: fields.len(),
            });
        }

        let transcript_id = fields[FEATURE];
        let gene_id = table.get_transcript_gene(&transcript_id).map_err(|_| {
            ParseError::TranscriptNotInGTEx {
                transcript: transcript_id.to_string(),
            }
        })?;

        let consequence_field = fields[CONSEQUENCE].splitn(2, '&').next().unwrap();

        let consequence_type = ConsequenceType::from_str(consequence_field).map_err(|_| {
            ParseError::UnsupportedConsequence {
                consequence: fields[CONSEQUENCE].to_string(),
            }
        })?;

        let protein_affecting = consequence_type.get_bool("protein_affecting").unwrap();

        Ok(Consequence {
            consequence_type: consequence_type.clone(),
            protein_affecting: protein_affecting,
            gene_id: gene_id.to_string(),
            transcript_id: transcript_id.to_string(),
            group_columns: vec![
                gene_id.to_string(),
                protein_affecting.to_string(),
                consequence_type.to_string(),
            ],
        })
    }

    pub fn parse_from_bcsq(csq: &str, table: &GTExTable) -> Result<Self, ParseError> {
        const CONSEQUENCE: usize = 0;
        const TRANSCRIPT_ID: usize = 2;

        let fields: Vec<&str> = csq.split('|').collect();
        const MIN_LEN: usize = 4;
        if fields.len() < MIN_LEN {
            return Err(ParseError::TooFewFields {
                csq: csq.to_string(),
                expected: MIN_LEN,
                found: fields.len(),
            });
        }

        let transcript_id = Consequence::strip_ensembl_version(fields[TRANSCRIPT_ID]);
        let gene_id = table.get_transcript_gene(&transcript_id).map_err(|_| {
            ParseError::TranscriptNotInGTEx {
                transcript: transcript_id.to_string(),
            }
        })?;

        let consequence_field = fields[CONSEQUENCE].splitn(2, '&').next().unwrap();

        let consequence_type = ConsequenceType::from_str(consequence_field).map_err(|_| {
            ParseError::UnsupportedConsequence {
                consequence: fields[CONSEQUENCE].to_string(),
            }
        })?;

        let protein_affecting = consequence_type.get_bool("protein_affecting").unwrap();

        Ok(Consequence {
            consequence_type: consequence_type.clone(),
            protein_affecting: protein_affecting,
            gene_id: gene_id.to_string(),
            transcript_id: transcript_id,
            group_columns: vec![
                gene_id.to_string(),
                protein_affecting.to_string(),
                consequence_type.to_string(),
            ],
        })
    }

    pub fn from_fields(
        consequence: String,
        gene_id: String,
        transcript_id: String,
        group_columns: Vec<String>,
    ) -> Result<Self, ParseError> {
        let consequence_field = consequence.splitn(2, '&').next().unwrap();

        let consequence_type = ConsequenceType::from_str(consequence_field).map_err(|_| {
            ParseError::UnsupportedConsequence {
                consequence: consequence,
            }
        })?;

        let protein_affecting = consequence_type.get_bool("protein_affecting").unwrap();

        let mut complete_group_columns = Vec::with_capacity(group_columns.len() + 3);
        complete_group_columns.push(gene_id.clone());
        complete_group_columns.push(protein_affecting.to_string());
        complete_group_columns.push(consequence_type.to_string());
        complete_group_columns.extend(group_columns);

        Ok(Consequence {
            consequence_type: consequence_type.clone(),
            protein_affecting: protein_affecting,
            gene_id: gene_id,
            transcript_id: transcript_id,
            group_columns: complete_group_columns,
        })
    }
}
