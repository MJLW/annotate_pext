use crate::gtex_table::GTExTable;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Consequence {
    pub transcript_id: String,
    pub gene_id: String,
    pub group_columns: Vec<String>,
    pub protein_coding: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    TooFewFields {
        csq: String,
        expected: usize,
        found: usize,
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
                "too few CSQ fields: expected at least {expected}, found {found} in \"{csq}\""
            ),
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
        const IMPACT: usize = 2;
        const FEATURE: usize = 6;
        const BIOTYPE: usize = 7;
        const LOFTEE: usize = 24;
        const GROUP_COLUMNS: [usize; 2] = [CONSEQUENCE, LOFTEE];

        let fields: Vec<&str> = csq.split('|').collect();
        const MIN_LEN: usize = 25;
        if fields.len() < MIN_LEN {
            return Err(ParseError::TooFewFields {
                csq: csq.to_string(),
                expected: MIN_LEN,
                found: fields.len(),
            });
        }

        let mut group_columns: Vec<String> = GROUP_COLUMNS
            .iter()
            .map(|&i| fields[i].to_string())
            .collect();

        let cds_affecting = fields[BIOTYPE] == "protein_coding" && fields[IMPACT] != "MODIFIER";
        group_columns.push(cds_affecting.to_string());

        let transcript_id = fields[FEATURE];
        let gene_id = table.get_transcript_gene(&transcript_id).map_err(|_| {
            ParseError::TranscriptNotInGTEx {
                transcript: transcript_id.to_string(),
            }
        })?;

        Ok(Consequence {
            gene_id: gene_id.to_string(),
            transcript_id: transcript_id.to_string(),
            group_columns: group_columns,
            protein_coding: cds_affecting,
        })
    }

    pub fn parse_from_bcsq(csq: &str, table: &GTExTable) -> Result<Self, ParseError> {
        const CONSEQUENCE: usize = 0;
        const TRANSCRIPT_ID: usize = 2;
        const BIOTYPE: usize = 3;
        const GROUP_COLUMNS: [usize; 1] = [CONSEQUENCE];

        let fields: Vec<&str> = csq.split('|').collect();
        const MIN_LEN: usize = 4;
        if fields.len() < MIN_LEN {
            return Err(ParseError::TooFewFields {
                csq: csq.to_string(),
                expected: MIN_LEN,
                found: fields.len(),
            });
        }

        let mut group_columns: Vec<String> = GROUP_COLUMNS
            .iter()
            .map(|&i| fields[i].to_string())
            .collect();

        let cds_affecting = fields[BIOTYPE] == "protein_coding" && fields[TRANSCRIPT_ID].len() > 0;
        group_columns.push(cds_affecting.to_string());

        let transcript_id = Consequence::strip_ensembl_version(fields[TRANSCRIPT_ID]);
        let gene_id = table.get_transcript_gene(&transcript_id).map_err(|_| {
            ParseError::TranscriptNotInGTEx {
                transcript: transcript_id.to_string(),
            }
        })?;

        Ok(Consequence {
            gene_id: gene_id.to_string(),
            transcript_id: transcript_id,
            group_columns: group_columns,
            protein_coding: cds_affecting,
        })
    }

    pub fn from_fields(
        gene_id: String,
        transcript_id: String,
        biotype: String,
        group_columns: Vec<String>,
    ) -> Self {
        let mut complete_group_columns = Vec::with_capacity(group_columns.len() + 3);
        complete_group_columns.push(gene_id.clone());
        complete_group_columns.push(biotype.clone());
        complete_group_columns.extend(group_columns);

        Consequence {
            gene_id: gene_id,
            transcript_id: transcript_id,
            group_columns: complete_group_columns,
            protein_coding: biotype == "protein_coding",
        }
    }
}
