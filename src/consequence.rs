use std::str::FromStr;

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
        }
    }
}

impl std::error::Error for ParseError {}

impl Consequence {
    pub fn parse_from_vep(csq: &str) -> Result<Self, ParseError> {
        // Gene_id (ensembl or symbol) is added to group columns as well to assist in deduplication logic
        const CONSEQUENCE: usize = 1;
        const IMPACT: usize = 2;
        const GENE_ID: usize = 4;
        const FEATURE: usize = 6;
        const BIOTYPE: usize = 7;
        const LOFTEE: usize = 24;
        const GROUP_COLUMNS: [usize; 3] = [GENE_ID, CONSEQUENCE, LOFTEE];

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

        Ok(Consequence {
            gene_id: fields[GENE_ID].to_string(),
            transcript_id: fields[FEATURE].to_string(),
            group_columns: group_columns,
            protein_coding: cds_affecting,
        })
    }

    pub fn parse_from_bcsq(csq: &str) -> Result<Self, ParseError> {
        const CONSEQUENCE: usize = 0;
        const GENE_ID: usize = 1;
        const TRANSCRIPT_ID: usize = 2;
        const BIOTYPE: usize = 3;
        const GROUP_COLUMNS: [usize; 2] = [GENE_ID, CONSEQUENCE];

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

        Ok(Consequence {
            gene_id: fields[GENE_ID].to_string(),
            transcript_id: fields[TRANSCRIPT_ID].to_string(),
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

impl FromStr for Consequence {
    type Err = ParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Consequence::parse_from_vep(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VEP_LINE: &str = "-|intron_variant|MODIFIER|FANCD2|ENSG00000144554|Transcript|ENST00000287647|protein_coding||3/42||||||||||1||||Homo_sapiens.GRCh37.87.gff3.chr.sorted.gz|||||";

    #[test]
    fn parses_all_fields() {
        let c = Consequence::parse_from_vep(VEP_LINE).unwrap();
        assert_eq!(c.transcript_id, "ENST00000287647");
        assert_eq!(c.protein_coding, true);
    }

    #[test]
    fn works_via_from_str() {
        let c: Consequence = VEP_LINE.parse().unwrap();
        assert_eq!(c.transcript_id, "ENST00000287647");
    }

    #[test]
    fn rejects_truncated_input() {
        assert!(Consequence::parse_from_vep("-|intron_variant|MODIFIER").is_err());
    }
}
