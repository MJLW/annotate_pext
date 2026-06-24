use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Consequence {
    pub consequence: String,
    pub transcript_id: String,
    pub gene_id: String,
    pub group_columns: Vec<String>,
    pub protein_coding: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    TooFewFields { expected: usize, found: usize },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::TooFewFields { expected, found } => write!(
                f,
                "too few CSQ fields: expected at least {expected}, found {found}"
            ),
        }
    }
}

impl std::error::Error for ParseError {}

impl Consequence {
    pub fn parse_from_vep(csq: &str) -> Result<Self, ParseError> {
        // Gene_id (ensembl or symbol) is added to group columns as well to assist in deduplication logic
        const CONSEQUENCE: usize = 1;
        const GENE_ID: usize = 4;
        const FEATURE: usize = 6;
        const GROUP_COLUMNS: [usize; 3] = [4, 1, 24];
        const BIOTYPE: usize = 7;

        let fields: Vec<&str> = csq.split('|').collect();
        let min_len = 25;
        if fields.len() < min_len {
            return Err(ParseError::TooFewFields {
                expected: min_len,
                found: fields.len(),
            });
        }

        // println!(
        //     "{:?}, {:?}",
        //     fields[BIOTYPE],
        //     fields[BIOTYPE] == "protein_coding"
        // );
        Ok(Consequence {
            consequence: fields[CONSEQUENCE].to_string(),
            gene_id: fields[GENE_ID].to_string(),
            transcript_id: fields[FEATURE].to_string(),
            group_columns: GROUP_COLUMNS
                .iter()
                .map(|&i| fields[i].to_string())
                .collect(),
            protein_coding: fields[BIOTYPE] == "protein_coding",
        })
    }

    // pub fn parse_from_bcsq(csq: &str) -> Result<Self, ParseError> {
    //
    // }

    pub fn consequence_terms(&self) -> impl Iterator<Item = &str> {
        self.consequence.split('&')
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
        assert_eq!(c.consequence, "intron_variant");
        assert_eq!(c.transcript_id, "ENST00000287647");
        assert_eq!(c.protein_coding, true);
    }

    #[test]
    fn works_via_from_str() {
        let c: Consequence = VEP_LINE.parse().unwrap();
        assert_eq!(c.transcript_id, "ENST00000287647");
    }

    #[test]
    fn multiple_consequence_terms() {
        let line = "-|missense_variant&splice_region_variant|MODERATE|G|ENSG1|Transcript|ENST1|protein_coding";
        let c = Consequence::parse_from_vep(line).unwrap();
        let terms: Vec<&str> = c.consequence_terms().collect();
        assert_eq!(terms, ["missense_variant", "splice_region_variant"]);
    }

    #[test]
    fn rejects_truncated_input() {
        assert!(Consequence::parse_from_vep("-|intron_variant|MODIFIER").is_err());
    }
}
