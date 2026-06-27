use derive_more::{Deref, DerefMut, From, IntoIterator};

use crate::consequence::Consequence;

#[derive(Debug, Clone, PartialEq, Eq, Deref, DerefMut, IntoIterator, From)]
#[deref(forward)]
pub struct Consequences(#[into_iterator(owned, ref, ref_mut)] Vec<Consequence>);

impl FromIterator<Consequence> for Consequences {
    fn from_iter<I: IntoIterator<Item = Consequence>>(iter: I) -> Self {
        Consequences(iter.into_iter().collect())
    }
}

impl Consequences {
    pub fn unique_genes(&self) -> Vec<&str> {
        let mut gene_ids: Vec<&str> = self.iter().map(|c| c.gene_id.as_str()).collect();

        gene_ids.sort();
        gene_ids.dedup();

        gene_ids
    }

    pub fn unique_transcripts(&self) -> Vec<&str> {
        let mut transcript_ids: Vec<&str> = self.iter().map(|c| c.transcript_id.as_str()).collect();

        println!("{:?}", transcript_ids);

        transcript_ids.sort();
        transcript_ids.dedup();

        println!("{:?}", transcript_ids);

        transcript_ids
    }

    pub fn unique_groups(&self) -> Vec<Vec<&str>> {
        let mut groups: Vec<Vec<&str>> = self
            .iter()
            .map(|c| c.group_columns.iter().map(|gc| gc.as_str()).collect())
            .collect();

        groups.sort();
        groups.dedup();

        groups
    }
}
