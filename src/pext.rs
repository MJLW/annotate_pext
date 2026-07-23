use std::{collections::HashMap, error::Error};

use crate::{consequences::Consequences, gtex_table::GTExTable};

trait ToStringVec {
    fn to_string_vec(&self) -> Vec<String>;
}

impl ToStringVec for [&str] {
    fn to_string_vec(&self) -> Vec<String> {
        self.iter().map(|s| s.to_string()).collect()
    }
}

fn column_sums(row_matrix: &[&Vec<f32>]) -> Vec<f32> {
    row_matrix.iter().fold(
        vec![0.0; row_matrix.first().map_or(0, |row| row.len())], // Get number of columns
        |mut acc, row| {
            acc.iter_mut().zip(*row).for_each(|(a, v)| *a += v);
            acc
        },
    )
}

pub fn calculate_pext(
    annotations: Consequences,
    table: &GTExTable,
) -> Result<Option<Vec<Option<f32>>>, Box<dyn Error>> {
    let protein_coding_annotations: Consequences = annotations
        .clone()
        .into_iter()
        .filter(|a| a.protein_coding)
        .collect();

    let annotation_groups = protein_coding_annotations.unique_groups();
    let genes = protein_coding_annotations.unique_genes();
    let transcript_csq_pairs: Vec<(&str, &str)> = protein_coding_annotations
        .iter()
        .map(|a| (a.transcript_id.as_str(), a.consequence.as_str()))
        .collect();

    let transcript_csq_counts =
        transcript_csq_pairs
            .into_iter()
            .fold(HashMap::<(&str, &str), usize>::new(), |mut m, x| {
                *m.entry(x).or_default() += 1;
                m
            });

    let has_incorrectly_split_multiallelic =
        transcript_csq_counts.iter().any(|(_, &count)| count > 1);

    if annotation_groups.len() == 0 || has_incorrectly_split_multiallelic {
        return Ok(None);
    }

    let mut annotation_scores: Vec<(Vec<String>, Option<f32>)> =
        Vec::with_capacity(annotation_groups.len());
    for gene in genes {
        let gene_annotations: Vec<&Vec<&str>> = annotation_groups
            .iter()
            .filter(|a| a.last().unwrap() == &gene)
            .collect();

        if gene_annotations.len() == 0 {
            // No protein coding consequences for this gene
            continue;
        }

        let gene_id = gene_annotations
            .first()
            .unwrap()
            .last()
            .unwrap()
            .to_string();

        let gene_tpms = table.get_gene_transcript_tpms(&gene_id)?;
        let summed_gene_tpms = column_sums(&gene_tpms);

        for annotation in gene_annotations {
            let owned_annotation = annotation.to_string_vec();

            // Get unique transcripts
            let mut transcripts: Vec<&str> = protein_coding_annotations
                .iter()
                .filter(|a| a.group_columns == owned_annotation)
                .map(|a| a.transcript_id.as_str())
                .collect();

            transcripts.sort();
            transcripts.dedup();

            if transcripts.len() == 0 {
                annotation_scores.push((owned_annotation, None));
                continue;
            }

            let annotation_tpms_result: Result<Vec<&Vec<f32>>, _> = transcripts
                .into_iter()
                .map(|transcript| table.get_transcript_tpms(transcript))
                .collect();

            if !annotation_tpms_result.is_ok() {
                annotation_scores.push((owned_annotation, None));
                continue;
            }

            let annotation_tpms = annotation_tpms_result?;

            let summed_annotation_tpms: Vec<f32> = column_sums(&annotation_tpms);

            let ratios: Vec<f32> = summed_annotation_tpms
                .iter()
                .zip(&summed_gene_tpms)
                .filter(|&(_, &g)| g != 0.0)
                .map(|(c, g)| c / g)
                .collect();

            let score: f32 = ratios.iter().sum::<f32>() / ratios.len() as f32;

            if score.is_nan() {
                annotation_scores.push((owned_annotation, None));
                continue;
            }
            annotation_scores.push((owned_annotation, Some(score)));
        }
    }

    // Match back scores to original annotation order
    let annotated_scores: Vec<Option<f32>> = annotations
        .iter()
        .map(|a| {
            annotation_scores
                .iter()
                .find(|s| (s.0 == a.group_columns) && a.protein_coding)
        })
        .map(|s| match Some(s) {
            Some(Some(annotation_score)) => annotation_score.1,
            Some(_) => None,
            None => None,
        })
        .collect();

    Ok(Some(annotated_scores))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_column_sums() {
        let row_matrix = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];

        assert_eq!(
            column_sums(&row_matrix.iter().collect::<Vec<&Vec<f32>>>()),
            vec![12.0, 15.0, 18.0]
        );
    }
}
