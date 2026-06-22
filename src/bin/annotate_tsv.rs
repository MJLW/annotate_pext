use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    path::{Path, PathBuf},
};

use annotate_pext::{gtex_table::GTExRow, gtex_table::GTExTable, utils::build_tsv_reader};
use clap::Parser;
use csv::StringRecord;

#[derive(Parser, Debug)]
#[command(author = "MJLW", version = "0.0.1", about = "", long_about = "")]
struct Args {
    #[arg(long)]
    variants: PathBuf,

    #[arg(long, value_delimiter = ',')]
    variant_columns: Vec<String>,

    #[arg(long)]
    transcript_id_column: String,

    #[arg(long)]
    gene_id_column: String,

    #[arg(long)]
    biotype_column: String,

    #[arg(long, value_delimiter = ',')]
    group_columns: Vec<String>,

    #[arg(long)]
    tpms: PathBuf,

    #[arg(long)]
    output: PathBuf,
}

struct GroupingColumns {
    variant_columns: Vec<String>,
    transcript_column: String,
    gene_column: String,
    biotype_column: String,
    group_columns: Vec<String>,
}

#[derive(Debug)]
struct GroupingData {
    variant_values: Vec<String>,
    transcript_id: String,
    gene_id: String,
    protein_coding: bool,
    group_values: Vec<String>,
}

fn find_header_index<S: AsRef<str>>(
    header: HashMap<&str, usize>,
    column: S,
) -> Result<usize, Box<dyn Error>> {
    let index: usize = *header
        .iter()
        .find(|(col_name, _)| column.as_ref() == **col_name)
        .map(|(_, index)| index)
        .ok_or(format!(
            "Could not find column '{}' in TSV.",
            column.as_ref()
        ))?;

    Ok(index)
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

// TODO: Rewrite this method, it is a mess.
fn annotate_tsv<P: AsRef<Path>>(
    input_path: P,
    grouping: &GroupingColumns,
    table: &GTExTable,
    output_path: P,
) -> Result<(), Box<dyn Error>> {
    let mut rdr = build_tsv_reader(input_path)?;

    let header = rdr.headers()?.clone();
    let header_index: HashMap<&str, usize> =
        header.iter().enumerate().map(|(i, h)| (h, i)).collect();

    // TODO: Figure out a way to do this without cloning the index every time.
    let variant_indices_result: Result<Vec<usize>, _> = grouping
        .variant_columns
        .iter()
        .map(|variant_column| find_header_index(header_index.clone(), variant_column))
        .collect();
    let variant_indices = variant_indices_result?;

    let gene_index: usize = find_header_index(header_index.clone(), &grouping.gene_column)?;

    let transcript_index: usize =
        find_header_index(header_index.clone(), &grouping.transcript_column)?;

    let biotype_index: usize = find_header_index(header_index.clone(), &grouping.biotype_column)?;

    let group_indices_result: Result<Vec<usize>, _> = grouping
        .group_columns
        .iter()
        .map(|group_column| find_header_index(header_index.clone(), group_column))
        .collect();
    let group_indices = group_indices_result?;

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(File::create(&output_path)?);

    for s in &header {
        wtr.write_field(s)?;
    }
    wtr.write_field("PEXT")?;
    wtr.write_record(None::<&[u8]>)?;
    wtr.flush()?;

    let mut history: Vec<(GroupingData, StringRecord)> = Vec::new();
    for result in rdr.records() {
        let record = result?;

        // Get variant information
        let variant_values_result: Result<Vec<&str>, _> = variant_indices
            .iter()
            .map(|i| {
                record
                    .get(*i)
                    .ok_or("Could not find a variant column value, your column may not exist or the TSV may be malformed")
            })
            .collect();

        let variant_values: Vec<String> = variant_values_result?
            .into_iter()
            .map(|s| s.to_string())
            .collect();

        // TODO: This should be done using peeking.
        if history.len() > 0 {
            let previous: &(GroupingData, StringRecord) = history.last().unwrap();
            if previous.0.variant_values != variant_values {
                // History contains all rows for a specific variant because we found a new one
                let protein_coding_annotations: Vec<&(GroupingData, StringRecord)> = history
                    .iter()
                    .filter(|(group, _)| group.protein_coding)
                    .collect();

                let mut group_values: Vec<Vec<String>> = protein_coding_annotations
                    .iter()
                    .map(|(group, _)| {
                        let mut group_values = group.group_values.clone();
                        group_values.push(group.gene_id.to_string());
                        group_values
                    })
                    .collect();

                group_values.sort();
                group_values.dedup();

                if group_values.len() == 0 {
                    for (_, record) in &history {
                        for s in record {
                            wtr.write_field(s)?;
                        }
                        wtr.write_field(".")?;
                        wtr.write_record(None::<&[u8]>)?;
                    }
                    wtr.flush()?;
                    history.clear();
                } else {
                    // Already tested that there is at least one group, so we can just unwrap
                    let gene_id = group_values.first().unwrap().last().unwrap();
                    let gene_transcripts: Vec<&GTExRow> = table.get_gene_transcripts(gene_id);

                    let gene_transcript_tpms: Vec<&Vec<f32>> =
                        gene_transcripts.iter().map(|row| &row.tpms).collect();

                    let gene_tpms: Vec<f32> = column_sums(&gene_transcript_tpms);

                    let grouped_transcript_ids: Vec<Vec<&str>> = group_values
                        .iter()
                        .map(|unique_group| {
                            let transcript_ids: Vec<&str> = protein_coding_annotations
                                .iter()
                                .map(|(group, _)| {
                                    let mut group_values = group.group_values.clone();
                                    group_values.push(group.gene_id.to_string());
                                    (group_values, group.transcript_id.as_str())
                                })
                                .filter(|(group, _)| unique_group == group)
                                .map(|(_, transcript_id)| transcript_id)
                                .collect();

                            transcript_ids
                        })
                        .collect();

                    let mut group_scores: Vec<Option<f32>> = Vec::with_capacity(group_values.len());
                    for transcript_ids in grouped_transcript_ids {
                        let group_transcript_tpms: Vec<&Vec<f32>> = gene_transcripts
                            .iter()
                            .filter(|row| transcript_ids.contains(&row.key.transcript_id.as_str()))
                            .map(|row| &row.tpms)
                            .collect();

                        let group_tpms: Vec<f32> = column_sums(&group_transcript_tpms);

                        let ratios: Vec<f32> = group_tpms
                            .iter()
                            .zip(&gene_tpms)
                            .filter(|&(_, &gene)| gene > 0.0)
                            .map(|(group, gene)| group / gene)
                            .collect();

                        let score: f32 = ratios.iter().sum::<f32>() / ratios.len() as f32;

                        if score.is_nan() {
                            group_scores.push(None);
                            continue;
                        }

                        group_scores.push(Some(score));
                    }

                    for (grouping, record) in &history {
                        let mut group = grouping.group_values.clone();
                        group.push(grouping.gene_id.clone());

                        let score: Option<f32> = group_values
                            .iter()
                            .zip(&group_scores)
                            .find(|(score_group, _)| group == **score_group)
                            .map(|(_, &score)| score)
                            .unwrap_or(None);

                        for s in record {
                            wtr.write_field(s)?;
                        }

                        if score.is_some() {
                            wtr.write_field(format!("{:.2}", score.unwrap()))?;
                        } else {
                            wtr.write_field(".")?;
                        }
                        wtr.write_record(None::<&[u8]>)?;
                    }
                    wtr.flush()?;

                    history.clear();
                }
            }
        }

        // Get annotations
        let gene: &str = record.get(gene_index).ok_or(
            "No gene value found, your column may not exist or your TSV may be malformed.",
        )?;

        let transcript: &str = record.get(transcript_index).ok_or(
            "No transcript value found, your column may not exist or your TSV may be malformed.",
        )?;

        let biotype: &str = record.get(biotype_index).ok_or(
            "No biotype value found, your column may not exist or your TSV may be malformed.",
        )?;

        let group_values_result: Result<Vec<&str>, _> = group_indices
            .iter()
            .map(|i| {
                record
                    .get(*i)
                    .ok_or("Could not find a group column value, your column may not exist or your TSV may be malformed")
            })
            .collect();
        let group_values: Vec<String> = group_values_result?
            .into_iter()
            .map(|s| s.to_string())
            .collect();

        let record_keys = GroupingData {
            variant_values: variant_values,
            gene_id: gene.to_string(),
            transcript_id: transcript.to_string(),
            protein_coding: biotype == "protein_coding",
            group_values: group_values,
        };

        history.push((record_keys, record.clone()));
    }

    // Process leftover in history
    let protein_coding_annotations: Vec<&(GroupingData, StringRecord)> = history
        .iter()
        .filter(|(group, _)| group.protein_coding)
        .collect();

    let mut group_values: Vec<Vec<String>> = protein_coding_annotations
        .iter()
        .map(|(group, _)| {
            let mut group_values = group.group_values.clone();
            group_values.push(group.gene_id.to_string());
            group_values
        })
        .collect();

    group_values.sort();
    group_values.dedup();

    if group_values.len() == 0 {
        for (_, record) in &history {
            for s in record {
                wtr.write_field(s)?;
            }
            wtr.write_field(".")?;
            wtr.write_record(None::<&[u8]>)?;
        }
        wtr.flush()?;
        return Ok(());
    }

    // Already tested that there is at least one group, so we can just unwrap
    let gene_id = group_values.first().unwrap().last().unwrap();
    let gene_transcripts: Vec<&GTExRow> = table.get_gene_transcripts(gene_id);

    let gene_transcript_tpms: Vec<&Vec<f32>> =
        gene_transcripts.iter().map(|row| &row.tpms).collect();

    let gene_tpms: Vec<f32> = column_sums(&gene_transcript_tpms);

    let grouped_transcript_ids: Vec<Vec<&str>> = group_values
        .iter()
        .map(|unique_group| {
            let transcript_ids: Vec<&str> = protein_coding_annotations
                .iter()
                .map(|(group, _)| {
                    let mut group_values = group.group_values.clone();
                    group_values.push(group.gene_id.to_string());
                    (group_values, group.transcript_id.as_str())
                })
                .filter(|(group, _)| unique_group == group)
                .map(|(_, transcript_id)| transcript_id)
                .collect();

            transcript_ids
        })
        .collect();

    let mut group_scores: Vec<Option<f32>> = Vec::with_capacity(group_values.len());
    for transcript_ids in grouped_transcript_ids {
        let group_transcript_tpms: Vec<&Vec<f32>> = gene_transcripts
            .iter()
            .filter(|row| transcript_ids.contains(&row.key.transcript_id.as_str()))
            .map(|row| &row.tpms)
            .collect();

        let group_tpms: Vec<f32> = column_sums(&group_transcript_tpms);

        let ratios: Vec<f32> = group_tpms
            .iter()
            .zip(&gene_tpms)
            .filter(|&(_, &gene)| gene > 0.0)
            .map(|(group, gene)| group / gene)
            .collect();

        let score: f32 = ratios.iter().sum::<f32>() / ratios.len() as f32;

        if score.is_nan() {
            group_scores.push(None);
            continue;
        }

        group_scores.push(Some(score));
    }

    for (grouping, record) in &history {
        let mut group = grouping.group_values.clone();
        group.push(grouping.gene_id.clone());

        let score: Option<f32> = group_values
            .iter()
            .zip(&group_scores)
            .find(|(score_group, _)| group == **score_group)
            .map(|(_, &score)| score)
            .unwrap_or(None);

        for s in record {
            wtr.write_field(s)?;
        }

        if score.is_some() {
            wtr.write_field(format!("{:.2}", score.unwrap()))?;
        } else {
            wtr.write_field(".")?;
        }
        wtr.write_record(None::<&[u8]>)?;
    }
    wtr.flush()?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let grouping = GroupingColumns {
        variant_columns: args.variant_columns,
        gene_column: args.gene_id_column,
        transcript_column: args.transcript_id_column,
        biotype_column: args.biotype_column,
        group_columns: args.group_columns,
    };

    let table = GTExTable::read(args.tpms)?;
    annotate_tsv(args.variants, &grouping, &table, args.output)?;

    Ok(())
}
