use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    path::{Path, PathBuf},
};

use annotate_pext::{
    consequence::Consequence, consequences::Consequences, gtex_table::GTExTable,
    pext::calculate_pext, utils::build_tsv_reader,
};
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
    biotype_column: String,

    #[arg(long, value_delimiter = ',')]
    group_columns: Vec<String>,

    #[arg(long)]
    tpms: PathBuf,

    #[arg(long, default_value = "PEXT")]
    output_tag: String,

    #[arg(long)]
    output: PathBuf,
}

struct GroupingColumns {
    variant_columns: Vec<String>,
    transcript_column: String,
    biotype_column: String,
    group_columns: Vec<String>,
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

// TODO: Rewrite this method, it is a mess.
fn annotate_tsv<P: AsRef<Path>, S: AsRef<str>>(
    input_path: P,
    grouping: &GroupingColumns,
    table: &GTExTable,
    output_tag: S,
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
    wtr.write_field(output_tag.as_ref())?;
    wtr.write_record(None::<&[u8]>)?;
    wtr.flush()?;

    let mut history: Vec<(Vec<String>, Consequence, StringRecord)> = Vec::new();
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
            let previous: &(Vec<String>, Consequence, StringRecord) = history.last().unwrap();
            if previous.0 != variant_values {
                // History contains all rows for a specific variant because we found a new one
                // Process existing history
                let annotations: Consequences = history
                    .iter()
                    .map(|(_, consequence, _)| consequence.to_owned())
                    .collect();

                let pext_scores = calculate_pext(annotations, table)?;

                if pext_scores.is_none() {
                    for (_, _, record) in &history {
                        for s in record {
                            wtr.write_field(s)?;
                        }
                        wtr.write_field(".")?;
                        wtr.write_record(None::<&[u8]>)?;
                    }
                    wtr.flush()?;
                    history.clear();
                } else {
                    for ((_, _, record), score) in history.iter().zip(pext_scores.unwrap()) {
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

        let gene = table.get_transcript_gene(transcript)?;

        let consequence = Consequence::from_fields(
            gene.to_string(),
            transcript.to_string(),
            biotype.to_string(),
            group_values,
        );

        history.push((variant_values, consequence, record.clone()));
    }

    let annotations: Consequences = history
        .iter()
        .map(|(_, consequence, _)| consequence.to_owned())
        .collect();

    let pext_scores = calculate_pext(annotations, table)?;

    if pext_scores.is_none() {
        for (_, _, record) in &history {
            for s in record {
                wtr.write_field(s)?;
            }
            wtr.write_field(".")?;
            wtr.write_record(None::<&[u8]>)?;
        }
        wtr.flush()?;
        history.clear();
    }

    for ((_, _, record), score) in history.iter().zip(pext_scores.unwrap()) {
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

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let grouping = GroupingColumns {
        variant_columns: args.variant_columns,
        transcript_column: args.transcript_id_column,
        biotype_column: args.biotype_column,
        group_columns: args.group_columns,
    };

    let table = GTExTable::read(args.tpms)?;
    annotate_tsv(
        args.variants,
        &grouping,
        &table,
        &args.output_tag,
        args.output,
    )?;

    Ok(())
}
