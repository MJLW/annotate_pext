use std::{
    error::Error,
    fs::File,
    io::BufWriter,
    path::{Path, PathBuf},
};

use annotate_pext::{consequence::Consequence, gtex_table::GTExTable};
use clap::{Parser, ValueEnum};
use noodles_vcf::{
    header::record::value::{
        Map,
        map::{
            Info,
            info::{Number, Type},
        },
    },
    variant::{
        RecordBuf,
        io::Write,
        record_buf::info::field::{Value, value::Array},
    },
};

#[derive(Debug, Clone, ValueEnum)]
#[value(rename_all = "verbatim")]
enum CSQCaller {
    VEP,
    BCSQ,
}

impl CSQCaller {
    fn default_csq_tag(&self) -> &'static str {
        match self {
            CSQCaller::VEP => "CSQ",
            CSQCaller::BCSQ => "BCSQ",
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(author = "MJLW", version = "0.0.1", about = "", long_about = "")]
struct Args {
    #[arg(long)]
    variants: PathBuf,

    #[arg(long)]
    csq_caller: CSQCaller,

    /// Consequence INFO tag. [Optional, defaults to: VEP=CSQ or BCSQ=BCSQ])
    #[arg(long)]
    csq_tag: Option<String>,

    #[arg(long)]
    tpms: PathBuf,

    #[arg(long, default_value = "PEXT")]
    output_tag: String,

    #[arg(long)]
    output: PathBuf,
}

impl Args {
    fn csq_tag(&self) -> String {
        self.csq_tag
            .clone()
            .unwrap_or_else(|| self.csq_caller.default_csq_tag().to_string())
    }
}

fn get_csq_values<'a, S: AsRef<str>>(
    record: &RecordBuf,
    csq_tag: S,
) -> Result<Vec<String>, Box<dyn Error>> {
    let tag = csq_tag.as_ref();

    let info = record.info();

    let value = info
        .get(tag)
        .ok_or_else(|| format!("Failed to find consequence tag '{tag}' for record."))?
        .ok_or_else(|| format!("No annotated value for consequence tag '{tag}' for record."))?;

    let values = match value {
        Value::Array(Array::String(values)) => values
            .iter()
            .map(|s| s.clone().unwrap().to_string())
            .collect::<Vec<String>>(),
        Value::String(s) => vec![s.to_string()],

        other => return Err(format!("CSQ tag '{tag}' had unexpected type: {other:?}").into()),
    };

    Ok(values)
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

fn annotate_vcf<P: AsRef<Path>, S: AsRef<str>>(
    table: &GTExTable,
    variants: P,
    caller: &CSQCaller,
    csq_tag: S,
    output_tag: S,
    output: P,
) -> Result<(), Box<dyn Error>> {
    let mut reader =
        noodles_vcf::io::reader::Builder::default().build_from_path(variants.as_ref())?;

    // Add PEXT score to header
    let mut header = reader.read_header()?;
    let info_def = Map::<Info>::new(Number::Unknown, Type::Float, "PEXT score");
    header
        .infos_mut()
        .insert(output_tag.as_ref().to_string(), info_def);

    let mut writer = File::create(output)
        .map(BufWriter::new)
        .map(noodles_vcf::io::Writer::new)?;
    writer.write_header(&header)?;

    for result in reader.record_bufs(&header) {
        let mut record = result?;

        let values = get_csq_values(&record, &csq_tag)?;
        let annotations_result: Result<Vec<Consequence>, _> = match caller {
            CSQCaller::VEP => values
                .iter()
                .map(|s| Consequence::parse_from_vep(s.as_str()))
                .collect(),
            CSQCaller::BCSQ => values
                .iter()
                .map(|s| Consequence::parse_from_vep(s.as_str()))
                .collect(),
        };

        // Get unique annotation combinations
        let annotations = annotations_result?;
        let pc_annotations: Vec<&Consequence> =
            annotations.iter().filter(|a| a.protein_coding).collect();
        let mut unique_pc_annotations: Vec<&[String]> = pc_annotations
            .iter()
            .map(|a| a.group_columns.as_slice())
            .collect();

        unique_pc_annotations.sort();
        unique_pc_annotations.dedup();

        let mut unique_genes: Vec<&str> = annotations
            .iter()
            .filter(|a| a.protein_coding)
            .map(|a| a.gene_id.as_str())
            .collect();

        unique_genes.sort();
        unique_genes.dedup();

        if unique_pc_annotations.len() == 0 {
            // No protein coding consequences, add empty annotation (".")
            writer.write_variant_record(&header, &record)?;
            continue;
        }

        let mut combination_scores: Vec<(Vec<String>, Option<f32>)> =
            Vec::with_capacity(unique_pc_annotations.len());
        for gene in unique_genes {
            let gene_combinations: Vec<&&[String]> = unique_pc_annotations
                .iter()
                .filter(|a| a.get(0).unwrap() == gene)
                .collect();

            if gene_combinations.len() == 0 {
                // No protein coding consequences for this gene
                continue;
            }

            // This one is a bit messy, but unfortunately necessary
            // The unwraps can be done without checks because the values cannot logically be missing
            // Because BCSQ doesn't have ENSEMBL gene ids, we have to go from transcript -> gene, so
            // we just select the first transcript and go to gene from that
            let first_transcript_id = pc_annotations
                .iter()
                .find(|a| a.gene_id == *gene_combinations.first().unwrap().get(0).unwrap())
                .unwrap()
                .transcript_id
                .to_string();

            let gene_id = table.get_transcript_gene(&first_transcript_id)?;

            let gene_tpms = table.get_gene_transcript_tpms(&gene_id)?;
            let summed_gene_tpms = column_sums(&gene_tpms);

            for combination in gene_combinations {
                let combination_annotations: Vec<&&Consequence> = pc_annotations
                    .iter()
                    .filter(|a| &a.group_columns == combination)
                    .collect();

                if combination_annotations.len() == 0 {
                    combination_scores.push((combination.to_vec(), None));
                    continue;
                }

                let combination_tpms: Result<Vec<&Vec<f32>>, _> = combination_annotations
                    .iter()
                    .map(|a| table.get_transcript_tpms(&a.transcript_id))
                    .collect();

                let summed_combination_tpms: Vec<f32> = column_sums(&combination_tpms?);

                let ratios: Vec<f32> = summed_combination_tpms
                    .iter()
                    .zip(&summed_gene_tpms)
                    .filter(|&(_, &g)| g > 0.0)
                    .map(|(c, g)| c / g)
                    .collect();

                let score: f32 = ratios.iter().sum::<f32>() / ratios.len() as f32;

                if score.is_nan() {
                    combination_scores.push((combination.to_vec(), None));
                    continue;
                }
                combination_scores.push((combination.to_vec(), Some(score)));
            }
        }

        // Match back scores to original annotation order
        let annotated_scores: Vec<Option<f32>> = annotations
            .iter()
            .map(|a| {
                combination_scores
                    .iter()
                    .find(|s| (s.0 == a.group_columns) && a.protein_coding)
            })
            .map(|s| match Some(s) {
                Some(Some(annotation_score)) => annotation_score.1,
                Some(_) => None,
                None => None,
            })
            .collect();

        let formatted_scores: Vec<Option<String>> = annotated_scores
            .iter()
            .map(|v| v.map(|x| format!("{x:.2}")))
            .collect();

        let value = noodles_vcf::variant::record_buf::info::field::Value::Array(
            noodles_vcf::variant::record_buf::info::field::value::Array::String(formatted_scores),
        );

        record
            .info_mut()
            .insert(String::from(output_tag.as_ref()), Some(value));
        writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let table = GTExTable::read(args.clone().tpms)?;
    annotate_vcf(
        &table,
        args.clone().variants,
        &args.csq_caller,
        &args.csq_tag(),
        &args.output_tag,
        args.output,
    )?;

    Ok(())
}
