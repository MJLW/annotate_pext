use std::{
    collections::HashMap,
    error::Error,
    fs::read_to_string,
    path::{Path, PathBuf},
};

use clap::Parser;
use serde::Deserialize;

use annotate_pext::gtex_table::GTExTable;
use annotate_pext::utils::build_tsv_reader;

#[derive(Parser, Debug)]
#[command(
    author = "MJLW",
    version = "0.0.1",
    about = "Creates a condensed TPM matrix from the GTEx TPMs to speed up PEXT annotations.",
    long_about = "Creates a condensed TPM matrix from the GTEx TPMs to speed up PEXT annotations. Takes median TPM over tissues, resulting in a single value per tissue."
)]
struct Args {
    #[arg(long)]
    gtex_tpms: PathBuf,

    #[arg(long)]
    gtex_sample_attributes: PathBuf,

    #[arg(long)]
    tissue_blacklist: PathBuf,

    #[arg(long)]
    transcript_whitelist: Option<PathBuf>,

    #[arg(long, default_value_t = default_min_samples())]
    min_samples_per_tissue: usize,

    // gff3: PathBuf,
    #[arg(long)]
    output: PathBuf,
}

fn default_min_samples() -> usize {
    100
}

#[derive(Debug, PartialEq, Eq, Deserialize)]
struct SampleAttribute {
    #[serde(rename = "SAMPID")]
    sample: String,

    #[serde(rename = "SMTSD")]
    tissue: String,
}

fn read_lines<P: AsRef<Path>>(path: P) -> Result<Vec<String>, Box<dyn Error>> {
    Ok(read_to_string(path)?
        .lines()
        .map(|s| s.to_string())
        .collect::<Vec<String>>())
}

fn read_sample_attributes<P: AsRef<Path>, S: AsRef<str>>(
    path: P,
    tissue_blacklist: &[S],
) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    let mut rdr = build_tsv_reader(path)?;

    let sample_attributes: Result<Vec<SampleAttribute>, _> =
        rdr.deserialize::<SampleAttribute>().collect();

    let filtered: Vec<SampleAttribute> = sample_attributes?
        .into_iter()
        .filter(|attribute| {
            !tissue_blacklist
                .iter()
                .any(|t| t.as_ref() == attribute.tissue.as_str())
        })
        .collect();

    let samples_per_tissue: HashMap<String, Vec<String>> =
        filtered.iter().fold(HashMap::new(), |mut acc, attribute| {
            let group = acc
                .entry(attribute.tissue.to_string())
                .or_insert(Vec::new());
            group.push(attribute.sample.to_string());
            acc
        });

    Ok(samples_per_tissue)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let blacklisted_tissues: Vec<String> = read_lines(args.tissue_blacklist)?;
    let samples_per_tissue =
        read_sample_attributes(args.gtex_sample_attributes, &blacklisted_tissues)?;

    // Should we filter to coding transcripts or not?
    let transcript_whitelist: Option<Vec<String>> =
        if let Some(file_path) = args.transcript_whitelist {
            Some(read_lines(file_path)?)
        } else {
            None
        };

    let table = GTExTable::create_from_gtex(
        args.gtex_tpms,
        &samples_per_tissue,
        transcript_whitelist,
        args.min_samples_per_tissue.try_into()?,
    )?;
    table.write(args.output)?;

    Ok(())
}
