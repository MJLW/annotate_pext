use std::{
    collections::HashMap,
    error::Error,
    fs::{File, read_to_string},
    io::{BufReader, Read},
    path::{Path, PathBuf},
};

use clap::Parser;
use flate2::read::GzDecoder;
use noodles_gff::io::Reader as GFF3Reader;
use serde::Deserialize;

use annotate_pext::gtex_table::GTExTable;
use annotate_pext::utils::build_tsv_reader;

const TRANSCRIPT_TYPE: &[u8] = "transcript".as_bytes();
const ATTRIBUTE_TRANSCRIPT_ID: &[u8] = "transcript_id".as_bytes();
const ATTRIIBUTE_TRANSCRIPT_TYPE: &[u8] = "transcript_type".as_bytes();

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
    coding_transcripts: PathBuf,
    // gff3: PathBuf,
    #[arg(long)]
    output: PathBuf,
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

fn read_coding_transcripts_from_gff3<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<String>, Box<dyn Error>> {
    let file = File::open(&path)?;
    let file_type_reader: Box<dyn Read> = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };

    let mut reader = GFF3Reader::new(BufReader::new(file_type_reader));

    let mut coding_transcripts: Vec<String> = Vec::new();
    for result in reader.record_bufs() {
        let record = result?;

        if record.ty() != TRANSCRIPT_TYPE {
            continue;
        }

        let transcript_type = record
            .attributes()
            .get(ATTRIIBUTE_TRANSCRIPT_TYPE)
            .ok_or_else(|| format!("Failed to parse transcript_type attribute: {:?}", record))?
            .as_string()
            .ok_or_else(|| {
                format!(
                    "Failed to parse transcript_type attribute as String: {:?}",
                    record
                )
            })?;

        if transcript_type != "protein_coding" {
            continue;
        }

        let transcript_id = record
            .attributes()
            .get(ATTRIBUTE_TRANSCRIPT_ID)
            .ok_or_else(|| format!("Failed to parse transcript_id attribute: {:?}", record))?
            .as_string()
            .ok_or_else(|| {
                format!(
                    "Failed to parse transcript_id attribute as String: {:?}",
                    record
                )
            })?;

        coding_transcripts.push(transcript_id.to_string());
    }

    Ok(coding_transcripts)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let blacklisted_tissues: Vec<String> = read_lines(args.tissue_blacklist)?;
    let samples_per_tissue =
        read_sample_attributes(args.gtex_sample_attributes, &blacklisted_tissues)?;
    // let coding_transcripts = read_coding_transcripts_from_gff3(args.gff3)?;
    let coding_transcripts = read_lines(args.coding_transcripts)?;

    let table =
        GTExTable::create_from_gtex(args.gtex_tpms, &samples_per_tissue, &coding_transcripts)?;
    table.write(args.output)?;

    Ok(())
}
