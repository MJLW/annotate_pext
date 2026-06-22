use std::{collections::HashMap, error::Error, fs::File, path::Path};

use csv::StringRecord;
use rustc_hash::FxHashSet;
use serde::Deserialize;

use crate::utils::build_tsv_reader;

#[derive(Debug, Deserialize)]
pub struct GTExRowKey {
    pub transcript_id: String,
    pub gene_id: String,
}

pub struct GTExRow {
    pub key: GTExRowKey,
    pub tpms: Vec<f32>,
}

pub struct GTExTable {
    rows: Vec<GTExRow>,
    // Tissue order is the same for each row (columns)
    tissues: Vec<String>,
}

fn strip_ensembl_version<S: AsRef<str>>(s: S) -> Result<String, Box<dyn Error>> {
    Ok(s.as_ref()
        .split_once(".")
        .map(|(id, _)| id)
        .ok_or_else(|| format!("Failed to remove version from ENSEMBL ID: {}", s.as_ref()))?
        .to_string())
}

impl GTExTable {
    pub fn create_from_gtex<P: AsRef<Path>, S: AsRef<str>>(
        path: P,
        samples_per_tissue: &HashMap<String, Vec<String>>,
        coding_transcripts: &[S],
    ) -> Result<Self, Box<dyn Error>> {
        let mut rdr = build_tsv_reader(path)?;

        // Get and validate header
        let headers: StringRecord = rdr.headers()?.clone();
        let header_index: HashMap<&str, usize> =
            headers.iter().enumerate().map(|(i, h)| (h, i)).collect();
        let tissue_indices: Vec<(&str, Vec<usize>)> = samples_per_tissue
            .iter()
            .map(|(tissue, samples)| {
                let idxs = samples
                    .iter()
                    .filter_map(|s| header_index.get(s.as_str()).copied())
                    .collect::<Vec<usize>>();
                (tissue.as_str(), idxs)
            })
            .filter(|(_, samples)| samples.len() > 0)
            .collect();

        // Create HashSet for coding_transcripts for quick lookup
        let hashed_coding_transcripts: FxHashSet<String> = coding_transcripts
            .into_iter()
            .map(|s| s.as_ref().to_string())
            .collect();

        // Get median TPM per tissue for each coding transcript
        let mut rows: Vec<GTExRow> = Vec::new();
        let mut record = csv::StringRecord::new();
        while rdr.read_record(&mut record)? {
            let key: GTExRowKey = record.deserialize(Some(&headers))?;

            // Not coding, skip!
            if !hashed_coding_transcripts
                .contains(strip_ensembl_version(key.transcript_id.clone())?.as_str())
            {
                continue;
            }

            let mut medians: Vec<f32> = Vec::with_capacity(tissue_indices.len());

            // Calculate median per tissue
            for (_, idxs) in &tissue_indices {
                let mut tpms: Vec<f32> = idxs
                    .iter()
                    .map(|&i| record[i].parse::<f32>())
                    .collect::<Result<_, _>>()?;
                tpms.sort_by(f32::total_cmp);

                let median = tpms[tpms.len() / 2];
                medians.push(median);
            }

            let unversioned_key: GTExRowKey = GTExRowKey {
                transcript_id: strip_ensembl_version(key.transcript_id)?,
                gene_id: strip_ensembl_version(key.gene_id)?,
            };

            rows.push(GTExRow {
                key: unversioned_key,
                tpms: medians,
            });
        }

        let tissues: Vec<String> = tissue_indices
            .iter()
            .map(|(tissue, _)| tissue.to_string())
            .collect();

        Ok(Self { rows, tissues })
    }

    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let mut rdr = build_tsv_reader(path)?;

        // Get and validate header
        let headers: StringRecord = rdr.headers()?.clone();
        let tissues: Vec<String> = headers.into_iter().map(|s| s.to_string()).skip(2).collect();

        let mut rows: Vec<GTExRow> = Vec::new();
        let mut record = csv::StringRecord::new();
        while rdr.read_record(&mut record)? {
            let key: GTExRowKey = record.deserialize(Some(&headers))?;

            let tpms: Vec<f32> = record
                .into_iter()
                .skip(2)
                .map(|s| s.parse::<f32>().unwrap())
                .collect();

            rows.push(GTExRow { key, tpms });
        }

        Ok(Self { rows, tissues })
    }

    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn Error>> {
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(File::create(&path)?);

        // Write header
        wtr.write_field("gene_id")?;
        wtr.write_field("transcript_id")?;

        for tissue in &self.tissues {
            wtr.write_field(tissue)?;
        }

        wtr.write_record(None::<&[u8]>)?;

        // Write data
        for row in &self.rows {
            wtr.write_field(row.key.gene_id.to_string())?;
            wtr.write_field(row.key.transcript_id.to_string())?;

            for tissue_tpm in &row.tpms {
                wtr.write_field(format!("{:.2}", tissue_tpm))?;
            }

            wtr.write_record(None::<&[u8]>)?;
        }

        wtr.flush()?;

        Ok(())
    }

    pub fn get_gene_transcripts<S: AsRef<str>>(&self, gene_id: S) -> Vec<&GTExRow> {
        self.rows
            .iter()
            .filter(|row| row.key.gene_id == gene_id.as_ref())
            .collect()
    }

    pub fn get_transcript<S: AsRef<str>>(
        &self,
        transcript_id: S,
    ) -> Result<&GTExRow, Box<dyn Error>> {
        let transcript_row: &GTExRow = self
            .rows
            .iter()
            .find(|row| row.key.transcript_id == transcript_id.as_ref())
            .ok_or_else(|| {
                format!(
                    "Could not find transcript `{}` in the condensed GTEx table.",
                    transcript_id.as_ref()
                )
            })?;

        Ok(transcript_row)
    }
}
