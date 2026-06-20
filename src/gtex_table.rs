use std::{collections::HashMap, error::Error, fs::File, path::Path};

use csv::StringRecord;
use rustc_hash::FxHashSet;
use serde::Deserialize;

use crate::utils::build_tsv_reader;

#[derive(Debug, Deserialize)]
struct GTExRowKey {
    transcript_id: String,
    gene_id: String,
}

struct GTExRow {
    key: GTExRowKey,
    tpms: Vec<f32>,
}

pub struct GTExTable {
    rows: Vec<GTExRow>,
    // Tissue order is the same for each row (columns)
    tissues: Vec<String>,
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
            if !hashed_coding_transcripts.contains(key.transcript_id.as_str()) {
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

            rows.push(GTExRow {
                key: key,
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
}
