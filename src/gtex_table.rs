use std::{collections::HashMap, error::Error, fs::File, path::Path};

use csv::StringRecord;
use rustc_hash::{FxHashMap, FxHashSet};
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
    transcript_tpms: FxHashMap<String, Vec<f32>>,
    // Tissue order is the same for each row (columns)
    pub tissues: Vec<String>,
    pub keys: Vec<GTExRowKey>,
    gene_transcript_map: FxHashMap<String, Vec<String>>,
    transcript_gene_map: FxHashMap<String, String>,
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
        transcript_whitelist: Option<Vec<S>>,
        min_samples_per_tissue: usize,
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
            .filter(|(_, samples)| samples.len() > min_samples_per_tissue)
            .collect();

        // Create HashSet for coding_transcripts for quick lookup
        let has_transcript_whitelist = transcript_whitelist.is_some();
        let hashed_transcript_whitelist: FxHashSet<String> =
            if let Some(result) = transcript_whitelist {
                result.into_iter().map(|s| s.as_ref().to_string()).collect()
            } else {
                FxHashSet::default()
            };

        // Get median TPM per tissue for each coding transcript
        let mut keys: Vec<GTExRowKey> = Vec::new();
        let mut transcript_tpms: FxHashMap<String, Vec<f32>> = FxHashMap::default();
        let mut record = csv::StringRecord::new();
        while rdr.read_record(&mut record)? {
            let key: GTExRowKey = record.deserialize(Some(&headers))?;

            // If coding transcripts passed and not coding, skip!
            if has_transcript_whitelist
                && hashed_transcript_whitelist
                    .contains(strip_ensembl_version(key.transcript_id.clone())?.as_str())
            {
                continue;
            }

            let mut medians: Vec<f32> = Vec::with_capacity(tissue_indices.len());

            // Calculate median per tissue
            for (_, idxs) in &tissue_indices {
                let mut tissue_tpms: Vec<f32> = idxs
                    .iter()
                    .map(|&i| record[i].parse::<f32>())
                    .collect::<Result<_, _>>()?;
                tissue_tpms.sort_by(f32::total_cmp);

                let tissue_median = tissue_tpms[tissue_tpms.len() / 2];
                medians.push(tissue_median);
            }

            let unversioned_key: GTExRowKey = GTExRowKey {
                transcript_id: strip_ensembl_version(key.transcript_id)?,
                gene_id: strip_ensembl_version(key.gene_id)?,
            };

            transcript_tpms.insert(unversioned_key.transcript_id.to_string(), medians);
            keys.push(unversioned_key);
        }

        let tissues: Vec<String> = tissue_indices
            .iter()
            .map(|(tissue, _)| tissue.to_string())
            .collect();

        let gene_transcript_map: FxHashMap<String, Vec<String>> =
            keys.iter().fold(FxHashMap::default(), |mut acc, key| {
                acc.entry(key.gene_id.to_string())
                    .or_default()
                    .push(key.transcript_id.to_string());
                acc
            });

        let transcript_gene_map: FxHashMap<String, String> = keys
            .iter()
            .map(|key| (key.transcript_id.to_string(), key.gene_id.to_string()))
            .collect();

        Ok(Self {
            transcript_tpms,
            tissues,
            keys,
            gene_transcript_map,
            transcript_gene_map,
        })
    }

    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let mut rdr = build_tsv_reader(path)?;

        // Get and validate header
        let headers: StringRecord = rdr.headers()?.clone();
        let tissues: Vec<String> = headers.into_iter().map(|s| s.to_string()).skip(2).collect();

        let mut keys: Vec<GTExRowKey> = Vec::new();
        let mut transcript_tpms: FxHashMap<String, Vec<f32>> = FxHashMap::default();
        let mut record = csv::StringRecord::new();
        while rdr.read_record(&mut record)? {
            let key: GTExRowKey = record.deserialize(Some(&headers))?;

            let tpms: Vec<f32> = record
                .into_iter()
                .skip(2)
                .map(|s| s.parse::<f32>().unwrap())
                .collect();

            transcript_tpms.insert(key.transcript_id.to_string(), tpms);
            keys.push(key);
        }

        let gene_transcript_map: FxHashMap<String, Vec<String>> =
            keys.iter().fold(FxHashMap::default(), |mut acc, key| {
                acc.entry(key.gene_id.to_string())
                    .or_default()
                    .push(key.transcript_id.to_string());
                acc
            });

        let transcript_gene_map: FxHashMap<String, String> = keys
            .iter()
            .map(|key| (key.transcript_id.to_string(), key.gene_id.to_string()))
            .collect();

        Ok(Self {
            transcript_tpms,
            tissues,
            keys,
            gene_transcript_map,
            transcript_gene_map,
        })
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
        for (key, tpms) in self
            .keys
            .iter()
            .map(|key| (key, self.transcript_tpms.get(&key.transcript_id).unwrap()))
        {
            wtr.write_field(key.gene_id.to_string())?;
            wtr.write_field(key.transcript_id.to_string())?;

            for tissue_tpm in tpms {
                wtr.write_field(format!("{:.2}", tissue_tpm))?;
            }

            wtr.write_record(None::<&[u8]>)?;
        }

        wtr.flush()?;

        Ok(())
    }

    pub fn get_transcript_tpms<S: AsRef<str>>(
        &self,
        transcript_id: S,
    ) -> Result<&Vec<f32>, Box<dyn Error>> {
        let transcript_tpms = self
            .transcript_tpms
            .get(transcript_id.as_ref())
            .ok_or_else(|| {
                format!(
                    "Could not find transcript `{}` in the condensed GTEx table.",
                    transcript_id.as_ref()
                )
            })?;

        Ok(transcript_tpms)
    }

    pub fn get_gene_transcripts(&self, gene_id: &String) -> Result<&Vec<String>, Box<dyn Error>> {
        let transcript_ids: &Vec<String> =
            self.gene_transcript_map.get(gene_id).ok_or_else(|| {
                format!(
                    "Could not find gene_id `{}` in the condensed GTEx table.",
                    gene_id
                )
            })?;

        Ok(transcript_ids)
    }

    pub fn get_gene_transcript_tpms(
        &self,
        gene_id: &String,
    ) -> Result<Vec<&Vec<f32>>, Box<dyn Error>> {
        let transcript_ids: &Vec<String> = self.get_gene_transcripts(gene_id)?;

        transcript_ids
            .iter()
            .map(|transcript_id| self.get_transcript_tpms(transcript_id))
            .collect()
    }

    pub fn get_transcript_gene<S: AsRef<str>>(
        &self,
        transcript_id: S,
    ) -> Result<&String, Box<dyn Error>> {
        let gene_id: &String = self
            .transcript_gene_map
            .get(transcript_id.as_ref())
            .ok_or_else(|| {
                format!(
                    "Could not find transcript `{}` in the condensed GTEx table.",
                    transcript_id.as_ref()
                )
            })?;

        Ok(gene_id)
    }
}
