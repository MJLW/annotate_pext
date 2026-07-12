use std::{error::Error, fs::File, io::Read, path::Path};

use csv::{Reader, ReaderBuilder};
use flate2::read::GzDecoder;

pub fn build_tsv_reader<P: AsRef<Path>>(path: P) -> Result<Reader<Box<dyn Read>>, Box<dyn Error>> {
    let file = File::open(&path)?;
    let file_type_reader: Box<dyn Read> = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("gz") => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };

    let rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(file_type_reader);

    Ok(rdr)
}
