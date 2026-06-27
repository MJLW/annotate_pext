use std::{
    error::Error,
    fs::File,
    io::BufWriter,
    io::Write,
    path::{Path, PathBuf},
};

use annotate_pext::{
    consequence::Consequence,
    consequences::Consequences,
    gtex_table::GTExTable,
    record::{Annotated, format_scores},
    utils::calculate_pext,
};
use clap::{Parser, ValueEnum};
use noodles_bgzf;
use noodles_vcf::{
    Header, Record,
    header::record::value::{
        Map,
        map::{
            Info,
            info::{Number, Type},
        },
    },
    variant::{
        // RecordBuf,
        io::Write as _,
        // record_buf::info::field::{Value, value::Array},
        record::info::field::{Value, value::Array},
    },
};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

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
    record: &Record,
    header: &Header,
    csq_tag: S,
) -> Result<Option<Vec<String>>, Box<dyn Error>> {
    let tag = csq_tag.as_ref();

    let info = record.info();

    let option_value = info.get(header, tag);

    if option_value.is_none() {
        return Ok(None);
    }

    let value = option_value
        .unwrap()?
        .ok_or_else(|| format!("No annotated value for consequence tag '{tag}' for record."))?;

    let values = match value {
        Value::Array(Array::String(values)) => values
            .iter()
            .map(|s| s.unwrap().unwrap().to_string())
            .collect::<Vec<String>>(),
        Value::String(s) => vec![s.to_string()],

        other => return Err(format!("CSQ tag '{tag}' had unexpected type: {other:?}").into()),
    };

    Ok(Some(values))
}

fn build_vcf_writer<P: AsRef<Path>>(
    path: P,
) -> Result<noodles_vcf::io::Writer<Box<dyn Write>>, Box<dyn Error>> {
    let file = File::create(&path)?;
    let buf = BufWriter::with_capacity(1 << 20, file);
    let inner: Box<dyn Write> = match path.as_ref().extension().and_then(|e| e.to_str()) {
        Some("gz" | "bgz") => Box::new(noodles_bgzf::io::Writer::new(buf)),
        _ => Box::new(buf),
    };

    Ok(noodles_vcf::io::Writer::new(inner))
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

    let mut writer = build_vcf_writer(output)?;

    let mut header = reader.read_header()?;
    let info_def = Map::<Info>::new(Number::Unknown, Type::Float, "PEXT score");
    header
        .infos_mut()
        .insert(output_tag.as_ref().to_string(), info_def);

    writer.write_header(&header)?;

    let mut record = Record::default();
    let mut joined = String::new();
    while reader.read_record(&mut record)? > 0 {
        // Parse consequences from tag
        let values = get_csq_values(&record, &header, &csq_tag)?;

        if values.is_none() {
            writer.write_variant_record(&header, &record)?;
            continue;
        }

        let annotations_result: Result<Consequences, _> = match caller {
            CSQCaller::VEP => values
                .unwrap()
                .iter()
                .map(|s| Consequence::parse_from_vep(s.as_str()))
                .collect(),
            CSQCaller::BCSQ => values
                .unwrap()
                .iter()
                .map(|s| Consequence::parse_from_bcsq(s.as_str()))
                .collect(),
        };

        // Calculate and write PEXT scores
        let pext_scores = calculate_pext(annotations_result?, table)?;

        match pext_scores {
            None => {
                writer.write_variant_record(&header, &record)?;
            }
            Some(scores) => {
                format_scores(&mut joined, scores.as_slice());

                // Wrapper for lazy record to allow for adding an INFO tag,
                // performs slightly better than eager RecordBuf
                let annotated = Annotated {
                    inner: &record,
                    tag: output_tag.as_ref(),
                    joined: std::mem::take(&mut joined),
                };
                writer.write_variant_record(&header, &annotated)?;
                joined = annotated.joined;
            }
        }
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
