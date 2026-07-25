//! Annotates a LOFTEE-like HC/LC call so that variant sets annotated with
//! bcftools csq can be grouped the same way as variant sets annotated with
//! VEP + LOFTEE.
//!
//! Two modes:
//!
//!   --mode vcf   reads INFO/BCSQ, writes INFO/<output-tag> with one value per
//!                BCSQ entry, in the same order, so downstream splitting stays
//!                aligned.
//!
//!   --mode tsv   reads a flattened variant-transcript table and appends a
//!                single column. Use this to put Roulette and the proband
//!                tables through an identical path, and to check these calls
//!                against a real vepLoF column.
//!
//! See `lof.rs` for exactly how much of LOFTEE this does and does not do.

use std::{
    error::Error,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use annotate_pext::{
    lof::{Classifier, Rule, load_gff},
    record::Annotated,
    utils::build_tsv_reader,
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
        io::Write as _,
        record::info::field::{Value, value::Array},
    },
};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Mode {
    Vcf,
    Tsv,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum RuleArg {
    /// LC if the variant is in the final coding exon.
    LastExon,
    /// LC if in the final coding exon or within --nmd-window bases of the last junction.
    NmdEscape,
    /// LC if in the terminal --end-trunc-fraction of the CDS. Closest to LOFTEE.
    EndTrunc,
}

impl From<RuleArg> for Rule {
    fn from(value: RuleArg) -> Self {
        match value {
            RuleArg::LastExon => Rule::LastExon,
            RuleArg::NmdEscape => Rule::NmdEscape,
            RuleArg::EndTrunc => Rule::EndTrunc,
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(
    author = "MJLW",
    version = "0.1.0",
    about = "LOFTEE-like HC/LC calls from a GFF3, for BCSQ VCFs and flattened TSVs."
)]
struct Args {
    #[arg(long)]
    variants: PathBuf,

    /// The same GFF3 that was handed to `bcftools csq`.
    #[arg(long)]
    gff: PathBuf,

    #[arg(long, value_enum, default_value = "vcf")]
    mode: Mode,

    #[arg(long, value_enum, default_value = "last-exon")]
    rule: RuleArg,

    #[arg(long, default_value_t = 0.05)]
    end_trunc_fraction: f64,

    #[arg(long, default_value_t = 50)]
    nmd_window: u64,

    #[arg(long, default_value = "LOF")]
    output_tag: String,

    #[arg(long)]
    output: PathBuf,

    // ---- VCF mode ----
    #[arg(long, default_value = "BCSQ")]
    csq_tag: String,

    // ---- TSV mode ----
    #[arg(long, default_value = "POS")]
    pos_column: String,

    #[arg(long, default_value = "TRANSCRIPT")]
    transcript_id_column: String,

    #[arg(long, default_value = "CONSEQUENCE")]
    consequence_column: String,
}

/// bcftools csq: consequence|gene|transcript|biotype|strand|aa|nt
const BCSQ_CONSEQUENCE: usize = 0;
const BCSQ_TRANSCRIPT: usize = 2;
const BCSQ_MIN_FIELDS: usize = 4;

fn get_csq_values(
    record: &Record,
    header: &Header,
    csq_tag: &str,
) -> Result<Option<Vec<String>>, Box<dyn Error>> {
    let info = record.info();
    let option_value = info.get(header, csq_tag);

    if option_value.is_none() {
        return Ok(None);
    }

    let value = match option_value.unwrap()? {
        Some(value) => value,
        None => return Ok(None),
    };

    let values = match value {
        Value::Array(Array::String(values)) => values
            .iter()
            .map(|s| s.unwrap().unwrap().to_string())
            .collect::<Vec<String>>(),
        Value::String(s) => vec![s.to_string()],
        other => return Err(format!("CSQ tag '{csq_tag}' had unexpected type: {other:?}").into()),
    };

    Ok(Some(values))
}

fn join_calls(buf: &mut String, calls: &[&str]) {
    buf.clear();
    for (i, call) in calls.iter().enumerate() {
        if i > 0 {
            buf.push(',');
        }
        buf.push_str(call);
    }
}

fn build_vcf_writer<P: AsRef<Path>>(
    path: P,
) -> Result<noodles_vcf::io::Writer<Box<dyn Write>>, Box<dyn Error>> {
    let file = File::create(&path)?;
    let buf = BufWriter::with_capacity(1 << 20, file);
    let inner: Box<dyn Write> = match path.as_ref().extension().and_then(|e| e.to_str()) {
        Some("gz") | Some("bgz") => Box::new(noodles_bgzf::io::Writer::new(buf)),
        _ => Box::new(buf),
    };
    Ok(noodles_vcf::io::Writer::new(inner))
}

fn annotate_vcf(args: &Args, classifier: &mut Classifier) -> Result<(), Box<dyn Error>> {
    let mut reader = noodles_vcf::io::reader::Builder::default().build_from_path(&args.variants)?;
    let mut writer = build_vcf_writer(&args.output)?;

    let mut header = reader.read_header()?;
    let info_def = Map::<Info>::new(
        Number::Unknown,
        Type::String,
        "LOFTEE-like HC/LC call, one value per BCSQ entry",
    );
    header.infos_mut().insert(args.output_tag.clone(), info_def);
    writer.write_header(&header)?;

    let mut record = Record::default();
    let mut joined = String::new();

    while reader.read_record(&mut record)? > 0 {
        let values = get_csq_values(&record, &header, &args.csq_tag)?;

        let values = match values {
            None => {
                writer.write_variant_record(&header, &record)?;
                continue;
            }
            Some(values) => values,
        };

        let position: u64 = match record.variant_start() {
            Some(start) => usize::from(start?) as u64,
            None => {
                writer.write_variant_record(&header, &record)?;
                continue;
            }
        };

        let calls: Vec<&str> = values
            .iter()
            .map(|csq| {
                // bcftools emits "@<pos>" back-references for entries that
                // belong to another record. They carry no transcript.
                if csq.starts_with('@') {
                    return ".";
                }
                let fields: Vec<&str> = csq.split('|').collect();
                if fields.len() < BCSQ_MIN_FIELDS {
                    return ".";
                }
                classifier
                    .call(
                        fields[BCSQ_TRANSCRIPT],
                        position,
                        fields[BCSQ_CONSEQUENCE],
                    )
                    .as_str()
            })
            .collect();

        join_calls(&mut joined, &calls);

        let annotated = Annotated {
            inner: &record,
            tag: args.output_tag.as_str(),
            joined: std::mem::take(&mut joined),
        };
        writer.write_variant_record(&header, &annotated)?;
        joined = annotated.joined;
    }

    Ok(())
}

fn find_header_index(header: &csv::StringRecord, column: &str) -> Result<usize, Box<dyn Error>> {
    header
        .iter()
        .position(|name| name == column)
        .ok_or_else(|| format!("Could not find column '{column}' in TSV.").into())
}

fn annotate_tsv(args: &Args, classifier: &mut Classifier) -> Result<(), Box<dyn Error>> {
    let mut reader = build_tsv_reader(&args.variants)?;
    let header = reader.headers()?.clone();

    let pos_index = find_header_index(&header, &args.pos_column)?;
    let transcript_index = find_header_index(&header, &args.transcript_id_column)?;
    let consequence_index = find_header_index(&header, &args.consequence_column)?;

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(File::create(&args.output)?);

    for field in &header {
        writer.write_field(field)?;
    }
    writer.write_field(&args.output_tag)?;
    writer.write_record(None::<&[u8]>)?;

    let mut record = csv::StringRecord::new();
    while reader.read_record(&mut record)? {
        let position: u64 = record
            .get(pos_index)
            .ok_or("Malformed TSV: missing position column")?
            .parse()?;
        let transcript = record
            .get(transcript_index)
            .ok_or("Malformed TSV: missing transcript column")?;
        let consequence = record
            .get(consequence_index)
            .ok_or("Malformed TSV: missing consequence column")?;

        let call = classifier.call(transcript, position, consequence);

        for field in &record {
            writer.write_field(field)?;
        }
        writer.write_field(call.as_str())?;
        writer.write_record(None::<&[u8]>)?;
    }

    writer.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let transcripts = load_gff(&args.gff)?;
    let mut classifier = Classifier::new(
        transcripts,
        args.rule.into(),
        args.end_trunc_fraction,
        args.nmd_window,
    );

    eprintln!(
        "[annotate_lof] loaded {} transcripts from {}",
        classifier.transcript_count(),
        args.gff.display()
    );

    match args.mode {
        Mode::Vcf => annotate_vcf(&args, &mut classifier)?,
        Mode::Tsv => annotate_tsv(&args, &mut classifier)?,
    }

    let counts = classifier.counts;
    eprintln!(
        "[annotate_lof] HC={} LC={} NA={} (of which transcript_absent={}, transcript_noncoding={})",
        counts.hc,
        counts.lc,
        counts.not_applicable,
        counts.transcript_absent,
        counts.transcript_noncoding
    );

    // A pLoF whose transcript is missing from the GFF is silently downgraded to
    // ".", which would then group with every non-pLoF consequence. Make that
    // loud rather than letting it quietly change pext groups.
    if counts.transcript_absent > 0 {
        eprintln!(
            "[annotate_lof] WARNING: {} pLoF consequences referenced a transcript that is not in \
             the GFF. Check that the GFF matches the annotation source.",
            counts.transcript_absent
        );
    }

    Ok(())
}
