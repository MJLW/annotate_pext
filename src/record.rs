//! A `variant::Record` wrapper that appends a single string-array INFO field
//! to a lazy `vcf::Record` on write, delegating every other field to the inner
//! record (so untouched fields are serialized from their stored bytes rather
//! than re-parsed and re-formatted).

use std::{borrow::Cow, fmt::Write as _, io};

use noodles_core::Position;
use noodles_vcf::{
    self as vcf, Header,
    variant::record::{
        AlternateBases, Filters, Ids, Info, ReferenceBases, Samples,
        info::field::{
            Value,
            value::{Array, array::Values},
        },
    },
};

pub struct Annotated<'r> {
    pub inner: &'r vcf::Record,
    pub tag: &'r str,
    pub joined: String,
}

impl vcf::variant::Record for Annotated<'_> {
    fn reference_sequence_name<'a, 'h: 'a>(&'a self, _header: &'h Header) -> io::Result<&'a str> {
        Ok(self.inner.reference_sequence_name())
    }

    fn variant_start(&self) -> Option<io::Result<Position>> {
        self.inner.variant_start()
    }

    fn ids(&self) -> Box<dyn Ids + '_> {
        Box::new(self.inner.ids())
    }

    fn reference_bases(&self) -> Box<dyn ReferenceBases + '_> {
        Box::new(self.inner.reference_bases())
    }

    fn alternate_bases(&self) -> Box<dyn AlternateBases + '_> {
        Box::new(self.inner.alternate_bases())
    }

    fn quality_score(&self) -> Option<io::Result<f32>> {
        self.inner.quality_score()
    }

    fn filters(&self) -> Box<dyn Filters + '_> {
        Box::new(self.inner.filters())
    }

    fn info(&self) -> Box<dyn Info + '_> {
        Box::new(AppendedInfo {
            inner: Box::new(self.inner.info()),
            tag: self.tag,
            joined: &self.joined,
        })
    }

    fn samples(&self) -> io::Result<Box<dyn Samples + '_>> {
        Ok(Box::new(self.inner.samples()))
    }
}

struct AppendedInfo<'a> {
    inner: Box<dyn Info + 'a>,
    tag: &'a str,
    joined: &'a str,
}

fn appended_value(joined: &str) -> Value<'_> {
    let values: Box<dyn Values<'_, Cow<'_, str>> + '_> = Box::new(joined);
    Value::Array(Array::String(values))
}

impl Info for AppendedInfo<'_> {
    fn is_empty(&self) -> bool {
        false // we always contribute at least our one field
    }

    fn len(&self) -> usize {
        self.inner.len() + 1
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        if key == self.tag {
            Some(Ok(Some(appended_value(self.joined))))
        } else {
            self.inner.get(header, key)
        }
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        let extra = std::iter::once(Ok((self.tag, Some(appended_value(self.joined)))));
        Box::new(self.inner.iter(header).chain(extra))
    }
}

pub fn format_scores(buf: &mut String, scores: &[Option<f32>]) {
    buf.clear();
    for (i, v) in scores.iter().enumerate() {
        if i > 0 {
            buf.push(',');
        }
        match v {
            Some(x) => {
                let _ = write!(buf, "{x:.2}");
            }
            None => buf.push('.'),
        }
    }
}

/* ---- usage sketch (replaces your read_record_buf loop) -------------------

let mut record = vcf::Record::default();   // lazy, raw-line-backed
let mut joined = String::new();            // reused scratch buffer

while reader.read_record(&mut record)? > 0 {
    // Adapt CSQ extraction to the lazy INFO API, e.g.:
    //   record.info().get(&header, &csq_tag)  -> Option<Result<Option<Value>>>
    let pext_scores: Option<Vec<Option<f64>>> = /* your pext pipeline */;

    match pext_scores {
        None => {
            // unchanged -> write the lazy record straight through, no rebuild
            writer.write_variant_record(&header, &record)?;
        }
        Some(scores) => {
            format_scores(&mut joined, &scores);
            let annotated = Annotated {
                inner: &record,
                tag: output_tag.as_ref(),
                joined: std::mem::take(&mut joined), // hand the buffer over
            };
            writer.write_variant_record(&header, &annotated)?;
            joined = annotated.joined;               // take it back to reuse
        }
    }
}
--------------------------------------------------------------------------- */
