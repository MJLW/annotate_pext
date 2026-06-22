# PEXT Annotater

Annotates PEXT scores for protein-coding transcripts on a TSV (2026-06-22: this differs from the current implementation of PEXT as that uses ALL transcripts). There are two steps involved in annotating the scores:
 - Condensing the GTEX transcript TPMs
 - Annotating PEXT scores on TSV using the condensed table

Example usage:
```
./target/release/condense --gtex-tpms <gtex_transcript_tpms> \
    --gtex-sample-attributes <sample_attributes> \
    --tissue-blacklist <tissue_blacklist> \
    --coding-transcripts <coding_transcripts> \
    --output <condensed_table>

./target/release/annotate_tsv --variants <variants_tsv> \
    --variant-columns CHROM,POS,REF,ALT \
    --gene-id-column vepGene \
    --transcript-id-column vepFeature \
    --biotype-column vepBIOTYPE \
    --group-columns vepConsequence,vepLoF \
    --tpms <condensed_table> \
    --output <annotated_variants_tsv>
```

Where:
 - `gtex_transcript_tpms` is the GTEx bulk expression tpms file (see https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression).
 - `sample_attributes` is the GTEx sample attributes file (.txt) associated with the bulk expression tpms (see https://www.gtexportal.org/home/downloads/adult-gtex/metadata).
 - `tissue_blacklist` is a file that should contain all tissues that should be excluded (e.g, cell lines, sexual organs), separated by newlines.
 - `coding_transcripts` is a file that should contain all coding transcripts without version (e.g., `ENST00000335137`, not `ENST00000335137.4`, separated by newlines. Any transcripts not in this file will not be in the condensed table.
 - `variants_tsv` is a TSV file containing lines of variant-transcript level annotations (one transcript per line).
 - `*-column(s)` describe the format of the TSV file. Should hopefully be self explanatory.

