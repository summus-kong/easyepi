# CUT&Tag and CUT&RUN

Both measure a Tn5 (CUT&Tag) or MNase (CUT&RUN) cut site, so the analysis follows the
ATAC path rather than the ChIP path: apply the shift, call narrow peaks, and expect a
high signal-to-noise ratio at low sequencing depth.

```bash
easyepi bam2bw cuttag.bam -o tracks/ --assay cuttag --jobs 8
easyepi peaks tracks/cuttag.bw -o peaks/ --assay cuttag
easyepi qc --assay chip --bam cuttag.bam --bigwig tracks/cuttag.bw \
    --peak peaks/cuttag.narrowPeak --out qc/cuttag.json
```

`--assay cuttag` and `--assay cut&run` are aliases of the Tn5 defaults (shift=tn5,
10 bp bins, CPM). Use the unquoted spelling in a shell: `&` is an operator, so
`--assay cut&tag` must be quoted or written as `cuttag`.

## QC notes

CUT&Tag libraries often have a much higher FRiP than ChIP-seq at far lower depth. A
FRiP of 0.3 from 5 M reads is normal here and would be unremarkable for ChIP at 30 M.
Compare within the assay, not across.

## Assumptions

- the BAM is paired-end or single-end; both give cut-site positions
- `easyepi qc --assay chip` is used for the metric set because CUT&Tag shares the
  FRiP and enrichment metrics; there is no separate `cuttag` QC table
