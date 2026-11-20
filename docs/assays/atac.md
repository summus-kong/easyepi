# ATAC-seq

ATAC-seq measures accessibility through Tn5 insertion. The signal of interest is the
insertion site, not the read body, so `bam2bw --assay atac` applies the standard
+4/−5 shift before accumulating.

## Commands

```bash
easyepi bam2bw atac.bam -o tracks/ --assay atac --jobs 8
easyepi peaks atac.bam -o peaks/ --assay atac
easyepi qc --assay atac --bam atac.bam --bigwig tracks/atac.bw \
    --peak peaks/atac/atac_peaks.narrowPeak --out qc/atac.json
easyepi footprint --bigwig tracks/atac.bw --target peaks/atac/atac_peaks.narrowPeak \
    --chrom chr1 --start 0 --end 5000000 --flank 100 --out fp/
```

## Metrics and what they mean

| Metric | Good range | If out of range |
| --- | --- | --- |
| FRiP (genome-wide) | 0.2–0.5 | low: over-digestion, or peaks called over a sub-region |
| Fragment periodicity | visible 10 bp ladder | lost: nuclei prep or digestion problem |
| Nucleosome ratio | mononucleosome dominant | dinucleosome-heavy: under-digestion |
| TSS enrichment | > 5 | low: dead cells or failed transposition |

`qc` reports FRiP as `interpretable: false` when the peak set spans only a small
fraction of the mapped genome — a regional numerator over a global denominator is not
a FRiP, and the message says so.

### Which unit the track counts

Every BigWig is written next to a `.bw.json` that now carries a `counting` block, so
a track can be compared with another one knowing what its values are:

```json
"counting": {
  "unit": "fragment",
  "fragments": 3371034,
  "unpaired_reads": 56539,
  "per_chromosome": {"chr1": {"fragment": 3371034, "read": 56539}}
}
```

`fragment` means one row per template: mate two is skipped and the span runs to the
mate's outer edge. `read` means one row per alignment, which is what an unpaired
record gets. **The two are not interchangeable** — counting both mates would double
every paired track, so a downstream comparison across a fragment-counted and a
read-counted track is not valid.

`unit` reports what was actually counted, not what the BAM header claims. A file
declared paired-end whose records carry no mate still reports `read`.

Counts are per chromosome and are never summed across chromosomes: chunks are
fetched with a pad and clipped, which keeps each chromosome right but would make a
genome total double count a read that crosses a boundary.

The cut-site path (`--assay atac`, the default Tn5 shift) accumulates read ends
rather than fragment bodies, so it writes **no** `counting` block. That absence is
intentional: there is no single fragment unit to name, and reporting one would be
wrong.

## Assumptions

- the BAM is coordinate-sorted and indexed
- reads are paired-end for fragment-length metrics; single-end still yields FRiP and
  signal enrichment
- `--assay atac` implies Tn5 shift, 10 bp bins and CPM; override with `--shift none`
  if the BAM was already shifted
