# Peak calling: narrow vs broad, and choosing the right mode per mark

This note explains when to use `narrow_peaks` and when to use `broad_peaks`, and
what input formats each accepts. Getting this wrong does not raise an error — it
silently produces the wrong shape of result, which is why it is documented here
rather than left to the docstring.

## The two modes

| Mode | Command | MACS3 flag | Output | Fits |
|---|---|---|---|---|
| **narrow** | `easyepi.specific.chip.narrow_peaks` | — | `*_peaks.narrowPeak` | Sharp, focal enrichment |
| **broad** | `easyepi.specific.chip.broad_peaks` | `--broad` | `*_peaks.broadPeak` (+ `.gappedPeak`) | Broad, domain-scale enrichment |

A narrowPeak record carries a **summit** (column 10) and a point-like shape. A
broadPeak record has **no summit**: the peak is a region, and the ranking columns
describe the region as a whole.

## Which mode for which mark

| Mark | Mode | Why |
|---|---|---|
| **H3K4me3** | narrow | Promoter-proximal, sharp, few hundred bp |
| **H3K27ac** | narrow | Enhancer/promoter, sharp |
| **H3K4me1** | narrow | Enhancer, broader than H3K27ac but still focal |
| **H3K27me3** | **broad** | Polycomb domains, tens to hundreds of kb |
| **H3K9me3** | **broad** | Heterochromatin domains |
| **H3K36me3** | **broad** | Gene-body marks, very wide |

Calling a broad mark in narrow mode is the common mistake: MACS3 then reports
only the local maxima inside each domain, so one large repressive domain appears
as a scatter of short fragments and the true domain size is lost. In the ENCODE
H3K27me3 benchmark the difference is visible directly — narrow mode yields many
short fragments, broad mode yields the expected wide domains (median 360 bp over
the tested regions, against 252 bp for H3K4me3 narrow).

## Accepted input formats

Peak calling auto-detects the format from the file suffix:

| Input | `-f` sent to MACS3 | Notes |
|---|---|---|
| `.bam` | `BAMPE` (default) or `BAM` with `--single-end` | Uses observed fragment boundaries |
| `.bedGraph` / `.bdg` | `BED` | Adds `--nomodel --extsize 147` |
| `.bed` | `BED` | Same as bedGraph |
| `.bw` / `.bigWig` | — | **Rejected with an explicit error** |

### Why bigWig is rejected

MACS3 cannot read bigWig. Rather than let it emit a confusing binary-parse
failure, easyepi stops early and tells you the conversion:

```bash
# convert once, then call peaks on the bedGraph
bigWigToBedGraph signal.bw signal.bedGraph
python -m easyepi.specific.chip.broad_peaks \
    --project histone.toml --genome-size hs --out results/broad
```

A coverage track has no fragment-length model, so `--nomodel --extsize 147` is
added automatically. Without it MACS3 guesses a fragment size from a signal
track and peak widths become arbitrary.

## Example: the four-mark histone panel

```bash
# sharp marks
python -m easyepi.specific.chip.narrow_peaks \
    --project H3K4me3.toml --genome-size hs --out results/H3K4me3
python -m easyepi.specific.chip.narrow_peaks \
    --project H3K27ac.toml --genome-size hs --out results/H3K27ac
python -m easyepi.specific.chip.narrow_peaks \
    --project H3K4me1.toml --genome-size hs --out results/H3K4me1

# repressive domain
python -m easyepi.specific.chip.broad_peaks \
    --project H3K27me3.toml --genome-size hs --out results/H3K27me3
```

Each run writes `<out>/<sample>/<sample>_peaks.narrowPeak` (or `.broadPeak`) and a
`*_peaks_manifest.json` recording the exact MACS3 command, so the call is
reproducible without reading the logs.

## Downstream note

`chip-qc --peak` and the signal-enrichment metrics accept either extension, so
broad and narrow results can be QC'd through the same command. When both marks
are in one project, keep the extension in the manifest so the report shows which
peak type each sample used.
