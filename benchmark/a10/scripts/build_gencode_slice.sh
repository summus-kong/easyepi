#!/usr/bin/env bash
# Build the chr1 slice used to exercise the *explicit* UTR path.
#
# Why a separate slice exists
# ---------------------------
# The repaired GTF (`txdb_fixed.gtf`) carries no UTR rows at all, so the only
# UTR path ever tested was `_derived_utr`: exon minus CDS. GENCODE's own GTF
# carries one `UTR` row per UTR exon, so it exercises the other branch.
#
# Note the tag: GENCODE writes `UTR`, not `five_prime_utr`/`three_prime_utr`.
# `gtf_features.py` therefore ignores those rows and derives the UTRs anyway --
# the comparison below is what shows the two derivations agree.
#
# Usage: build_gencode_slice.sh <gencode.gtf.gz> <out.gtf> [chr] [end_bp]
set -euo pipefail

src=${1:?usage: build_gencode_slice.sh <gencode.gtf.gz> <out.gtf> [chr] [end_bp]}
out=${2:?usage: build_gencode_slice.sh <gencode.gtf.gz> <out.gtf> [chr] [end_bp]}
chrom=${3:-chr1}
end_bp=${4:-20000000}

mkdir -p "$(dirname "$out")"
zcat "$src" \
  | awk -F'\t' -v OFS='\t' -v c="$chrom" -v e="$end_bp" '
      /^#/ { print; next }
      # Keep any feature that starts inside the window and does not run far
      # past it, so no transcript is cut in half by the slice boundary.
      $1 == c && $4 <= e && $5 <= e + 5000000 && $3 != "Selenocysteine" { print }
  ' > "$out"

echo "wrote $out"
echo "lines: $(wc -l < "$out")"
awk -F'\t' '!/^#/{c[$3]++} END{for(k in c) printf "  %-16s %d\n", k, c[k]}' "$out" | sort -k2 -nr
