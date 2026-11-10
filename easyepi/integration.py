"""Gene-centric integration of normalized BigWig tracks, peaks and RNA result tables. Nearest-TSS links are candidates, not causal regulation."""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


class IntegrationError(ValueError):
    """Invalid integration input contract."""


@dataclass(frozen=True)
class GeneRegion:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str = "+"

    @property
    def tss(self) -> int:
        """Transcription start site coordinate implied by this gene record."""
        return self.start if self.strand != "-" else self.end


def read_gene_bed(path: str | Path) -> list[GeneRegion]:
    """Read BED4+ gene intervals with unique IDs and valid coordinates."""
    records: list[GeneRegion] = []
    seen: set[str] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise IntegrationError(f"gene BED line {line_no} needs at least four columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise IntegrationError(f"gene BED line {line_no} has invalid coordinates") from exc
            gene_id = fields[3].strip()
            if not gene_id or end <= start:
                raise IntegrationError(f"gene BED line {line_no} has an invalid ID or interval")
            if gene_id in seen:
                raise IntegrationError(f"duplicate gene ID in BED: {gene_id}")
            seen.add(gene_id)
            records.append(GeneRegion(gene_id, fields[0], start, end, fields[5] if len(fields) >= 6 else "+"))
    if not records:
        raise IntegrationError("gene BED contains no valid records")
    return records


def _bin_mean(bigwig, chrom: str, start: int, end: int, bins: int) -> np.ndarray:
    if bins < 1 or end <= start:
        return np.full(bins, np.nan, dtype=float)
    edges = np.linspace(start, end, bins + 1).astype(int)
    values = []
    for left, right in zip(edges[:-1], edges[1:]):
        if right <= left:
            values.append(np.nan)
            continue
        try:
            value = bigwig.stats(chrom, int(max(0, left)), int(right), type="mean")[0]
        except RuntimeError:
            value = None
        # For an existing chromosome, uncovered BigWig intervals represent zero signal.
        # Missing chromosomes remain NaN and must not be treated as uncovered intervals.
        values.append(0.0 if value is None else float(value))
    return np.asarray(values, dtype=float)


def gene_signal_matrix(
    bigwig_path: str | Path,
    genes: Sequence[GeneRegion],
    upstream: int = 2000,
    downstream: int = 2000,
    bins_upstream: int = 40,
    bins_body: int = 100,
    bins_downstream: int = 40,
) -> pd.DataFrame:
    """Build strand-oriented upstream, scaled gene-body and downstream signal matrices."""
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - dependency contract
        raise ImportError("gene integration requires pyBigWig") from exc
    bigwig = pyBigWig.open(str(bigwig_path))
    if bigwig is None:
        raise IntegrationError(f"cannot open BigWig: {bigwig_path}")
    try:
        chrom_sizes = bigwig.chroms()
        rows = []
        for gene in genes:
            if gene.chrom not in chrom_sizes:
                rows.append(np.full(bins_upstream + bins_body + bins_downstream, np.nan))
                continue
            if gene.strand == "-":
                pieces = (
                    _bin_mean(bigwig, gene.chrom, gene.end, gene.end + upstream, bins_upstream)[::-1],
                    _bin_mean(bigwig, gene.chrom, gene.start, gene.end, bins_body)[::-1],
                    _bin_mean(bigwig, gene.chrom, gene.start - downstream, gene.start, bins_downstream)[::-1],
                )
            else:
                pieces = (
                    _bin_mean(bigwig, gene.chrom, gene.start - upstream, gene.start, bins_upstream),
                    _bin_mean(bigwig, gene.chrom, gene.start, gene.end, bins_body),
                    _bin_mean(bigwig, gene.chrom, gene.end, gene.end + downstream, bins_downstream),
                )
            rows.append(np.concatenate(pieces))
    finally:
        bigwig.close()
    columns = [f"up_{i + 1}" for i in range(bins_upstream)]
    columns += [f"body_{i + 1}" for i in range(bins_body)]
    columns += [f"down_{i + 1}" for i in range(bins_downstream)]
    return pd.DataFrame(rows, index=[gene.gene_id for gene in genes], columns=columns)


def _row_zscore(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.to_numpy(dtype=float)
    means = np.nanmean(values, axis=1, keepdims=True)
    scales = np.nanstd(values, axis=1, keepdims=True)
    scales[scales == 0] = 1.0
    return pd.DataFrame((values - means) / scales, index=frame.index, columns=frame.columns)


def deterministic_kmeans(values: np.ndarray, clusters: int, seed: int = 1, n_init: int = 20, max_iter: int = 300) -> tuple[np.ndarray, float]:
    """Run scikit-learn KMeans with an explicit random seed; return one-based labels and inertia."""
    if values.ndim != 2 or len(values) < clusters or clusters < 2:
        raise IntegrationError("clusters must be >=2 and no larger than the number of complete genes")
    if not np.isfinite(values).all() or len(np.unique(values,axis=0)) < clusters:
        raise IntegrationError("Requested clusters exceed the number of distinct complete profiles")
    from sklearn.cluster import KMeans
    model = KMeans(n_clusters=clusters, random_state=seed, n_init=n_init, max_iter=max_iter)
    labels = model.fit_predict(values)
    return labels + 1, float(model.inertia_)


def cluster_tracks(matrices: Mapping[str, pd.DataFrame], clusters: int, seed: int = 1) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Concatenate standardized gene profiles while retaining only complete measurements."""
    if not matrices:
        raise IntegrationError("at least one signal track is required")
    common = None
    for frame in matrices.values():
        common = frame.index if common is None else common.intersection(frame.index)
    assert common is not None
    common = common.sort_values()
    complete = common
    for frame in matrices.values():
        complete = complete[frame.loc[complete].notna().all(axis=1)]
    if len(complete) == 0:
        raise IntegrationError("no genes have complete signal across all requested tracks")
    features = []
    for name, frame in matrices.items():
        z = _row_zscore(frame.loc[complete])
        z.columns = [f"{name}|{column}" for column in z.columns]
        features.append(z)
    feature_frame = pd.concat(features, axis=1).replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    labels, inertia = deterministic_kmeans(feature_frame.to_numpy(), clusters, seed=seed)
    labels_frame = pd.DataFrame({"gene_id": feature_frame.index, "cluster": labels}).sort_values(["cluster", "gene_id"])
    return labels_frame, feature_frame, inertia


def link_peaks_to_genes(peaks: str | Path, genes: Sequence[GeneRegion], max_distance: int = 100_000) -> pd.DataFrame:
    """Assign candidate nearest-TSS links within an explicit maximum distance."""
    by_chrom: dict[str, list[GeneRegion]] = {}
    for gene in genes:
        by_chrom.setdefault(gene.chrom, []).append(gene)
    rows = []
    with Path(peaks).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise IntegrationError(f"peak BED line {line_no} has fewer than three columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise IntegrationError(f"peak BED line {line_no} has invalid coordinates") from exc
            center = (start + end) // 2
            candidates = by_chrom.get(fields[0], [])
            if not candidates:
                continue
            gene = min(candidates, key=lambda item: abs(center - item.tss))
            distance = center - gene.tss
            if abs(distance) <= max_distance:
                rows.append({"peak_id": f"peak_{line_no}", "chrom": fields[0], "start": start, "end": end,
                             "gene_id": gene.gene_id, "signed_distance_to_tss": distance,
                             "link_method": "nearest_tss", "max_distance": max_distance})
    return pd.DataFrame(rows, columns=["peak_id", "chrom", "start", "end", "gene_id", "signed_distance_to_tss", "link_method", "max_distance"])


def read_rna_results(path: str | Path, gene_column: str = "gene_id") -> pd.DataFrame:
    """Read RNA results and remove only Ensembl-style version suffixes."""
    frame = pd.read_csv(path, sep="\t")
    if gene_column not in frame:
        raise IntegrationError(f"RNA result lacks requested gene column: {gene_column}")
    frame = frame.copy()
    frame["gene_id"] = frame[gene_column].astype(str).str.replace(r"\.\d+$", "", regex=True)
    if frame["gene_id"].duplicated().any():
        raise IntegrationError("RNA result contains duplicate gene IDs after version normalization")
    return frame


def run_integration(
    gene_bed: str | Path,
    tracks: Mapping[str, str | Path],
    out_dir: str | Path,
    clusters: int = 4,
    seed: int = 1,
    peaks: str | Path | None = None,
    rna_results: str | Path | None = None,
    rna_gene_column: str = "gene_id",
    max_link_distance: int = 100_000,
) -> dict[str, object]:
    """Write auditable gene-centric integration tables and a manifest."""
    if not tracks:
        raise IntegrationError("at least one named BigWig track is required")
    genes = read_gene_bed(gene_bed)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    matrices = {name: gene_signal_matrix(path, genes) for name, path in tracks.items()}
    for name, matrix in matrices.items():
        matrix.to_csv(out / f"gene_matrix.{name}.tsv.gz", sep="\t", compression="gzip")
    labels, features, inertia = cluster_tracks(matrices, clusters=clusters, seed=seed)
    labels.to_csv(out / "gene_clusters.tsv", sep="\t", index=False)
    features.to_csv(out / "cluster_features.tsv.gz", sep="\t", compression="gzip")
    result: dict[str, object] = {"genes_input": len(genes), "genes_clustered": len(labels), "tracks": list(tracks),
                                 "clusters": clusters, "seed": seed, "inertia": inertia,
                                 "outputs": {"clusters": str(out / "gene_clusters.tsv"), "features": str(out / "cluster_features.tsv.gz")}}
    if rna_results:
        rna = read_rna_results(rna_results, rna_gene_column)
        joint = labels.merge(rna, on="gene_id", how="left", validate="one_to_one")
        joint.to_csv(out / "gene_clusters_rna.tsv", sep="\t", index=False)
        result["outputs"]["rna_integration"] = str(out / "gene_clusters_rna.tsv")
        result["rna_genes_matched"] = int(joint["gene_id"].isin(rna["gene_id"]).sum())
    if peaks:
        links = link_peaks_to_genes(peaks, genes, max_distance=max_link_distance)
        links.to_csv(out / "peak_gene_links.tsv", sep="\t", index=False)
        result["outputs"]["peak_gene_links"] = str(out / "peak_gene_links.tsv")
        result["peaks_linked"] = len(links)
    (out / "integration_manifest.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    result["outputs"]["manifest"] = str(out / "integration_manifest.json")
    return result
