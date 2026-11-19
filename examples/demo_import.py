"""Easyepi import-usage demo.

Run in the ``easyepi`` conda environment:

    python examples/demo_import.py

This shows the Python-library API of Easyepi (as opposed to the CLI): every
resolution here is a plain function call returning numpy arrays / dicts.
"""

import numpy as np

import easyepi as epi


def section(title):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)


def main():
    section("1. Version & library API surface")
    print("easyepi", epi.__version__)
    public = [n for n in epi.__all__ if not n.startswith("__")]
    print("public API (%d):" % len(public))
    print(", ".join(public))

    section("2. Core data model (PeakSet / EpigenomeData)")
    peaks = epi.PeakSet(
        chrom=np.array(["chr1", "chr1", "chr2"]),
        start=np.array([100, 900, 500]),
        end=np.array([300, 1100, 700]),
        signal=np.array([5.0, 3.0, 8.0]),
    )
    print("n peaks:", len(peaks), "| widths:", peaks.width.tolist())
    ed = epi.EpigenomeData(name="demo", peaks=peaks)
    ed.add_sample("S1", {"group": "ctl"})
    ed.add_sample("S2", {"group": "ctl"})
    ed.add_sample("S3", {"group": "trt"})
    print("sample groups:", ed.sample_groups())

    section("3. Descriptive statistics & tests")
    a = np.array([1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    b = np.array([3.0, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    print("describe(a):", epi.describe(a))
    print("log2fc(a,b):", round(epi.log2fc(a, b), 4))
    print("cohens_d(a,b):", round(epi.cohens_d(a, b), 4))
    print("pearson:", round(epi.correlation(a, b), 4))
    print("t-test p:", epi.t_test(a, b)["pvalue"])
    print("BH:", epi.benjamini_hochberg([0.01, 0.04, 0.031, 0.6, 0.2]).round(4))

    section("4. Differential analysis (numpy matrix -> structured result)")
    matrix = np.array([
        [12, 14, 13, 55, 60, 58],
        [5, 4, 6, 3, 2, 4],
        [20, 21, 18, 2, 1, 3],
        [8, 9, 7, 8, 9, 10],
    ], dtype=float)
    mask_a = np.array([True, True, True, False, False, False])
    res = epi.differential_features(matrix, mask_a, method="wilcoxon")
    print("== per-region differential ==")
    print("log2fc  :", res["log2fc"].round(3).tolist())
    print("cohens_d:", res["cohens_d"].round(3).tolist())
    print("pvalue  :", res["pvalue"].round(4).tolist())
    print("padj    :", res["padj"].round(4).tolist())

    section("5. Normalization")
    m = np.array([[1000, 2000], [500, 900]], dtype=float)
    print("cpm:", epi.cpm_normalize(m, np.array([1e6, 2e6])).round(2).tolist())
    factors = epi.compute_normalization_factors(
        ["S1", "S2"], lib_sizes={"S1": 1e6, "S2": 2e6}, spikeins={"S1": 1.0, "S2": 0.5}
    )
    print("factors:", factors.as_dict())

    section("6. Nucleosome occupancy & periodicity")
    sig = np.random.RandomState(0).randn(200).cumsum()
    print("occupancy z-score head:", epi.occupancy_profile(sig)[:4].round(3).tolist())
    per = epi.periodicity(sig)
    print("periodicity:", per)

    print("\nAll import-API demos completed OK.")


if __name__ == "__main__":
    main()
