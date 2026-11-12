"""Independent Python tools for bulk epigenomics data and signal summaries.

Import functions or select a module under universal/specific. Established
optional tools provide peak calling, chromatin states and regulator activity.
See the reliability documentation for exploratory method limitations.
"""

__version__ = "0.5.2"

# --- Core data model -------------------------------------------------------
from easyepi.core.models import (  # noqa: F401
    CountMatrix,
    EpigenomeData,
    NormalizationFactors,
    PeakSet,
    read_peak_table,
)

# --- Statistics primitives -------------------------------------------------
from easyepi.core import stats as _stat
from easyepi.core.stats import (  # noqa: F401
    benjamini_hochberg,
    cohens_d,
    correlation,
    describe,
    fisher_exact,
    hypergeom_survival,
    log2fc,
    mann_whitney_u,
    t_test,
    wilcoxon,
)

# --- Analysis surface (semantic functions for import use) --------------------
from easyepi.normalize import (  # noqa: F401
    compute_normalization_factors,
    cpm_normalize,
    rpgc_normalize,
    spikein_normalize,
)

from easyepi import washu  # noqa: F401
from easyepi import footprint  # noqa: F401
from easyepi.footprint import (  # noqa: F401
    run_footprint,
    read_site_bed,
)
from easyepi import pubfig  # noqa: F401

from easyepi.differential import (  # noqa: F401
    differential_features,
    differential_test,
    read_count_matrix,
)

from easyepi.nucleosome import (  # noqa: F401
    bigwig_signal,
    nucleosome_scores,
    occupancy_profile,
    periodicity,
)

from easyepi.qc import bam_qc_summary  # noqa: F401

from easyepi.viz import (  # noqa: F401
    make_heatmap, make_profile, plot_matrix_heatmap, plot_profile,
)
from easyepi.plot import (  # noqa: F401
    plot_bar, plot_distribution, plot_heatmap, plot_volcano, save_figure,
    render_module_figures,
)

from easyepi.peaks import (  # noqa: F401
    bam_to_binned,
    bedgraph_to_signal,
    bigwig_to_signal,
    call_peaks_from_bam,
    peaks_to_bed,
    poisson_peaks,
)

from easyepi.annotate import annotate_peaks, read_intervals  # noqa: F401

from easyepi.integration import (  # noqa: F401
    IntegrationError, GeneRegion, cluster_tracks, gene_signal_matrix,
    link_peaks_to_genes, read_gene_bed, read_rna_results, run_integration,
)
from easyepi.project import (  # noqa: F401
    ProjectConfig, ProjectError, SampleRecord, apply_overrides, load_project,
    validate_project, write_inputs_tsv,
)
from easyepi.logging_utils import (  # noqa: F401
    RunRecorder, run_external, setup_logging,
)
from easyepi.assay_commands import run_assay_qc  # noqa: F401
from easyepi.activity import estimate_tf_activity  # noqa: F401
from easyepi.methylation import (  # noqa: F401
    MethylationError, differential_methylation, profile_methylation,
    read_methylation_table, summarize_methylation,
)

from easyepi.motifs import (  # noqa: F401
    discover_motifs,
    read_pwm,
    scan_sequences,
    formats,
    homer_like,
    homer_workflow,
    external_tools,
)

# Convenience re-export of some stats helpers
describe = _stat.describe  # noqa: F811

__all__ = [
    # version
    "__version__", "estimate_tf_activity",
    # data model
    "EpigenomeData", "PeakSet", "CountMatrix", "NormalizationFactors",
    "read_peak_table",
    # stats
    "describe", "correlation", "t_test", "wilcoxon", "mann_whitney_u",
    "fisher_exact", "hypergeom_survival", "cohens_d", "log2fc",
    "benjamini_hochberg",
    # normalize
    "compute_normalization_factors", "cpm_normalize", "rpgc_normalize",
    "spikein_normalize",
    # differential
    "differential_features", "differential_test", "read_count_matrix",
    # nucleosome
    "occupancy_profile", "periodicity", "nucleosome_scores", "bigwig_signal",
    # qc
    "bam_qc_summary",
    # plotting: table-based (plot_bar/plot_heatmap/...) and array-based helpers
    "plot_profile", "make_profile", "plot_matrix_heatmap", "make_heatmap",
    "plot_bar", "plot_distribution", "plot_volcano", "plot_heatmap",
    "save_figure", "render_module_figures",
    # peaks
    "poisson_peaks", "call_peaks_from_bam", "bam_to_binned", "bedgraph_to_signal",
    "bigwig_to_signal", "peaks_to_bed",
    # annotate
    "annotate_peaks", "read_intervals",
    # motifs
    "read_pwm", "scan_sequences", "discover_motifs",
    "formats", "homer_like", "homer_workflow", "external_tools",
    "washu", "footprint", "pubfig",
    "run_footprint", "read_site_bed",
    # unified project contract and run logging
    "ProjectConfig", "ProjectError", "SampleRecord", "load_project", "validate_project",
    "apply_overrides", "write_inputs_tsv", "setup_logging", "RunRecorder", "run_external",
    "run_assay_qc",
    "MethylationError", "read_methylation_table", "summarize_methylation",
    "profile_methylation", "differential_methylation",
    # Gene-centric multi-omics integration
    "IntegrationError", "GeneRegion", "read_gene_bed", "gene_signal_matrix", "cluster_tracks",
    "link_peaks_to_genes", "read_rna_results", "run_integration",
]
