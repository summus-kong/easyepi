"""Write a project file from a folder of BAM/BigWig/peak files.

The project TOML is the only state this module produces; it records the samples,
the assay guessed for each, and the reference paths, so a later command can take
``--project`` instead of a long argument list.

There is deliberately no "run everything" entry point: the suite is step-by-step,
and an orchestrator hid which step was running and made a partial failure hard to
read. ``easyepi init`` writes the file; the analysis commands consume it.
"""

from __future__ import annotations

from pathlib import Path

from easyepi.quickstart import build_project, project_to_toml


def run_init(args) -> int:
    """CLI driver for ``easyepi init``: write a project file from a folder."""
    project, notes = build_project(
        args.source,
        project_id=getattr(args, "project_id", None),
        genome_build=getattr(args, "genome_build", None),
        annotation=getattr(args, "annotation", None),
        motif_database=getattr(args, "motif_database", None))
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(project_to_toml(project), encoding="utf-8")

    print(f"[easyepi] init: {len(project['samples'])} sample(s) -> {target}")
    for row in project["samples"]:
        kinds = ", ".join(k for k in ("bam", "bigwig", "peak", "methylation")
                          if k in row)
        print(f"    {row['sample_id']:20s} {row.get('assay', '?'):10s} {kinds}")
    for note in notes:
        print(f"  note: {note}")
    print("  next: review the assay labels and paths, then `easyepi quick "
          f"{target} --out results`")
    return 0
