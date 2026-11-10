"""Run Linux verification with the selected optional dependencies installed."""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description="Verify the Easyepi Linux runtime")
    parser.add_argument("--out",default=".verification/linux")
    args=parser.parse_args()
    if not sys.platform.startswith("linux"):
        parser.error("This verification target requires Linux")
    root=Path(__file__).resolve().parents[1]
    output=(root/args.out).resolve()
    output.mkdir(parents=True,exist_ok=True)
    environment=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",MPLBACKEND="Agg")
    results=[]
    commands=[[sys.executable,"-m","pytest","-q","-p","no:cacheprovider"]]
    modules=["easyepi.cli"]
    for folder in (root/"easyepi/universal",root/"easyepi/specific"):
        modules.extend(".".join(path.relative_to(root).with_suffix("").parts)
                       for path in sorted(folder.rglob("*.py")) if path.name!="__init__.py")
    commands.extend([sys.executable,"-m",module,"--help"] for module in modules)
    # The new analysis entries are subcommands of easyepi.cli, so they are not
    # discovered by the module scan above. Check each one explicitly; a missing
    # subcommand is a release blocker, not a cosmetic issue.
    for subcommand in ("tf-annotate","nucleosome-positioning","gene-panel",
                       "differential-peaks","histone-states","footprint","qc",
                       "integrate"):
        commands.append([sys.executable,"-m","easyepi.cli",subcommand,"--help"])
    for index,command in enumerate(commands):
        completed=subprocess.run(command,cwd=root,env=environment,text=True,capture_output=True)
        (output/f"check_{index:02d}.log").write_text(completed.stdout+completed.stderr,encoding="utf-8")
        results.append({"command":command,"returncode":completed.returncode})
        print(f"[Easyepi][VERIFY][{'PASS' if completed.returncode==0 else 'FAIL'}] {' '.join(command)}")
    versions={}
    for package in ("easyepi","numpy","scipy","pandas","pysam","pyBigWig","biopython","scikit-learn","decoupler"):
        try:
            versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package]="not installed"
    report={"platform":platform.platform(),"python":platform.python_version(),"versions":versions,"checks":results}
    (output/"verification.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    return int(any(result["returncode"] for result in results))


if __name__=="__main__":
    raise SystemExit(main())
