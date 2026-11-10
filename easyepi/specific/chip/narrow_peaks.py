"""CHIP narrow peak calling with MACS3."""
from easyepi._peak_calling import main as _main


def main(argv=None):
    """Command-line entry point for this module."""
    return _main(argv, assay='CHIP', broad=False)


if __name__ == "__main__":
    raise SystemExit(main())
