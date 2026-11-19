"""Standalone ChromHMM train entry."""
import sys
from easyepi.specific.histone.chromhmm import main as _main


def main(argv=None):
    """Command-line entry point for this module."""
    return _main([*(sys.argv[1:] if argv is None else argv), '--run-learn'])


if __name__ == "__main__":
    raise SystemExit(main())
