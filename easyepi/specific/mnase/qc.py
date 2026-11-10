"""MNase-seq QC."""
from easyepi._dispatch import dispatch


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    return dispatch("mnase-qc", argv)


if __name__ == "__main__":
    raise SystemExit(main())

