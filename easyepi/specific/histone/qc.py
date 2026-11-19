"""Independent assay input QC."""
from easyepi._dispatch import dispatch


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    return dispatch("histone-qc", argv)


if __name__ == "__main__":
    raise SystemExit(main())

