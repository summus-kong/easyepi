"""Summarize input QC without generating downstream biological conclusions."""
from easyepi._dispatch import dispatch


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    return dispatch("make-qc-report", argv)


if __name__ == "__main__":
    raise SystemExit(main())

