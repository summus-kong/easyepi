"""Annotate peak genomic features and distance to transcription start sites."""
from easyepi._dispatch import dispatch


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    return dispatch("annotate", argv)


if __name__ == "__main__":
    raise SystemExit(main())

