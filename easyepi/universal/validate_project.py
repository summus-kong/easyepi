"""Validate project inputs and record resolved references."""
from easyepi._dispatch import dispatch


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    return dispatch("validate-project", argv)


if __name__ == "__main__":
    raise SystemExit(main())

