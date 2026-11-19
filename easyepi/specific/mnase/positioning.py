"""Reserved entry for an established MNase nucleosome-positioning backend."""
import argparse


def main(argv=None):
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="MNase positioning is not implemented; use an established caller and import its BED output")
    parser.parse_args(argv)
    parser.error("No validated MNase positioning backend is connected. Occupancy summaries are available in specific.mnase.profile.")


if __name__ == "__main__":
    raise SystemExit(main())
