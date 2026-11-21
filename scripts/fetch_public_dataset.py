"""Populate the verified local public data cache without running statistics."""
import argparse
from pathlib import Path

from aginggenomics.ingest.geo import acquire
from aginggenomics.public_studies import ACCESSION, ANNOTATION, MATRIX


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--study", choices=[ACCESSION, "GSE40279"], default=ACCESSION)
    args = parser.parse_args()
    if args.study == "GSE40279":
        from aginggenomics.methylation_study import acquire_study

        try:
            paths, records = acquire_study(args.cache, offline=args.offline)
        except OSError as exc:
            parser.error(str(exc))
        for path, record in zip(paths, records):
            print(f"{record['sha256']}  {path}")
        return
    for artifact in (MATRIX, ANNOTATION):
        try:
            path, record = acquire(artifact, args.cache / ACCESSION, offline=args.offline)
        except OSError as exc:
            parser.error(str(exc))
        print(f"{record['sha256']}  {path}")


if __name__ == "__main__":
    main()
