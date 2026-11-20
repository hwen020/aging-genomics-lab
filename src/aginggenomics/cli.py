"""Command-line interface; python -m aginggenomics.cli --help."""
from __future__ import annotations

import argparse
from pathlib import Path

from .reporting import export_analysis
from .simulation import write_demo


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aginggenomics", description="Exploratory aging genomics analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo", help="Regenerate the clearly labeled synthetic demonstration")
    demo.add_argument("--out", type=Path, default=Path("data/demo"))
    analyze = sub.add_parser("analyze", help="Analyze two CSV inputs on the normalized log2 scale")
    analyze.add_argument("--expression", type=Path, required=True)
    analyze.add_argument("--metadata", type=Path, required=True)
    analyze.add_argument("--out", type=Path, default=Path("output/analysis"))
    analyze.add_argument("--reference", default=None)
    analyze.add_argument("--target", default=None)
    analyze.add_argument("--simulated", action="store_true", help="Flag files as synthetic in outputs")
    analyze.add_argument("--ack-normalized-log2", action="store_true",
                         help="Required for non-demo imports: confirm values are pre-normalized log2 expression")
    analyze.add_argument("--source", default="user-imported")
    public = sub.add_parser("public", help="Analyze the pinned GSE5086 public microarray study")
    public.add_argument("--study", choices=["GSE5086"], default="GSE5086")
    public.add_argument("--cache", type=Path, default=Path("data/cache"))
    public.add_argument("--out", type=Path, default=Path("output/GSE5086"))
    public.add_argument("--offline", action="store_true", help="Require a verified local cache")
    analyze.add_argument("--data-type", choices=["normalized-log2", "rna-seq-counts"], default="normalized-log2")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "demo":
        write_demo(args.out)
        print(f"Synthetic dataset written to {args.out}")
    elif args.command == "public":
        from .public_analysis import analyze_study, build_export
        from .public_studies import load_study

        try:
            study = load_study(args.cache, offline=args.offline)
            artifacts, provenance = build_export(study, analyze_study(study))
            args.out.mkdir(parents=True, exist_ok=True)
            for filename, artifact in artifacts.items():
                (args.out / filename).write_text(artifact, encoding="utf-8")
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(f"GSE5086: {provenance['qc']['n_samples']} samples analyzed. Report written to {args.out}")
    else:
        if args.data_type == "rna-seq-counts":
            parser.error("Raw RNA-seq counts require a count-aware model such as DESeq2.")
        if not args.simulated and not args.ack_normalized_log2:
            parser.error("For real/user data, pass --ack-normalized-log2 after verifying input scale")
        provenance = export_analysis(args.expression, args.metadata, args.out,
                                     reference=args.reference, target=args.target,
                                     simulated=args.simulated, source=args.source)
        print(f"Report written to {args.out}; analyzed {provenance['qc']['n_samples']} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
