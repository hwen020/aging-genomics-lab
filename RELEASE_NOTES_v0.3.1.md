# Aging Genomics Lab v0.3.1

- Added a read-only Portfolio Demo for GSE5086 and GSE40279 with aggregate p-value charts, example feature estimates, methods and generation timestamps. Synthetic Demo remains interactive.
- Added independent numerical checks for HC3 effects and uncertainty, M-value clipping, plate coding and BH-FDR. The published-study report distinguishes unique paper records, unique GPL570 probes and matched pairs.
- Replaced the EWAS parser's pandas tokenizer with a bounded row reader after reproducing a Windows memory-allocation failure. Malformed rows now produce explicit errors.
- Added a reproducible snapshot builder, citation and contributing guides, deployment configuration, and local and CI health checks.
- The complete GEO acquisition, computation and exports are available in the local research workflow.

See the [scientific audit](docs/RELEASE_AUDIT.md), [deployment limits](docs/DEPLOYMENT.md) and [release checklist](RELEASE_CHECKLIST.md).
