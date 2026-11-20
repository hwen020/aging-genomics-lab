# Public study cards

## GSE5086 · analyzed in P1

Zahn et al., *Transcriptional Profiling of Aging in Human Muscle Reveals a Common Aging Signature*, PLoS Genetics 2(7):e115 (2006), PMID 16789832. The submitters profiled 81 surgically obtained human skeletal muscle samples, ages 16–89, on GPL570. The analyzed subset is the 62 rectus abdominis samples, ages 24–83. GEO accessions identify samples in the output. GEO's sample characteristics supply age and sex, and the sample source field supplies anatomical site. The matrix has no batch labels. [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5086) · [Original study](https://doi.org/10.1371/journal.pgen.0020115)

The submitters normalized probe-set intensities with dChip's perfect-match-only model, without a lower cutoff. The deposited measurements are on a linear intensity scale. The project validates the source metadata and pinned bytes, then removes probes with nonfinite or zero signal in the selected cohort and transforms the remaining signals with log2. GPL570 identifiers remain the statistical units. Gene symbols and Entrez IDs are joined as annotations, including one-to-many mappings. The paper averaged probe sets by gene before transforming. [GEO sample and value description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM114716) · [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) · [Paper methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC1513263/)

| Official artifact | Compressed bytes | SHA-256 |
|---|---:|---|
| [GSE5086 Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5086/matrix/GSE5086_series_matrix.txt.gz) | 15,793,968 | `2b2d82b8a3839e8516edc1b55bca87ee84093762741ee5ed215d7838f7c1f016` |
| [GPL570 annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz) | 8,471,521 | `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394` |

These are project-pinned hashes of the bytes fetched over HTTPS. Each successful acquisition records its UTC timestamp, URL, length and checksum in the ignored cache. Each load checks the bytes again. A changed upstream artifact requires review and an explicit pin update. Downloads use three attempts, 15-second socket timeout, a 90-second transfer budget per attempt and a 25 MiB compressed limit per file. Interrupted transfers retry from the beginning. The cache and reports containing GEO sample records stay local.

GEO permits public downloading and notes that submitters may assert rights in deposited material. The project found no separate dataset license on these GEO series records. The repository's MIT license covers its code. [GEO data-use disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)

## GSE164471 · reviewed candidate

Tumasian et al., *Skeletal muscle transcriptome in healthy aging*, Nature Communications 12:2014 (2021), PMID 33795677. The paper describes vastus lateralis biopsies from 53 GESTALT participants ages 22–83. GEO includes ages and sex in sample records, an annotated ENSG count file and a separate annotated ENST TPM file. The GEO processing record specifies Illumina HiSeq 2500, STAR alignment to hg19 with Ensembl v82 and featureCounts gene quantification. [GEO series and supplementary formats](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164471) · [Sample processing record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5011664) · [Original study and access statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC8016876/)

GSE5086 was selected because its complete processed matrix, sample metadata and probe annotation can be verified and analyzed on Windows with the project's Python dependencies. GSE164471 remains a catalog entry. Its gene-count file requires a validated count-aware design with donor and covariate checks and an independently managed DESeq2 or edgeR environment. TPM is a distinct quantification product.
