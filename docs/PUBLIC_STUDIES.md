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

## GSE40279 · analyzed in P2

Hannum et al., *Genome-wide methylation profiles reveal quantitative views of human aging rates*, Molecular Cell 49(2):359–367 (2013), PMID 23177740. GEO describes 656 whole-blood samples measured on the Illumina HumanMethylation450 array. The study uses the author-processed average beta table and the recorded age, sex and plate. GPL13534 supplies hg19 coordinates, gene symbols, CpG island context and nearby-SNP annotations. GEO makes the files public under its deposited-data terms. [GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) · [Original paper](https://doi.org/10.1016/j.molcel.2012.10.016) · [GEO terms](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)

| Pinned source | Downloaded bytes | SHA-256 |
|---|---:|---|
| [Author beta matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_average_beta.txt.gz) | 1,223,258,643 | `d980065c653381ce29b5bfa54f462a1dd99853bb3796cccd0ec15d78043eb953` |
| [GEO brief sample records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279&targ=gsm&form=text&view=brief) | 1,100,547 | `b007e0448000fb177c405cceb704dd4165a996ee9408dd22bdfdb0877b910fa8` |
| [Submitter sample key](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_sample_key.txt.gz) | 4,563 | `dbb7c510da6a90bd2ff203564288b3b3400298c9d50217377a347f2ef2a9e7f9` |
| [GPL13534 SOFT platform prefix](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/soft/GPL13534_family.soft.gz) | 68,157,440 | `c8ca5ad3c0ac18fcf2870daa63c127bcd5f6fea128f292d769ce335b3027822e` |

The platform prefix contains the complete platform table. The client requires HTTP 206 and the exact Content-Range for this asset before reading it. The pin covers bytes 0 through 68,157,439 inclusive. The much larger platform family download includes sample data beyond the annotation table. Partial downloads retain a local `.part` file and resume only when the response range and size agree. The completed bytes must pass SHA-256 before promotion to the cache.
