# Data

The raw data are **not stored in this repository** (the count matrix is about 72 MB). The notebook
(`notebooks/PBMC_scRNAseq_Scanpy_GSM5227130.ipynb`) and `scripts/01_setup_and_data.py` download them
automatically from NCBI GEO into `data/GSM5227130/` and verify them before any analysis.

## Source

| Field | Value |
|---|---|
| GEO series | [GSE171555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171555) |
| GEO sample | [GSM5227130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5227130) — `PBMC_healthy [52BA23-00]` |
| Download folder | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5227nnn/GSM5227130/suppl/ |
| Format | Cell Ranger 3.1.0 filtered feature-barcode matrix (10x Genomics Single Cell 5′ v1, GRCh38) |
| Size | 13,801 cells × 36,601 features; 18,588,006 non-zero entries |

## Files and MD5 checksums

| File | MD5 |
|---|---|
| `GSM5227130_52BA23-00_barcodes.tsv.gz` | `542a03cbcc30842a55b7bc462cbb3f7f` |
| `GSM5227130_52BA23-00_features.tsv.gz` | `8361c1748151e98503e656c4e2d89a48` |
| `GSM5227130_52BA23-00_matrix.mtx.gz` | `9b4cadcc807f8b3c68d7cf31a9f630ca` |

## Manual download

If the automatic download does not work on your network, download the three files above from the GEO
sample page and place them in `data/GSM5227130/`, or put a `.tar` archive containing them in `data/`.
The code checks the MD5 checksums either way.

## Citation

Yu C, Littleton S, Giroux NS, *et al.* Mucosal-associated invariant T cell responses differ by sex in
COVID-19. *Med.* 2021;2(6):755–772.e5. doi:10.1016/j.medj.2021.04.008
