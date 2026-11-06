# Scanpy re-analysis of a healthy human PBMC sample (GSM5227130, GSE171555)

## Data

- Sample: GSM5227130, PBMC_healthy [52BA23-00], healthy control, day 0 (GEO series GSE171555).
- Chemistry: 10x Genomics Single Cell 5' v1; Cell Ranger 3.1.0 filtered matrix, GRCh38.
- Files: the three processed GEO files for this one sample, verified by MD5 checksum.
  This is a one-sample subset; the full series has 48 samples from 24 participants.
- Source study: Yu C, Littleton S, Giroux NS, et al. Mucosal-associated invariant T cell responses differ by sex in COVID-19. Med. 2021;2(6):755-772.e5. doi:10.1016/j.medj.2021.04.008

## Methods

Cells with 300–3,500 detected genes and <15% mitochondrial reads were kept, and genes detected
in fewer than 3 cells were removed. These thresholds were chosen for this re-analysis and differ
from Yu et al. (2021), whose criteria would retain 12,148 cells instead of 12,809.
Doublet scores were computed with Scrublet on raw counts. Counts were normalised to 10,000 per
cell and log-transformed. 2,000 highly variable genes (Seurat flavour) were selected after excluding
TCR/BCR V genes, then scaled; PCA (first 30 of 50 components used), a 15-nearest-neighbour graph, UMAP and Leiden clustering
(igraph, resolution 0.5) followed. Cluster markers were found with a Wilcoxon test on the full
log-normalised matrix and cell types assigned from marker evidence.
Software: Python 3.12.3, scanpy 1.11.5, anndata 0.12.19, numpy 2.0.2, pandas 2.2.3, scipy 1.16.3, igraph 1.0.0.

## Results

Of 13,801 cells, 12,809 passed quality control
(92.8%). Leiden clustering identified 15
clusters, annotated into 15 cell-type labels.

After excluding 654 cells flagged as non-mononuclear, low quality
or likely doublets, 12,155 mononuclear immune cells remained. In this sample, T cells
made up 64.7% of these cells, monocytes 18.9%, NK cells 9.1% and B cells 6.3%.

### Mononuclear immune cell composition

| Cell type | Cells | % of immune cells |
|---|---|---|
| Naive CD4+ T cells | 3,590 | 29.54% |
| CD4+ memory T cells | 2,415 | 19.87% |
| CD14+ classical monocytes | 2,037 | 16.76% |
| CD8+ effector memory T cells (GZMK+) | 1,133 | 9.32% |
| NK cells | 1,107 | 9.11% |
| B cells | 772 | 6.35% |
| Naive CD8+ T cells | 486 | 4.00% |
| FCGR3A+ non-classical monocytes | 259 | 2.13% |
| Regulatory T cells (Tregs) | 243 | 2.00% |
| Dendritic cells (cDC2 and pDC) | 87 | 0.72% |
| Haematopoietic progenitors (HSPCs) | 26 | 0.21% |

### Flagged populations (excluded from immune composition)

| Population | Cells | Share |
|---|---|---|
| Platelets | 185 | 1.44% of all cells |
| Erythrocytes (RBC contamination) | 170 | 1.33% of all cells |
| Unresolved (likely doublets) | 153 | 1.19% of all cells |
| Low-quality cells (mito/MALAT1-high) | 146 | 1.14% of all cells |

## Limitations

- One healthy donor only: proportions describe this sample, not a population, and no
  healthy-versus-disease comparison is made.
- Wilcoxon marker statistics are cell-level and descriptive; they are not donor-level
  differential expression.
- Annotation is marker-based. Closely related subsets (naive versus central memory T cells,
  cDC2 versus pDC within the 87-cell dendritic-cell cluster) would benefit from
  reference-based annotation (for example CellTypist or Azimuth).
- MAIT and gamma-delta T cells, the focus of the source study, could not be resolved as separate
  clusters: TRAV1-2 and SLC4A10 transcripts were detected in very few cells. They are likely
  contained within the memory T-cell and NK clusters.
- Scrublet's score distribution was not bimodal on this sample (automatic threshold called
  1 cell(s)), so no cells were removed on Scrublet calls alone. Instead, one cluster with
  mixed-lineage markers, few genes per cell and the highest mean doublet score was flagged
  as likely doublets and excluded from the composition.
- Ambient-RNA correction (for example SoupX) was not performed; low-level HBB signal across
  clusters indicates ambient red-cell RNA.

## Conclusion

This project demonstrates a reproducible Python/Scanpy workflow on public 10x Genomics PBMC data:
pinned software, checksum-verified input, quality control, doublet scoring, normalisation,
highly variable gene selection, PCA, UMAP, Leiden clustering, marker-gene analysis and
evidence-based cell-type annotation with explicit handling of low-quality and contaminating cells.
