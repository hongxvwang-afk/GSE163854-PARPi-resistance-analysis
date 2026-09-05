# Analysis methods

## Data source and sample selection

Gene-level RNA-sequencing count data for the PH039 HGSOC PDX model were
obtained from GEO dataset GSE163854. Five PARP inhibitor-sensitive tumors,
including three source-sensitive tumors (P2A1-P2A3) and two niraparib-treated
sensitive tumors (F1A3 and F1A4), were compared with five niraparib-resistant
tumors from resistant passage 2 (F351 and F352) and resistant passage 3
(F451-F453). Mixed-response tumors F1A5 and F251 were excluded.

## Heatmap processing

Raw counts were normalized using a median-of-ratios size-factor approach. A
gene-specific geometric mean was calculated from positive counts, and the size
factor for each sample was defined as the median ratio of its counts to the
corresponding gene-specific geometric means. Size factors were rescaled to have
a geometric mean of one. Normalized counts were transformed as
log2(normalized count + 1).

A curated panel of 70 genes involved in BLM-BTRR/core HR, ATR-checkpoint
signaling, replication-fork protection/Fanconi-anemia repair, DNA-end resection,
R-loop processing, cell-cycle/replication-stress regulation, NHEJ/alternative
end joining, and PI3K-AKT signaling was displayed. Expression was standardized
separately for each gene across the ten included tumors using:

```text
Z_ij = (x_ij - mean_i) / SD_i
```

where `x_ij` is the log2-transformed expression of gene `i` in tumor `j` and
`SD_i` is the sample standard deviation across the ten tumors. The color scale
was clipped to -2.4 to +2.4. Genes and samples were displayed in a predefined
biological order without unsupervised clustering. The heatmap is descriptive;
no statistical tests were performed on the row Z-scores.

## Core-HR differential-expression analysis

Differential expression was analyzed from raw counts using PyDESeq2 version
0.5.2. Duplicate gene symbols, if present, were collapsed by summing their
counts. Genes with zero total counts across the ten included tumors were
removed. All remaining genes were fitted using a negative-binomial generalized
linear model with the design `~ condition`, where condition was classified as
sensitive or resistant. Resistant-versus-sensitive log2 fold changes were
estimated without additional fold-change shrinkage.

Statistical significance was assessed using the DESeq2 Wald test. Cook's
distance filtering and independent filtering were enabled. P values were
adjusted across all genes entering the differential-expression analysis using
the Benjamini-Hochberg procedure.

A defined panel of 15 BLM-BTRR/core HR genes (BLM, TOP3A, RMI1, RMI2, BRCA1,
BRCA2, PALB2, BARD1, BRIP1, RAD51, RAD51B, RAD51C, RAD51D, XRCC2 and XRCC3) was
extracted from the genome-wide results. Genes were ranked by absolute log2 fold
change, and the eight largest changes were displayed. Statistical labels were:
`ns`, FDR >= 0.05; `*`, FDR < 0.05; `**`, FDR < 0.01; and `***`, FDR < 0.001.
The x-axis break accommodates the substantially larger RAD51C fold change and
does not affect the analysis.

## Statistical interpretation

All analyzed tumors derive from one serially treated PH039 PDX lineage. They
should not be treated as independent patient-level biological replicates.
Accordingly, the DESeq2 results quantify variation within this PDX model and
should be interpreted as exploratory or descriptive evidence rather than
population-level inference.

