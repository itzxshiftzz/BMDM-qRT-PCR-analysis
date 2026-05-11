# RT-qPCR  Pipeline
## What the script does

1. **Loads** multiple sheets from an Excel file (one sheet per gene).
2. ** QC** — filters by housekeeping gene Cq and removes outliers via IQR.
3. **Univariate analysis** — 4 plots per gene + CSVs + text report.
4. **Multivariate analysis** — Heatmap, PCA, Volcano plots, Clustergram and Correlogram.

---

## You will need to edit SECTION 1

```r
input_file  <- 'path/to/your/file.xlsx'
output_dir  <- 'path/to/your/output/folder/'

gene_panel <- c(
  "ExcelSheetName" = "InternalLabel",
  ...
)
```

Everything else in the script adapts from what you define here.

---

## Configurable parameters in SECTION 1

| Parameter | Description | Default |
|---|---|---|
| `HK_CQ_THRESHOLD` | Maximum acceptable Cq for the housekeeping gene | `25` |
| `IQR_MULTIPLIER` | IQR multiplier for ΔCq outlier detection | `3` |
| `group_baseline` | Exact name of the control group | `"Control"` |
| `group_m1` | M1 group name (reference for volcano/efficacy) | `"LPS"` |
| `group_m2` | M2 group name | `"IL-4/IL_13"` |
| `groups_treatment` | Vector of treatment group names | `c("DEXA","MCC950","TDV19")` |
| `HEATMAP_MODE` | Heatmap colour scale: raw or standardised values | `"log2fc"` |
| `CLUSTER_GENES` | Cluster genes in the heatmap | `TRUE` |
| `CLUSTER_CONDITIONS` | Cluster conditions in the heatmap | `FALSE` |
| `DIST_METHOD` | Distance method for clustering | `"pearson"` |
| `CLUST_METHOD` | Linkage method for clustering | `"complete"` |

---

## Hardcoded parts

These elements are fixed in the script logic. If you need to change them,
you will have to edit the code.

- **Excel column names** — the script expects exactly `...1`, `...2`,
  `AVCp Housekeeping`, `∆Cq` and `Fold change`. If your template uses different
  headers, edit the `load_sheet()` function in SECTION 3.

- **The 4 univariate comparisons** — Efficacy , Polarization ,
  MOA and Benchmarking are built from the groups defined in SECTION 1,
  but the logic of *which groups go into each plot* is hardcoded in the SECTION 5
  loop. For example, Benchmarking always compares `MCC950` vs `TDV19`.

- **Colour palette** — uses `ggsci::pal_npg("nrc")`. To change it, edit the
  `cb_palette` line in SECTION 2.

- **Multiple testing correction in volcano plots** — always BH
  (Benjamini–Hochberg). Editable inside `s9_run_volcano()` in SECTION 6e.

- **Volcano label ** — genes are labelled when `p < 0.05 & |log2FC| > 1`,
  plus the gene with the highest |log2FC| per comparison. Editable in the
  `label_genes` variable in SECTION 6e.

- **Clustergram visual parameters** — `cellwidth = 46`, `cellheight = 14`,
  colour breaks fixed at `[-2.5, 2.5]`. Edit  in SECTION 6f.

---

## Requirements

**R packages**:

```r
install.packages(c(
  "readxl", "tidyverse", "rstatix", "ggprism", "ggpubr",
  "ggsci", "car", "ggrepel", "RColorBrewer", "ggdendro",
  "cowplot", "pheatmap", "ggcorrplot"
))
```

**Excel format** — one sheet per gene, following the structure of
`plantilla_qPCR_analysis.xlsx`. Sheet names must match the keys of `gene_panel`
exactly (case-sensitive).

---

## Output structure

```
output_dir/
├── univariate/
│   ├── IL1b/
│   │   ├── IL1b_Plot1_Efficacy.png / .pdf
│   │   ├── IL1b_Plot2_Polarization.png / .pdf
│   │   ├── IL1b_Plot3_MOA.png / .pdf
│   │   ├── IL1b_Plot4_Benchmarking.png / .pdf
│   │   ├── IL1b_Descriptive.csv
│   │   ├── IL1b_PostHoc_All.csv
│   │   ├── IL1b_Normality.csv
│   │   ├── IL1b_QC_HK_Removed.csv
│   │   ├── IL1b_QC_IQR_Removed.csv
│   │   └── IL1b_Report.txt
│   ├── TNF/ ...
│   └── (one subdirectory per gene)
└── multivariate/
    ├── PlotA_Heatmap.png / .pdf
    ├── PlotB_PCA.png / .pdf
    ├── PlotC_Volcano.png / .pdf
    ├── PlotD_SampleClustergram.png / .pdf
    ├── PlotE_Correlogram.png / .pdf
    ├── PCA_Scores.csv
    ├── PCA_Loadings.csv
    ├── PCA_Variance.csv
    ├── Heatmap_Log2FC.csv
    ├── Heatmap_GeneClusterOrder.csv  ← only if CLUSTER_GENES = TRUE
    ├── Volcano_Stats.csv
    ├── Clustergram_SampleOrder.csv
    ├── Correlation_r.csv
    └── Correlation_p.csv
```

---

##  Statistical logic (univariate)

The script selects the appropriate test automatically based on the data:

```
Are all groups normally distributed? (Shapiro-Wilk)
    NO  →  Kruskal-Wallis + Dunn post-hoc (BH correction)
    YES →  Equal variances? (Levene's test)
               YES  →  One-way ANOVA + Tukey HSD
               NO   →  Welch ANOVA + Games-Howell

(For 2-group comparisons)
    NOT normal  →  Wilcoxon rank-sum test
    Normal      →  Welch t-test
```

---

## Limitations

- **Small sample sizes** — Shapiro-Wilk has very low power with n < 5.
  With groups of small replicates, the test will almost always assume normality
  even when it cannot be properly verified. Results should be interpreted
  with caution.

- **Median ** — if an entire condition has no data for a given gene
  (e.g. all samples removed during QC), the global median is imputed.
  The script raises a `warning` in that case but does not stop execution.

- **Volcano plots with small n** — with low replicates per group, the t-test
  has very low statistical power. Volcano plots should be treated as
  exploratory.

- **Correlogram** — requires a minimum of 3 total samples. With few samples,
  Pearson p-values are unreliable.

- **Sequential execution** — if one gene fails inside the univariate loop,
  the entire script stops. Wrap the loop body in `tryCatch()` if you want
  execution to continue with the remaining genes after an isolated error.

---

## Quick notes

- All plots are saved as both `.png` (300 dpi) and `.pdf` simultaneously.
- Log₂FC in the heatmap and volcano is computed as `−ΔΔCq`
  (ΔCq_group − ΔCq_control, negated), so positive values indicate
  upregulation relative to the control.
- In the PCA, ΔCq values are negated (`× −1`) so that loading arrows point
  in the biologically intuitive direction (higher expression = positive arrow).
