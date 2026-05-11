# =============================================================================
# RT-qPCR  PIPELINE  | Macrophage-Fibroblast Crosstalk
# Univariate (all genes) + Multivariate (Heatmap · PCA · Volcano ·
#                                         Clustergram · Correlogram)
# -----------------------------------------------------------------------------
# TO ADAPT TO NEW EXPERIMENTS: edit only SECTION 1.
# The pipeline auto-loops over all genes for univariate analysis and feeds
# the same QC-passed data into the multivariate module.
# =============================================================================


# =============================================================================
# SECTION 1: CONFIGURATION  ←  SECTION YOU MUST EDIT
# =============================================================================

input_file  <- ''
output_dir  <- ''

# --- Gene panel --------------------------------------------------------------
# Each entry: "ExcelSheetName" = "InternalLabel"
# Sheet names must match EXACTLY (case-sensitive).
# InternalLabel is used in filenames and plot titles.
gene_panel <- c(
  "IL-1b" = "IL1b",
  "TNF"   = "TNF",
  "COX2"  = "COX2",
  "NFKB"  = "NFKB",
  "DUSP1" = "DUSP1"
)

# --- Quality control ---------------------------------------------------------
HK_CQ_THRESHOLD <- 25   # HK Cq above this = amplification failure → remove
IQR_MULTIPLIER  <- 3    # IQR multiplier for per-condition ΔCq outlier removal

# --- Experimental design -----------------------------------------------------
# Values must match EXACT spelling in your Excel condition column.
group_baseline   <- "Control"
group_m1         <- "LPS"
group_m2         <- "IL-4/IL_13"
groups_treatment <- c("DEXA", "MCC950", "TDV19")
condition_levels <- c(group_baseline, group_m1, group_m2, groups_treatment)

# --- Multivariate settings ---------------------------------------------------
PCA_ELLIPSE_LEVEL  <- 0.95
HEATMAP_MODE       <- "log2fc"   # "log2fc" | "zscore"
CLUSTER_GENES      <- TRUE
CLUSTER_CONDITIONS <- FALSE
DIST_METHOD        <- "pearson"  # "pearson" | "euclidean"
CLUST_METHOD       <- "complete" # "complete" | "average" | "ward.D2"
SHOW_DENDROGRAMS   <- TRUE
DEND_GENE_WIDTH    <- 0.18
DEND_COND_HEIGHT   <- 0.20


# =============================================================================
# SECTION 2: PACKAGES
# =============================================================================
# Run once to install:
# install.packages(c("readxl","tidyverse","rstatix","ggprism","ggpubr",
#                    "ggsci","car","ggrepel","RColorBrewer","ggdendro",
#                    "cowplot","pheatmap","ggcorrplot"))

suppressPackageStartupMessages({
  library(readxl);     library(tidyverse);  library(rstatix)
  library(ggprism);    library(ggpubr);     library(ggsci)
  library(car);        library(ggrepel);    library(RColorBrewer)
  library(ggdendro);   library(cowplot);    library(pheatmap)
  library(ggcorrplot)
})

# Output sub-directories
uni_dir   <- file.path(output_dir, "univariate")
multi_dir <- file.path(output_dir, "multivariate")
for (d in c(uni_dir, multi_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cb_palette <- setNames(ggsci::pal_npg("nrc")(length(condition_levels)),
                       condition_levels)

# Convenience: sheet names and internal labels as separate vectors
sheet_names <- names(gene_panel)

gene_labels      <- unname(gene_panel)          # plain character vector, no names
sheet_to_label   <- gene_panel                  # named vector for lookup inside loop


# =============================================================================
# SECTION 3: SHARED HELPER FUNCTIONS
# =============================================================================

# --- QC helpers --------------------------------------------------------------
iqr_flag <- function(x, k) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iq <- q3 - q1
  if (iq == 0) return(rep(TRUE, length(x)))
  x >= (q1 - k * iq) & x <= (q3 + k * iq)
}

load_sheet <- function(sheet, gene_name) {
  read_excel(input_file, sheet = sheet) |>
    rename(Condition  = `...1`, Sample = `...2`,
           HK_Cq = `AVCp Housekeeping`, DeltaCq = `∆Cq`,
           FoldChange = `Fold change`) |>
    tidyr::fill(Condition, .direction = "down") |>
    select(Condition, Sample, HK_Cq, DeltaCq, FoldChange) |>
    filter(!is.na(Sample), !is.na(DeltaCq),
           !Sample %in% c("Unnamed: 1", "...2", "Sample")) |>
    mutate(
      across(c(HK_Cq, DeltaCq, FoldChange),
             ~ as.numeric(str_replace_all(as.character(.x), ",", "."))),
      Condition = factor(str_trim(Condition), levels = condition_levels),
      Gene      = gene_name
    ) |>
    filter(!is.na(Condition))
}

apply_qc <- function(df) {
  df |>
    filter(!is.na(HK_Cq), HK_Cq <= HK_CQ_THRESHOLD) |>
    group_by(Condition) |>
    mutate(keep = iqr_flag(DeltaCq, k = IQR_MULTIPLIER)) |>
    ungroup() |>
    filter(keep) |>
    select(-keep) |>
    mutate(Condition = droplevels(Condition))
}

# --- Statistical engine ------------------------------------------------------
add_sig <- function(df, col = "p.adj") {
  df |> mutate(sig_label = case_when(
    .data[[col]] < 0.001 ~ "***",
    .data[[col]] < 0.01  ~ "**",
    .data[[col]] < 0.05  ~ "*",
    TRUE                 ~ "ns"))
}

check_normality <- function(df, grp = "Condition", val = "DeltaCq") {
  res   <- df |> group_by(across(all_of(grp))) |>
    rstatix::shapiro_test(!!sym(val)) |> ungroup()
  min_n <- min(count(df, across(all_of(grp)))$n)
  list(results    = res,
       all_normal = all(res$p >= 0.05, na.rm = TRUE) && min_n >= 3,
       min_n      = min_n)
}

check_variance <- function(df, grp = "Condition", val = "DeltaCq") {
  lev <- car::leveneTest(as.formula(paste(val, "~", grp)), data = df, center = median)
  list(equal_var = lev$`Pr(>F)`[1] >= 0.05)
}

run_multi_stats <- function(df, grp = "Condition", val = "DeltaCq", label = "") {
  message("    [Stats] ", label)
  if (n_distinct(df[[grp]]) < 2) return(NULL)
  frm  <- as.formula(paste(val, "~", grp))
  norm <- check_normality(df, grp, val)
  if (!norm$all_normal) {
    return(list(test_name = "Kruskal-Wallis", posthoc_name = "Dunn (BH)",
                test_result    = df |> rstatix::kruskal_test(frm),
                posthoc_result = df |> rstatix::dunn_test(frm, p.adjust.method = "BH") |> add_sig(),
                normality = norm$results))
  }
  lev <- check_variance(df, grp, val)
  if (lev$equal_var) {
    return(list(test_name = "One-Way ANOVA", posthoc_name = "Tukey HSD",
                test_result    = df |> rstatix::anova_test(frm),
                posthoc_result = df |> rstatix::tukey_hsd(frm) |> add_sig(),
                normality = norm$results))
  }
  list(test_name = "Welch ANOVA", posthoc_name = "Games-Howell",
       test_result    = df |> rstatix::welch_anova_test(frm),
       posthoc_result = df |> rstatix::games_howell_test(frm) |> add_sig(),
       normality = norm$results)
}

run_two_stats <- function(df, grp = "Condition", val = "DeltaCq", label = "") {
  message("    [Stats] ", label)
  if (n_distinct(df[[grp]]) != 2) return(NULL)
  frm  <- as.formula(paste(val, "~", grp))
  norm <- check_normality(df, grp, val)
  if (!norm$all_normal) {
    res <- df |> rstatix::wilcox_test(frm, paired = FALSE) |>
      mutate(p.adj = p) |> add_sig("p.adj")
    return(list(test_name = "Wilcoxon", posthoc_name = "N/A",
                test_result = res, posthoc_result = res, normality = norm$results))
  }
  res <- df |> rstatix::t_test(frm, paired = FALSE, var.equal = FALSE) |>
    mutate(p.adj = p) |> add_sig("p.adj")
  list(test_name = "Welch t-test", posthoc_name = "N/A",
       test_result = res, posthoc_result = res,
       normality = norm$results)
}

# --- Plot helpers ------------------------------------------------------------
build_annot <- function(ph, ref = NULL, x_pos, max_fc, y_step) {
  if (is.null(ph) || nrow(ph) == 0) return(tibble())
  out <- if (!is.null(ref)) filter(ph, group1 == ref | group2 == ref) else ph
  out |> filter(p.adj < 0.05) |>
    mutate(xmin = x_pos[as.character(group1)],
           xmax = x_pos[as.character(group2)],
           label = sig_label) |>
    filter(!is.na(xmin), !is.na(xmax)) |>
    arrange(abs(xmin - xmax)) |>
    mutate(y.position = max_fc + y_step * row_number())
}

make_uni_plot <- function(df, stats_obj, title_str, subtitle_str, ref = NULL) {
  if (is.null(stats_obj)) return(NULL)
  lvls   <- levels(df$Condition)
  x_pos  <- setNames(seq_along(lvls), lvls)
  max_fc <- max(df$FoldChange, na.rm = TRUE)
  annot  <- build_annot(stats_obj$posthoc_result, ref, x_pos, max_fc, max_fc * 0.13)
  
  ggplot(df, aes(Condition, FoldChange, color = Condition)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.70, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.45, linewidth = 0.65, color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.22, linewidth = 0.65, color = "black") +
    { if (nrow(annot) > 0)
      ggpubr::stat_pvalue_manual(annot, label = "label",
                                 xmin = "xmin", xmax = "xmax",
                                 y.position = "y.position",
                                 tip.length = 0.015, size = 4,
                                 inherit.aes = FALSE) } +
    scale_color_manual(values = cb_palette[lvls]) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = NULL,
         y = expression(paste("Relative Expression  (2"^{-Delta*Delta*Cq}, ")")),
         title = title_str, subtitle = subtitle_str) +
    theme_prism(base_size = 12) +
    theme(plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle      = element_text(size = 8, hjust = 0.5, color = "grey45"),
          axis.text.x        = element_text(angle = 30, hjust = 1, size = 10),
          panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4),
          plot.margin        = margin(12, 20, 10, 10))
}

save_gg <- function(p, path, w = 7.5, h = 5.5) {
  ggsave(paste0(path, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(path, ".pdf"), p, width = w, height = h)
}

save_csv <- function(df, path) {
  write.csv(as.data.frame(df), paste0(path, ".csv"), row.names = FALSE)
}

fmt_sub <- function(s) if (!is.null(s)) paste0(s$test_name, " + ", s$posthoc_name) else "N/A"


# =============================================================================
# SECTION 4: VALIDATE SHEETS EXIST
# =============================================================================

message("\n[0] Validating Excel sheets")
avail <- readxl::excel_sheets(input_file)
miss  <- setdiff(sheet_names, avail)
if (length(miss) > 0) {
  stop("Sheet(s) not found: ", paste(miss, collapse = ", "),
       "\nAvailable: ", paste(avail, collapse = ", "))
}
message("  All ", length(sheet_names), " sheets found: ",
        paste(sheet_names, collapse = ", "))


# =============================================================================
# SECTION 5: UNIVARIATE LOOP — all genes
# =============================================================================
# For each gene: load → QC → 4 plots → CSVs → text report
# All outputs go to output_dir/univariate/<GENE>/

message("\n", strrep("=", 60))
message("  UNIVARIATE MODULE — ", length(gene_panel), " genes")
message(strrep("=", 60))

# Storage for multivariate module (built while looping)
all_long_clean <- list()

for (sheet in sheet_names) {
  
 
  glabel <- sheet_to_label[sheet]
  gdir   <- file.path(uni_dir, glabel)
  dir.create(gdir, showWarnings = FALSE)
  
  message("\n--- ", glabel, " (sheet: '", sheet, "') ---")
  
  # Load & QC
  raw   <- load_sheet(sheet, glabel)
  clean <- apply_qc(raw)
  all_long_clean[[glabel]] <- clean   # keep for multivariate
  
  removed_hk  <- raw |> filter(is.na(HK_Cq) | HK_Cq > HK_CQ_THRESHOLD)
  removed_iqr <- raw |>
    filter(!is.na(HK_Cq), HK_Cq <= HK_CQ_THRESHOLD) |>
    group_by(Condition) |>
    mutate(keep = iqr_flag(DeltaCq, k = IQR_MULTIPLIER)) |>
    ungroup() |>
    filter(!keep)
  
  message("  QC: loaded=", nrow(raw),
          " | HK removed=", nrow(removed_hk),
          " | IQR removed=", nrow(removed_iqr),
          " | final=", nrow(clean))
  
  # Descriptive stats
  desc <- clean |>
    group_by(Condition) |>
    summarise(n        = n(),
              Mean_FC  = mean(FoldChange),
              SEM_FC   = sd(FoldChange) / sqrt(n()),
              Mean_dCq = mean(DeltaCq),
              SEM_dCq  = sd(DeltaCq) / sqrt(n()),
              .groups  = "drop")
  
  # Statistical comparisons
  data_polar <- clean |>
    filter(Condition %in% c(group_baseline, group_m1, group_m2)) |>
    mutate(Condition = droplevels(Condition))
  data_moa <- clean |>
    filter(Condition %in% groups_treatment) |>
    mutate(Condition = droplevels(Condition))
  data_bench <- clean |>
    filter(Condition %in% c("MCC950", "TDV19")) |>
    mutate(Condition = droplevels(Condition))
  
  s_all   <- run_multi_stats(clean,       label = "All conditions")
  s_polar <- run_multi_stats(data_polar,  label = "Polarization")
  s_moa   <- run_multi_stats(data_moa,    label = "MOA")
  s_bench <- run_two_stats(data_bench,    label = "Benchmarking")
  
  # 4 Plots
  p1 <- make_uni_plot(clean, s_all,
                      paste0(glabel, " | Efficacy Landscape"),
                      paste0("All conditions vs ", group_m1, " | ", fmt_sub(s_all)),
                      ref = group_m1)
  p2 <- make_uni_plot(data_polar, s_polar,
                      paste0(glabel, " | Polarization Spectrum"),
                      paste0("M0 / M1 / M2 | ", fmt_sub(s_polar)))
  p3 <- make_uni_plot(data_moa, s_moa,
                      paste0(glabel, " | Mechanism of Action"),
                      paste0("DEXA vs MCC950 vs TDV19 | ", fmt_sub(s_moa)))
  p4 <- make_uni_plot(data_bench, s_bench,
                      paste0(glabel, " | Clinical Benchmarking"),
                      paste0("MCC950 vs TDV19 | ",
                             if (!is.null(s_bench)) s_bench$test_name else "N/A"))
  
  for (pair in list(list(p1, "Plot1_Efficacy"),    list(p2, "Plot2_Polarization"),
                    list(p3, "Plot3_MOA"),          list(p4, "Plot4_Benchmarking"))) {
    if (!is.null(pair[[1]])) {
      save_gg(pair[[1]], file.path(gdir, paste0(glabel, "_", pair[[2]])))
    }
  }
  
  # CSV exports
  save_csv(desc,         file.path(gdir, paste0(glabel, "_Descriptive")))
  save_csv(removed_hk,  file.path(gdir, paste0(glabel, "_QC_HK_Removed")))
  save_csv(removed_iqr |> select(Condition, Sample, HK_Cq, DeltaCq, FoldChange),
           file.path(gdir, paste0(glabel, "_QC_IQR_Removed")))
  
  for (pair in list(list(s_all,   "PostHoc_All"),
                    list(s_polar, "PostHoc_Polar"),
                    list(s_moa,   "PostHoc_MOA"),
                    list(s_bench, "PostHoc_Bench"),
                    list(s_all,   "Normality"))) {
    obj <- pair[[1]]
    tag <- pair[[2]]
    if (!is.null(obj)) {
      dat <- if (tag == "Normality") obj$normality else obj$posthoc_result
      if (!is.null(dat)) {
        save_csv(dat, file.path(gdir, paste0(glabel, "_", tag)))
      }
    }
  }
  
  # Text report
  rpt <- file.path(gdir, paste0(glabel, "_Report.txt"))
  capture.output({
    cat("RT-qPCR REPORT —", glabel, "\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat(strrep("=", 55), "\n\n[QC]\n")
    cat("  Loaded:", nrow(raw), "| HK removed:", nrow(removed_hk),
        "| IQR removed:", nrow(removed_iqr), "| Final:", nrow(clean), "\n\n")
    cat("[DESCRIPTIVE — Fold Change]\n")
    print(select(desc, Condition, n, Mean_FC, SEM_FC))
    
    report_block <- function(label, s, ref = NULL) {
      cat("\n[", label, "]\n")
      if (is.null(s)) { cat("  Skipped\n"); return() }
      cat("  Test:", s$test_name, "+", s$posthoc_name, "\n")
      ph  <- if (!is.null(ref)) {
        filter(s$posthoc_result, group1 == ref | group2 == ref)
      } else {
        s$posthoc_result
      }
      sig <- filter(ph, p.adj < 0.05)
      if (nrow(sig) == 0) {
        cat("  No significant pairs\n")
      } else {
        for (i in seq_len(nrow(sig))) {
          cat(sprintf("  %s vs %s: p.adj=%.4f [%s]\n",
                      sig$group1[i], sig$group2[i],
                      sig$p.adj[i],  sig$sig_label[i]))
        }
      }
    }
    
    report_block("Efficacy Landscape",   s_all,   ref = group_m1)
    report_block("Polarization",         s_polar)
    report_block("MOA",                  s_moa)
    report_block("Benchmarking",         s_bench)
    cat("\n", strrep("=", 55), "\nEND\n")
  }, file = rpt)
  
  message("  Outputs → ", gdir)
}

message("\n  Univariate complete — ", length(gene_panel), " genes processed.")


# =============================================================================
# SECTION 6: MULTIVARIATE MODULE
# =============================================================================
# Uses the QC-passed long data collected during the univariate loop.
# gene_labels is already a plain (unnamed) character vector — safe for
# all_of() inside across() and select().

message("\n", strrep("=", 60))
message("  MULTIVARIATE MODULE")
message(strrep("=", 60))

# --- 6a: Build wide matrices -------------------------------------------------
long_all <- bind_rows(all_long_clean)   # Gene column already set per iteration

impute_median_wide <- function(wide_df, gene_cols) {
  out <- wide_df |>
    group_by(Condition) |>
    mutate(across(all_of(gene_cols),
                  ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) |>
    ungroup()
  out |> mutate(across(all_of(gene_cols), ~ {
    n_na <- sum(is.na(.x))
    if (n_na > 0) {
      warning("[impute] '", cur_column(), "': ", n_na,
              " NA(s) remain — global median used.", call. = FALSE)
    }
    ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)
  }))
}

wide_dcq <- long_all |>
  select(Condition, Sample, Gene, DeltaCq) |>
  pivot_wider(names_from = Gene, values_from = DeltaCq, values_fn = mean) |>
  impute_median_wide(gene_labels)

wide_fc <- long_all |>
  select(Condition, Sample, Gene, FoldChange) |>
  pivot_wider(names_from = Gene, values_from = FoldChange, values_fn = mean) |>
  impute_median_wide(gene_labels)

message("  Wide matrix: ", nrow(wide_dcq), " samples x ",
        length(gene_labels), " genes")
print(count(wide_dcq, Condition))


# --- 6b: Heatmap helpers -----------------------------------------------------
compute_dist <- function(mat, method) {
  if (method == "pearson") {
    zv  <- apply(mat, 1, var, na.rm = TRUE) < 1e-10
    if (any(zv)) mat <- mat[!zv, , drop = FALSE]
    as.dist(1 - cor(t(mat), method = "pearson", use = "pairwise.complete.obs"))
  } else {
    dist(mat, method = method)
  }
}

make_dend <- function(hclust_obj, n, flip = FALSE) {
  dd <- ggdendro::dendro_data(as.dendrogram(hclust_obj), type = "rectangle")
  p  <- ggplot(ggdendro::segment(dd)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey35", linewidth = 0.45) +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    theme_void()
  
  if (flip) {
    p + coord_flip() + theme(plot.margin = margin(12, 0, 8, 5))
  } else {
    p + theme(plot.margin = margin(5, 5, 0, 5))
  }
}


# --- 6c: Plot A — Heatmap ----------------------------------------------------
message("[M1] Heatmap")

mean_dcq <- wide_dcq |>
  group_by(Condition) |>
  summarise(across(all_of(gene_labels), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

ctrl_row <- mean_dcq |> filter(Condition == group_baseline) |>
  select(all_of(gene_labels))
if (nrow(ctrl_row) == 0) stop("Baseline '", group_baseline, "' not found.")

log2fc_mat <- mean_dcq |>
  mutate(across(all_of(gene_labels),
                ~ -(.x - as.numeric(ctrl_row[, cur_column()])))) |>
  column_to_rownames("Condition")
log2fc_mat <- log2fc_mat[condition_levels[condition_levels %in%
                                            rownames(log2fc_mat)],
                         gene_labels, drop = FALSE]

display_mat          <- t(log2fc_mat)
log2fc_for_threshold <- display_mat

if (HEATMAP_MODE == "zscore") {
  bad <- names(which(apply(display_mat, 1, sd, na.rm = TRUE) < 1e-10))
  if (length(bad) > 0) {
    display_mat          <- display_mat[!rownames(display_mat) %in% bad, , drop = FALSE]
    log2fc_for_threshold <- log2fc_for_threshold[
      !rownames(log2fc_for_threshold) %in% bad, , drop = FALSE]
  }
  display_mat <- t(scale(t(display_mat)))
}

do_clust_g <- CLUSTER_GENES      && nrow(display_mat) >= 2
do_clust_c <- CLUSTER_CONDITIONS && ncol(display_mat) >= 2

gene_factor_levels <- if (do_clust_g) {
  gene_hclust <- hclust(compute_dist(log2fc_for_threshold, DIST_METHOD), CLUST_METHOD)
  gene_order  <- gene_hclust$labels[gene_hclust$order]
  message("  Gene cluster order: ", paste(gene_order, collapse = " → "))
  gene_order
} else {
  rev(rownames(display_mat))
}

cond_order <- if (do_clust_c) {
  cond_hclust <- hclust(compute_dist(t(log2fc_for_threshold), DIST_METHOD), CLUST_METHOD)
  cond_hclust$labels[cond_hclust$order]
} else {
  condition_levels[condition_levels %in% colnames(display_mat)]
}

txt_thr <- max(log2fc_for_threshold, na.rm = TRUE) * 0.60

heatmap_long <- as.data.frame(display_mat) |>
  rownames_to_column("Gene") |>
  pivot_longer(-Gene, names_to = "Condition", values_to = "DisplayVal") |>
  left_join(as.data.frame(log2fc_for_threshold) |>
              rownames_to_column("Gene") |>
              pivot_longer(-Gene, names_to = "Condition", values_to = "Log2FC"),
            by = c("Gene", "Condition")) |>
  mutate(Gene      = factor(Gene,      levels = gene_factor_levels),
         Condition = factor(Condition, levels = cond_order))

fill_scale <- if (HEATMAP_MODE == "log2fc") {
  scale_fill_gradientn(
    colours = c("white","#FEE5D9","#FCAE91","#FB6A4A","#CB181D","#67000D"),
    limits  = c(0, max(heatmap_long$DisplayVal, na.rm = TRUE)),
    name    = "Log\u2082FC",
    guide   = guide_colourbar(barwidth = 0.85, barheight = 5.5,
                              ticks = FALSE, title.position = "top",
                              title.hjust = 0.5))
} else {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits  = c(-2.5, 2.5),
    name    = "Z-score",
    guide   = guide_colourbar(barwidth = 0.85, barheight = 5.5,
                              ticks = FALSE, title.position = "top",
                              title.hjust = 0.5))
}

hm_title <- if (HEATMAP_MODE == "log2fc") {
  "Macrophage Phenotypic Heatmap\nColour = Log\u2082FC | Numbers = Log\u2082FC"
} else {
  "Macrophage Phenotypic Heatmap\nColour = Row Z-score | Numbers = Log\u2082FC"
}

p_heatmap <- ggplot(heatmap_long, aes(Condition, Gene, fill = DisplayVal)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label  = sprintf("%.1f", Log2FC),
                colour = ifelse(Log2FC > txt_thr, "white", "grey15")), size = 3.5) +
  scale_colour_identity() +
  fill_scale +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = hm_title, x = NULL, y = NULL) +
  theme_prism(base_size = 12) +
  theme(plot.title      = element_text(face = "bold", size = 12, hjust = 0.5,
                                       margin = margin(b = 10)),
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y     = element_text(size = 12),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        panel.grid      = element_blank(),
        panel.border    = element_blank(),
        legend.position = "right",
        plot.margin     = margin(12, 5, 8, 5))

p_dend_gene <- if (do_clust_g && SHOW_DENDROGRAMS) {
  make_dend(gene_hclust, nrow(display_mat), flip = TRUE)
} else {
  NULL
}
p_dend_cond <- if (do_clust_c && SHOW_DENDROGRAMS) {
  make_dend(cond_hclust, ncol(display_mat), flip = FALSE)
} else {
  NULL
}

hg <- !is.null(p_dend_gene)
hc <- !is.null(p_dend_cond)


final_heatmap <- if (!hg && !hc) {
  p_heatmap
} else if (hg && !hc) {
  cowplot::plot_grid(p_dend_gene, p_heatmap, ncol = 2,
                     rel_widths = c(DEND_GENE_WIDTH, 1 - DEND_GENE_WIDTH),
                     align = "h", axis = "tb")
} else if (!hg && hc) {
  cowplot::plot_grid(p_dend_cond, p_heatmap, nrow = 2,
                     rel_heights = c(DEND_COND_HEIGHT, 1 - DEND_COND_HEIGHT),
                     align = "v", axis = "lr")
} else {
  top <- cowplot::plot_grid(ggplot() + theme_void(), p_dend_cond, ncol = 2,
                            rel_widths = c(DEND_GENE_WIDTH, 1 - DEND_GENE_WIDTH))
  bot <- cowplot::plot_grid(p_dend_gene, p_heatmap, ncol = 2,
                            rel_widths = c(DEND_GENE_WIDTH, 1 - DEND_GENE_WIDTH),
                            align = "h", axis = "tb")
  cowplot::plot_grid(top, bot, nrow = 2,
                     rel_heights = c(DEND_COND_HEIGHT, 1 - DEND_COND_HEIGHT),
                     align = "v", axis = "lr")
}

save_gg(final_heatmap, file.path(multi_dir, "PlotA_Heatmap"), w = 8.5, h = 5.5)


# --- 6d: Plot B — PCA --------------------------------------------------------
message("[M2] PCA")

pca_mat <- -1 * (wide_dcq |> select(all_of(gene_labels)) |> as.matrix())
rownames(pca_mat) <- paste(wide_dcq$Condition, wide_dcq$Sample, sep = "_")
pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
pct_var <- round(summary(pca_res)$importance[2, ] * 100, 1)

pca_scores <- as_tibble(pca_res$x[, 1:2]) |>
  mutate(Condition = factor(wide_dcq$Condition, levels = condition_levels),
         Sample    = wide_dcq$Sample)

centroids <- pca_scores |>
  group_by(Condition) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), n = n(), .groups = "drop") |>
  mutate(label = paste0(Condition, "\n(n=", n, ")"))

p_pca_scores <- ggplot(pca_scores, aes(PC1, PC2, colour = Condition, fill = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey85", linewidth = 0.3) +
  stat_ellipse(type = "t", level = PCA_ELLIPSE_LEVEL, geom = "polygon",
               alpha = 0.05, show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.8, show.legend = FALSE) +
  geom_label_repel(data = centroids, aes(label = label), fontface = "bold",
                   size = 3, fill = "white", alpha = 0.9, box.padding = 0.5,
                   show.legend = FALSE) +
  scale_colour_manual(values = cb_palette) +
  scale_fill_manual(values   = cb_palette) +
  labs(x = paste0("PC1 (", pct_var["PC1"], "%)"),
       y = paste0("PC2 (", pct_var["PC2"], "%)"),
       title = "A. Macrophage Phenotypic Clusters") +
  theme_prism(base_size = 11)

loadings_df <- as_tibble(pca_res$rotation[, 1:2], rownames = "Gene")

p_pca_load <- ggplot(loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey85", linewidth = 0.3) +
  annotate("path", colour = "grey90", linewidth = 0.3,
           x = cos(seq(0, 2 * pi, length.out = 100)),
           y = sin(seq(0, 2 * pi, length.out = 100))) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")),
               colour = "firebrick", linewidth = 0.8) +
  geom_text_repel(aes(x = PC1, y = PC2, label = Gene),
                  fontface = "bold.italic", size = 4) +
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  labs(x = "Weight on PC1", y = "Weight on PC2",
       title = "B. Gene Contributions (Loadings)") +
  theme_prism(base_size = 11)

save_gg(cowplot::plot_grid(p_pca_scores, p_pca_load, ncol = 2, rel_widths = c(1.1, 0.9)),
        file.path(multi_dir, "PlotB_PCA"), w = 10, h = 5)


# --- 6e: Plot C — Volcano plots ----------------------------------------------
message("[M3] Volcano plots")

s9_run_volcano <- function(dcq_wide, ref_cond, test_conds, genes) {
  purrr::map_dfr(test_conds, function(tc) {
    ref  <- dcq_wide |> filter(Condition == ref_cond)
    tst  <- dcq_wide |> filter(Condition == tc)
    if (nrow(ref) < 2 || nrow(tst) < 2) return(NULL)
    purrr::map_dfr(genes, function(g) {
      xr <- na.omit(ref[[g]]); xt <- na.omit(tst[[g]])
      log2fc <- NA_real_; pval <- NA_real_
      if (length(xr) >= 2 && length(xt) >= 2) {
        if (var(xr) < 1e-10 && var(xt) < 1e-10) {
          log2fc <- mean(xr) - mean(xt)
        } else {
          tt <- tryCatch(t.test(xr, xt, var.equal = FALSE), error = function(e) NULL)
          if (!is.null(tt)) { log2fc <- mean(xr) - mean(xt); pval <- tt$p.value }
        }
      }
      tibble(comparison = paste0(tc, "_vs_", ref_cond),
             treatment  = tc, gene = g, log2fc, p_value = pval)
    })
  }) |>
    group_by(comparison) |>
    mutate(p_adj      = p.adjust(p_value, method = "BH"),
           neg_log10p = -log10(p_value),
           sig = case_when(
             is.na(p_value)                    ~ "No test",
             p_value < 0.05 & abs(log2fc) > 1 ~ "p<0.05 & |log2FC|>1",
             p_value < 0.05                    ~ "p<0.05 only",
             abs(log2fc) > 1                   ~ "|log2FC|>1 only",
             TRUE                              ~ "NS")) |>
    ungroup()
}

vstat <- s9_run_volcano(wide_dcq, group_m1, groups_treatment, gene_labels)

vcols <- c("p<0.05 & |log2FC|>1" = "#D62728", "p<0.05 only" = "#FF7F0E",
           "|log2FC|>1 only" = "#1F77B4", "NS" = "grey70", "No test" = "grey40")
vcols <- vcols[names(vcols) %in% unique(vstat$sig)]

ymax  <- max(quantile(vstat$neg_log10p, 0.99, na.rm = TRUE) * 1.1,
             -log10(0.05) * 1.5, na.rm = TRUE)
vstat <- vstat |> mutate(y_plot = pmin(neg_log10p, ymax, na.rm = TRUE))

label_genes <- union(
  vstat |> filter(sig == "p<0.05 & |log2FC|>1") |> pull(gene) |> unique(),
  vstat |> group_by(comparison) |>
    slice_max(abs(log2fc), n = 1, with_ties = FALSE) |> pull(gene) |> unique()
)

p_volcano <- ggplot(vstat, aes(log2fc, y_plot, colour = sig, fill = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey30", linewidth = 0.45) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "grey30", linewidth = 0.45) +
  geom_point(size = 3, alpha = 0.85, stroke = 0.25, shape = 21, colour = "white") +
  ggrepel::geom_text_repel(
    data = vstat |> filter(gene %in% label_genes),
    aes(label = gene), size = 3.2, fontface = "bold.italic",
    colour = "grey10", box.padding = 0.45, point.padding = 0.30,
    max.overlaps = 20, show.legend = FALSE,
    segment.colour = "grey60", segment.size = 0.3) +
  annotate("text", x = Inf, y = -log10(0.05) + 0.05 * ymax,
           label = "p = 0.05", hjust = 1.05, vjust = 0,
           size = 2.8, colour = "grey40", fontface = "italic") +
  scale_fill_manual(values   = vcols, name = "Significance") +
  scale_colour_manual(values = vcols, name = "Significance") +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.06))) +
  facet_wrap(~ comparison, ncol = length(groups_treatment),
             labeller = labeller(comparison = setNames(
               paste0(groups_treatment, " vs ", group_m1),
               paste0(groups_treatment, "_vs_", group_m1)))) +
  labs(title    = paste0("Volcano Plots — Reference: ", group_m1),
       subtitle = paste0("log\u2082FC = mean(\u0394Cq_LPS) \u2212 mean(\u0394Cq_treatment)",
                         " | Welch t-test, BH corrected"),
       x = expression(log[2] ~ "Fold Change"),
       y = expression(-log[10] ~ "(p-value)")) +
  theme_prism(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle    = element_text(size = 8, hjust = 0.5, colour = "grey40"),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        panel.spacing    = unit(0.6, "cm"),
        legend.position  = "bottom")

save_gg(p_volcano, file.path(multi_dir, "PlotC_Volcano"),
        w = 5.5 * length(groups_treatment), h = 6.5)
save_csv(vstat |>
           select(comparison, treatment, gene, log2fc, p_value, p_adj, sig) |>
           arrange(comparison, p_value),
         file.path(multi_dir, "Volcano_Stats"))


# --- 6f: Plot D — Sample-level clustergram -----------------------------------
message("[M4] Clustergram")

s_mat <- -1 * (wide_dcq |>
                 mutate(rl = paste0(Condition, "_", Sample)) |>
                 column_to_rownames("rl") |>
                 select(all_of(gene_labels)) |>
                 as.matrix())

zv <- names(which(apply(s_mat, 2, var, na.rm = TRUE) < 1e-10))
if (length(zv) > 0) s_mat <- s_mat[, !colnames(s_mat) %in% zv, drop = FALSE]
s_mat_z <- scale(s_mat)

pearson_dist <- function(m) {
  cm <- cor(t(m), method = "pearson", use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  as.dist(1 - cm)
}

s_annot <- data.frame(
  Condition = factor(as.character(wide_dcq$Condition),
                     levels = condition_levels[condition_levels %in%
                                                 as.character(wide_dcq$Condition)]),
  row.names = paste0(wide_dcq$Condition, "_", wide_dcq$Sample))

s_cg <- pheatmap::pheatmap(
  mat                      = s_mat_z,
  clustering_distance_rows = pearson_dist(s_mat_z),
  clustering_distance_cols = pearson_dist(t(s_mat_z)),
  clustering_method        = "complete",
  scale                    = "none",
  color      = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
  breaks     = seq(-2.5, 2.5, length.out = 101),
  annotation_row    = s_annot,
  annotation_colors = list(Condition = cb_palette[levels(s_annot$Condition)]),
  border_color = "white", cellwidth = 46, cellheight = 14,
  fontsize = 10, fontsize_row = 8, fontsize_col = 11, angle_col = 45,
  main = "Sample-Level Clustergram\nZ-score (\u2212\u0394Cq) | Pearson, complete",
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", " 0", "+1", "+2"),
  silent = TRUE)

cg_h <- max(5.5, 0.22 * nrow(s_mat_z) + 2.5)
for (ext in c("png", "pdf")) {
  path <- file.path(multi_dir, paste0("PlotD_SampleClustergram.", ext))
  if (ext == "png") {
    png(path, width = 8.5, height = cg_h, units = "in", res = 300)
  } else {
    pdf(path, width = 8.5, height = cg_h)
  }
  grid::grid.newpage()
  grid::grid.draw(s_cg$gtable)
  dev.off()
}
save_csv(data.frame(Order  = seq_along(s_cg$tree_row$order),
                    Sample = rownames(s_mat_z)[s_cg$tree_row$order]),
         file.path(multi_dir, "Clustergram_SampleOrder"))
message("  Saved: PlotD_SampleClustergram")


# --- 6g: Plot E — Gene-gene correlogram --------------------------------------
message("[M5] Correlogram")

expr_mat <- -1 * (wide_dcq |> select(all_of(gene_labels)) |> as.matrix())
n_samp   <- nrow(expr_mat)

if (n_samp >= 3) {
  cor_r <- cor(expr_mat, method = "pearson", use = "pairwise.complete.obs")
  cor_p <- tryCatch(ggcorrplot::cor_pmat(expr_mat, method = "pearson"),
                    error = function(e) {
                      matrix(NA_real_, nrow = ncol(expr_mat), ncol = ncol(expr_mat),
                             dimnames = list(gene_labels, gene_labels))
                    })
  
  p_corr <- ggcorrplot::ggcorrplot(
    corr = cor_r, p.mat = cor_p, type = "lower", method = "square",
    lab = TRUE, lab_size = 4.5, sig.level = 0.05,
    insig = "pch", pch = 4, pch.col = "grey20", pch.cex = 6,
    tl.cex = 11, tl.col = "grey10", tl.srt = 45, outline.col = "white",
    ggtheme     = ggplot2::theme_minimal(base_size = 11),
    colors      = c("#2166AC", "white", "#D6604D"),
    title       = paste0("Gene-Gene Pearson Correlation (n = ", n_samp, ")"),
    legend.title = "Pearson r") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text       = ggplot2::element_text(size = 11, face = "bold.italic"),
      legend.key.size = ggplot2::unit(0.55, "cm"),
      plot.margin     = ggplot2::margin(10, 15, 10, 10)) +
    ggplot2::labs(caption = "\u00D7 = p > 0.05 (unadjusted Pearson t-test)")
  
  save_gg(p_corr, file.path(multi_dir, "PlotE_Correlogram"), w = 5.5, h = 5.2)
  save_csv(as.data.frame(cor_r) |> tibble::rownames_to_column("Gene"),
           file.path(multi_dir, "Correlation_r"))
  save_csv(as.data.frame(cor_p) |> tibble::rownames_to_column("Gene"),
           file.path(multi_dir, "Correlation_p"))
  message("  Saved: PlotE_Correlogram")
}


# --- 6h: Multivariate CSV exports --------------------------------------------
message("[M6] CSV exports")

write_out <- function(df, name) {
  save_csv(as.data.frame(df), file.path(multi_dir, name))
  message("  ", name)
}
write_out(pca_scores |> select(Condition, Sample, PC1, PC2), "PCA_Scores")
write_out(as.data.frame(pca_res$rotation) |> rownames_to_column("Gene"), "PCA_Loadings")
write_out(as.data.frame(summary(pca_res)$importance) |>
            rownames_to_column("Metric"), "PCA_Variance")
write_out(as.data.frame(log2fc_for_threshold) |>
            rownames_to_column("Gene"), "Heatmap_Log2FC")
if (do_clust_g) {
  write_out(data.frame(Gene = gene_order, Order = seq_along(gene_order)),
            "Heatmap_GeneClusterOrder")
}


# =============================================================================
# SECTION 7: FINAL SUMMARY
# =============================================================================

message("\n", strrep("=", 60))
message("  UNIFIED PIPELINE COMPLETE")
message(strrep("=", 60))
message("  Genes: ", paste(gene_labels, collapse = " | "))
message("  QC:    HK_Cq <= ", HK_CQ_THRESHOLD,
        " | IQR multiplier = ", IQR_MULTIPLIER)
message("  PCA:   PC1=", pct_var["PC1"], "% | PC2=", pct_var["PC2"], "%")
if (do_clust_g) {
  message("  Gene cluster order: ", paste(gene_order, collapse = " → "))
}
message("\n  Univariate outputs : ", uni_dir)
message("  Multivariate outputs: ", multi_dir)
message(strrep("=", 60))

# =============================================================================
# END OF SCRIPT
# =============================================================================