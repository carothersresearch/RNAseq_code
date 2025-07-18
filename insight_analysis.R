library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(knitr)
library(pheatmap)
library(tibble)

# ==== Directories ====
enrich_dir <- "enrichment_analysis"
output_dir <- "insight_analysis"
dir.create(output_dir, showWarnings = FALSE)

# ==== Helpers ====
parse_list_string <- function(s) {
  if (is.na(s) || s == "[]" || s == "") return(character(0))
  strsplit(gsub("^\\[|'|\\]|'$", "", s), "', '")[[1]]
}

to_html_table <- function(df, n = 10) {
  if (is.null(df) || nrow(df) == 0) return("<p><i>None found.</i></p>")
  df %>%
    head(n) %>%
    kable(
      format = "html",
      table.attr = "border='1' cellpadding='5' cellspacing='0' style='table-layout:auto; width:100%; border-collapse:collapse; text-align:left;'"
    ) %>%
    as.character()
}

remove_locus_columns <- function(df) {
  if (is.null(df)) return(NULL)
  patterns <- c("geneID", "core_enrichment")
  locus_cols <- names(df)[sapply(names(df), function(x) any(grepl(paste(patterns, collapse="|"), x, ignore.case=TRUE)))]
  if (length(locus_cols) > 0) df <- df %>% select(-all_of(locus_cols))
  df
}

# Helper to convert list columns to comma-separated strings for HTML tables
format_lists <- function(df, cols) {
  for (col in cols) {
    if (col %in% names(df)) {
      df[[col]] <- purrr::map_chr(df[[col]], function(x) {
        if (is.null(x) || length(x) == 0) {
          return("")
        }
        if (is.character(x)) {
          return(paste(x, collapse = ", "))
        }
        return(as.character(paste(x, collapse = ", ")))
      })
    }
  }
  return(df)
}

# === Heatmap plotting helper ===
plot_heatmap_top_genes <- function(genes_df, expr_mat, title_prefix, suffix, output_dir) {
  top_genes <- head(genes_df$Locus_Tag, 10)
  if(length(top_genes) == 0) return(NULL)  # no genes
  
  # Filter expr matrix rows for these genes
  expr_sub <- expr_mat[top_genes, , drop = FALSE]
  # Get protein names for labels
  top_proteins <- genes_df %>% filter(Locus_Tag %in% top_genes) %>% pull(`Protein Name`)
  
  # Replace rownames with Protein Names for nicer labels
  rownames(expr_sub) <- top_proteins
  
  # Scale rows (genes) to z-scores for heatmap visualization
  expr_scaled <- t(scale(t(expr_sub)))
  
  # Output filename
  out_svg <- file.path(output_dir, paste0(title_prefix, "_heatmap_condition_", suffix, ".svg"))
  
  # Save as SVG
  svg(filename = out_svg, width = 12, height = 4)
  
  # Generate heatmap
  pheatmap(
    expr_scaled,
    main = paste(title_prefix, "Top 10 Genes - Condition", suffix),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    fontsize_row = 12,
    fontsize_col = 12,
    angle_col = 45,
    legend = TRUE,
    annotation_legend = TRUE,
    heatmap_legend_param = list(
      title = "Scaled\nExpression")
  )
  
  # Close device
  dev.off()
  
  return(out_svg)
}

# === Load TPM expression matrix ===
# Path to your TPM counts CSV (rows=Locus_Tag, cols=samples)
expr_path <- "deg_analysis/TPM_matrix.csv"
expr_mat <- read_csv(expr_path) %>% 
  column_to_rownames(var = "Locus_Tag") %>% 
  as.matrix()

message("Loaded TPM expression matrix with dimensions: ", paste(dim(expr_mat), collapse = " x "))

# ==== Read Data ====
deg_df <- read_csv("deg_results_all.csv")

annot <- read_csv("chatgpt_annotations.csv") %>%
  mutate(
    KO_IDs_list = map(KO_IDs, parse_list_string),
    GO_IDs_list = map(GO_IDs, parse_list_string),
    GO_Names_list = map(GO_Names, parse_list_string),
    Tags_list = map(Tags, parse_list_string),
    KEGG_Pathways_list = map(
      KEGG_Pathway_IDs ,
      ~ {
        pathways <- parse_list_string(.x)       # returns character vector
        str_replace(pathways, "^ko", "map")     # replace any "ko..." with "map..."
      }
    )
  )

message("deg_df columns: ", paste(names(deg_df), collapse=", "))
message("annot columns: ", paste(names(annot), collapse=", "))

# Ensure deg_df has Locus_Tag column
if (!"Locus_Tag" %in% names(deg_df)) {
  if ("gene" %in% names(deg_df)) {
    deg_df <- deg_df %>% rename(Locus_Tag = gene)
    message("Renamed 'gene' column in deg_df to 'Locus_Tag'")
  } else {
    stop("deg_results_all.csv must have a 'Locus_Tag' or 'gene' column!")
  }
}

long_deg <- deg_df %>%
  pivot_longer(
    cols = matches("^log2FC_|^padj_|^sig_"),
    names_to = c(".value", "condition"),
    names_pattern = "(log2FC|padj|sig)_(.*)"
  ) %>%
  left_join(annot, by = "Locus_Tag") %>%
  filter(sig == "yes") %>%
  mutate(condition = str_remove(condition, "\\.csv$"))

message("long_deg sample after removing .csv:")
print(head(long_deg[, c("Locus_Tag","condition","log2FC","sig")]))

# ==== Count DEGs per condition & metabolism pathway ====
library(KEGGREST)
kegg_hierarchy <- keggList("pathway") %>% enframe(name = "pathway_id", value = "pathway_name")

metabolism_ids <- kegg_hierarchy %>%
  filter(str_detect(tolower(pathway_name), "metabolism")) %>%
  pull(pathway_id)

long_deg_meta <- long_deg %>%
  filter(KEGG_Pathways_list %in% metabolism_ids)

# Expand list-column into long format
long_deg_unnested <- long_deg %>%
  unnest(KEGG_Pathways_list)

# Debug
message("Dimensions after unnesting: ", paste(dim(long_deg_unnested), collapse = " x "))
head(long_deg_unnested[, c("Locus_Tag", "KEGG_Pathways_list")])

# Now filter for metabolism pathways
long_deg_meta <- long_deg_unnested %>%
  filter(KEGG_Pathways_list %in% metabolism_ids)

message("Number of DEGs for metabolism pathways: ", nrow(long_deg_meta))

count_df <- long_deg_meta %>%
  group_by(condition, KEGG_Pathways_list) %>%
  summarise(n_deg = n_distinct(Locus_Tag), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n_deg, values_fill = list(n_deg = 0)) %>%
  column_to_rownames(var = "KEGG_Pathways_list")

# Save CSV
readr::write_csv(as.data.frame(count_df), file.path(output_dir, "metabolism_deg_counts_per_condition.csv"))

# Save heatmap
svg(file.path(output_dir, "metabolism_deg_counts_heatmap.svg"), width = 10, height = 12)
pheatmap(
  count_df,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  main = "DEGs per Metabolism Pathway",
  cluster_rows = FALSE, cluster_cols = FALSE,
  fontsize_row = 8, fontsize_col = 10
)
dev.off()

# ==== Read metadata and parse experimental design ====
metadata_df <- read_csv("metadata.csv")

conditions <- grep("^condition_", names(metadata_df), value = TRUE)
condition_info <- lapply(conditions, function(cond_col) {
  data <- metadata_df %>%
    filter(!is.na(.data[[cond_col]]) & .data[[cond_col]] != "") %>%
    group_by(label = .data[[cond_col]]) %>%
    summarise(
      reps = paste(unique(`biological replicates`), collapse = ", "),
      .groups = "drop"
    )
  list(condition = cond_col, groups = data)
})
names(condition_info) <- conditions

condition_cols <- grep("^condition_", names(metadata_df), value = TRUE)

html_lines <- c(
  "<html>",
  "<head>",
  "<title>Enrichment & DEG Insights</title>",
  "<style>",
  "body { font-family: Arial, sans-serif; padding: 1em; background: #fafafa; }",
  "h2 { color: #1a5276; border-bottom: 3px solid #1a5276; padding-bottom: 0.2em; }",
  "h3 { color: #117a65; margin-top: 1.5em; }",
  "table { border-collapse: collapse; width: 100%; table-layout: fixed; word-wrap: break-word; }",
  "th, td { border: 1px solid #ccc; padding: 0.5em; text-align: left; font-size: 0.9em; }",
  "th { background-color: #e8f8f5; color: #0e6655; }",
  "hr { border: none; border-top: 2px solid #d5d8dc; margin: 2em 0; }",
  ".section-wrapper { padding: 1em; background: white; border: 1px solid #ddd; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); margin-bottom: 2em; }",
  ".metric-summary { background-color: #f0f9ff; padding: 0.5em; border-radius: 4px; color: #34495e; font-size: 0.9em; }",
  "</style>",
  "</head>",
  "<body>",
  "<h1 style='color:#2e86c1;'>Enrichment & DEG Insights</h1>"
)

for (cond_full in condition_cols) {
  message(sprintf("Processing condition: %s", cond_full))
  suffix <- sub("^condition_", "", cond_full)
  
  design_summary <- if (cond_full %in% names(condition_info)) {
    design <- condition_info[[cond_full]]$groups
    paste(
      apply(design, 1, function(row) sprintf("%s: %s", row[["label"]], row[["reps"]])),
      collapse = "<br/>"
    )
  } else {
    "No design information available."
  }
  
  sub_df <- long_deg %>% filter(condition == cond_full)
  message(sprintf("DEGs for condition %s: %d rows", cond_full, nrow(sub_df)))
  
  up_genes <- sub_df %>% filter(log2FC > 0) %>% arrange(-log2FC)
  down_genes <- sub_df %>% filter(log2FC < 0) %>% arrange(log2FC)
  
  # Transcription factors and unknowns
  tf_genes <- sub_df %>%
    filter(str_detect(tolower(`Functional_Group_Name`), "transcription"))
  
  unknown_genes <- sub_df %>%
    filter(
      is.na(`Functional_Group_Name`) |
        `Functional_Group_Name` == "" |
        str_detect(tolower(`Functional_Group_Name`), "unknown|hypothetical")
    )
  
  # Get neighboring locus for unknowns
  annot_ordered <- annot %>% arrange(Locus_Tag)
  locus_vec <- annot_ordered$Locus_Tag
  ko_vec <- sapply(annot_ordered$Tags, function(x) if (length(x) == 0) NA_character_ else x[[1]])
  
  unknown_genes <- unknown_genes %>%
    rowwise() %>%
    mutate(
      idx = which(locus_vec == Locus_Tag),
      upstream_Locus_Tag = ifelse(length(idx) > 0 && idx > 1, locus_vec[idx - 1], NA_character_),
      upstream_tags = ifelse(length(idx) > 0 && idx > 1, ko_vec[idx - 1], NA_character_),
      downstream_Locus_Tag = ifelse(length(idx) > 0 && idx < length(locus_vec), locus_vec[idx + 1], NA_character_),
      downstream_tags = ifelse(length(idx) > 0 && idx < length(locus_vec), ko_vec[idx + 1], NA_character_)
    ) %>%
    ungroup() %>%
    select(-idx)
  
  # Format list columns
  list_cols_deg <- c("KO_IDs_list", "KEGG_Pathways_list", "GO_Names_list", "Tags_list")
  up_genes <- format_lists(up_genes, list_cols_deg)
  down_genes <- format_lists(down_genes, list_cols_deg)
  tf_genes <- format_lists(tf_genes, list_cols_deg)
  unknown_genes <- format_lists(unknown_genes, c("Tags_list"))
  
  # Generate and save heatmaps for top 10 up and down genes
  heatmap_up_file <- plot_heatmap_top_genes(up_genes, expr_mat, "upregulated", suffix, output_dir)
  heatmap_down_file <- plot_heatmap_top_genes(down_genes, expr_mat, "downregulated", suffix, output_dir)
  heatmap_tf_file <- plot_heatmap_top_genes(tf_genes, expr_mat, "tf", suffix, output_dir)
  heatmap_unknown_file <- plot_heatmap_top_genes(unknown_genes, expr_mat, "unknown", suffix, output_dir)
  
  
  # Save CSVs — all rows per table
  write_csv(
    up_genes %>% select(
      Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list,
      Tags_list, Subcategories, Confidence_Score, Notes
    ),
    file.path(output_dir, sprintf("all_upregulated_%s.csv", suffix))
  )
  write_csv(
    down_genes %>% select(
      Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list,
      Tags_list, Subcategories, Confidence_Score, Notes
    ),
    file.path(output_dir, sprintf("all_downregulated_%s.csv", suffix))
  )
  write_csv(
    tf_genes %>% select(
      Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list,
      Tags_list, Subcategories, Confidence_Score, Notes
    ),
    file.path(output_dir, sprintf("all_tf_%s.csv", suffix))
  )
  write_csv(
    unknown_genes %>% select(
      Locus_Tag, `Protein Name`, log2FC, Tags_list,
      Subcategories, Confidence_Score, Notes,
      upstream_Locus_Tag, upstream_tags,
      downstream_Locus_Tag, downstream_tags
    ),
    file.path(output_dir, sprintf("all_unknown_%s.csv", suffix))
  )
  
  # Read enrichments and plots
  go_file <- file.path(enrich_dir, sprintf("GO_enrichment_condition_%s.csv", suffix))
  kegg_file <- file.path(enrich_dir, sprintf("KEGG_enrichment_condition_%s.csv", suffix))
  gsea_file <- file.path(enrich_dir, sprintf("GSEA_condition_%s.csv", suffix))
  
  go_df <- if (file.exists(go_file)) read_csv(go_file, show_col_types = FALSE) else NULL
  kegg_df <- if (file.exists(kegg_file)) read_csv(kegg_file, show_col_types = FALSE) else NULL
  gsea_df <- if (file.exists(gsea_file)) read_csv(gsea_file, show_col_types = FALSE) else NULL
  
  go_file_disk <- file.path(enrich_dir, sprintf("GO_barplot_condition_%s.svg", suffix))
  kegg_file_disk <- file.path(enrich_dir, sprintf("KEGG_dotplot_condition_%s.svg", suffix))
  gsea_file_disk <- file.path(enrich_dir, sprintf("GSEA_plot_condition_%s.svg", suffix))
  
  go_img_path <- file.path("..", "enrichment_analysis", basename(go_file_disk))
  kegg_img_path <- file.path("..", "enrichment_analysis", basename(kegg_file_disk))
  gsea_img_path <- file.path("..", "enrichment_analysis", basename(gsea_file_disk))
  
  # Append HTML report section
  html_lines <- c(
    html_lines,
    "<details class='section-wrapper'>",
    sprintf("<summary><h2>Condition %s</h2></summary>", suffix),
    "<div class='section-content'>",
    sprintf("<p class='metric-summary'><i>Experimental design:</i><br/>%s</p>", design_summary),
    sprintf(
      "<div class='metric-summary'>Upregulated: %d | Downregulated: %d | Total DEGs: %d</div>",
      nrow(up_genes), nrow(down_genes), nrow(sub_df)
    ),
    "<h3>Top Upregulated DEGs</h3>",
    to_html_table(
      head(up_genes, 10) %>% select(
        Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list, Tags_list,
        Subcategories, Confidence_Score, Notes
      ), 10
    ),
    if (!is.null(heatmap_up_file)) sprintf("<img src='%s' alt='Heatmap Upregulated Genes - Condition %s' style='max-width:100%%; height:auto;'>", file.path("..", output_dir, basename(heatmap_up_file)), suffix) else "",
    "<h3>Top Downregulated DEGs</h3>",
    to_html_table(
      head(down_genes, 10) %>% select(
        Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list, Tags_list,
        Subcategories, Confidence_Score, Notes
      ), 10
    ),
    if (!is.null(heatmap_down_file)) sprintf("<img src='%s' alt='Heatmap Downregulated Genes - Condition %s' style='max-width:100%%; height:auto;'>", file.path("..", output_dir, basename(heatmap_down_file)), suffix) else "",
    "<h3>Transcription Factors</h3>",
    to_html_table(
      head(tf_genes, 10) %>% select(
        Locus_Tag, `Protein Name`, KO_IDs_list, log2FC, KEGG_Pathways_list, Tags_list,
        Subcategories, Confidence_Score, Notes
      ), 10
    ),
    if (!is.null(heatmap_tf_file)) sprintf("<img src='%s' alt='Heatmap TF Genes - Condition %s' style='max-width:100%%; height:auto;'>", file.path("..", output_dir, basename(heatmap_tf_file)), suffix) else "",
    "<h3>Genes of Unknown Function</h3>",
    to_html_table(
      head(unknown_genes, 10) %>% select(
        Locus_Tag, `Protein Name`, log2FC, Tags_list,
        Subcategories, Confidence_Score, Notes,
        upstream_Locus_Tag, upstream_tags,
        downstream_Locus_Tag, downstream_tags
      ), 10
    ),
    if (!is.null(heatmap_unknown_file)) sprintf("<img src='%s' alt='Heatmap Unknown Genes - Condition %s' style='max-width:100%%; height:auto;'>", file.path("..", output_dir, basename(heatmap_unknown_file)), suffix) else "",
    "<h3>GO Enrichment</h3>",
    to_html_table(remove_locus_columns(go_df), 10),
    if (file.exists(go_file_disk)) sprintf("<img src='%s' alt='GO plot'>", go_img_path) else "<p><i>No GO plot available.</i></p>",
    "<h3>KEGG Pathway Enrichment</h3>",
    to_html_table(remove_locus_columns(kegg_df), 10),
    if (file.exists(kegg_file_disk)) sprintf("<img src='%s' alt='KEGG plot'>", kegg_img_path) else "<p><i>No KEGG plot available.</i></p>",
    "<h3>Gene Set Enrichment Analysis (GSEA)</h3>",
    to_html_table(remove_locus_columns(gsea_df), 10),
    if (file.exists(gsea_file_disk)) sprintf("<img src='%s' alt='GSEA plot'>", gsea_img_path) else "<p><i>No GSEA plot available.</i></p>",
    "</div>",
    "</details>"
  )
}


html_lines <- c(html_lines, "</body></html>")
writeLines(html_lines, file.path(output_dir, "insight_summary_report.html"))
message("✅ Report saved successfully to ", file.path(output_dir, "insight_summary_report.html"))
