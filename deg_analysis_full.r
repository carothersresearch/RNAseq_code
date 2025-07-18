# run_deseq2.R
library(DESeq2)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(EnhancedVolcano)
library(tibble)
library(stringi)

dir.create("deg_analysis")

# Load sample meta
meta <- read.csv("/data/metadata.csv", stringsAsFactors = FALSE)

# Loop through columns, convert those with "condition" in their name to factor
for (colname in colnames(meta)) {
  if (grepl("condition", colname, ignore.case = TRUE)) {
    meta[[colname]] <- factor(meta[[colname]])
    message("Converted ", colname, " to factor.")
  }
}

for (file_name in meta$fileName) {
  file_path <- file.path("/data", file_name)
  
  # Detect encoding using stringi
  enc_info <- stri_enc_detect(readBin(file_path, what = "raw", n = 10000))[[1]]
  detected_enc <- enc_info$Encoding[1]
  
  if (detected_enc %in% c("UTF-16LE", "UTF-16")) {
    # Open and read file in UTF-16LE
    con <- file(file_path, open = "r", encoding = "UTF-16LE")
    text_utf8 <- readLines(con)
    close(con)
    
    # Overwrite file as UTF-8
    writeLines(text_utf8, file_path, useBytes = TRUE)
    
    message("Converted to UTF-8: ", file_path)
  } else {
    message("Skipped (not UTF-16): ", file_path)
  }
}

# STEP 2: Read and merge HTSeq count files
count_list <- lapply(meta$fileName, function(f) {
  read.table(f, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene", f))
})

# Merge by gene
merged_counts <- Reduce(function(x, y) merge(x, y, by = "gene"), count_list)
rownames(merged_counts) <- merged_counts$gene
# Vector of unwanted feature names from HTSeq summary
htseq_summary_features <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
# After merging all count files but BEFORE running DESeq2, remove these rows:
merged_counts <- merged_counts[!merged_counts$gene %in% htseq_summary_features, ]
count_matrix <- merged_counts[, -1]  # Remove gene column
colnames(count_matrix) <- meta$sampleName

# Combine technical replicates by summing columns
combined_counts <- sapply(unique(meta$`technical.replicates`), function(g) {
  cols <- meta$sampleName[meta$`technical.replicates` == g]
  rowSums(count_matrix[, cols, drop = FALSE])
})

# STEP 4: Prepare final count matrix
combined_counts <- as.data.frame(combined_counts)
rownames(combined_counts) <- rownames(count_matrix)

# Get metadata for technical replicates only, one row each
tech_meta <- meta %>%
  distinct(technical.replicates, .keep_all = TRUE) %>%
  rename(sample = technical.replicates)

rownames(tech_meta) <- tech_meta$sample

# Reorder rows of tech_meta to match columns of combined_counts
tech_meta <- tech_meta[match(colnames(combined_counts), tech_meta$sample), ]

# 7. Identify all condition columns dynamically
condition_cols <- grep("^condition_", colnames(tech_meta), value = TRUE)

# 8. Run DESeq2 for each condition column
results_list <- list()

# Clean all condition columns: convert "" or blanks to NA
tech_meta[condition_cols] <- lapply(tech_meta[condition_cols], function(x) {
  x <- trimws(x)                  # remove whitespace
  x[x == ""] <- NA               # replace empty strings with NA
  factor(x)                     # convert to factor
})

# Then filter out samples with NA for each condition before running DESeq2
results_list <- list()

for(cond in condition_cols) {
  # Subset samples without NA in this condition
  keep <- !is.na(tech_meta[[cond]])
  sub_counts <- combined_counts[, keep, drop = FALSE]
  sub_meta <- tech_meta[keep, , drop = FALSE]
  
  # Check number of levels after filtering
  cond_levels <- unique(sub_meta[[cond]])
  if(length(cond_levels) != 2) {
    warning(paste("Skipping", cond, "- requires exactly 2 factor levels, found:", paste(cond_levels, collapse = ", ")))
    next
  }
  
  # Make sure factor levels are set
  sub_meta[[cond]] <- factor(sub_meta[[cond]], levels = c("control", "treated"))
  
  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = sub_counts,
                                colData = sub_meta,
                                design = as.formula(paste("~", cond)))
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c(cond, "treated", "control"))
  
  results_list[[cond]] <- res
  
  # Convert results to a data.frame and add gene names as a column
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Add a significance column
  res_df$significant <- with(res_df, 
                             ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
                                    "yes", "no"))
  
  # Write to CSV
  write.csv(res_df, file = file.path("deg_analysis", paste0("results_", cond, ".csv")),
            row.names = FALSE)
  
  # Calculate summary stats
  res_df <- as.data.frame(res) %>% rownames_to_column("gene")
  n_up <- sum(res_df$log2FoldChange > 1 & res_df$padj < 0.05, na.rm = TRUE)
  n_down <- sum(res_df$log2FoldChange < -1 & res_df$padj < 0.05, na.rm = TRUE)
  summary_text <- paste0("Upregulated: ", n_up, "\nDownregulated: ", n_down)
  
  # Open SVG device
  svg(filename = paste0("deg_analysis/volcano_", cond, ".svg"), width = 8, height = 6)
  
  # Plot volcano with annotation
  volcano_plot <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    title = paste("Volcano Plot -", cond),
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0
  )
  
  # Create the annotation layer separately
  annotation_layer <- annotate(
    "text",
    x = mean(range(res$log2FoldChange, na.rm=TRUE)),
    y = mean(range(-log10(res$pvalue), na.rm=TRUE)),
    label = paste0("Upregulated: ", n_up, "\nDownregulated: ", n_down),
    hjust = 0, vjust = 1,
    size = 5,
    color = "black",
    fontface = "bold"
  )
  
  print(volcano_plot + annotation_layer)
  
  # Close SVG device and save file
  dev.off()
  
  
  ### Select significant DEGs
  res <- results(dds)
  res <- lfcShrink(dds, coef=2, type="apeglm")
  sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))
  
  sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  # Split into up and down regulated
  up_genes <- sig_res[sig_res$log2FoldChange > 0, ]
  down_genes <- sig_res[sig_res$log2FoldChange < 0, ]
  
  # Order by adjusted p-value ascending (more significant first)
  up_genes <- up_genes[order(up_genes$padj), ]
  down_genes <- down_genes[order(down_genes$padj), ]
  
  # Take top 50 from each
  top_up <- head(rownames(up_genes), 100)
  top_down <- head(rownames(down_genes), 100)
  
  # Combine top genes
  top_genes <- c(top_up, top_down)
  
  # Get normalized counts for top genes
  norm_counts <- counts(dds, normalized=TRUE)
  norm_counts_subset <- norm_counts[top_genes, ]
  
  # Scale by gene (z-score)
  scaled_counts <- t(scale(t(norm_counts_subset)))
  
  # Sample annotation
  sample_annotation <- as.data.frame(colData(dds)[, "biological.replicates", drop=FALSE])
  
  # Plot heatmap and save as SVG
  svg(filename = paste0("deg_analysis/heatmap_", cond, ".svg"), width = 16, height = ceiling(nrow(scaled_counts) * 0.2))
  pheatmap(
    scaled_counts,
    annotation_col = sample_annotation,
    show_rownames = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
    cellwidth = 12,
    cellheight = 10,
    main = "Top 100 Up and Downregulated DEGs"
  )
  dev.off()
  
}


### Raw reads and TPM
# 1. Load raw counts matrix (genes x samples)
write.csv(combined_counts, "deg_analysis/counts_matrix.csv")

# 2. Load gene lengths
gene_lengths <- read.delim("gene_lengths.tsv", stringsAsFactors = FALSE)

# 3. Subset genes present in both counts and gene lengths
common_genes <- intersect(rownames(combined_counts), gene_lengths$locus_tag)
raw_counts_sub <- combined_counts[common_genes, ]
gene_lengths_sub <- gene_lengths[match(common_genes, gene_lengths$locus_tag), ]

# 4. Calculate TPM
length_kb <- gene_lengths_sub$length / 1000
rpk <- sweep(raw_counts_sub, 1, length_kb, FUN = "/")
scaling_factor <- colSums(rpk)
tpm <- sweep(rpk, 2, scaling_factor, FUN = "/") * 1e6

write.csv(tpm, "deg_analysis/TPM_matrix.csv")



### Sample-to-Sample Distance
dds <- DESeqDataSetFromMatrix(countData = combined_counts,
                              colData = tech_meta,
                              design = ~ 1) 

# Transform counts with variance stabilizing transformation (VST)
vsd <- vst(dds, blind = TRUE)

# Extract the transformed expression matrix
mat <- assay(vsd)

# Compute sample-to-sample distances
sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)

# Optional: add sample names
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- colnames(mat)

# Add condition labels
sample_annotation <- as.data.frame(colData(dds)[, "biological.replicates", drop=FALSE])

svg(filename = "deg_analysis/sample_distance_heatmap.svg", width = 8, height = 6)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         annotation_col = sample_annotation,
         main = "Sample-to-Sample Distance Heatmap")

dev.off()


library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# Define directory containing the CSVs
csv_dir <- "deg_analysis"

# List all result CSV files
csv_files <- list.files(csv_dir, pattern = "^results_.*\\.csv$", full.names = TRUE)

# Read and tag each file with the condition name
all_dfs <- map(csv_files, function(file) {
  df <- read_csv(file)
  # Extract condition name from filename (e.g., results_condition_3.csv → condition_3)
  cond <- str_remove(basename(file), "^results_|\\.csv$")
  df$condition <- cond
  df
})

# Combine into long format
long_df <- bind_rows(all_dfs)

# Optional: Ensure gene column exists and clean
long_df <- long_df %>% 
  rename(gene = gene) %>%
  mutate(condition = factor(condition))

# Make wide format for log2FoldChange
lfc_wide <- long_df %>%
  select(gene, condition, log2FoldChange) %>%
  pivot_wider(names_from = condition, values_from = log2FoldChange, names_prefix = "log2FC_")

# Make wide format for padj
padj_wide <- long_df %>%
  select(gene, condition, padj) %>%
  pivot_wider(names_from = condition, values_from = padj, names_prefix = "padj_")

# Make wide format for significance
sig_wide <- long_df %>%
  select(gene, condition, significant) %>%
  pivot_wider(names_from = condition, values_from = significant, names_prefix = "sig_")

# Merge all into one wide table
wide_df <- reduce(list(lfc_wide, padj_wide, sig_wide), full_join, by = "gene")

# Preview the merged wide table
print(head(wide_df))

# Save to CSV
write_csv(wide_df, file = "deg_results_all.csv")




## Use ggplot2::theme_bw()
library('ggplot2')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("regionReport")

library(regionReport)
## Create the HTML report
report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
                       intgroup = c('condition', 'type'), outdir = 'deg_analysis',
                       output = 'index', theme = theme_bw())








