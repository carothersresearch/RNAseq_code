# Load required libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(enrichplot)
library(knitr)
library(GO.db)

dir.create("enrichment_analysis")

# --- User inputs: file paths ---
deg_file <- "deg_results_all.csv"
annot_file <- "chatgpt_annotations.csv"
num_conditions <- 6  # number of conditions in DEG file

# --- Read DEG data ---
deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)

# --- Read annotation data ---
annot <- read.csv(annot_file, stringsAsFactors = FALSE)

# --- Helper function to parse KO_IDs and GO_IDs columns ---
parse_list_string <- function(s) {
  if (is.na(s) || s == "[]") return(character(0))
  s_clean <- gsub("^\\['|\\']$", "", s)
  if (nchar(s_clean) == 0) return(character(0))
  strsplit(s_clean, "', '")[[1]]
}

annot$KO_IDs_list <- lapply(annot$KO_IDs, parse_list_string)
annot$GO_IDs_list <- lapply(annot$GO_IDs, parse_list_string)

# --- Build TERM2GENE data frames ---

# KO TERM2GENE
ko_term2gene <- annot %>%
  dplyr::select(Locus_Tag, KO_IDs_list) %>%
  tidyr::unnest(cols = c(KO_IDs_list)) %>%
  dplyr::filter(KO_IDs_list != "") %>%
  dplyr::rename(gene = Locus_Tag, KO = KO_IDs_list)

# GO TERM2GENE
go_term2gene <- annot %>%
  dplyr::select(Locus_Tag, GO_IDs_list) %>%
  tidyr::unnest(GO_IDs_list) %>%
  dplyr::rename(gene = Locus_Tag, GO = GO_IDs_list) %>%
  dplyr::filter(!is.na(GO), GO != "") %>%
  dplyr::select(GO, gene) %>%
  dplyr::rename(term = GO)

unique_terms <- unique(go_term2gene$term)
go_term2name <- data.frame(
  term = unique_terms,
  name = sapply(unique_terms, function(x) {
    term_obj <- GO.db::GOTERM[[x]]
    if (is.null(term_obj)) {
      return(NA_character_)  # fallback if not found
    } else {
      return(term_obj@Term)
    }
  }),
  stringsAsFactors = FALSE
)

go_term2name <- go_term2name %>% dplyr::filter(!is.na(name))
go_term2gene <- go_term2gene %>% dplyr::filter(term %in% go_term2name$term)

# --- Function to run enrichment per condition ---

run_enrichment_for_condition <- function(condition_num, deg_data, ko_term2gene, go_term2gene, mode = c("ORA", "GSEA", "both")) {
  
  mode <- match.arg(mode)
  
  sig_col <- paste0("sig_condition_", condition_num, ".csv")
  logfc_col <- paste0("log2FC_condition_", condition_num, ".csv")
  
  # Get significant genes for condition
  sig_genes <- deg_data$gene[deg_data[[sig_col]] == "yes"]
  message(sprintf("Condition %d: %d significant genes", condition_num, length(sig_genes)))
  
  # Print first 10 sig genes for debugging
  message(sprintf("Condition %d: Sample sig genes: %s", condition_num, paste(head(sig_genes, 10), collapse = ", ")))
  
  # Map sig locus tags to KO IDs for KEGG enrichment
  sig_kos <- unique(ko_term2gene$KO[ko_term2gene$gene %in% sig_genes])
  message(sprintf("Condition %d: %d mapped KO terms from sig genes", condition_num, length(sig_kos)))
  
  # KEGG enrichment
  kk <- enrichKEGG(
    gene = sig_kos,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  message(sprintf("Condition %d: KEGG enrichment results - %d terms", condition_num, ifelse(is.null(kk), 0, nrow(kk@result))))
  
  go_enrich <- NULL
  gsea_res <- NULL
  
  # GO universe for ORA
  go_universe <- unique(go_term2gene$gene)
  message(sprintf("Condition %d: GO universe size = %d", condition_num, length(go_universe)))
  
  # Print first 10 universe genes
  message(sprintf("Condition %d: Sample GO universe genes: %s", condition_num, paste(head(go_universe, 10), collapse = ", ")))
  
  # GO enrichment (over-representation)
  if (mode %in% c("ORA", "both")) {
    # Check overlap of sig_genes and TERM2GENE genes
    overlap_genes <- sum(sig_genes %in% go_universe)
    message(sprintf("Condition %d: Overlap between sig genes and GO TERM2GENE genes = %d", condition_num, overlap_genes))
    
    if (overlap_genes > 0) {
      # Run enricher with explicit universe and relaxed cutoffs for debugging
      go_enrich <- enricher(
        gene = sig_genes,
        TERM2GENE = go_term2gene,
        TERM2NAME = go_term2name,
        universe = go_universe,
        pvalueCutoff = 0.1,      # Relaxed cutoff to check results
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 500
      )
      n_go_terms <- ifelse(is.null(go_enrich) || nrow(go_enrich) == 0, 0, nrow(go_enrich))
      message(sprintf("Condition %d: GO enrichment (ORA) results - %d terms", condition_num, n_go_terms))
      
      if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
        go_enrich@result <- dplyr::left_join(
          go_enrich@result,
          go_term2name,
          by = c("ID" = "term")
        )
        go_enrich@result$Description <- go_enrich@result$name
      }
      
      # Test ORA on subset of sig genes
      if(n_go_terms == 0 && length(sig_genes) > 20) {
        message(sprintf("Condition %d: Testing ORA on first 20 sig genes...", condition_num))
        test_enrich <- enricher(
          gene = head(sig_genes, 20),
          TERM2GENE = go_term2gene,
          universe = go_universe,
          pvalueCutoff = 0.2,
          pAdjustMethod = "BH",
          minGSSize = 5,
          maxGSSize = 300
        )
        if(!is.null(test_enrich) && nrow(test_enrich) > 0) {
          message(sprintf("Condition %d: Subset ORA found %d terms", condition_num, nrow(test_enrich)))
        } else {
          message(sprintf("Condition %d: Subset ORA found no terms", condition_num))
        }
      }
      
    } else {
      message(sprintf("Condition %d: No overlap with GO TERM2GENE genes. Skipping ORA GO enrichment.", condition_num))
    }
  }
  
  # Prepare ranked gene list for GSEA (log2FC)
  ranked_df <- deg_data %>%
    dplyr::select(gene, !!sym(logfc_col)) %>%
    dplyr::filter(!is.na(!!sym(logfc_col))) %>%
    arrange(desc(!!sym(logfc_col)))
  
  ranked_vec <- ranked_df[[logfc_col]]
  names(ranked_vec) <- ranked_df$gene
  
  # GSEA with GO TERM2GENE
  if (mode %in% c("GSEA", "both")) {
    if (length(ranked_vec) > 0) {
      gsea_res <- GSEA(
        geneList = ranked_vec,
        TERM2GENE = go_term2gene,
        pAdjustMethod = "BH",
        verbose = FALSE,
        eps = 0
      )
      n_gsea_terms <- ifelse(is.null(gsea_res) || is.null(gsea_res@result), 0, nrow(gsea_res@result))
      message(sprintf("Condition %d: GSEA results - %d terms", condition_num, n_gsea_terms))
      
      if (!is.null(gsea_res) && !is.null(gsea_res@result) && nrow(gsea_res@result) > 0) {
        gsea_res@result <- dplyr::left_join(
          gsea_res@result,
          go_term2name,
          by = c("ID" = "term")
        )
        gsea_res@result$Description <- gsea_res@result$name
      }
      
    } else {
      message(sprintf("Condition %d: Ranked gene list empty. Skipping GSEA.", condition_num))
    }
  }
  
  # Save CSV results
  write.csv(as.data.frame(kk), paste0("enrichment_analysis/KEGG_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  
  if (!is.null(go_enrich)) {
    write.csv(as.data.frame(go_enrich), paste0("enrichment_analysis/GO_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  } else {
    message(sprintf("Condition %d: No GO enrichment results to save.", condition_num))
  }
  
  if (!is.null(gsea_res)) {
    write.csv(as.data.frame(gsea_res), paste0("enrichment_analysis/GSEA_condition_", condition_num, ".csv"), row.names = FALSE)
  } else {
    message(sprintf("Condition %d: No GSEA results to save.", condition_num))
  }
  
  # Save plots as SVG
  svg(paste0("enrichment_analysis/KEGG_dotplot_condition_", condition_num, ".svg"), width = 8, height = 6)
  print(dotplot(kk, showCategory = 10) + ggtitle(paste("KEGG Pathway Enrichment - Condition", condition_num)))
  dev.off()
  
  if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
    svg(paste0("enrichment_analysis/GO_barplot_condition_", condition_num, ".svg"), width = 8, height = 6)
    print(barplot(go_enrich, showCategory = 5) + ggtitle(paste("GO Enrichment - Condition", condition_num)))
    dev.off()
  } else {
    message(sprintf("Condition %d: No GO enrichment plot to save.", condition_num))
  }
  
  svg(paste0("enrichment_analysis/GSEA_plot_condition_", condition_num, ".svg"), width = 8, height = 6)
  if (!is.null(gsea_res) && !is.null(gsea_res@result) && nrow(gsea_res@result) > 0) {
    print(gseaplot2(gsea_res, geneSetID = gsea_res@result$ID[1], title = paste("GSEA Plot - Condition", condition_num)))
    dev.off()
  } else {
    message(sprintf("Condition %d: No GSEA plot to save.", condition_num))
  }
  
  message(paste("Condition", condition_num, "enrichment analysis complete."))
  
  # Return enrichment results to build summary report later
  list(
    condition = condition_num,
    kegg = kk,
    go = go_enrich,
    gsea = gsea_res
  )
}

# --- Run enrichment analyses for all conditions ---

all_results <- list()
for (cond in 1:num_conditions) {
  all_results[[cond]] <- run_enrichment_for_condition(cond, deg_data, ko_term2gene, go_term2gene,mode = "both")
}


# Helper function to format tables as HTML
format_top_hits <- function(df, n = 10) {
  rownames(df) <- NULL
  if (nrow(df) == 0) return("<p>No enriched terms found.</p>")
  
  cols_to_try <- c("ID","Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID", "Count",
                   "enrichmentScore", "NES", "setSize")
  existing_cols <- intersect(cols_to_try, colnames(df))
  
  kable(df %>%
          arrange(p.adjust) %>%
          slice_head(n = n) %>%
          dplyr::select(all_of(existing_cols)),
        format = "html", table.attr = "border='1' cellpadding='5' cellspacing='0'") %>%
    as.character()
}

# Start HTML report content
html_report <- c(
  "<html>",
  "<head>",
  "<title>Enrichment Analysis Summary Report</title>",
  "<style>",
  "body { font-family: Arial, sans-serif; max-width: 900px; margin: auto; padding: 1em; }",
  "h1, h2, h3 { color: #2c3e50; }",
  "table { border-collapse: collapse; width: 100%; margin-bottom: 2em; }",
  "th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }",
  "th { background-color: #f4f4f4; }",
  "img { max-width: 100%; height: auto; margin-bottom: 2em; }",
  "</style>",
  "</head>",
  "<body>",
  "<h1>Enrichment Analysis Summary Report</h1>"
)

# Append per-condition results
for (res in all_results) {
  cnd <- res$condition
  kk_df <- as.data.frame(res$kegg)
  if (ncol(kk_df) > 0 && (is.null(colnames(kk_df)[1]) || colnames(kk_df)[1] == "")) {
    colnames(kk_df)[1] <- "ID"
  }
  go_df <- as.data.frame(res$go)
  gsea_df <- as.data.frame(res$gsea)
  
  html_report <- c(html_report,
                   sprintf("<h2>Condition %d</h2>", cnd),
                   
                   "<h3>KEGG Pathway Enrichment</h3>",
                   format_top_hits(kk_df, 10),
                   sprintf("<img src='enrichment_analysis/KEGG_dotplot_condition_%d.svg' alt='KEGG dotplot Condition %d'>", cnd, cnd),
                   
                   "<h3>GO Term Enrichment</h3>",
                   format_top_hits(go_df, 10),
                   sprintf("<img src='enrichment_analysis/GO_barplot_condition_%d.svg' alt='GO barplot Condition %d'>", cnd, cnd),
                   
                   "<h3>GSEA Results</h3>",
                   format_top_hits(gsea_df, 10),
                   sprintf("<img src='enrichment_analysis/GSEA_plot_condition_%d.svg' alt='GSEA plot Condition %d'>", cnd, cnd)
  )
}

# Close HTML tags
html_report <- c(html_report,
                 "</body>",
                 "</html>")

# Write to file
writeLines(html_report, "enrichment_summary_report.html")

message("HTML report 'enrichment_summary_report.html' generated successfully.")