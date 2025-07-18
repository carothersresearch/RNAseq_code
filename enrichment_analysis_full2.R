# Load required libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(enrichplot)
library(knitr)
library(GO.db)

dir.create("enrichment_analysis2")

# --- User inputs: file paths ---
deg_file <- "deg_results_all.csv"
annot_file <- "chatgpt_annotations.csv"
num_conditions <- 6  # number of conditions in DEG file

# --- Read data ---
deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)
annot <- read.csv(annot_file, stringsAsFactors = FALSE)

# --- Helper to parse list-like columns ---
parse_list_string <- function(s) {
  if (is.na(s) || s == "[]") return(character(0))
  s_clean <- gsub("^\\['|\\']$", "", s)
  if (nchar(s_clean) == 0) return(character(0))
  strsplit(s_clean, "', '")[[1]]
}

# Parse KO/GO
annot$KO_IDs_list <- lapply(annot$KO_IDs, parse_list_string)
annot$GO_IDs_list <- lapply(annot$GO_IDs, parse_list_string)

# Parse Tags & Subcategories ### NEW
parse_list_string_simple <- function(s) {
  if (is.na(s) || s == "[]") return(character(0))
  gsub("\\[|\\]|'", "", s) %>% strsplit(",\\s*") %>% unlist()
}
tags_term2gene <- annot %>%
  dplyr::select(Locus_Tag, Tags) %>%
  dplyr::mutate(Tags = lapply(Tags, parse_list_string_simple)) %>%
  tidyr::unnest(Tags) %>%
  dplyr::filter(Tags != "") %>%
  dplyr::rename(gene = Locus_Tag, term = Tags)

subcat_term2gene <- annot %>%
  dplyr::select(Locus_Tag, Subcategories) %>%
  dplyr::mutate(Subcategories = lapply(Subcategories, parse_list_string_simple)) %>%
  tidyr::unnest(Subcategories) %>%
  dplyr::filter(Subcategories != "") %>%
  dplyr::rename(gene = Locus_Tag, term = Subcategories)

# Build KO TERM2GENE
ko_term2gene <- annot %>%
  dplyr::select(Locus_Tag, KO_IDs_list) %>%
  tidyr::unnest(cols = c(KO_IDs_list)) %>%
  dplyr::filter(KO_IDs_list != "") %>%
  dplyr::rename(gene = Locus_Tag, KO = KO_IDs_list)

# Build GO TERM2GENE
go_term2gene <- annot %>%
  dplyr::select(Locus_Tag, GO_IDs_list) %>%
  tidyr::unnest(GO_IDs_list) %>%
  dplyr::rename(gene = Locus_Tag, term = GO_IDs_list) %>%
  dplyr::filter(!is.na(term) & term != "")

unique_terms <- unique(go_term2gene$term)
go_term2name <- data.frame(
  term = unique_terms,
  name = sapply(unique_terms, function(x) {
    term_obj <- GO.db::GOTERM[[x]]
    if (is.null(term_obj)) return(NA_character_) else return(term_obj@Term)
  }),
  stringsAsFactors = FALSE
) %>% dplyr::filter(!is.na(name))

go_term2gene <- go_term2gene %>% dplyr::filter(term %in% go_term2name$term)

# Helper for Fisher's enrichment ### NEW
enrich_categorical <- function(term2gene_df, sig_genes, all_genes) {
  terms <- unique(term2gene_df$term)
  results <- lapply(terms, function(t) {
    genes_in_term <- unique(term2gene_df$gene[term2gene_df$term == t])
    a <- sum(sig_genes %in% genes_in_term)
    b <- sum(!sig_genes %in% genes_in_term)
    c <- sum(all_genes %in% genes_in_term) - a
    d <- length(all_genes) - (a+b+c)
    test <- fisher.test(matrix(c(a, b, c, d), nrow=2), alternative = "greater")
    c(
      Term = t,
      Sig_in_term = a,
      Bg_in_term = a+c,
      pvalue = test$p.value
    )
  })
  results_df <- do.call(rbind, results) %>% as.data.frame(stringsAsFactors=FALSE)
  results_df$p.adjust <- p.adjust(as.numeric(results_df$pvalue), method="BH")
  results_df <- results_df %>% arrange(p.adjust)
  return(results_df)
}

# Enrichment per condition
run_enrichment_for_condition <- function(condition_num, deg_data, ko_term2gene, go_term2gene, mode = c("ORA","GSEA","both")) {
  mode <- match.arg(mode)
  sig_col <- paste0("sig_condition_", condition_num, ".csv")
  logfc_col <- paste0("log2FC_condition_", condition_num, ".csv")
  
  sig_genes <- deg_data$gene[deg_data[[sig_col]] == "yes"]
  all_genes <- unique(deg_data$gene)
  
  # KEGG enrichment
  sig_kos <- unique(ko_term2gene$KO[ko_term2gene$gene %in% sig_genes])
  kk <- enrichKEGG(
    gene = sig_kos,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  # GO enrichment
  go_enrich <- NULL
  go_universe <- unique(go_term2gene$gene)
  overlap_genes <- sum(sig_genes %in% go_universe)
  if ((mode %in% c("ORA","both")) && overlap_genes > 0) {
    go_enrich <- enricher(
      gene = sig_genes,
      TERM2GENE = go_term2gene,
      TERM2NAME = go_term2name,
      universe = go_universe,
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
      go_enrich@result <- dplyr::left_join(
        go_enrich@result,
        go_term2name,
        by = c("ID" = "term")
      )
      go_enrich@result$Description <- go_enrich@result$name
    }
  }
  
  # GSEA
  gsea_res <- NULL
  if (mode %in% c("GSEA","both")) {
    ranked_df <- deg_data %>%
      dplyr::select(gene, !!sym(logfc_col)) %>%
      dplyr::filter(!is.na(!!sym(logfc_col))) %>%
      arrange(desc(!!sym(logfc_col)))
    ranked_vec <- ranked_df[[logfc_col]]
    names(ranked_vec) <- ranked_df$gene
    gsea_res <- GSEA(
      geneList = ranked_vec,
      TERM2GENE = go_term2gene,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
      gsea_res@result <- dplyr::left_join(
        gsea_res@result,
        go_term2name,
        by = c("ID" = "term")
      )
      gsea_res@result$Description <- gsea_res@result$name
    }
  }
  
  # NEW: Enrichment for Tags and Subcategories
  tags_enrich_df <- enrich_categorical(tags_term2gene, sig_genes, all_genes)
  write.csv(tags_enrich_df, file = paste0("enrichment_analysis2/Tags_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  subcat_enrich_df <- enrich_categorical(subcat_term2gene, sig_genes, all_genes)
  write.csv(subcat_enrich_df, file = paste0("enrichment_analysis2/Subcategories_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  
  # Save plots and CSVs
  write.csv(as.data.frame(kk), file = paste0("enrichment_analysis2/KEGG_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  if (!is.null(go_enrich)) {
    write.csv(as.data.frame(go_enrich), file = paste0("enrichment_analysis2/GO_enrichment_condition_", condition_num, ".csv"), row.names = FALSE)
  }
  if (!is.null(gsea_res)) {
    write.csv(as.data.frame(gsea_res), file = paste0("enrichment_analysis2/GSEA_condition_", condition_num, ".csv"), row.names = FALSE)
  }
  
  # Plots
  svg(paste0("enrichment_analysis2/KEGG_dotplot_condition_", condition_num, ".svg"), width=8, height=6)
  print(dotplot(kk, showCategory=10) + ggtitle(paste("KEGG - Condition", condition_num)))
  dev.off()
  
  if (!is.null(go_enrich)) {
    svg(paste0("enrichment_analysis2/GO_barplot_condition_", condition_num, ".svg"), width=8, height=6)
    print(barplot(go_enrich, showCategory=10) + ggtitle(paste("GO - Condition", condition_num)))
    dev.off()
  }
  
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    svg(paste0("enrichment_analysis2/GSEA_plot_condition_", condition_num, ".svg"), width=8, height=6)
    print(gseaplot2(gsea_res, geneSetID=gsea_res@result$ID[1], title=paste("GSEA - Condition", condition_num)))
    dev.off()
  }
  
  # NEW: Plot top 10 Tags
  svg(paste0("enrichment_analysis2/Tags_enrichment_condition_", condition_num, ".svg"), width=8, height=6)
  ggplot(tags_enrich_df %>% slice_head(n=10), aes(x=reorder(Term,-as.numeric(pvalue)), y=-log10(as.numeric(pvalue)))) +
    geom_bar(stat="identity") + coord_flip() + theme_minimal() + ggtitle(paste("Tags - Condition", condition_num)) +
    xlab(NULL) + ylab("-log10(p-value)")
  dev.off()
  
  # NEW: Plot top 10 Subcategories
  svg(paste0("enrichment_analysis2/Subcategories_enrichment_condition_", condition_num, ".svg"), width=8, height=6)
  ggplot(subcat_enrich_df %>% slice_head(n=10), aes(x=reorder(Term,-as.numeric(pvalue)), y=-log10(as.numeric(pvalue)))) +
    geom_bar(stat="identity") + coord_flip() + theme_minimal() + ggtitle(paste("Subcategories - Condition", condition_num)) +
    xlab(NULL) + ylab("-log10(p-value)")
  dev.off()
  
  list(
    condition = condition_num,
    kegg = kk,
    go = go_enrich,
    gsea = gsea_res,
    tags = tags_enrich_df,
    subcategories = subcat_enrich_df
  )
}

# Run all conditions
all_results <- list()
for (cond in 1:num_conditions) {
  all_results[[cond]] <- run_enrichment_for_condition(cond, deg_data, ko_term2gene, go_term2gene, mode="both")
}

# Prepare HTML report
format_top_hits <- function(df, n=10) {
  rownames(df) <- NULL
  if (nrow(df) == 0) return("<p>No enriched terms found.</p>")
  cols_to_try <- c("Term","Description","GeneRatio","BgRatio","pvalue","p.adjust","geneID","Count","NES","enrichmentScore","setSize")
  existing_cols <- intersect(cols_to_try, colnames(df))
  kable(df %>% arrange(p.adjust) %>% slice_head(n=n) %>% select(all_of(existing_cols)), format="html", table.attr="border='1' cellpadding='5' cellspacing='0'") %>% as.character()
}

html_report <- c(
  "<html><head><title>Enrichment Analysis Summary</title>",
  "<style>body{font-family:sans-serif;max-width:900px;margin:auto;}h1,h2,h3{color:#2c3e50;}table{border-collapse:collapse;width:100%;margin-bottom:2em;}th,td{border:1px solid #ccc;padding:8px;}th{background:#f4f4f4;}img{max-width:100%;margin-bottom:2em;}</style>",
  "</head><body><h1>Enrichment Analysis Summary Report</h1>"
)

for (res in all_results) {
  cnd <- res$condition
  kk_df <- as.data.frame(res$kegg)
  go_df <- as.data.frame(res$go)
  gsea_df <- as.data.frame(res$gsea)
  tags_df <- as.data.frame(res$tags)
  subcat_df <- as.data.frame(res$subcategories)
  
  html_report <- c(html_report,
                   sprintf("<h2>Condition %d</h2>", cnd),
                   "<h3>KEGG</h3>",
                   format_top_hits(kk_df),
                   sprintf("<img src='enrichment_analysis2/KEGG_dotplot_condition_%d.svg'>", cnd),
                   "<h3>GO</h3>",
                   format_top_hits(go_df),
                   sprintf("<img src='enrichment_analysis2/GO_barplot_condition_%d.svg'>", cnd),
                   "<h3>GSEA</h3>",
                   format_top_hits(gsea_df),
                   sprintf("<img src='enrichment_analysis2/GSEA_plot_condition_%d.svg'>", cnd),
                   "<h3>Tags Enrichment</h3>",
                   format_top_hits(tags_df),
                   sprintf("<img src='enrichment_analysis2/Tags_enrichment_condition_%d.svg'>", cnd),
                   "<h3>Subcategories Enrichment</h3>",
                   format_top_hits(subcat_df),
                   sprintf("<img src='enrichment_analysis2/Subcategories_enrichment_condition_%d.svg'>", cnd)
  )
}

html_report <- c(html_report, "</body></html>")
writeLines(html_report, "enrichment_summary_report.html")

message("✅ Analysis complete. Check 'enrichment_summary_report.html' for the summary!")
