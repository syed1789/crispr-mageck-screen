#!/usr/bin/env Rscript
# Functional enrichment (GO Biological Process + KEGG) on the top
# negatively-selected ("essential") genes from the MAGeCK RRA hit-calling
# in ../results/hl60_kbm7.gene_summary.txt.
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
gene_summary_path <- ifelse(length(args) >= 1, args[1], "/work/results/hl60_kbm7.gene_summary.txt")
out_dir <- ifelse(length(args) >= 2, args[2], "/work/enrichment")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gs <- fread(gene_summary_path)

# Rank by MAGeCK's negative-selection FDR: genes whose sgRNAs dropped out
# between T0 and endpoint, i.e. essential-gene candidates in this screen.
top_neg <- gs[order(`neg|fdr`)][1:30]
cat(sprintf("Top %d negatively-selected genes (lowest neg|fdr):\n", nrow(top_neg)))
print(top_neg[, .(id, `neg|score`, `neg|p-value`, `neg|fdr`, `neg|rank`, `neg|lfc`)])

gene_list <- top_neg$id

ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2
)

ego_df <- as.data.frame(ego)
fwrite(ego_df, file.path(out_dir, "go_bp_enrichment.tsv"), sep = "\t")
cat(sprintf("\nGO:BP enrichment: %d terms passed pvalueCutoff=0.2 (written to go_bp_enrichment.tsv)\n", nrow(ego_df)))

entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ekegg <- enrichKEGG(
  gene          = entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2
)
ekegg_df <- as.data.frame(ekegg)
fwrite(ekegg_df, file.path(out_dir, "kegg_enrichment.tsv"), sep = "\t")
cat(sprintf("KEGG enrichment: %d pathways passed pvalueCutoff=0.2 (written to kegg_enrichment.tsv)\n", nrow(ekegg_df)))

if (nrow(ego_df) > 0) {
  p <- dotplot(ego, showCategory = min(15, nrow(ego_df))) +
    ggtitle("GO:BP enrichment — top MAGeCK-depleted genes (HL60/KBM7 screen)")
  ggsave(file.path(out_dir, "go_bp_dotplot.png"), p, width = 8, height = 6, dpi = 150)
  cat("Saved go_bp_dotplot.png\n")
} else {
  cat("No GO:BP terms passed the cutoff — skipping dotplot.\n")
}

cat("\nDone.\n")
