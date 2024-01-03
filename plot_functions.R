##########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Romane Cathelin
# 
# Date: 2024-12-20
# 
# Script Name: Plot functions
# 
# Notes:
# 


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################

packages.list <- c("ggplot2", "DESeq2", "EnhancedVolcano")
for (package in packages.list) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
sapply(packages.list, library, character.only = TRUE)


###########################################################################
#
#                                 CODE
#
###########################################################################

# Function to plot z-scores on a heatmap, and cluster by the condition used in DESeq2 
# Arguments:
#   dds: DESeq2 object
#   log2FC_cutoff: log2FC threshold 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   heatmap
plot_heatmap_zscores_categories <- function(dds, log2FC_cutoff = 2, padj_cutoff = 0.1){
  res <- results(dds)
  res.df <- as.data.frame(res)
  # Filter data to only keep significant according to cutoff 
  res.df <- res.df[(abs(res.df$log2FoldChange) > log2FC_cutoff) & ((!is.na(res.df$padj) & res.df$padj < padj_cutoff)),]
  print(dim(res.df))
  # Extract normalized counts for significant genes
  mat <- counts(dds, normalized = TRUE)[rownames(res.df),]
  mat.z <- t(scale(t(mat)))
  # Extract categories
  categories <- resultsNames(dds)[[2]]
  categories <- unlist(strsplit(categories, split = "_"))
   p <- pheatmap(mat.z,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         labels_col = colnames(mat.z),
         annotation_col = data.frame(categories = dds[[1]]),
         annotation_colors = list(categories = setNames(c("#3498db", "#e74c3c"), categories[c(2, 4)])),
         main = "Gene Expression Heatmap",
         border_color = "white",
         show_colnames = TRUE,
         show_rownames = TRUE
  )
  return(p)

}

# Function to plot log2FC on a heatmap, and cluster by the condition used in DESeq2 
# Arguments:
#   dds: DESeq2 object
#   log2FC_cutoff: log2FC threshold 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   heatmap
plot_heatmap_log2FC_categories <- function(dds, log2FC_cutoff = 2, padj_cutoff = 0.1){
  res <- results(dds)
  res.df <- as.data.frame(res)
  res.df <- res.df[(abs(res.df$log2FoldChange) > log2FC_cutoff) & ((!is.na(res.df$padj) & res.df$padj < padj_cutoff)),]

  
  # Extract normalized counts for significant genes
  mat <- counts(dds, normalized = TRUE)[rownames(res.df),]
  categories <- resultsNames(dds)[[2]]
  categories <- unlist(strsplit(categories, split = "_"))
  mat_log <- log2(mat + 1)  # Adding 1 to avoid log(0)
  # Annotate columns (patients) with sex information for the Heatmap function
  pheatmap(mat_log,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           labels_col = colnames(mat_log),
           annotation_col = data.frame(categories = dds[[1]]),
           annotation_colors = list(categories = setNames(c("#3498db", "#e74c3c"), categories[c(2, 4)])),
           main = "Heatmap",
           fontsize_row = 8,
           fontsize_col = 8
  )
}


# Function to plot a volcano-Plot  
# Arguments:
#   dds: DESeq2 object
#   title_custom: title of the plot
#   log2FC_cutoff: log2FC threshold 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   plot
volcanoplot <- function(dds, title_custom, FC_cutoff = 0.5, padj_cutoff = 0.1, ylim_custom = 'default'){
  res <- results(dds)
  
  p <- EnhancedVolcano(toptable = res,
                  x = "log2FoldChange",
                  y = "padj",
                  lab = rownames(res),
                  # xlim = c(-5, +5),
                  ylim = ylim_custom,
                  pCutoff = padj_cutoff,
                  FCcutoff = FC_cutoff,
                  # pointSize = 2.0,
                  title = paste0(title_custom, " \n Fold change cutoff = ", FC_cutoff, ", p-adj cutoff = ", padj_cutoff),
                  legendLabels=c(
                    'Not significant',
                    'Pass Log2 fold-change cutoff',
                    'Pass p-adj cutoff',
                    'Pass both p-adj & Log2 fold change')
  )
  return(p)
}