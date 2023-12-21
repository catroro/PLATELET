##########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Romane Cathelin
# 
# Date: 2024-12-21
# 
# Script Name: Processing bulkRNA functions
# 
# Notes:
# 


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################


library(dplyr)
library(ggplot2)
library(DESeq2)
library(reader)



###########################################################################
#
#                                 CODE
#
###########################################################################


# Function to perform quality control filtering based on read count
# Arguments:
#   raw.count.df:
#   min.read.cutoff: 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   list: 
QC.filter.readcount <- function(raw.count.df, min.read.cutoff,out.path) {
  
  # Create a data frame for read count summary
  read.count.summary <- data.frame(sample.read.sum = colSums(raw.count.df))
  # Plot histogram for pre-filtering read counts
  pre.filter.plot <- (ggplot(data = read.count.summary, aes(x=sample.read.sum)) 
                      + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
                      + ggtitle("Counts for total reads per Sample prefilter")
                      + labs(x = "Number of reads in Sample", y = "Count")
                      + theme_bw())
  ggsave(cat.path(out.path, "prefiltered_samp_readcount_hist.pdf"), pre.filter.plot)
  
  # Select samples with read counts greater than min.read.cutoff
  # Save the filtered count table and the read count summary as metadata
  count.tab.filt <- raw.count.df[, colSums(raw.count.df) > min.read.cutoff]

  # Plot histogram for post-filtering read counts
  post.filter.plot <- (ggplot(data = data.frame(sample.read.sum = colSums(count.tab.filt)), aes(x=sample.read.sum))
                    + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
                    + ggtitle("Counts for total reads per Sample postfilter")
                    + labs(x = "Number of reads in Sample", y = "Count")
                    + theme_bw())
  ggsave(cat.path(out.path, "postfiltered_samp_readcount_hist.pdf"), pre.filter.plot)
  
  # Print the dimensions of the filtered count table
  print(dim(count.tab.filt))
  
  # Identify omitted samples with counts below min.read.cutoff
  omitted.samples <- colnames(raw.count.df[, colSums(raw.count.df) <= min.read.cutoff, drop = FALSE])
  
  # Return a list containing pre-filtering and post-filtering histograms,
  # filtered count table, and a data frame with outlier read counts
  return(list(prefilter.hist = pre.filter.plot, postfilter.hist = post.filter.plot, 
              count.tab.filt = count.tab.filt, 
              outlier.counts = data.frame(read.count = colSums(raw.count.df[, omitted.samples, drop = FALSE]))))
}



# Function to perform quality control filtering based on gene count
# Arguments:
#   raw.count.df:
#   min.read.cutoff: 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   list: 
QC.filter.genecount <- function(count.tab.filt2, min.count.cutoff, min.samples.cutoff,out.path) {
  
  # Create a data frame for gene means
  gene.means.data <- data.frame(genemeans = log2(rowMeans(count.tab.filt2)))
  gene.means.data <- gene.means.data[is.finite(rowSums(gene.means.data)),,drop=FALSE]
  
  # Plot histogram for pre-filtering gene means
  pre.filter.plot <- (ggplot(data = gene.means.data, aes(x=genemeans)) 
    + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
    + ggtitle("Counts for rowMeans per Gene prefilter")
    + labs(x = "log2(RowMean)", y = "Count")
    + theme_bw())
  # Applying pre-filtering: Minimum value of min.count.cutoff in at least min.samples.cutoff of samples
  count.tab.filt <- count.tab.filt2[rowSums(count.tab.filt2 >= min.count.cutoff) >= min.samples.cutoff,]
  ggsave(cat.path(out.path, "prefiltered_gene_rowmeans_hist.pdf"), pre.filter.plot)
  
  # Create a data frame for post-filtering gene means
  post.filter.data <- data.frame(genemeans = log2(rowMeans(count.tab.filt)))
  post.filter.data <- post.filter.data[is.finite(rowSums(post.filter.data)),,drop=FALSE]
  
  # Plot histogram for post-filtering gene means
  post.filter.plot <- (ggplot(data = post.filter.data, aes(x=genemeans)) 
    + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
    + ggtitle("Counts for rowMeans per Gene postfilter")
    + labs(x = "log2(RowMean)", y = "Count")
    + theme_bw())
  ggsave(cat.path(out.path, "postfiltered_gene_rowmeans_hist.pdf"), post.filter.plot)
  
  # Print the dimensions of the filtered count table
  print(dim(count.tab.filt))
  
  # Return a list containing pre-filtering and post-filtering histograms and the filtered count table
  return(list(prefilter.hist = pre.filter.plot, postfilter.hist = post.filter.plot, count.tab.filt = count.tab.filt))
}



# Function for DESeq normalization
# Arguments:
#   raw.count.df:
#   min.read.cutoff: 
#   padj_cutoff: p-adjusted threshold
# Outputs:
#   list:
DEseq.Normalization <- function(count.table, meta.table, out.path) {
  
  # Create DESeqDataSet and perform normalization
  dds <- DESeqDataSetFromMatrix(countData = count.table, colData = meta.table, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  norm.count.tab <- counts(dds, normalized = TRUE)
  
  # Visualizations
  raw.v.norm.tab <- cbind(rawcount = colSums(count.table), normcount = colSums(norm.count.tab), 
                          norm.to.raw.ratio = colSums(norm.count.tab) / colSums(count.table))
  raw.v.norm.tab <- raw.v.norm.tab[order(raw.v.norm.tab[, 2], decreasing = TRUE), ]

  # Histogram for raw counts data
  p <- (ggplot(data =  as.data.frame(raw.v.norm.tab), aes(x=rawcount)) 
    + geom_histogram(color = "black", bins = 10, alpha = 0.3, position="identity") 
    + ggtitle("Raw counts per Sample")
    + labs(x = "Reads in sample", y = "Count")
    + theme_bw())
  ggsave(cat.path(out.path, "raw_count_per_sample_grid1.pdf"), p)
  
  # Histogram for normalized data
  p2 <- (ggplot(data =  as.data.frame(raw.v.norm.tab), aes(x=normcount)) 
    + geom_histogram(color = "black", bins = 10, alpha = 0.3, position="identity") 
    + ggtitle("Normalized counts per Sample")
    + labs(x = "Reads in sample", y = "Count")
    + theme_bw())
  ggsave(cat.path(out.path, "norm_count_per_sample_grid2.pdf"), p2)
  
  
  
  scaled.readcount.tab <- scale(data.frame(colSums(norm.count.tab)))
  blowup.samples <- scaled.readcount.tab[abs(scaled.readcount.tab) > 4, 1, drop = FALSE]
  colnames(blowup.samples) <- "norm_diff"
  
  # Identify outliers based on normalized count differences
  readcount.tab <- data.frame(colSums(norm.count.tab))
  scaled.readcount.tab <- scale(readcount.tab)
  blowup.samples <- scaled.readcount.tab[abs(scaled.readcount.tab) > 4, 1, drop = FALSE]
  colnames(blowup.samples) <- "norm_diff"
  # Scatter plot comparison between raw and normalized counts
  indata <- as.data.frame(raw.v.norm.tab[, c(1, 2)])

  datalabels = rownames(blowup.samples)
  indata$labels <- ""
  indata[datalabels,"labels"] <- datalabels


  p3 <- (ggplot(indata, aes(x=rawcount, y=normcount, label = labels)) 
        + geom_point(alpha = 0.7, size=1) 
        + theme_bw()    
        + geom_label_repel(mapping = aes(size = NULL, label = labels, fill = NULL),
                         box.padding = unit(0.5, "lines"), point.padding = unit(0.5, "lines"),
                         size = 2, max.overlaps = 10000)
        + ggtitle("Raw to norm count comparison")
        + xlab("Raw")
        + ylab("Norm")
  )
  
  ggsave(cat.path(out.path, "norm_ratio_per_sample_grid3.pdf"), p3)
  

  # Write normalized count table to file
  write.table(norm.count.tab, file = paste(out.path, "normcounttab.txt", sep = ""), sep = "\t",
              col.names = NA, row.names = TRUE, quote = FALSE)
  
  # Return a list containing the normalized count table and identified outliers
  return(list(normcounttab = norm.count.tab, deseqnorm_outliers = blowup.samples))
}








calculate.gene.metrics <- function(gene.list, wb.mtx, plt.mtx, path) {
  gene.out.df <- data.frame()
  cohort <- tail(unlist(strsplit(path, split = "/")), n = 1)
  print(cohort)
  for (gene in gene.list) {
    wb.gene.expr <- t(wb.mtx[gene, , drop = FALSE])
    plt.gene.expr <- t(plt.mtx[gene, , drop = FALSE])
    corr.plt.wb <- corr.test(wb.gene.expr, plt.gene.expr, method = "spearman")
    
    gene.out.df <- rbind(gene.out.df, c(gene, mean(wb.gene.expr), mean(plt.gene.expr), corr.plt.wb$r[1, 1], corr.plt.wb$p[1, 1]))
  }
  print(dim(gene.out.df))
  colnames(gene.out.df) <- c("gene", "wb.avg", "plt.avg", "Rval", "Pval")
  gene.out.df[, c("wb.avg", "plt.avg", "Rval", "Pval")] <- apply(gene.out.df[, c("wb.avg", "plt.avg", "Rval", "Pval")], 2, as.numeric)
  gene.out.df[, c("wb.avg.log2", "plt.avg.log2")] <- log2(gene.out.df[, c("wb.avg", "plt.avg")] + 1)
  print(colnames(gene.out.df))
  
  title <- paste0(cohort, ": Whole blood expression correlate to platelet expression")
  p <- ggplot(gene.out.df, aes(x = wb.avg.log2, y = plt.avg.log2)) +
    geom_point(alpha = 0.3, size = 1) +
    theme_bw() +
    sm_statCorr(corr_method = 'spearman', color = '#8DD3C7') +
    ggtitle(title) +
    xlab("Whole blood expression (log2)") +
    ylab("Platelet expression (log2)")
  print(p)
  print(cat.path(path, "shared.genes.post.process.csv"))
  write.csv(gene.out.df, cat.path(path, "shared.genes.post.process.csv"), row.names = FALSE)
  print(cat.path(path, "scatterplot.plt.wb.pdf"))
  ggsave(cat.path(path, "scatterplot.plt.wb.pdf"), p)
  return(gene.out.df)
}







# Function to plot genes subset
plot.genes.sub <- function(top.genes, tissue, path) {
  # Define column name dynamically based on tissue
  col <- paste0(tissue, ".avg.log2")
  
  # Select top 20 genes
  genes.sub <- top.genes[, c("gene", "plt_rank", "wb_rank", substitute(col))] %>%  
    slice_head(n = 20) 
  
  # Set up color palette
  nb.cols <- 20
  pal <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
  
  # Set up labels and plot based on tissue
  if (tissue == "wb") { 
    lab <- c("WB ranking", "PLT ranking")
    plot <- ggplot(genes.sub, aes(axis1 = wb_rank , axis2 = plt_rank, y = genes.sub[, length(genes.sub)]))
  } else {
    lab <- c("PLT ranking", "WB ranking")
    plot <- ggplot(genes.sub, aes(axis1 = plt_rank , axis2 = wb_rank, y = genes.sub[, length(genes.sub)]))
  }
  
  # Create the plot
  p <- plot +
    geom_alluvium(aes(fill = gene)) +
    geom_stratum(width = 0.2, aes(fill = gene)) +
    theme_bw() +
    labs(x = " ", y = " ") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    ggtitle("Most highly expressed genes in platelet and whole blood") +
    scale_x_continuous(breaks = c(1, 2), labels = lab) +
    scale_fill_manual(values = pal)
  
  # Save the plot
  ggsave(cat.path(path, paste0(tissue, "_alluvium_plot_ranking.pdf")), p)
  
  # top.genes.plt.pace %>%
  #   dplyr::select(gene,plt_rank, plt_avg,wb_rank, wb_avg) %>%
  #   kbl("html", col.names = c("Gene", "Rank", "Avg expr", "Rank", "Avg expr")) %>%
  #   kable_styling(bootstrap_options = "striped", full_width = F) %>%
  #   add_header_above(c(" " = 1, "Platelet" = 2, "Whole blood" = 2))%>%
  #   kable_minimal()%>%
  #   save_kable(paste0(out.path, "pace_analysis/", "PLTtable_topgenes.html"))
}

# Function to analyze top genes
analyze.top.genes <- function(gene.out.df, path) {
  # Rank genes based on tissue-specific averages
  gene.out.df.rank <- gene.out.df %>%
    mutate(plt_rank = rank(-plt.avg),
           wb_rank = rank(-wb.avg))
  
  # Select top 100 genes for each tissue
  top.genes.wb <- gene.out.df.rank %>%
    arrange(wb_rank) %>%
    slice_head(n = 100) 
  
  top.genes.plt <- gene.out.df.rank %>%
    arrange(plt_rank) %>%  
    slice_head(n = 100) 
  
  # Print top genes for whole blood
  print(top.genes.wb)
  
  # Plot genes for platelet and whole blood
  plot.genes.sub(top.genes.plt, "plt", path)
  plot.genes.sub(top.genes.wb, "wb", path)
  return(list("wb"= top.genes.wb, "plt"= top.genes.plt))
}



top.gene.enrichment <- function(top.genes.list,path){
  
  deg_wb = top.genes.list$wb$gene
  deg_plt = top.genes.list$plt$gene
  
  deg_wb_entrez <- (mapIds(org.Hs.eg.db, keys = deg_wb, column = "ENTREZID", keytype = "SYMBOL"))
  deg_plt_entrez <- (mapIds(org.Hs.eg.db, keys = deg_plt, column = "ENTREZID", keytype = "SYMBOL"))
  
  go.wb <- enrichGO(deg_wb_entrez, OrgDb=org.Hs.eg.db,pvalueCutoff  = 0.05)
  go.plt <- enrichGO(deg_plt_entrez, OrgDb=org.Hs.eg.db,pvalueCutoff  = 0.05)
  
  kegg.wb <- enrichKEGG(deg_wb_entrez,pvalueCutoff  = 0.05)
  kegg.plt <- enrichKEGG(deg_plt_entrez,pvalueCutoff  = 0.05)
  
  # gene_ids_list <- list(wb = deg_wb_entrez, plt = deg_plt_entrez)
  # ck <- compareCluster(geneCluster = gene_ids_list, fun = enrichKEGG )   #,  OrgDb='org.Hs.eg.db')
  # ck2 <- compareCluster(geneCluster = gene_ids_list, fun = enrichGO,  OrgDb='org.Hs.eg.db')
  
  pdf(cat.path(path, "enrichment_analysis.pdf"))
  print(dotplot(kegg.wb, showCategory=20, title = "KEGG - whole blood"))
  print(dotplot(kegg.plt, showCategory=20, title = "KEGG - platelet"))
  print(dotplot(go.wb, showCategory=20, title = "GO - whole blood"))
  print(dotplot(go.plt, showCategory=10, title = "GO - platelet") )
  dev.off()
}
