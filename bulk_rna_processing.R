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
QC.filter.readcount <- function(raw.count.df, min.read.cutoff,out.file.path) {
  
  # Create a data frame for read count summary
  read.count.summary <- data.frame(sample.read.sum = colSums(raw.count.df))
  # Plot histogram for pre-filtering read counts
  pre.filter.plot <- (ggplot(data = read.count.summary, aes(x=sample.read.sum)) 
                      + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
                      + ggtitle("Counts for total reads per Sample prefilter")
                      + labs(x = "Number of reads in Sample", y = "Count")
                      + theme_bw())
  ggsave(file.path(out.file.path, "prefiltered_samp_readcount_hist.pdf"), pre.filter.plot)
  
  # Select samples with read counts greater than min.read.cutoff
  # Save the filtered count table and the read count summary as metadata
  count.tab.filt <- raw.count.df[, colSums(raw.count.df) > min.read.cutoff]

  # Plot histogram for post-filtering read counts
  post.filter.plot <- (ggplot(data = data.frame(sample.read.sum = colSums(count.tab.filt)), aes(x=sample.read.sum))
                    + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
                    + ggtitle("Counts for total reads per Sample postfilter")
                    + labs(x = "Number of reads in Sample", y = "Count")
                    + theme_bw())
  ggsave(file.path(out.file.path, "postfiltered_samp_readcount_hist.pdf"), pre.filter.plot)
  
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
QC.filter.genecount <- function(count.tab.filt2, min.count.cutoff, min.samples.cutoff,out.file.path) {
  
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
  ggsave(file.path(out.file.path, "prefiltered_gene_rowmeans_hist.pdf"), pre.filter.plot)
  
  # Create a data frame for post-filtering gene means
  post.filter.data <- data.frame(genemeans = log2(rowMeans(count.tab.filt)))
  post.filter.data <- post.filter.data[is.finite(rowSums(post.filter.data)),,drop=FALSE]
  
  # Plot histogram for post-filtering gene means
  post.filter.plot <- (ggplot(data = post.filter.data, aes(x=genemeans)) 
    + geom_histogram(color = "black", bins = 50, alpha = 0.3, position="identity") 
    + ggtitle("Counts for rowMeans per Gene postfilter")
    + labs(x = "log2(RowMean)", y = "Count")
    + theme_bw())
  ggsave(file.path(out.file.path, "postfiltered_gene_rowmeans_hist.pdf"), post.filter.plot)
  
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
DEseq.Normalization <- function(count.table, meta.table, out.file.path) {
  
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
  ggsave(file.path(out.file.path, "raw_count_per_sample_grid1.pdf"), p)
  
  # Histogram for normalized data
  p2 <- (ggplot(data =  as.data.frame(raw.v.norm.tab), aes(x=normcount)) 
    + geom_histogram(color = "black", bins = 10, alpha = 0.3, position="identity") 
    + ggtitle("Normalized counts per Sample")
    + labs(x = "Reads in sample", y = "Count")
    + theme_bw())
  ggsave(file.path(out.file.path, "norm_count_per_sample_grid2.pdf"), p2)
  
  
  
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
  
  ggsave(file.path(out.file.path, "norm_ratio_per_sample_grid3.pdf"), p3)
  

  # Write normalized count table to file
  write.table(norm.count.tab, file = paste(out.file.path, "normcounttab.txt", sep = ""), sep = "\t",
              col.names = NA, row.names = TRUE, quote = FALSE)
  
  # Return a list containing the normalized count table and identified outliers
  return(list(normcounttab = norm.count.tab, deseqnorm_outliers = blowup.samples))
}
