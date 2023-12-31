library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
#BiocManager::install("fgsea")
library('fgsea')
library('GSEABase')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('data/verse_counts.tsv', 'data/sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts <- read.table(counts_csv, sep = "\t", header = TRUE, row.names = 1)
  meta <- read_csv(metafile_csv)
  #row_data <- counts["gene"]
  meta_subset <- meta[meta$timepoint %in% selected_times,]
  counts <- counts[,meta_subset$samplename]
  #rownames(counts) <- row_data$gene
  col_data <- DataFrame(samplename = meta_subset$samplename, timepoint = factor(meta_subset$timepoint, levels = selected_times))
  #meta_subset$timepoint <- relevel(as.factor(meta_subset$timepoint),ref = 'vP0')
  se <- SummarizedExperiment(
    assays=list(counts = as.matrix(counts)),
    colData=col_data,
    #rowData=row_data
  )
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = ~timepoint)
  # compute normalization factor
  #dds <- estimateSizeFactors(dds)
  # extract the normalized counts
  #deseq_norm_counts <- as_tibble(counts(dds,normalized=TRUE)) %>%
  #  mutate(gene=se@elementMetadata@listData$gene) %>%
  #  relocate(gene)
  dds <- DESeq(dds)
  #resultsNames(dds)
  # set vP0 as reference
  res <- results(dds,  contrast = c("timepoint", "vAd", "vP0"))
  #p0_vs_Ad_de <- as_tibble(res) %>%
  #  mutate(genes=se@elementMetadata@listData$gene) %>%
  #  relocate(genes) %>%
  #  arrange(pvalue)
  my_list <- list(dds, as.data.frame(res))
  return(my_list)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(results[[2]], .10)
label_res <- function(deseq2_res, padj_threshold) {
  #get the result table
  #deseq2_res <- deseq2_res[[2]]
  UP <- which(deseq2_res$padj<padj_threshold & deseq2_res$log2FoldChange>0)
  DOWN <- which(deseq2_res$padj<padj_threshold & deseq2_res$log2FoldChange<0)
  deseq2_res <- deseq2_res %>% mutate(volc_plot_status = rep("NS",nrow(deseq2_res))) %>% relocate(volc_plot_status,.before =  log2FoldChange)
  deseq2_res$volc_plot_status[UP] = "UP"
  deseq2_res$volc_plot_status[DOWN] = "DOWN"
  deseq2_res <- cbind(genes = rownames(deseq2_res), deseq2_res)
  return(as_tibble(deseq2_res))
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  p <- dplyr::select(labeled_results, pvalue) %>% ggplot(aes(x = pvalue)) + 
    geom_histogram(binwidth = 0.01, fill = "lightblue2", color = "black") +
    labs(title = "Histogram of raw p-values from DE analysis (vP0 vs. vAd)",
         x = "pvalue",
         y = "count") +
    theme_minimal()
  return(p)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  p <- dplyr::filter(labeled_results, labeled_results$padj<0.1) %>% dplyr::select(log2FoldChange) %>%
    ggplot(aes(x = log2FoldChange)) +
    geom_histogram(binwidth = 0.2, fill = "lightblue2", color = "black") +
    labs(title = "Histogram of log2FoldChange for DE genes (vP0 vs. vAd)",
         x = "log2FoldChange",
         y = "count") +
    theme_minimal()
  return(p)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  #rank data
  rank <- rank(labeled_results$padj)
  top10_gene <- labeled_results[rank<=num_genes,]
  #get normalized counts with gene name
  dds <- results[[1]] 
  counts <- as.data.frame(counts(dds, normalized = TRUE)[top10_gene$genes,])
  counts <- cbind(genes = rownames(counts), counts)
  #transform the data
  counts_long <- counts %>%
    pivot_longer(cols = -genes, names_to = "SampleName", values_to = "Count")
  counts_long$SampleName <- as.factor(counts_long$SampleName)
  counts_long$Count <- log10(counts_long$Count + 1)
  #plot
  p <- ggplot(data = counts_long, aes(x = genes, y = Count, color = SampleName)) +
    geom_point() +
    labs(x = "Gene", y = "Counts") +
    theme_minimal() +
    scale_y_continuous(trans = "log10") +
    theme(legend.position = "right",axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  p <- mutate(labeled_results, `FDR<0.05`=padj<0.05, `-log10(padj)` = -log10(padj)) %>% 
    ggplot(aes(x = log2FoldChange, y = `-log10(padj)`, color=volc_plot_status)) +
    geom_point()+
    geom_hline(yintercept = 0,linetype = "dashed", colour = "black")
  labs(title = "Volcano plot of DESeq2 differential expression results (vP0 vs vAd)")
  return(p)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  #read id2 file
  id <- read_delim(id2gene_path,col_names = FALSE) %>% dplyr::select(genes = X1, symbol = X2)
  #add the gene name as a new column
  join <- left_join(labeled_results, id, by = "genes") # by function could be c('gene' = 'id')
  #remove NA from df
  join_rmna <- na.omit(join)
  #sort the data in descending order
  desc_data <- join_rmna[order(-join_rmna$log2FoldChange),] #or arrange(desc(join_rmna$log2FoldChange))
  vector <- setNames(desc_data$log2FoldChange,nm = desc_data$symbol)
  return(vector)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  gene_sets <-  getGmt(con= gmt_file_path)
  hallmark_pathways_GSEABase <- geneIds(gene_sets)
  #remove NAs from rnk_list
  #rnk_list <- rnk_list[!names(rnk_list) == "NA"]
  #rnk_list <- rnk_list[!sapply(rnk_list, is.na)]
  results <- fgsea(hallmark_pathways_GSEABase,rnk_list, minSize = min_size, maxSize = max_size)
  results <- as_tibble(results)
  return(results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  #filter by FDR threshold and subset significant gene sets by the direction of their NES
  top_positive_nes <- fgsea_results %>%
    filter(padj < .25 & NES > 0) %>%
    slice_max(NES, n=num_paths)
  top_negative_nes <- fgsea_results %>%
    filter(padj < .25 & NES < 0) %>%
    slice_min(NES, n=num_paths)
  top_fgsea_results <- rbind(top_positive_nes,top_negative_nes)
  #plot
  p <- top_fgsea_results %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
    ggplot() +
    geom_bar(aes(x=pathway, y=NES, fill = ifelse(NES > 0, 'Positive', 'Negative')), stat='identity', width = 0.7) +
    scale_fill_manual(values = c('Positive' = 'red', 'Negative' = 'blue')) + 
    theme_minimal() +
    ggtitle('fgsea results for Hallmark MSigDB gene sets') +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') +
    coord_flip() + 
    theme(axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 5),
          plot.title = element_text(size = 8), legend.position = "none")
  
  return(p)
}

