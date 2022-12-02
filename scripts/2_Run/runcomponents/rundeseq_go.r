###############
# Script to Run DESeq2 and Gene Ontology Searches
# Author: Justin Currie, Pranav Mahableshwarkar
# Date: 09/10/2022
# Note: Part of RNA-Seq pipeline designed to process counts tables and perform
#       gene ontology searches (GProfiler2) on significant results. 
###############

# load libraries
library(optparse)   
library(gprofiler2) 
library(DESeq2) 
library(tidyverse)  
library(pheatmap)  
library(vsn)
library(hexbin)  
library(EnhancedVolcano)  
library(patchwork) 
library(ashr) 


option_list <- list(
  make_option(c("-r", "--rawcounts"), type = "character", default = NULL,
              help = "raw counts table filepath (.txt)"),
  make_option(c("-c", "--numcontrol"), type = "numeric", default = NULL,
              help = "number of control replicates"),
  make_option(c("-e", "--numexperimental"), type = "numeric", default = NULL,
              help = "number of experimental replicates"),
  make_option(c("-g", "--genename"), type = "character", default = NULL,
              help = "The name of the gene of interest"),
  make_option(c("-p", "--perturbation"), type = "character", default = NULL,
              help = "The condition of the experimental group. Generally a gene
              name followed by some kind of method of perturbing gene expression
              or protein activity. ex. 'Deaf1RNAi', 'phoKO'."),
  make_option(c("-d", "--direction"), type = "character", default = NULL,
              help = "the direction of the expression or activity perturbation on 
              the protein in question. Either 'up' or 'down'. If the expression 
              or activity of the protein was increased, 'up'. If decreased, 'down'."),
  make_option(c("-f", "--flybaseconverter"), type = "character", default = NULL,
              help = "the filepath to a csv file containing 3 columns. First
              column is FlyBaseIDs outputted by featureCounts. Second is corresponding
              updated FlyBaseIDs. Third is corresponding gene symbols"),
  make_option(c("-o", "--out"), type = "character", default = "",
              help = "the path to the directory that results csvs and images will
              be written to.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

RunDESeq2GO = function(counts, numControl, numExpr, conditionString, symbols, outputDir){
  # allow for flexibility in specification of output directory (adds "/" if necessary)
  if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != "/") {
    outputDir <- paste(outputDir, "/", sep = "")
  }
  
  # get only columns with FBgn IDs and counts
  counts <- read.table(counts)
  colnames(counts) <- counts[c(1),]
  counts <- counts[-c(1), c(1, 7:length(colnames(counts)))]
  
  # order columns so that controls are first
  counts <- counts[, order(colnames(counts))]
  
  # change columnames of counts table
  colnames <- c("Gene")
  for (i in 1:numControl) {
    colnames <- append(colnames, paste("control", as.character(i), sep = ""))
  }
  for (i in 1:numExpr) {
    colnames <- append(colnames, paste(conditionString, as.character(i), sep = ""))
  }
  colnames(counts) <- colnames

  # get rid of extra lncRNA
  counts <- counts[-c(3063),]
  
  # Make all counts numeric prior to running. 
  for (i in 2:length(colnames(counts))) {
    counts[[i]] <- as.numeric(counts[[i]])
  }
  print(class(counts[,c(2)])) 
  
  # get rid of genes with 0 counts for every sample 
  counts <- counts[rowSums(counts[,2:length(colnames(counts))]) != 0,]
  
  # create metadata table
  condition <- c()
  for (i in 1:numControl) {
    condition <- append(condition, "control")
  }
  for (i in 1:numExpr) {
    condition <- append(condition, "experimental")
  }
  replicate <- c()
  for (i in 1:numControl) {
    replicate <- append(replicate, i)
  }
  for (i in 1:numExpr) {
    replicate <- append(replicate, i)
  }
  metadata <- data.frame(id = colnames(counts)[-1],
                         condition = condition,
                         replicate = replicate)
  write.csv(metadata, paste(outputDir, "metadata.csv", sep = ""))
  write.csv(counts, paste(outputDir, "counts.csv", sep = ""))
  
  
  # differential expression 
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~condition,
                                tidy = TRUE) 
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "experimental", "control"), tidy = TRUE)
  res <- dplyr::tbl_df(res)
  res <- dplyr::arrange(res, padj)
  
  # make table of significant results
  res_sig <- res[res$padj<0.05,]
  res_sig <- na.omit(res_sig)

  
  # separate upregulated and downregulated genes and write CSVs for each
  upreg <- res_sig %>%
    filter(log2FoldChange > 0) %>%
    arrange(pvalue)
  write.csv(upreg, paste(outputDir, conditionString, "_upreg_genes.csv", sep = ""))
  summary <- data.frame(condition = "upregulated genes", numOfGenes = length(upreg$row))
  downreg <- res_sig %>%
    filter(log2FoldChange < 0) %>%
    arrange(pvalue)
  write.csv(downreg, paste(outputDir, conditionString, "_downreg_genes.csv", sep = ""))
  summary <- rbind(summary, c("downregulated genes", length(downreg$row)))
  
  
  ### make quality check plots ###
  # apply a variance stabilizing transformation
  vsd <- vst(dds, blind = FALSE)
  # apply a regularized log transformation
  rld <- rlog(dds, blind = FALSE)
  # transform with log2(n + 1)
  ntd <- normTransform(dds)
  
  # make png file of vsd sdPLot
  png(filename = paste(outputDir, conditionString, "_vsd_meanSD.png", sep = ""))
  plot1 <- meanSdPlot(assay(vsd), xlab = "ranks(mean) vst")
  print(plot1)
  dev.off()
  
  # make pca plots
  pcaData1 <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
  percentVar1 <- round(100 * attr(pcaData1, "percentVar"))
  pcaData2 <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
  percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
  pcaData3 <- plotPCA(ntd, intgroup = c("condition"), returnData = TRUE)
  percentVar3 <- round(100 * attr(pcaData3, "percentVar"))
  
  png(filename = paste(outputDir, conditionString, "_pca.png", sep = ""), height = 600, width = 1800)
  plot1 <- ggplot(pcaData1, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text(
      label=gsub("[[:punct:][:blank:]]|aligned|sorted|bam","",rownames(pcaData1)),
      nudge_x = 0, nudge_y = 0.5, check_overlap = F, hjust = "inward") +
    xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
    ggtitle(paste(conditionString, " PCA (vst)", sep = "")) 
  plot2 <- ggplot(pcaData2, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text(
      label=gsub("[[:punct:][:blank:]]|aligned|sorted|bam","",rownames(pcaData2)),
      nudge_x = 0, nudge_y = 0.5, check_overlap = F, hjust = "inward") +
    xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
    ggtitle(paste(conditionString, " PCA (rlog)", sep = "")) 
  plot3 <- ggplot(pcaData3, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text(
      label=gsub("[[:punct:][:blank:]]|aligned|sorted|bam","",rownames(pcaData3)),
      nudge_x = 0, nudge_y = 0.5, check_overlap = F, hjust = "inward") +
    xlab(paste0("PC1: ",percentVar3[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar3[2],"% variance")) + 
    ggtitle(paste(conditionString, " PCA (normTransform)", sep = "")) 
  print(plot1 + plot2 + plot3)
  dev.off()
  
  # make ma plots
  plotMA(dds, ylim = c(-10, 10), main = paste(conditionString, " MA plot", sep = ""))
  
  resLFC <- lfcShrink(dds, coef = "condition_experimental_vs_control", type = "ashr")
  plotMA(resLFC, ylim = c(-10, 10), main = paste(conditionString, " LFC shrink MA plot"), sep = "")
  
  png(filename = paste(outputDir, conditionString, "_lfcshrink_MA.png", sep = ""))
  plot <- plotMA(resLFC, ylim = c(-10, 10), main = paste(conditionString, " LFC shrink MA plot"), sep = "")
  print(plot)
  dev.off()
  
  # make volcano plot
  png(filename = paste(outputDir, conditionString, "_volcano.png", sep = ""), height = 1080, width = 1100)
  plot <- EnhancedVolcano(res,
                          lab = res$row,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          title = conditionString,
                          pointSize = 1.5,
                          labSize = 6.0,
                          pCutoff = 0.05,
                          FCcutoff = 0.5,
                          colAlpha = 1,
                          legendLabels = c('NotSig', 'Log2FC', 'pvalue',
                                           'pvalue & Log2FC'),
                          legendPosition = 'right',
                          legendLabSize = 16,
                          legendIconSize = 6.0) +
    ggplot2::coord_cartesian(xlim = c(-10, 10))
  print(plot)
  dev.off()


  # Perform Gene Ontology on All Significant Results
  write.csv(res_sig, paste(outputDir, conditionString, "_all_sig.csv", sep = ""))
  significant <- read.csv(paste(outputDir, conditionString, "_all_sig.csv", sep = ""))
  go_one_all <- as.character(significant[, c(2)])
  
  gostres_one <- gost(query = go_one_all, organism = "dmelanogaster", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  
  p <- gostplot(gostres_one, capped = FALSE, interactive = FALSE)
  publish_gostplot(p, highlight_terms = gostres_one$result[c(1:3),], 
                   width = NA, height = NA, 
                   filename = paste(outputDir, conditionString, "_allGeneOnt.png", sep = ""))

  publish_gosttable(gostres_one, highlight_terms = gostres_one$result[c(1:150),],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"),
                    filename = paste(outputDir, conditionString, "_allGeneOntTable.pdf", sep = ""))
  
  # Perform Gene Ontology on All UPREGULATED Significant Results
  upsig <- read.csv(paste(outputDir, conditionString, "_upreg_genes.csv", sep = ""))
  go_two_up <- as.character(upsig[, c(2)])
  
  gostres_two <- gost(query = go_two_up, organism = "dmelanogaster", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  w <- gostplot(gostres_two, capped = FALSE, interactive = FALSE)
  publish_gostplot(w, highlight_terms = gostres_two$result[c(1:3),], 
                   width = NA, height = NA, 
                   filename = paste(outputDir, conditionString, "_UpReg_GO.png", sep = ""))
  
  publish_gosttable(gostres_two, highlight_terms = gostres_two$result[c(1:150),],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"),
                    filename = paste(outputDir, conditionString, "_UpReg_GO.pdf", sep = ""))
  
  # Perform Gene Ontology on DOWNREGULATED Significant Results
  downsig <- read.csv(paste(outputDir, conditionString, "_downreg_genes.csv", sep = ""))
  go_three_down <- as.character(downsig[, c(2)])
  
  gostres_three <- gost(query = go_three_down, organism = "dmelanogaster", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  t <- gostplot(gostres_three, capped = FALSE, interactive = FALSE)
  publish_gostplot(t, highlight_terms = gostres_three$result[c(1:3),], 
                   width = NA, height = NA, 
                   filename = paste(outputDir, conditionString, "_DownReg_GO.png", sep = ""))
  
  publish_gosttable(gostres_three, highlight_terms = gostres_three$result[c(1:150),],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"),
                    filename = paste(outputDir, conditionString, "_DownReg_GO.pdf", sep = ""))

  }

RunDESeq2GO(opt$rawcounts, opt$numcontrol, opt$numexperimental, opt$perturbation, opt$flybaseconverter, opt$out)

