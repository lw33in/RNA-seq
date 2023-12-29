# RNA-seq downstream analysis using DESeq2 v1.42.0

# Note: DESeq2 requires raw integer read counts for performing accurate DGE analysis. 
# For gene quantification from Salmon, Sailfish, Kallisto, or RSEM, use the tximport package to import the count data.
# The normalized read counts should not be used in DESeq2 analysis. 
# DESeq2 internally normalizes the count data correcting for differences in the library sizes as sequencing depth influence the read counts (sample-specific effect). 
# DESeq2 does not consider gene length for normalization as gene length is constant for all samples (it may not have significant effect on DGE analysis). 

# Load packages ----------------------------------------------------------------
# install.packages(c("tximport", "readr", "dplyr", "DESeq2", "ggplot2", "GenomicFeatures",
                   # "org.Hs.eg.db","knitr","kableExtra","xlsx","pheatmap", "RColorBrewer", "genefilter"))
library(tximport)
library(readr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(knitr)
library(kableExtra)
library(xlsx)
library(pheatmap)
library(RColorBrewer)
library(genefilter)

# Set output options -----------------------------------------------------------
opts_chunk$set(
  fig.width = 5,
  fig.height = 4,
  warning = FALSE,
  message = FALSE
)

# Data preparation -------------------------------------------------------------
# Generate sample metadata
samples_data <- data.frame(run = c("ctrl_rep1", "ctrl_rep2", "ctrl_rep3", 
                                   "exp_rep1", "exp_rep2", "exp_rep3"), 
                           condition = rep(c("control", "exp"), each=3), 
                           stringsAsFactors = TRUE)
row.names(samples_data) <- samples_data$run 
samples_data

# Create paths to count files
data_dir <- "star_salmon/output"
count_file_paths <- file.path(data_dir, samples_data$run, "count")
names(count_file_paths) <- samples_data$run
count_file_paths

# Load TxDb object -------------------------------------------------------------
# Build TxDb object from reference genome GTF if it doesn't exist
gtf_file <- "reference_genome_metadata/gencode.v32.primary_assembly.annotation.gtf"
txdb_file <- "reference_genome_metadata/gencode.v32.annotation.sqlite"
if (file.exists(txdb_file)) {
  txdb <- loadDb(txdb_file)
} else {
  txdb <- makeTxDbFromGFF(gtf_file, 
                          organism = "Homo sapiens", 
                          dataSource = "https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38_2020A")
  saveDb(txdb, txdb_file)
}
# View TxDb object
txdb
genes(txdb)

# Load metadata linking transcripts to genes
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Import counts from each sample into txi & match geneIDs with tx names
txi <- tximport(count_file_paths, type = "salmon", tx2gene = tx2gene)
head(txi$counts)

# Count QC ---------------------------------------------------------------------
# Create a function that prints various statistics and plots for count data.
explore_count_data <- function(mat) {
  # View data and stats
  writeLines("Dimensions\n")
  print(dim(mat))
  writeLines("\nSummary\n")
  print(summary(mat))
  writeLines("\nSample Data\n")
  print(head(mat))
  
  # Plot Standard Deviation vs. Mean for count matrix.
  # Add 1 to mean and standard deviation to account for 0 mean or SD as both axes are log-scaled.
  data.frame(count_mean = rowMeans(mat) + 1, count_sd = rowSds(mat, useNames=FALSE) + 1) %>% 
    ggplot(aes(x = count_mean, y = count_sd)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    theme (
      plot.title = element_text(hjust = 0.5),
    ) +
    labs (
      x = "Mean",
      y = "Standard Deviation",
      title = "Count: SD vs. Mean"
    )
}

raw_counts_mat <- assay(dds)
explore_count_data(raw_counts_mat)

# Build DESeq2 object ----------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples_data,
                                design = ~ condition)
dds
# Set ctrl as the reference level/condition
dds$condition <- relevel(dds$condition, ref = "ctrl")

# Pre-filtering ----------------------------------------------------------------
# Remove genes with less than 5 reads per sample
keep <- rowSums(counts(dds) >= 5) >= dim(dds)[2]
table(keep)
dds <- dds[keep,]

explore_count_data(assay(dds))

# Differential Expression Analysis (DEA) ---------------------------------------
# Note: Two options to run DEA 
# 1) All-in-one-function: run `DESeq()`, Everything from normalization to linear modeling was carried out by the use of a single function.
# 2) Run individual functions in a step-wise manner, rather than a single call. 
#    - Step 1: Estimate size factors
#    - Step 2: Estimate gene-wise dispersion
#            - fit curve to gene-wise estimates
#            - shrink gene-wise dispersion estimates
#    - Step 3: GLM fit for each gene

# Option 1: Use `DESeq()`
dds <- DESeq(dds)

# # Option 2: Use individual functions
# # Step 1: Estimate size factors
# # Generates size factors for median of ratios normalization
# dds <- estimateSizeFactors(dds)
# # Check the size factors
# sizeFactors(dds)
# # Check the total number of reads for each sample
# colSums(counts(dds))
# # Check the total depth after normalization
# colSums(counts(dds, normalized=T))
# 
# # Step 2: Estimate gene-wise dispersion
# # The DESeq2 dispersion estimates are inversely related to the mean and directly related to variance. 
# # Based on this relationship, the dispersion is higher for small mean counts and lower for large mean counts. 
# # The dispersion estimates for genes with the same mean will differ only based on their variance. 
# # Therefore, the dispersion estimates reflect the variance in gene expression for a given mean value.
# dds <- estimateDispersionsGeneEst(dds) 
# dispersions(dds) <- mcols(dds)$dispGeneEst
# 
# # Calculates signficance of log2FoldChange
# dds <- nbinomWaldTest(dds, maxit=1000)

# Extract results --------------------------------------------------------------
# Results desciption: Result table from a DESeq analysis.
#   baseMean: Average of the count values.
#   log2FoldChange: Measurement of the extent to which a gene’s expression changed due to treatment.
#   lfcSE: The standard error for the log2FoldChange.
#   stat: Z-statistic = log2FoldChange/lfcSE

fdr_cutoff <- 0.05
res <- results(dds, alpha=fdr_cutoff, pAdjustMethod="fdr")
res

# Check outliers
summary(res)
# Check for NA adjusted p-values
table(is.na(res$padj))
# Check genes with adjusted p-value < my fdr cutoff 
sum(res$padj < fdr_cutoff, na.rm=TRUE)

# DEA results exploration ------------------------------------------------------
# Dispersion curve
plotDispEsts(dds)

# MA plot (unshrunken)
plotMA(res, ylim=c(-2,2))
# MA plot (shrunken)
#   Shrunken Plot - Reduce noise associated with LFC from low count genes.
#   Shrinkage of LFC estimates for visualization and ranking
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="exp_vs_ctrl", type="apeglm")
plotMA(res_LFC, ylim=c(-2,2))
# For weakly expressed genes, we have no chance of seeing differential expression, 
#    because the low read counts suffer from so high Poisson noise that any biological effect 
#    is drowned in the uncertainties from the read counting

# Histogram of pvalue
# Plot all
hist(res$pvalue, breaks=20, col="grey50", border="white")
# Plot without low count genes
hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")

# Gene annotation --------------------------------------------------------------
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=rownames(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, 
                     keys=row.names(res), 
                     column="ENTREZID", 
                     keytype="ENSEMBL",
                     multiVals="first")

# Remove genes with no gene symbol (i.e., NA).
res <- na.omit(res)
res

# Normalization ----------------------------------------------------------------
# Normalize count data via regularized logarithm
#    which transforms data on the log2 scale (with variance stabilization) with respect to library size.
normalized_dds <- rlog(dds, blind=TRUE)

# Export results
write.csv(as.data.frame(normalized_dds), file="normalized_dds.csv")


# Sample distance --------------------------------------------------------------
head(assay(normalized_dds))

# Option 1: 
# Manually provide the sampleDists to the clustering_distance argument of the pheatmap function. 
#  Otherwise the pheatmap function would assume that the matrix contains the data values themselves, 
#  and would calculate distances between the rows/columns of the distance matrix
sampleDists <- dist( t( assay(normalized_dds) ) )
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( normalized_dds$condition, normalized_dds$run, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Option 2:
# Use the Poisson Distance, implemented in the CRAN package PoiClaClu
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( normalized_dds$condition, normalized_dds$run, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)

# Option 3:
# Use plotPCA which comes with DESeq2
plotPCA(normalized_dds, intgroup = c("condition", "run"))

# Option 4:
# Build the PCA plot using ggplot2, leveraging the plotPCA function to return the data used for plotting
pca_data <- plotPCA(normalized_dds, intgroup = c("condition", "run"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
qplot(PC1, PC2, color=dex, shape=cell, data=pca_data) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# Get the normalized count matrix ----------------------------------------------
normalized_counts_mat <- assay(normalized_dds)
# Merge res & normalized_counts_mat by their row names (ENSEMBL gene IDs).
counts_res_df <- merge(res, # FDR cutoff applied
                       normalized_counts_mat, 
                       by.x='row.names', 
                       by.y="row.names", 
                       all.x = TRUE) %>% 
  `rownames<-` (.[,"Row.names"]) %>%
  dplyr::select(-Row.names)

# Gene clustering --------------------------------------------------------------
# Select the 20 genes with the highest variance across samples, with the rlog transformed counts
topVarGenes <- head(order(-rowVars(assay(normalized_counts_mat))),20)
# Center each genes' values across samples, and plot a heatmap
heatmap_mat <- assay(normalized_counts_mat)[ topVarGenes, ]
heatmap_mat <- heatmap_mat - rowMeans(heatmap_mat) # Calculate variance
heatmap_df <- as.data.frame(colData(normalized_counts_mat)[,c("run","condition")])
pheatmap(heatmap_mat, annotation_col=heatmap_df)

# Output sorting ---------------------------------------------------------------
# Sort by log2FoldChange in descending order
counts_res_df <- counts_res_df[-order(counts_res_df$log2FoldChange),]
head(counts_res_df)

# Count visualization for gene of interest -------------------------------------
# Check the gene with the lowest FDR
topGene <- rownames(res)[which.min(res$padj)] 

# Option 1:
plotCounts(dds, gene=topGene, intgroup=c("dex"))

# Option 2: 
topGene_data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
ggplot(topGene_data, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")

# DEG identification -----------------------------------------------------------
# Create a funciton `split_and_filter` that filters log2FoldChange in merged counts/results dataframe, counts_res_df.
# Filters (by LFC) and stores by up- and down-regulated genes in own dataframes.
split_and_filter <- function(my_log2_fold_change = -1) {
  # If LFC valid, filter by LFC
  if (my_log2_fold_change > 0) {
    counts_res_df <- counts_res_df[abs(counts_res_df$log2FoldChange) >= my_log2_fold_change, ]
  }
  
  # Up-regulated genes
  up_reg_df <- counts_res_df %>% 
    filter(log2FoldChange > 0)
  # Down-regulated genes
  down_reg_df <- counts_res_df %>%
    arrange(log2FoldChange) %>%
    filter(log2FoldChange < 0)
  
  # Print stats on de-regulated genes.
  print(paste0("Number of de-regulated genes: ", nrow(counts_res_df)))
  print(paste0("Number of up-regulated genes: ", nrow(up_reg_df)))
  print(paste0("Number of down-regulated genes: ", nrow(down_reg_df)))
  
  # Print additional stats for LFC range for up- and down-regulated genes.
  writeLines(paste0("\nUp-regulated Genes - LFC Range Sorted by Magnitude:\n", 
                    "\tMinimum LFC: ", min(up_reg_df$log2FoldChange), 
                    "\n\tMaximum LFC: ", max(up_reg_df$log2FoldChange)))
  writeLines(paste0("Down-regulated Genes - LFC Range Sorted by Magnitude:\n", 
                    "\tMinimum LFC: ", max(down_reg_df$log2FoldChange), 
                    "\n\tMaximum LFC: ", min(down_reg_df$log2FoldChange)))
  
  # Return full counts dataframe, 
  #   up-regulated genes dataframe and down-regulated genes dataframe as list.
  list(counts_res_df, up_reg_df, down_reg_df)
}

# Get up- and down-regulated genes from unfiltered counts_res_df
df_list <- split_and_filter()
up_reg_df <- df_list[[2]]
down_reg_df <- df_list[[3]]

# Print Up-regulated genes
up_reg_df %>% 
  dplyr::select(symbol, log2FoldChange)

# Print Down-regulated genes
down_reg_df %>% 
  dplyr::select(symbol, log2FoldChange)

# Filter for Fold Change of 1.5
df_list <- split_and_filter(log2(1.5))

counts_res_fc_1.5_df <- df_list[[1]]
up_reg_fc_1.5_df <- df_list[[2]]
down_reg_fc_1.5_df <- df_list[[3]]

# Filter for Fold Change of 2.0
df_list <- split_and_filter(log2(2))

counts_res_fc_2_df <- df_list[[1]]
up_reg_fc_2_df <- df_list[[2]]
down_reg_fc_2_df <- df_list[[3]]

# Export DEGs ------------------------------------------------------------------
# Save Unfiltered, Up-regulated & Down-regulated genes to their own sheets and
#   Unfiltered - Up-regulated & Down-regulated to their own sheets within single excel file for client.
xlsx_file_name <- "DESeq2_DEG.xlsx"
write.xlsx(c("SHEET KEY",
             "Sheet2 -> Up-regulated genes sorted in descending order by magnitude of log2FoldChange", 
             "Sheet3 -> Down-regulated genes sorted in descending order by magnitude of log2FoldChange"), 
           file=xlsx_file_name, 
           sheetName="README", 
           row.names=FALSE,
           col.names=FALSE)
write.xlsx(up_reg_df %>% 
             tibble::rownames_to_column("ENSEMBL_Gene_ID") %>% 
             dplyr::select(ENSEMBL_Gene_ID, symbol, log2FoldChange), 
           file=xlsx_file_name, 
           sheetName="Sheet2 = Up-regulated Genes", 
           row.names=FALSE,
           append = TRUE,)
write.xlsx(down_reg_df %>% 
             tibble::rownames_to_column("ENSEMBL_Gene_ID") %>% 
             dplyr::select(ENSEMBL_Gene_ID, symbol, log2FoldChange), 
           file=xlsx_file_name, 
           sheetName="Sheet3 = Down-regulated Genes", 
           row.names=FALSE,
           append=TRUE)

# [OPTIONAL] Gene of interest exploration --------------------------------------
genes_of_interest <- "^ALK"
# genes_of_interest <- "RAS|RAF|MEK|ERK|SOS|GRB|HSP|PAK|SRC|EGFR|TGF|MYC|JUN|FOS|ETS|P13K|TIAM|AKT|
#                             RAC|RHO|BAD|BCL|MTOR|ERBB|FGFR|PDGFR|GPCR|AURKA|BRAF|MAPK|MAP2K|MET|MYC|SMARCA4|
#                             IDH1|TP53|TSC|SPRY|SPRED|SHC|RPS6KA|ROS1|ROCK|RHOA|RAL|PTEN|PRKAG|PRKAB|PIK3R|PIK3C|
#                             PDPK1|NF1|MLST8|MDM2|KSR|JUN|FOS|ERBB2|EIF4EBP1|DUSP|CDK|CCND1|ARHGAP35|ALK|GAP|RPS6"

counts_res_df %>%
  filter(grepl(genes_of_interest, symbol, ignore.case = TRUE)) %>%
  arrange(desc(log2FoldChange)) %>%
  knitr::kable(align = 'c') %>%
  kable_styling("striped") %>% 
  column_spec (1:(ncol(counts_res_df) + 1), border_left = TRUE, border_right = TRUE) %>%
  scroll_box(width = "100%", height = "500px")



sessionInfo()

