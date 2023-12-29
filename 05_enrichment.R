# RNA-seq downstream analysis - Enrichment analysis 

# Load packages ----------------------------------------------------------------
# install.packages(c("clusterProfiler", "ggnewscale", "ggridges", "DOSE", 
#  "ggplot2", "enrichplot","reactome.db", "AnnotationDbi","knitr","kableExtra")
library(clusterProfiler)
library(ggnewscale)
library(ggridges)
library(DOSE)
library(ggplot2)
library(enrichplot)
library(reactome.db)
library(AnnotationDbi)
library(knitr)
library(kableExtra)

# GSEA  ------------------------------------------------------------------------
# Load data --------------------------------------------------------------------
# Continuation of 03_DESeq2.R 
# # Gene annotation section (line 196))
# columns(org.Hs.eg.db)
# res$symbol <- mapIds(org.Hs.eg.db,
#                      keys=rownames(res),
#                      column="SYMBOL",
#                      keytype="ENSEMBL",
#                      multiVals="first")
# res$entrez <- mapIds(org.Hs.eg.db, 
#                      keys=row.names(res), 
#                      column="ENTREZID", 
#                      keytype="ENSEMBL",
#                      multiVals="first")
# 
# # Remove genes with no gene symbol (i.e., NA).
# res <- na.omit(res)
res

# Incorporate Reactome database ------------------------------------------------
# Subset the results table, normalized_dds, to only those genes for which the Reactome database has data 
#  whose Entrez ID we find in the respective key column of reactome.db and for which the DESeq2 test gave an adjusted p value that was not NA
res_forGSEA <- res[ res$entrez %in% keys( reactome.db, "ENTREZID" ) & !is.na( res$padj ) , ]
head(res2)
# Create a table with the mapping from Entrez IDs to Reactome Path IDs
reactomeTable <- AnnotationDbi::select( reactome.db, 
                                        keys=as.character(res_forGSEA$entrez), 
                                        keytype="ENTREZID", 
                                        columns=c("ENTREZID","REACTOMEID") )
head(reactomeTable)

# Transform reactomeTable into an incidence matrix (incm) ----------------------
# This is a Boolean matrix with one row for each Reactome Path and one column for each unique gene in res2, 
#   which tells us which genes are members of which Reactome Paths.

# Select only unique entrez IDs
ft2mat <- function(x, y) {
  ux = unique(x)
  uy = unique(y)
  im = matrix(FALSE, nrow=length(ux), ncol=length(uy), dimnames=list(ux, uy))
  im[ cbind(x, y) ] = TRUE
return(im)
}

incm <- with(reactomeTable, ft2mat(REACTOMEID, ENTREZID))

# Create incm
incm <- do.call( rbind, with(reactomeTable, tapply( 
  ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x ) ))
colnames(incm) <- res2$entrez
str(incm)

# Remove rows corresponding to Reactome Paths with < 20 or > 80 assigned genes
within <- function(x, lower, upper) (x>=lower & x<=upper)
incm <- incm[ within(rowSums(incm), lower=20, upper=80), ]

# Test whether the genes in a Reactome Path behave in a special way 
# Create a function to calculate a number of statistics, including a t-statistic to see 
#   whether the average of the genes’ log2 fold change values in the gene set is different from zero.
testCategory <- function( reactomeID ) {
  isMember <- incm[ reactomeID, ]
  data.frame( 
    reactomeID  = reactomeID,
    numGenes    = sum( isMember ),
    avgLFC      = mean( res_forGSEA$log2FoldChange[isMember] ),
    sdLFC       = sd( res_forGSEA$log2FoldChange[isMember] ),
    zValue      = mean( res_forGSEA$log2FoldChange[isMember] ) /sd( res_forGSEA$log2FoldChange[isMember] ),
    strength    = sum( res_forGSEA$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
    pvalue      = t.test( res_forGSEA$log2FoldChange[ isMember ] )$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]],
    stringsAsFactors = FALSE ) }

# Test the `testCategory` function with a Reactome Path ID
testCategory("446193")

# Use the `testCategory` function for all Reactome IDs and put results in a dataframe
reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )
reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )

# Create a list of a list of Reactome Paths which are significantly differentially expressed in the comparison of exp vs control, 
#  sorted according to sign and strength of the signal
reactomeResult_signif <- reactomeResult[ reactomeResult$padjust < 0.05, ]
head( reactomeResult_signif[ order(-reactomeResult_signif$strength), ] )

# GO  --------------------------------------------------------------------------
# Load data --------------------------------------------------------------------
# Continuation of 03_DESeq2.R 
# # Output sorting section (line 287)
# # Sort by log2FoldChange in descending order
# counts_res_df <- counts_res_df[-order(counts_res_df$log2FoldChange),]
head(counts_res_df)

# DEG identification section (line 351)
# # Filter for Fold Change of 1.5
# df_list <- split_and_filter(log2(1.5))
# 
# counts_res_fc_1.5_df <- df_list[[1]]
# up_reg_fc_1.5_df <- df_list[[2]]
# down_reg_fc_1.5_df <- df_list[[3]]
up_reg_fc_1.5_df
down_reg_fc_1.5_df

# Prepare data for GO ----------------------------------------------------------
# Create a function `get_LFC_vector` to retrieve log2FoldChange, row names are ENSEMBL gene IDs
#  create LFC vector
get_LFC_vector <- function (count_df) {
  gene_vec <- count_df %>% 
    dplyr::select(log2FoldChange) %>% 
    pull()
  names(gene_vec) <- rownames(count_df)
  # Sort in descending order for clusterProfiler
  gene_vec <- sort(gene_vec, decreasing = TRUE)
}
  
# GO ---------------------------------------------------------------------------
# Create a function `perform_go` to retrieve log2FoldChange, row names are ENSEMBL gene IDs
perform_go <- function (count_df) {
  gene_vec <- get_LFC_vector(count_df)
  set.seed(7788)
  gse_go <- gseGO(geneList = gene_vec, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, # Minimal size of each geneSet for analyzing
                  maxGSSize = 800, # Maximal size of genes annotated for testing
                  pvalueCutoff = 1, # Get all possible results by setting pvalueCutoff to 1
                  eps = 0, # This parameter sets the boundary for calculating the p value.
                  verbose = TRUE, 
                  seed = TRUE,
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "fdr")
  # Convert ENSEMBL IDs to gene symbols and sort by descending Normalized Enrichment Score (NES)
  # NES = Enrichment score after normalization across all gene sets.
  gse_go <- setReadable(gse_go, 'org.Hs.eg.db') %>% 
    arrange(desc(NES))
  # Print stats
  print(paste0("Number of GO Terms: ", nrow(gse_go)))
  gse_go
}

# Up-regulated genes, 1.5 FC
gse_go_upDEG <- perform_go(up_reg_fc_1.5_df)

# Down-regulated genes, 1.5 FC
gse_go_downDEG <- perform_go(down_reg_fc_1.5_df)

# Extract GO terms  ------------------------------------------------------------
# # All Terms
# gse_go %>%
#   knitr::kable(align = 'c', caption = "All Terms") %>%
#   kable_styling("striped") %>% 
#   column_spec (1:13, border_left = TRUE, border_right = TRUE) %>%
#   column_spec(13, width_min = "20em") %>%
#   scroll_box(width = "100%", height = "500px")

# Create a function `table_top50_go_terms` to generate tables for the top 50 up- and down-regulated terms by NES.
table_top50_go_terms <- function(my_gse_go, num_entries = 50) {
  # Top50 up-regulated terms
  table_1 <- as.data.frame(my_gse_go) %>%
    slice_head(n = num_entries) %>%
    knitr::kable(align = 'c', caption = "Top 50 Up-regulated Terms by NES") %>%
    kable_styling("striped") %>% 
    column_spec (1:(ncol(my_gse_go) + 1), border_left = TRUE, border_right = TRUE) %>%
    column_spec(ncol(my_gse_go) + 1, width_min = "20em") %>%
    scroll_box(width = "100%", height = "500px")
  
  # Top50 down-regulated terms
  table_2 <- as.data.frame(my_gse_go) %>%
    slice_tail(n = num_entries) %>%
    arrange(NES) %>%
    knitr::kable(align = 'c', caption = "Top 50 Down-regulated Terms by NES") %>%
    kable_styling("striped") %>% 
    column_spec (1:(ncol(my_gse_go) + 1), border_left = TRUE, border_right = TRUE) %>%
    column_spec(ncol(my_gse_go) + 1, width_min = "20em") %>%
    scroll_box(width = "100%", height = "500px")
  
  list(table_1, table_2)
}

# Up-regulated genes, 1.5 FC --> top GO terms
top_go_terms_upDEG <- table_top50_go_terms(gse_go_upDEG)
top_go_terms_upDEG[[1]] # top50 up terms
top_go_terms_upDEG[[2]] # top50 down terms

# Down-regulated genes, 1.5 FC --> top GO terms
top_go_terms_downDEG <- table_top50_go_terms(gse_go_downDEG)
top_go_terms_downDEG[[1]] # top50 up terms
top_go_terms_downDEG[[2]] # top50 down terms


sessionInfo()

