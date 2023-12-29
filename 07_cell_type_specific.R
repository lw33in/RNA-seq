# RNA-seq downstream analysis - Identify and explore cell-type specific genes

# Note: This RNA-seq data in rats. The objective is to identify lymph node specific genes.

# Load packages ----------------------------------------------------------------
# install.packages(c("tidyverse", "data.table", "DESeq2", "gplots", 
#   "RColorBrewer", "pheatmap")
library(tidyverse)
library(data.table)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(pheatmap)

# Load data --------------------------------------------------------------------
setwd("/rnaseq/data/")
# Load FPKM
gene_exp = read.table("gene_expression.tsv",sep="\t",header=T) 
# Data cleaning
rownames(gene_exp) = gene_exp$gene_symbol
fpkm = gene_exp[grepl("fpkm_",colnames(gene_exp))]
dim(fpkm) #  19824   144
colnames(fpkm) = gsub("fpkm_", "", colnames(fpkm)) # remove 'fpkm_'from sample names

# Load design_matrix
design_matrix <- read.table("design_matrix.tsv",sep="\t",header=T)
design_matrix <- design_matrix[grepl("ILN", design_matrix$tissue),]
design_matrix <- design_matrix[grepl(pattern = 'G1|G2|G3|G4', design_matrix$Group_number),]
# Reorder samples so that sample names match design_matrix
fpkm <- fpkm[,match(design_matrix$Sample_name,colnames(fpkm))] 
head(fpkm)

# Convert HPA to rat gene names ------------------------------------------------
setwd("/genelist/")
lymph_hpa = read.table('HPA_LymphoidTissue_RNA_TissueEnriched_201genes.tsv',sep="\t",header=T) 
rat_human_mapping = read.table('rat_human_mapping.id_ac.tsv',sep="\t",header=T,na.strings=c("","NA")) 
# Match 
lymph_rat = as.data.frame(rat_human_mapping$RAT_GENE_SYMBOL[ rat_human_mapping$HUMAN_GENE_SYMBOL %in% lymph_hpa$Gene])
lymph_rat = as.data.frame(lymph_rat[!apply(lymph_rat == "", 1, all),])
dim(lymph_rat) 

# Load tissue specific HPA converted rat genes 
setwd("/genelist/RatGenes/")
lymph_rat = read.table("LymphoidTissue_Rat.tsv" ,sep="\t",header=F) 

# Cell-type specific gene identification ---------------------------------------
# Look at lymph node data
fpkm_ILN <- fpkm[,grepl("D7_ILNC", colnames(fpkm))]
colnames(fpkm_ILN)  
# Identify significant rat genes in lymph nodes
lymph_genes_inILN <- fpkm_ILN[match(lymph_rat$V1, rownames(fpkm_ILN)),]
lymph_genes_inILN = na.omit(lymph_genes_inILN)
# Check number of matched rat genes  in 24 G1-G4 ILN samples (before removing low read genes)
dim(lymph_genes_inILN) 
# Create matrix for heatmap.2
lymph_genes_inILN = as.matrix(lymph_genes_inILN) 
rownames(lymph_genes_inILN)

# Heatmap ----------------------------------------------------------------------
setwd("/plots/")
pdf("lymph_node_genes.pdf",paper='special')
heatmap.2(lymph_genes_inILN,
          col= brewer.pal(11,"RdBu"),
          scale="row", trace="none", Colv =NA, Rowv = NA,
          cexRow=0.5, cexCol=0.5, key=TRUE, keysize=1.5,
          ColSideColors = c(rep("firebrick", 6),rep("lime green", 6),rep("dodger blue", 6),rep("purple",6))
)
legend(y=1.1, x=.25, xpd=TRUE, 
       legend = c("ILN Group 1", "ILN Group 2", "ILN Group 3", "ILN Group 4"), 
       col = c("firebrick", "lime green", "dodger blue", "purple"),  lty= 1,lwd = 3 ,cex=.5           
)

dev.off()


sessionInfo()

