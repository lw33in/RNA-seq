# RNA-seq data: Outlier removal and VarianceStabilizingTransformation (vst)

# Step 1: Remove outliers previously identified via PCA plot
# Step 2: Perform vst on filtered data, after removing outliers 
# Step 3: Further explore the data after vst via PCA plot, hclust, and heatmap (compare to previous PCA plot)

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(data.table)
library(DESeq2)
library(ggfortify)      # heatmap option 1 - unsupervised 
library(RColorBrewer)   # heatmap option 2 - hierarchical clustering   
library(pheatmap)       # heatmap option 2 - hierarchical clustering   

# Load data --------------------------------------------------------------------
datadir <- "/opt/projects/../RNA_seq/data"
plotdir <- "/opt/projects/../RNA_seq/vst_PCA/"
set(datadie)
# Load design_matrix
design_matrix <- read.table("design_matrix.tsv",sep="\t",header=T)
design_matrix <- design_matrix[grepl("ILN", design_matrix$tissue),] %>%
  design_matrix[grepl(pattern = 'G1|G2|G3|G4', design_matrix$Group_number),] 
dim(design_matrix) 

# Load raw_gene_count 
raw_gene_count = read.table("raw_gene_count.tsv",sep="\t",header=T)
colnames(raw_gene_count) = gsub("read_count_","",colnames(raw_gene_count))

# Load gene_count_matrix
gene_count_matrix = round(as.matrix(raw_gene_count[,4:ncol(raw_gene_count)]),0)
row.names(gene_count_matrix) = (as.vector(raw_gene_count$gene_symbol))
colnames(gene_count_matrix)[which(colnames(gene_count_matrix) %in% design_matrix$Sample_name)]
gene_count_matrix = gene_count_matrix[,match(design_matrix$Sample_name,colnames(gene_count_matrix))] 

# Remove outliers --------------------------------------------------------------
# Remove sample 1 and 10 outliers, sample name: ZP_1_D7_ILNC , V_10_D7_ILNC (identified previously by PCA plots)
design_matrix = design_matrix[!(design_matrix$Sample_name) %in% c("ZP_1_D7_ILNC","V_10_D7_ILNC"),]

# Create dds object ------------------------------------------------------------
dds=DESeqDataSetFromMatrix(countData = gene_count_matrix,
                           colData = design_matrix,
                           design= ~ Group)

# VarianceStabilizingTransformation (vst) --------------------------------------
vsd = varianceStabilizingTransformation(dds)
dists = dist(t(assay(vsd)))

# PCA --------------------------------------------------------------------------
setwd(plotdir)
pdf("ILN_vst_PCA.pdf",paper='special')
plotPCA(vsd, intgroup=c("Group")) + geom_text(
  aes(label=colnames(vsd)),vjust=2,size = 2)

dev.off()

# hclust (Hierarchical Clustering) ---------------------------------------------
pdf("ILN_vst_hclust.pdf",width=20,height=15,paper='special')
plot(hclust(dists))

dev.off()

# heatmap - option 1 - unsupervised --------------------------------------------
distsMatrix = as.matrix(dists)
pdf("ILN_vst_heatmap.pdf",width=10,height=10,paper='special')
autoplot(distsMatrix,  colour = 'Species')

dev.off()

# heatmap - option 2 - hierarchical clustering ---------------------------------
distsMatrix = as.matrix(dists)
rownames(distsMatrix) = colnames(vsd)
colnames(distsMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
groupnumber= data.frame(vsd$Group_number)
rownames(groupnumber) = colnames(vsd)
colnames(groupnumber) = "Group"

pdf("ILN_vst_heatmap_hierarchical clustering_removingS1S10.pdf")
pheatmap(distsMatrix,
         clustering_distance_rows=dists,
         clustering_distance_cols=dists,
         annotation_row = groupnumber,
         col=colors)

dev.off()

sessionInfo()
