# RNA-seq data: outlier removal and VarianceStabilizingTransformation (vst)

# Step 1: Remove outliers previously identified via PCA plot
# Step 2: Perform vst on filtered data, after removing outliers 
# Step 3: Further explore the data after vst via PCA plot, hclust, and heatmap (compare to previous PCA plot)

#=====================================================================================
# Directories 
#=====================================================================================
datadir = "/opt/projects/../RNA_seq/data"
plotdir = "/opt/projects/../RNA_seq/vst_PCA/"

#=====================================================================================
# Load packages
#=====================================================================================
library(tidyverse)
library(data.table)
library(DESeq2)
library(ggfortify)      # heatmap - option 1 - unsupervised 
library("RColorBrewer") # heatmap option 2 - hierarchical clustering   
library(pheatmap)       # heatmap option 2 - hierarchical clustering   

#=====================================================================================
# Load data
#=====================================================================================
setwd(datadir)

# design_matrix
design_matrix = read.table("design_matrix.tsv",sep="\t",header=T)
design_matrix = design_matrix[grepl("ILN", design_matrix$tissue),] # extract only ILN sample info
design_matrix = design_matrix[grepl(pattern = 'G1|G2|G3|G4', design_matrix$Group_number),] # extract only ILN samples from G1-G4
dim(design_matrix) # confirm there are 24 samples in ILN G1-G4   

# raw_gene_count 
raw_gene_count = read.table("raw_gene_count.tsv",sep="\t",header=T)
colnames(raw_gene_count) = gsub("read_count_","",colnames(raw_gene_count))

# gene_count_matrix
gene_count_matrix = round(as.matrix(raw_gene_count[,4:ncol(raw_gene_count)]),0)
row.names(gene_count_matrix) = (as.vector(raw_gene_count$gene_symbol))
colnames(gene_count_matrix)[which(colnames(gene_count_matrix) %in% design_matrix$Sample_name)]
gene_count_matrix = gene_count_matrix[,match(design_matrix$Sample_name,colnames(gene_count_matrix))] 

#=====================================================================================
# Remove outliers
#=====================================================================================
# remove sample 1 and 10 outliers, sample name: ZP_1_D7_ILNC , V_10_D7_ILNC (identified previously by PCA plots)
design_matrix = design_matrix[!(design_matrix$Sample_name) %in% c("ZP_1_D7_ILNC","V_10_D7_ILNC"),]
dim(design_matrix) # confirm there are 22 samples left in ILN G1-G4, after removing sample 1 and 10 outliers

#=====================================================================================
# Create dds object
#=====================================================================================
dds=DESeqDataSetFromMatrix(countData = gene_count_matrix,
                           colData = design_matrix,
                           design= ~ Group)

#=====================================================================================
# vst (VarianceStabilizingTransformation)
#=====================================================================================
vsd = varianceStabilizingTransformation(dds)
dists = dist(t(assay(vsd)))

#=====================================================================================
# PCA  
#=====================================================================================
setwd(plotdir)

pdf("ILN_vst_PCA_removingS1S10.pdf",paper='special')
plotPCA(vsd, intgroup=c("Group")) + geom_text(
  aes(label=colnames(vsd)),vjust=2,size = 2)
dev.off()

#=====================================================================================
# hclust (Hierarchical Clustering)
#=====================================================================================
setwd(plotdir)

pdf("ILN_vst_hclust_removingS10.pdf",width=20,height=15,paper='special')
plot(hclust(dists))
dev.off()

#=====================================================================================
# heatmap - option 1 - unsupervised   
#=====================================================================================
setwd(plotdir)

distsMatrix = as.matrix(dists)
pdf("ILN_vst_heatmap_removingS10.pdf",width=10,height=10,paper='special')
autoplot(distsMatrix,  colour = 'Species')
dev.off()

#=====================================================================================
# heatmap - option 2 - hierarchical clustering   
#=====================================================================================
setwd(plotdir)

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



