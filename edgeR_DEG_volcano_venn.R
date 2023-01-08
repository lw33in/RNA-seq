# RNA-seq DEG analysis using edgeR Bioconductor
# Data visualization using ggplot2, ggfortify, ggrepel, ggpubr, and VennDiagram packages

# Public GEO data: GSE145160 (PMID: 32406922)
# Experimental Samples: ZFX mutant construct transfected DKO cells (293T cells lacking both ZFX and ZNF711 via CRISPR knockout)
# Controls: 1) EGV transfected DKO cells, 2) DKO cells

# Data platform: Illumina NovaSeq 6000 (Homo sapiens), paired-end 150 bp sequencing 
# RNA-seq results were aligned to GENCODE v19 and reads were counted using STAR
# DEG cutoffs: absolute fold change > 1.5, p.value < 0.05

#=====================================================================================
# Directories 
#=====================================================================================
datadir = "/Users/../mutantRNA/"
resultdir = "/Users/../mutantRNA/results/"

#=====================================================================================
# Load packages
#=====================================================================================
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR") 
library(edgeR)
biocLite("baySeq")
library(baySeq)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggfortify")
library(ggfortify)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("ggpubr")
library(ggpubr)
# install.packages('VennDiagram')
library(VennDiagram)

#=====================================================================================
# Select gene counts for the experimental and control groups
#=====================================================================================
setwd(datadir)

# quantify to annotation model (E/M) - gene counts
gene_counts = read.delim("/Users/../gene_counts_GeneAnno.txt",header =T)
dim(gene_counts) # 27875 genes   27 samples
gene_counts[1:9,] 
rownames (gene_counts) = gene_counts$Gene 
gene_counts = gene_counts [,-1]
colnames(gene_counts)
gene_counts = gene_counts [,c(1,2,3,      # ZF_FL_rep1, rep2, rep3
                              4,5,6,      # ZF_1_8_rep1, rep2, rep3
                              7,8,9,      # ZF_9_13_rep1, rep2, rep3
                              10,11,12,   # ZF_11_13_rep1, rep2, rep3
                              13,14,15,   # ZF_zero_rep1, rep2, rep3
                              17,18,19,   # EGV_rep1, rep2, rep3
                              21,22,23) ] # RNA_EGV_rep1, rep2, rep3

#=====================================================================================
# Save gene counts for GEO submission 
#=====================================================================================
setwd(resultdir)
ZFX_ZF_FL_rep1 =  data.frame ( gene_counts [ , 1] )
rownames(USC384_ZFX_ZF_FL_rep1) = rownames(gene_counts)
colnames(USC384_ZFX_ZF_FL_rep1) = "ZFX_ZF_FL_rep1"
write.table(USC384_ZFX_ZF_FL_rep1, "./forGEO/USC_ZFX_ZF_FL_rep1_ReadsPerGene.txt",quote = FALSE,sep="\t",col.names = T,row.names = T)

#=====================================================================================
# Create a DGEList object
#=====================================================================================
zf_groups = c(rep("ZF_FL", 3), rep("ZF_1_8", 3), rep("ZF_9_13", 3), rep("ZF_9_11", 3), rep("ZF_zero", 3), rep("EGV", 3))
dgList = DGEList(counts=gene_counts, genes=gene_counts$Gene,group=factor(zf_groups))
dgList$samples # confirm 27875 genes   21 samples

#=====================================================================================
# Filter low expressed genes
#=====================================================================================
countsPerMillion = cpm(dgList)
summary(countsPerMillion) # explore the numeric data

# retain only those genes that are represented at least 1cpm reads in at least two samples (cpm=counts per million)
countCheck = countsPerMillion > 1 # at least 1cpm reads
countsPerMillion [ rownames(countsPerMillion) == "AMER3" , ] # use an example gene for sanity check
countCheck[1:20,]

# retain genes with at least two samples countsPerMillion > 1
keep = which(rowSums(countCheck) >= 2) # at least 2 samples have countsPerMillion > 1
dgList = dgList[keep,]
summary(cpm(dgList)) # compare this to the original summary
dim(dgList) # 16850 expressed genes, after filtering
dgList
length ( rownames ( dgList$counts ) )

# save 16850 genes post filtering
setwd(resultdir)
write.table(rownames (dgList$counts), "GeneCounts_postFildering_16850geneNames.txt", quote = FALSE, sep="\t", col.names=F, row.names=F)
write.table(dgList$counts ,"GeneCounts_postFildering_16850geneCounts.txt", quote = FALSE, sep="\t")

#=====================================================================================
# PCA plot
#=====================================================================================
gene_counts_forPCA = counts
rownames(gene_counts_forPCA) = NULL
df = gene_counts_forPCA[c(1, 2, 3, 4)]
autoplot(prcomp(gene_counts_forPCA))

#=====================================================================================
# Normalization
#=====================================================================================
# edgeR implements the trimmed mean of M-values (TMM) method. It is also based on the hypothesis that most genes are not DE.
?calcNormFactors # Calculate normalization factors to scale the raw library sizes.
dgList = calcNormFactors(dgList, method="TMM")  # method=c("TMM","RLE","upperquartile","none") 
# test different normalization methods 

#=====================================================================================
# MDS plot
#=====================================================================================
# examine inter-sample relationships by producing a plot based on mutlidimensional scaling.
par(mar=c(5,5,5,5))
plotMDS(dgList) # plotMDS(dgList_JL, method="bcv",col=as.numeric(dgList$samples$group))
pdf("../MDSplot.pdf",width=4.5,height=4,paper='special')
# biological coefficient of variation (BCV) is the square root of the dispersion parameter in the negative binomial model.
plotMDS(dgList, method="bcv",col=as.numeric(dgList$samples$group))
dev.off()

#=====================================================================================
# Fit a regression model in edgeR
#=====================================================================================
# design matrix, which describes the setup of the experiment 
design.mat = model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) = levels(dgList$samples$group)

# Fitting a model in edgeR takes several steps:
# 1) fit the common dispersion. 
# 2) fit a trended model (if you do not fit a trend, the default is to use the common dispersion as a trend). 
# 3) fit the tagwise dispersion which is a function of this model.

d1 = estimateCommonDisp(dgList, verbose=T) # Disp = 0.02729 , BCV = 0.1652  (TMM)
names(d1)
d1 =  estimateTagwiseDisp(d1)
names(d1)
# plotBCV() plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(d1)
# above plot shows a single estimate for the coefficient of variation is a bad model since tagwise dispersion does not follow the model 
# but decreases as the counts per million (CPM) increases.

# GLM estimates of dispersion
# In addition to the common and tagwise disperson, we can also estimate a generalized linear model (glm) fit using edgeR.
design.mat = model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) = levels(dgList$samples$group)
d2 = estimateGLMCommonDisp(dgList,design.mat)
d2 = estimateGLMTrendedDisp(d2,design.mat, method="auto") # can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.

d2 = estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
# biological coefficient of variation (BCV) is the square root of the dispersion parameter in the negative binomial model.
# plotBCV() plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.

# GLM testing for differential expression
design.mat
fit = glmFit(d2, design.mat)

#=====================================================================================
# Identify DEGs
#=====================================================================================
# (   1 ,        1,             1,             1,           1,          1,             1     )  
# ( EGV ,     ZF_1_8 ,     ZF_11_13 ,       ZF_9_11 ,    ZF_9_13 ,    ZF_FL ,       ZF_zero)
lrt_ZF_FL_vs_EGV     = glmLRT(fit, contrast=c(-1,0,0,0,0,1,0)) # "-1" is control  
lrt_ZF_1_8_vs_EGV    = glmLRT(fit, contrast=c(-1,1,0,0,0,0,0)) 
lrt_ZF_9_13_vs_EGV   = glmLRT(fit, contrast=c(-1,0,0,0,1,0,0))  
lrt_ZF_9_11_vs_EGV   = glmLRT(fit, contrast=c(-1,0,0,1,0,0,0))  
lrt_ZF_11_13_vs_EGV  = glmLRT(fit, contrast=c(-1,0,1,0,0,0,0)) 
lrt_ZF_zero_vs_EGV   = glmLRT(fit, contrast=c(-1,0,0,0,0,0,1)) 
lrt_ZF_FL_vs_ZF_zero =glmLRT(fit, contrast=c(0,0,0,0,0,1,-1)) 

# lfc = numeric value giving the desired absolute minimum log-fold-change = 0.176 (FC1.5)
de2 = decideTestsDGE(lrt_ZF_FL_vs_EGV, adjust.method="BH", p.value = 0.05) 
# BH = Benjamini and Hochberg's algorithm. By default, Benjamini and Hochberg's algorithm is used to control the false discovery rate (FDR).
summary(de2)
de2tags12 = rownames(d2)[as.logical(de2)]
# The function plotSmear generates a plot of the tagwise log-fold-changes against log-cpm (analogous to an MA-plot for microarray data). 
# DE tags are highlighted on the plot:
plotSmear(lrt_ZF_FL_vs_EGV, de.tags=de2tags12)
abline(h = c(-0.5849625, 0.5849625), col = "blue") # log2(1.5) = 0.5849625
dev.off()

#=====================================================================================
# Explore DEGs
#=====================================================================================
# test results for the n most significant tags are conveniently displayed by the topTags() function 
ZF_FL_vs_EGV = topTags(lrt_ZF_FL_vs_EGV, n=Inf, adjust.method="BH", sort.by="logFC") 

# Up-regulated genes
ZF_FL_vs_EGV_upTMM = subset (ZF_FL_vs_EGV, ( ZF_FL_vs_EGV$logFC > log2(1.5) ) & (ZF_FL_vs_EGV$FDR < 0.05) ) # FDR is adj-p-val
dim(ZF_FL_vs_EGV_upTMM)  
ZF_FL_vs_EGV_upTMM$genes = rownames(ZF_FL_vs_EGV_upTMM)
"TRAK2" %in% rownames(ZF_FL_vs_EGV_upTMM) # check certain genes 

# Down-regulated genes
ZF_FL_vs_EGV_downTMM = subset (ZF_FL_vs_EGV, ( ZF_FL_vs_EGV$logFC < - log2(1.5) ) & (ZF_FL_vs_EGV$FDR < 0.05)) # FDR is adj-p-val
dim(ZF_FL_vs_EGV_downTMM)  
ZF_FL_vs_EGV_downTMM$genes = rownames(ZF_FL_vs_EGV_downTMM)
"FBLN2" %in% rownames(ZF_FL_vs_EGV_downTMM) # check certain genes 

#=====================================================================================
# Volcano plot
#=====================================================================================
# Highlight genes that have an absolute fold change > 1.5 and a p-value < Bonferroni cut-off
ZF_11_13_vs_EGV$genes = rownames(ZF_11_13_vs_EGV)
ZF_11_13_vs_EGV$threshold = as.factor(abs(ZF_11_13_vs_EGV$logFC) > log2(1.5) & ZF_11_13_vs_EGV$FDR < 0.05)
# add names to the top 10 genes, label those dots with the gene name on the Volcano plot using geom_text_repel()
ZF_11_13_vs_EGV_volcano = ZF_11_13_vs_EGV[order(ZF_11_13_vs_EGV$logFC), ] # sort by ordered padj
# create a column to indicate which genes to label
ZF_11_13_vs_EGV_volcano$geneLabel = ""
ZF_11_13_vs_EGV_volcano$geneLabel[ZFX_ZF_11_13_vs_ZFX_EGV_volcano$genes =="LONRF2"] = "LONRF2"
ZF_11_13_vs_EGV_volcano[ZFX_ZF_11_13_vs_ZFX_EGV_volcano$genes == "LONRF2" , ]
View(ZF_11_13_vs_EGV_volcano)

# volcano plot
pdf("./ZF_11_13_vs_EGV_volcano.pdf",width=4.5,height=4,paper='special')
par(mar=c(3,3,2,1), mgp=c(1.2,0,0), tck=-.01)
cbPalette <- c("#999999", "#1e1aed") 
p=ggplot(data=ZFX_ZF_11_13_vs_ZFX_EGV_volcano, 
         aes(x=logFC, y=-log10(FDR), colour=threshold, 
             label = geneLabel 
         )) +
  geom_point(alpha=0.8, size=0.5) +
  labs(legend.position = "none") +
  xlim(c(-5.5,5.5) )+ 
  ylim(c(0, 52) )+
  xlab("log2(Fold Change)") + ylab("-log10(Adj-pVal)")  +
  theme_bw() +
  scale_colour_manual(values=cbPalette) +
  geom_text_repel(color = 'black',size = 5,point.padding = unit(0.2, "lines") , box.padding = unit(2, "lines"),)

p+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("xy.text", size = 12)

dev.off()

#=====================================================================================
# Venn diagram 
#=====================================================================================
venn.plot=venn.diagram(x=list("FL_up"= ZF_FL_vs_EGV_upTMM$genes,"zero_up"= ZF_zero_vs_EGV_upTMM$genes),
filename=c("./VennDiagram.tiff"),
euler.d=T,
#col="transparent",
fill=c(
  "#D8BFD8", # purple
  "#90EE90" # green 
  ),   
lty = "blank", # lines of circles
fontface = "bold",
fontfamily = "sans",
cat.col=c( 
  "#9932CC", # purple
  "#228B22" #  green
), 
cat.fontface = "bold", 
cex             = 1,  # size of the area label / sizes of numbers
cat.fontfamily = "sans",
cat.cex         = 1, 
cat.default.pos = "outer",
margin= 0.5, # size of the circles
label.col = "black", # the numbers 
scaled = TRUE,
ext.text = FALSE
)


