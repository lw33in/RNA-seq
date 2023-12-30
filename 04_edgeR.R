# RNA-seq downstream analysis using edgeR

# Input: count quantification

# Load packages ----------------------------------------------------------------
# install.packages("edgeR", "baySeq", "VennDiagram", "ggplot2", "ggrepel") 
library(edgeR)
library(baySeq)
library(VennDiagram)
library(ggplot2)
library(ggrepel)

# Set output options -----------------------------------------------------------
# opts_chunk$set(
#   fig.width = 5,
#   fig.height = 4,
#   warning = FALSE,
#   message = FALSE
# )

# Data preparation -------------------------------------------------------------
# Generate sample metadata
dko_groups <- rep(c("exp1", "exp2","exp3","ctrl"), each=3)
dko_counts <- read.delim("/output/geneCounts.txt",header =T)
dgList <- DGEList(counts=dko_counts, genes=dko_counts$Gene,group=factor(dko_groups))
dgList

# Gene-level filtering ---------------------------------------------------------
# Retain only those genes that are represented at least 1cpm reads in at least two samples.
# Note: cpm = count per million 
countsPerMillion <- cpm(dgList)
summary(countsPerMillion) 
countCheck <- countsPerMillion > 1 # at least 1cpm reads
keep <- which(rowSums(countCheck) >= 2) # at least two samples 
dgList <- dgList[keep,]
summary(cpm(dgList))
length(rownames(dgList$counts))

# Save filtered data------------------------------------------------------------
write.table(rownames(dgList$counts),"/rnaseq_downstream/post_filtering_exp_gene_names.txt",
            quote = FALSE, sep="\t", col.names=F, row.names=F)
write.table(dgList$counts,"/rnaseq_downstream/post_filtering_exp_gene_counts.txt",
            quote = FALSE, sep="\t")

# [Project Specific] Overlap with TF binding events-----------------------------
# How many ZFX bound genes are expressed?
exp_gene <- rownames(dgList$counts)
colname(exp_gene) <- "gene"
ZFXsummit_gene <- read.delim("/ZFX/summit_gene",header = F)
ZFX_exp_gene <- exp_gene[exp_gene$gene %in% ZFXsummit_gene$gene, ] %>%
  length()

# Normalization ----------------------------------------------------------------
dgList <- calcNormFactors(dgList, method="TMM")  
# Trimmed mean of M-values (TMM) method: hypothesis that most genes are not DE.
# `calcNormFactors` calculates normalization factors to scale the raw library sizes.
# Other normalization methods: c("TMM","RLE","upperquartile","none")

# Data exploration -------------------------------------------------------------
plotMDS(dgList)
plotMDS(dgList, method="bcv",col=as.numeric(dgList$samples$group))
# Plot mutlidimensional scaling (MDS) to examine inter-sample relationships.

# Model setup ------------------------------------------------------------------
design.mat <- model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)

# Dispersion estimation --------------------------------------------------------
# Step 1: fit the common dispersion --> estimate tagwise dispersions
d1 <- estimateCommonDisp(dgList, verbose=F) %>% 
  estimateTagwiseDisp()
plotBCV(d1)
# Step 2: estimate a generalized linear model (glm) fit
d2 <- estimateGLMCommonDisp(dgList, design.mat)
d2 <- estimateGLMTrendedDisp(d2, design.mat, method="auto")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
# Biological coefficient of variation (BCV) is the square root of the dispersion parameter in the negative binomial model.

# Differential expression analysis ---------------------------------------------
design.mat
fit <- glmFit(d2, design.mat)
# (1 ,1 , 1, 1)  = (ctrl, exp1, exp2, exp3)
lrt_exp1 <- glmLRT(fit, contrast=c(-1,1,0,0)) 
lrt_exp2 <- glmLRT(fit, contrast=c(-1,0,1,0)) 
lrt_exp3 <- glmLRT(fit, contrast=c(-1,0,0,1))  

de2 <- decideTestsDGE(lrt_exp1, adjust.method="BH", p.value = 0.05, lfc = 0.5849625) 
# Cutoffs: FDR < 0.05, abs|lfc|> 0.5849625 (FC 1.5)
# BH = Benjamini and Hochberg's algorithm is used to control the false discovery rate (FDR).
summary(de2)

de2tags12 <- rownames(d2)[as.logical(de2)] # so that DE tags are highlighted on the plot
plotSmear(lrt_Seven_vs_EGV, de.tags=de2tags12)
abline(h = c(-0.5849625, 0.5849625), col = "blue") 
# `plotSmear` generates a plot of the tagwise log-fold-changes against log-cpm (analogous to an MA-plot for microarray data). 
dev.off()

# Identify significant DEGs
exp1_DEG <- topTags(lrt_exp1, n=Inf, adjust.method="BH", sort.by="logFC")
exp1_DEG_up <- subset (exp1_DEG, exp1_DEG$logFC > log2(1.5) & exp1_DEG$FDR < 0.05)
exp1_DEG_down <- subset (exp1_DEG, exp1_DEG$logFC < -log2(1.5) & exp1_DEG$FDR < 0.05)

# Volcano plot -----------------------------------------------------------------
exp1_volcano <- exp1_DEG
# Highlight significant DEGs
exp1_volcano$threshold <- as.factor(abs(exp1_DEG$logFC) > log2(1.5) & exp1_DEG$FDR < 0.05)
# Add names to the top 10 DEGs using geom_text_repel()
exp1_volcano <- exp1_volcano[order(exp1_volcano$logFC), ] # order by FDR
exp1_volcano$geneLabel <- ""
exp1_volcano$geneLabel[1:10] <- "TRUE"
View(exp1_volcano)
# Plot object
exp1_volcano <- exp1_volcano[order(exp1_volcano$logFC), ] 
# Setup
pdf("/plots/exp1_volcano.pdf",width=3.6,height=3.1,paper='special')
par(mar=c(3,3,2,1), mgp=c(1.2,0,0), tck=-.01)
# Volcano plot
cbPalette <- c("#999999", "#CC0000")
ggplot(data=exp1_volcano, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  #geom_text_repel(aes(x = logFC, y = -log10(PValue), label = ifelse(geneLabel == T, rownames(pdko25_volcano),""))) +
  geom_point(alpha=0.8, size=0.3) +
  labs(legend.position = "none") +
  xlim(c(-9, 9)) + 
  ylim(c(0, 300)) +
  xlab("log2(Fold Change)") + ylab("-log10(Adj-pVal)") +
  theme_bw() +
  scale_colour_manual(values=cbPalette)

dev.off()

# [Project Specific] DEG exploration -------------------------------------------
# Check gene of interest 
NSD1 <- exp1_DEG [exp1_DEG$genes =="NSD1",]
# Venn diagram
venn.plot <- venn.diagram(x=list("exp1_DEG_up"= exp1_DEG_up$genes,
                                 "exp2_DEG_up"= exp2_DEG_up$genes,
                                 "exp2_DEG_up" = exp2_DEG_up$genes),
                          filename=c("/plots/venn.tiff"),
                          euler.d=T,
                          fill=c("#D8BFD8","#F5DEB3","#87CEFA"),   
                          lty = "blank", # lines of circles
                          #lwd = 25, # width of the circle's circumference 
                          fontface = "bold",
                          fontfamily = "sans",
                          cat.col = c("#9932CC","#CD853F","#4169E1"), # name color
                          cat.fontface = "bold", 
                          cex = 1,  # size of the area label / sizes of numbers
                          cat.fontfamily = "sans",
                          cat.cex = 1, 
                          cat.default.pos = "outer",
                          margin= 0.5, # size of the circles
                          label.col = "black", # the numbers 
                          scaled = TRUE,
                          ext.text = FALSE)


sessionInfo()
