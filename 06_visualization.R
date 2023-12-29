# RNA-seq downstream analysis - Data visualization 

# Venn diagram -----------------------------------------------------------------
# Input: list
# install.packages('VennDiagram')
library(VennDiagram)

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

dev.off()

# Volcano plot -----------------------------------------------------------------
# Input: dataframe
# install.packages(c("ggplot2", "ggfortify", "ggrepel"))
library(ggplot2)
library(ggrepel)
# library(ggfortify)
# library(ggpubr)

exp1_volcano <- exp1_DEG
# Highlight significant DEGs
exp1_volcano$threshold <- as.factor(abs(exp1_DEG$logFC) > log2(1.5) & exp1_DEG$FDR < 0.05)
# Add names to the top 10 DEGs using geom_text_repel()
exp1_volcano <- exp1_volcano[order(exp1_volcano$logFC), ] # order by FDR
exp1_volcano$geneLabel <- ""
exp1_volcano$geneLabel[1:10] <- "TRUE"
View(exp1_volcano)
# Make object
exp1_volcano <- exp1_volcano[order(exp1_volcano$logFC), ] 
# Setup
pdf("/plots/exp1_volcano.pdf",width=3.6,height=3.1,paper='special')
par(mar=c(3,3,2,1), mgp=c(1.2,0,0), tck=-.01)
# Plot
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

# Chow ruskey -----------------------------------------------------------------
# Input: list 
# install.packages(c("devtools", "Vennerable"))
library(devtools)
library(Vennerable)

vignette("Venn") # The function based documentation is sketchy. A better guide is the Venn vignette which you can see with

# Create input list
x=list("list1"=exp1$genes,
       "list2"=exp2$genes,
       "Oncogenes"=onc$genes  )
# Plot
ChowRuskey = Venn(x)
pdf("/plots/ChowRuskey.pdf",width=4.5,height=4,paper='special')
plot(ChowRuskey,
     type = "ChowRuskey",
     show = list(SetLabels = F, 
                 DarkMatter = F, 
                 FaceText = ""
     ))

dev.off()

# Heatmap -----------------------------------------------------------------
# Input: matrix
# install.packages(c("limma", "matlab", "Biobase","matlab","RColorBrewer"))
library(limma)
library(matlab)
library(Biobase)
library(matlab)
library(RColorBrewer)

# Create color palette:
Colors=rev(brewer.pal(3,"Spectral"))
Colors=colorRampPalette(Colors)(length(palette.breaks) - 1) 
# Set up 
source("/source/heatmap.3.R")
pdf("/plots/heatmap.pdf",width=3.6,height=3.1,paper='special')
par(cex.main=1)
# Plot
heatmap.3(matrix_forHeatmap, 
          #dendrogram = "col",
          col=jet.colors(1000000), 
          margins=c(10,10),
          sepcolor = "white",
          sepwidth = c(0.05, 0.05),
          #distfun = dist,
          Rowv=F, Colv=F,
          symbreaks=FALSE, key=T, 
          symkey=F, scale="none",
          density.info="none", trace="none", 
          #breaks=palette.breaks, 
          #ColSideColors=Colors,
          main="exp1", 
          cexRow=1, cexCol=1.8,
          keysize=1, KeyValueName="4 groups")

dev.off()

sessionInfo()
