# RNA-seq gene filtering: remove genes with low read counts

# Input data: gene count matrix
# Filtering criteria
# Filter 1: >50 % of samples should express gene (read count > 0)
# Filter 2: >25 % of samples should have read count >= 5
# Filter 3: >=50 % of samples in each group should express gene (read count > 0)

# Note: These are stringent filtering criteria. An alternative gene filtering 
#       criteria can be: remmove genes with less than 5 reads per sample

# Load data --------------------------------------------------------------------
dds < -DESeqDataSetFromMatrix(countData = gene_count_matrix, colData = design_matrix, design= ~ Group)
# dds dim: 19824 genes  144 samples

# Filter genes with low read counts --------------------------------------------
# Filter 1: >50 % of samples should express gene (read count > 0)
keep1 = rowSums(counts(dds) > 0) > 72
dds = dds[keep1,] #dim: 16078 144

# Filter 2: >25 % of samples should have read count >= 5
keep2 = rowSums(counts(dds) >= 5) > 36
dds = dds[keep2,] # dim: 15275 144

# Filter 3: >=50 % of samples in each group should express gene (read count > 0)
dds_matrix = as.matrix(gene_count_matrix[rownames(dds),],
                       colData = design_matrix,
                       design= ~ Group)
dim(dds_matrix)
keep3 = which( apply( t(rowsum(t(dds_matrix>0)+0, group=dds$Group_number)), MARGIN=1, FUN=function(x) { any(x<3) }) == FALSE) 
dds = dds[keep3,] # dim: 15275 144 

sessionInfo()