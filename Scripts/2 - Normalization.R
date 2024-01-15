rm(list = ls())
library(scater)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(edgeR)
library(scran)
options(stringsAsFactors=FALSE)

setwd("~/Desktop/Medulloblastoma")
outDir <- file.path("~/Desktop/Medulloblastoma/Output/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1 load the data
selected_sce <- readRDS(file.path(outDir, "selected_sce.rds"))
cell_types <- unique(selected_sce$tumor_subpopulation)

## choose the genes with higher expression and low-dropout rate 
dropout_cutoff <- 0.5
gene_select_mat <- matrix(FALSE,nrow=nrow(selected_sce),
                          ncol=length(cell_types),
                          dimnames = list(rownames(selected_sce),cell_types))
for(c in cell_types){
  each_sce <- selected_sce[,selected_sce$tumor_subpopulation == c]
  each_exp <- assay(each_sce,"expr")
  dropout_rate <- apply(each_exp,1, function(x) sum(x>0)/ncol(each_exp))
  select <- dropout_rate >= dropout_cutoff
  gene_select_mat[select,c] <- TRUE
}
print("the number of genes selected:")
print(sum(rowSums(gene_select_mat) >= length(cell_types)))
low_dropout_genes <- rownames(gene_select_mat)[rowSums(gene_select_mat) >= length(cell_types)]

#prepare the gene length file
#convert the TPM to count scale
all_gene_lengths <- read.table("~/Desktop/Medulloblastoma/Data/gene_length.txt",sep="\t",header=F,row.names=1)
tmp <- intersect(rownames(all_gene_lengths),rownames(selected_sce))
if (length(tmp) != nrow(selected_sce)){
  warning("check the length file")
}
genelen <- all_gene_lengths[rownames(selected_sce),]
genelen <- as.numeric(as.vector(genelen))

selected_impute_tpm <- tpm(selected_sce)
selected_impute_counts <- sweep(selected_impute_tpm, 1, genelen, FUN = "*")

# scran, deconvolution
scran.sf <- scran::calculateSumFactors(selected_impute_counts[low_dropout_genes,],clusters=selected_sce$tumor_subpopulation)
summary(scran.sf)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / scran.sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm+1)

#save
saveRDS(selected_impute_tpm_norm, file.path(outDir,"Deconvolution_tpm.rds"))


###Evaluation of the normalization method: check ratio distribution
all_cell_types <- as.vector(selected_sce$tumor_subpopulation)
norm_tpm <- selected_impute_tpm_norm
low_dropout_genes_tpm <- norm_tpm[low_dropout_genes, ]
low_dropout_genes_tpm_mean <- apply(low_dropout_genes_tpm, 1, function(x) by(x, all_cell_types, mean))
low_dropout_genes_tpm_ratio <- t(low_dropout_genes_tpm_mean) / colMeans(low_dropout_genes_tpm_mean)
dat <- reshape2::melt(low_dropout_genes_tpm_ratio)
p <- 
  ggplot(dat,aes(x=Var2,y=value)) +
  geom_boxplot(outlier.alpha=0.1)+ theme_classic() + 
  ylab("expression ratio") + xlab("") +    
  theme(axis.text.x = element_text(angle=45,hjust=1))
p
  
ggsave(file.path(outDir,paste0("Deconvolution_ratio_distribution.pdf")),p,width=3.5,height=2.5)
