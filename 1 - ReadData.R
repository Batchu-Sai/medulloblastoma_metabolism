rm(list = ls())
setwd("~/Desktop/Medulloblastoma")
source("Scripts/Utils.R")
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
outdir <- "./Output"
if(!dir.exists(outdir)) dir.create(outdir)

# Read the data

raw_count_file <- "Data/exprMatrix.tsv.gz"
metadata_file <- "Data/GSE155446_human_cell_metadata.csv"
metadata <- read.table(metadata_file,head=T,sep=",",quote=NULL,stringsAsFactors=F)
metadata <- metadata[!metadata$tumor_subpopulation %in% c('GP4-X1', 'GP4-X2' ,'GP4-X3'),]
metadata$otherPop <- paste0(metadata$subgroup,"-",metadata$coarse_cell_type)

# keep GP4 samples consistent with what original authors analyzed in their paper (Figure 3): 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8804892/ 

count_data <- read.table(raw_count_file,head=T,sep="\t",quote=NULL,stringsAsFactors=F)
colnames(count_data) <- gsub("X", "", colnames(count_data))
rownames(count_data) <- count_data$gene
count_data$gene <- NULL
count_data <- count_data[,colnames(count_data)%in% metadata$cell]

# Convert raw counts to tpm 

all_gene_lengths <- read.delim("Data/gene_length.txt", col.names= c('Gene', 'Length'), sep = "\t")

genelen <- all_gene_lengths[all_gene_lengths$Gene %in% rownames(count_data),]
count_data <- count_data[rownames(count_data) %in% genelen$Gene, ]

genelen  <- genelen %>%
  dplyr::slice(match(rownames(count_data), Gene))

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

tpms <- apply(count_data, 2, function(x) tpm(x, genelen$Length))
colSums(tpms[,1:5]) # check

# mark the metabolic genes

pathways <- gmtPathways("Data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(count_data)),row.names = rownames(count_data))

row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE

table(row_data$metabolic)

# FALSE  TRUE 
# 15845  1498 

# build scater object
metadata <- metadata[metadata$cell %in% colnames(count_data),]
sce <- SingleCellExperiment(
  assays = list(tpm = data.matrix(tpms),
                expr= data.matrix(log2(count_data+1))),
  colData = metadata,
  rowData = row_data
)


# filtering the cells.

#malignant cells
tumor_sce <- sce[,sce$coarse_cell_type == "malignant"]

#select tumor cells with atleast 50
tumor_sample_stats <- table(tumor_sce$tumor_subpopulation)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats>=50]
selected_tumor_sce <- tumor_sce[,tumor_sce$tumor_subpopulation %in% tumor_sample_select]

# Save as R data files
saveRDS(selected_tumor_sce,file.path("Output/","selected_sce.rds"))
