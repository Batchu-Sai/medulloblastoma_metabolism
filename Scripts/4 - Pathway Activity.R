rm(list=ls())
library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)


setwd("~/Desktop/Medulloblastoma")
outDir <- file.path("~/Desktop/Medulloblastoma/Output/")
pathway_file <- "~/Desktop/Medulloblastoma/Data/KEGG_metabolism.gmt"
source("~/Desktop/Medulloblastoma/Scripts/utils.R")
norm_rds_file <- file.path(outDir, "Deconvolution_tpm.rds")

#1. Loading the data
selected_sce <- readRDS(file.path(outDir, "selected_sce.rds"))

pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)

all_subgroups <- as.vector(selected_sce$subgroup)
subgroups <- unique(all_subgroups)

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_sce)[rowData(selected_sce)$metabolic])
norm_tpm <- readRDS(norm_rds_file)

set.seed(123)

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = list(pathway_names,subgroups))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = list(pathway_names,subgroups))

###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = (list(pathway_names, subgroups)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  
  #remove genes which are zeros in any subgroup to avoid extreme ratio value
  keep <- colnames(mean_exp_eachsubgroup)[colAlls(mean_exp_eachsubgroup>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  ratio_exp_eachsubgroup <- t(mean_exp_eachsubgroup) / colMeans(mean_exp_eachsubgroup)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachsubgroup,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachsubgroup,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  ratio_exp_eachsubgroup <- t(mean_exp_eachsubgroup) / colMeans(mean_exp_eachsubgroup)
  mean_exp_pathway <- apply(ratio_exp_eachsubgroup,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[subgroups]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[subgroups]
  
  ##shuffle 1000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(subgroups,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_subgroups_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachsubgroup_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:1000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_subgroups_list <- lapply(times,function(x) sample(all_subgroups)) 
  names(shuffle_subgroups_list) <- times
  mean_exp_eachsubgroup_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachsubgroup_list <- lapply(times,function(x) mean_exp_eachsubgroup_list[[x]] / rowMeans(mean_exp_eachsubgroup_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(subgroups),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- subgroups
  for(c in subgroups){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 1000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 1000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}

all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]

# heatmap
dat <- mean_expression_shuffle
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}

sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)

pheatmap(dat[sort_row,sort_column],
         cluster_cols = F,
         cluster_rows = F,
         na_col = "white", 
         fontsize = 9)

write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")


#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,1.5),breaks=0:1.5,labels=0:1.5)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=3,color="green")+
  scale_fill_manual(values = c("#000000",'darkred',"#009292","#db6d00","#ffb6db"))+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 16,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 16),
        axis.line=element_line(size=0.9,color="black"),
        axis.ticks = element_line(colour = "black",size=1),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(2.5, "mm"))


# ==============================================================================
all_subgroups <- as.vector(selected_sce$tumor_subpopulation)
subgroups <- unique(all_subgroups)

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle_subpopulation <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = list(pathway_names,subgroups))
mean_expression_noshuffle_subpopulation <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = list(pathway_names,subgroups))

###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(subgroups),dimnames = (list(pathway_names, subgroups)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next
  
  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  
  #remove genes which are zeros in any subgroup to avoid extreme ratio value
  keep <- colnames(mean_exp_eachsubgroup)[colAlls(mean_exp_eachsubgroup>0.001)]
  
  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
  
  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  ratio_exp_eachsubgroup <- t(mean_exp_eachsubgroup) / colMeans(mean_exp_eachsubgroup)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachsubgroup,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachsubgroup,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachsubgroup <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_subgroups, mean))
  ratio_exp_eachsubgroup <- t(mean_exp_eachsubgroup) / colMeans(mean_exp_eachsubgroup)
  mean_exp_pathway <- apply(ratio_exp_eachsubgroup,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle_subpopulation[p, ] <-  mean_exp_pathway[subgroups]
  mean_expression_noshuffle_subpopulation[p, ] <-  mean_exp_pathway[subgroups]
  
  ##shuffle 1000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(subgroups,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_subgroups_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachsubgroup_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:1000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_subgroups_list <- lapply(times,function(x) sample(all_subgroups)) 
  names(shuffle_subgroups_list) <- times
  mean_exp_eachsubgroup_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachsubgroup_list <- lapply(times,function(x) mean_exp_eachsubgroup_list[[x]] / rowMeans(mean_exp_eachsubgroup_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(subgroups),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- subgroups
  for(c in subgroups){
    if(is.na(mean_expression_shuffle_subpopulation[p,c])) next
    if(mean_expression_shuffle_subpopulation[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle_subpopulation[p,c]) / 1000 
    }else if(mean_expression_shuffle_subpopulation[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle_subpopulation[p,c]) / 1000
    }
    if(pval>0.01) mean_expression_shuffle_subpopulation[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}

all_NA <- rowAlls(is.na(mean_expression_shuffle_subpopulation))
mean_expression_shuffle_subpopulation <- mean_expression_shuffle_subpopulation[!all_NA,]

# heatmap
dat <- mean_expression_shuffle_subpopulation
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}

sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
mybreaks <- c(
  seq(0.1, 1, length.out=33),
  seq(1.01, 1.1, length.out=33),
  seq(1.11, 2,length.out=34)
) 

pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,na_col = "white", breaks = mybreaks, angle_col = 45, fontsize = 10,)

write.table(mean_expression_shuffle_subpopulation,file=file.path(outDir,"KEGGpathway_activity_shuffle_subpop.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue_subpop.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle_subpopulation)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,2),breaks=0:1.5,labels=0:1.5)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=3,color="black")+
  #scale_fill_manual(values = c("#000000",'darkred',"#009292","#db6d00","#ffb6db"))+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 16,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 16),
        axis.line=element_line(size=0.9,color="black"),
        axis.ticks = element_line(colour = "black",size=1),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(2.5, "mm"))



