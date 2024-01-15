rm(list=ls())
library(scater)
library(stringr)
library(pheatmap)
library(gtools)
library(scran)

setwd("~/Desktop/Medulloblastoma")
outDir <- file.path("~/Desktop/Medulloblastoma/Output/")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
pathway_file <- "~/Desktop/Medulloblastoma/Data/KEGG_metabolism.gmt"
source("~/Desktop/Medulloblastoma/Scripts/utils.R")


# Loading the data
selected_sce <- readRDS(file.path(outDir, "selected_sce.rds"))
selected_sce$tumor_subpopulation <- gsub("-", "_", selected_sce$tumor_subpopulation)
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]
subgroups <- unique(selected_sce$tumor_subpopulation)

#2.subgroup cells
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)

pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),subgroup=character())

for (t in subgroups){
  each_metabolic_sce <- selected_metabolic_sce[,selected_metabolic_sce$tumor_subpopulation==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"tpm")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- each_metabolic_tpm
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             subgroup=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  
  ####
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  
  runGSEA_preRank(pre_rank_matrix,pathway_file,t)
  #get the result
  result_dir <- list.files(path="~/Desktop/Medulloblastoma/Output/preRankResults",pattern = paste0("^",t,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}

#remove pvalue <0.05 pathways
min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
select_pathways <- names(min_pval)[(min_pval<=0.05)]
select_enrich_data_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]

#converto pvalue to -log10
pvals <- select_enrich_data_df$PVAL
pvals[pvals<=0] = 1e-10
select_enrich_data_df$PVAL <- -log10(pvals)

#sort
pathway_pv_sum <- by(select_enrich_data_df$PVAL,select_enrich_data_df$y,FUN=sum)
pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum,decreasing = T)]

# top 10
pathway_order <- pathway_order[1:10]
select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$y %in% pathway_order,]

select_enrich_data_df$x <- factor(select_enrich_data_df$x, levels = mixedsort(subgroups))
select_enrich_data_df$y <- factor(select_enrich_data_df$y,levels = pathway_order)
select_enrich_data_df$x <- gsub("_", "-", select_enrich_data_df$x)

write.table(select_enrich_data_df,file=file.path(outDir,"enrich_path_subpop.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#select_enrich_data_df <- read.csv("Output/enrich_path_subpop.txt",row.names=1, sep="\t", header = T)

##buble plot
p <- ggplot(select_enrich_data_df, aes(x = x, y = y, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  #ggtitle("pathway heterogeneity") +
  labs(x = NULL, y = NULL,
       size = "-log10 pvalue", color = "NES") +
  scale_size(range = c(0, 10)) +
  scale_color_gradient( low = "white", high = "red") +
  #scale_color_gradient2(low="red",mid="white",high="blue",midpoint = 1) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(.6, "cm"),
        legend.text = element_text(colour="black",size=12),
        axis.line = element_line(size=.8, colour = "black"),
        #panel.grid.major = element_line(colour = "#d3d3d3"),
        #panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 1.9),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 14,angle=45,hjust=1,vjust=0.9),
        axis.text.y=element_text(colour="black", size = 14)) +
  theme(plot.margin = unit(rep(1,4),"lines"))
p
##plot variance
pc_plotdata$subgroup <- factor(pc_plotdata$subgroup,levels=mixedsort(subgroups))
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~subgroup,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))
p
unlink(file.path(outDir,"preRankResults"),recursive=T)
unlink(file.path(outDir,"prerank.rnk"))

date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)

# For subgroup ==========================================================
selected_metabolic_sce$subgroup <- gsub("/", "_", selected_metabolic_sce$subgroup)
subgroups <- unique(selected_metabolic_sce$subgroup)

#2.subgroup cells
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)

pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),subgroup=character())

for (t in subgroups){
  each_metabolic_sce <- selected_metabolic_sce[,selected_metabolic_sce$subgroup==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"tpm")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- each_metabolic_tpm
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             subgroup=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  
  ####
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  
  runGSEA_preRank(pre_rank_matrix,pathway_file,t)
  #get the result
  result_dir <- list.files(path="~/Desktop/Medulloblastoma/Output/preRankResults",pattern = paste0("^",t,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}

#remove pvalue <0.05 pathways
min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
select_pathways <- names(min_pval)[(min_pval<=0.05)]
select_enrich_data_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]

#converto pvalue to -log10
pvals <- select_enrich_data_df$PVAL
pvals[pvals<=0] = 1e-10
select_enrich_data_df$PVAL <- -log10(pvals)

#sort
pathway_pv_sum <- by(select_enrich_data_df$PVAL,select_enrich_data_df$y,FUN=sum)
pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum,decreasing = T)]

# top 10
pathway_order <- pathway_order[1:10]
select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$y %in% pathway_order,]

select_enrich_data_df$x <- factor(select_enrich_data_df$x, levels = mixedsort(subgroups))
select_enrich_data_df$y <- factor(select_enrich_data_df$y,levels = pathway_order)
select_enrich_data_df$x <- gsub("_", "-", select_enrich_data_df$x)

write.table(select_enrich_data_df,file=file.path(outDir,"enrich_path_subgroup.txt"),row.names=T,col.names=T,quote=F,sep="\t")
#select_enrich_data_df <- read.csv("enrich_path_subpop.txt",row.names=1, sep="\t", header = T)

#bubble plot
ggplot(select_enrich_data_df, aes(x = x, y = y, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  labs(x = NULL, y = NULL,
       size = "-log10 pvalue", color = "NES") +
  scale_size(range = c(0, 8)) +
  scale_color_gradient( low = "white", high = "red") +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(.6, "cm"),
        legend.text = element_text(colour="black",size=8),
        axis.line = element_line(size=.8, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 1.5),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10,angle=45,hjust=1,vjust=0.9),
        axis.text.y=element_text(colour="black", size = 10)) +
  theme(plot.margin = unit(rep(1,4),"lines"))

##plot variance
pc_plotdata$subgroup <- factor(pc_plotdata$subgroup,levels=mixedsort(subgroups))

ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~subgroup,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))

unlink(file.path(outDir,"preRankResults"),recursive=T)
unlink(file.path(outDir,"prerank.rnk"))

date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)

 











