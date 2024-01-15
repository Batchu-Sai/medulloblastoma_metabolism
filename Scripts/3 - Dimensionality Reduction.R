rm(list = ls())
library(scater)
library(reshape2)
library(umap)
library(parallel)
library(cowplot)
library(Rtsne)
options(stringsAsFactors=FALSE)

setwd("~/Desktop/Medulloblastoma")
outDir <- file.path("~/Desktop/Medulloblastoma/Output/")

# Loading the tumor data -------------------------------------------------------
selected_sce <- readRDS("Output/selected_sce.rds")
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]

# UMAP / t-SNE ========================================================================
detectCores()
set.seed(12345)

umap_metabolic <- umap::umap(t(assay(selected_metabolic_sce,"expr")),
                        n_threads = 8,
                        verbose = T)

tsne_metabolic <- Rtsne::Rtsne(t(assay(selected_metabolic_sce,"expr")),
                             num_threads = 8,
                             verbose = T)

# Get output ------------------------------------------------------------------
umap_metabolic_out <- data.frame(x=umap_metabolic$layout[,1],y=umap_metabolic$layout[,2],
                                 Sample = colData(selected_metabolic_sce)$geo_sample_id,
                                 Subpopulation = colData(selected_metabolic_sce)$tumor_subpopulation,
                                 Subgroup = colData(selected_metabolic_sce)$subgroup)

tsne_metabolic_out <- data.frame(x=tsne_metabolic$Y[,1],y=tsne_metabolic$Y[,2],
                                 Sample = colData(selected_metabolic_sce)$geo_sample_id,
                                 Subpopulation = colData(selected_metabolic_sce)$tumor_subpopulation,
                                 Subgroup = colData(selected_metabolic_sce)$subgroup)

# Write output to file ---------------------------------------------------------

write.table(umap_metabolic_out,file=file.path(outDir, "umap_metabolic_output.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(tsne_metabolic_out,file=file.path(outDir, "tsne_metabolic_output.txt"),row.names=T,col.names=T,quote=F,sep="\t")


umap_metabolic_out <- read.csv(file.path(outDir, "umap_metabolic_output.txt"),row.names=1, sep="\t", header = T)
tsne_metabolic_out <- read.csv(file.path(outDir, "tsne_metabolic_output.txt"),row.names=1, sep="\t", header = T)

# Plot ------------------------------------------------------------------------

p1 <-
  ggplot(umap_metabolic_out) + 
  geom_point(aes(x, y, colour = Subpopulation), size = 0.001) +
  labs(x = "UMAP 1",y = "UMAP 2") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(6, 6, 6, 6))+
  labs(colour = "Subpopulation")

p2 <-
  ggplot(umap_metabolic_out) + 
  geom_point(aes(x, y, colour = Subgroup), size = 0.001) +
  labs(x = "UMAP 1",y = "UMAP 2") +
  scale_color_manual(values = c("#000000",'darkred',"#009292","#db6d00","#ffb6db"))+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(6, 6, 6, 6))+
  labs(colour = "Subgroup")+ 
  guides(fill = guide_legend(ncol = 1))


cowplot::plot_grid(p2, p1, ncol = 2)


p3 <-
  ggplot(tsne_metabolic_out) + 
  geom_point(aes(x, y, colour = Subpopulation), size = 0.001) +
  labs(x = "tSNE 1",y = "tSNE 2") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(6, 6, 6, 6))+
  labs(colour = "Subpopulation") +
  guides(fill = guide_legend(ncol = 1), 
       color = guide_legend(override.aes = list(size = 5)))
p4 <-
  ggplot(tsne_metabolic_out) + 
  geom_point(aes(x, y, colour = Subgroup), size = 0.001) +
  labs(x = "tSNE 1",y = "tSNE 2") +
  scale_color_manual(values = c("#000000",'darkred',"#009292","#db6d00","#ffb6db"))+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(6, 6, 6, 6))+
  labs(colour = "Subgroup")+ 
  guides(fill = guide_legend(ncol = 1), 
         color = guide_legend(override.aes = list(size = 5)))

cowplot::plot_grid(p4, p3, ncol = 2)

 
