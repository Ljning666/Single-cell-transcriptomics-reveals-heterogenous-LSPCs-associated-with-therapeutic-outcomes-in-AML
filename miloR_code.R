rm(list = ls())
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(stringr)
library(Seurat)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)
library(Matrix)
#1. 数据准备
setwd("F:/AML单细胞/new_result/result2/malignant/lspc")
lspc_scrna<-readRDS("lspc_scrna_harmony.rds")
lspc_scrna <- FindClusters(lspc_scrna, resolution = 0.1)
lspc_scrna<- JoinLayers(lspc_scrna)
table(lspc_scrna@meta.data$seurat_clusters)
meta<-lspc_scrna@meta.data
meta$new_type<-meta$type
meta[which(meta$new_type=="refractory"),26]<-"Refractory/Relapse"
meta[which(meta$new_type=="relapse"),26]<-"Refractory/Relapse"
lspc_scrna<-AddMetaData(lspc_scrna,metadata = meta)
lspc_scrna@meta.data$cluster_type = paste0("C",lspc_scrna@meta.data$seurat_clusters)
table(meta$cluster_type)
sce <- as.SingleCellExperiment(lspc_scrna)
#Visualize the data
plotReducedDim(sce, colour_by="seurat_clusters", dimred = "UMAP")
#2. 构建领域，对领域中的细胞计数
#d：用于 KNN 细化的降维数，类似于聚类降维时选择的dims。
#k：分布峰值在 50 到 100 之间，平均邻域大小超过 5 x N_samples。
#Create a Milo object
scmilo <- Milo(sce)
#Construct KNN graph
scmilo <- buildGraph(scmilo, k = 30, d = 30,reduced.dim = "PCA")
#Defining representative neighbourhoods on the KNN graph
set.seed(10)
scmilo <- makeNhoods(scmilo, prop = 0.05,
                     k = 30, d = 30, 
                     refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(scmilo)
#Counting cells in neighbourhoods
scmilo <- countCells(scmilo, 
                     meta.data = as.data.frame(colData(scmilo)), 
                     sample = "patient")
head(nhoodCounts(scmilo))
#3. 定义分组，差异测试
#Defining experimental design
sc_design <- data.frame(colData(scmilo))[,c("patient", "new_type")]

sc_design <- distinct(sc_design)
rownames(sc_design) <- sc_design$patient

#Computing neighbourhood connectivity
scmilo <- calcNhoodDistance(scmilo, d = 30, reduced.dim = "PCA")

#Testing
results <- testNhoods(scmilo, design = ~ new_type, 
                      design.df = sc_design, reduced.dim="PCA")
head(results)

results %>%
  arrange(SpatialFDR) %>%
  head() 
##4.可视化
#Inspecting DA testing results
ggplot(results, aes(PValue)) + geom_histogram(bins=50)
ggplot(results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) 

scmilo <- buildNhoodGraph(scmilo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(scmilo, dimred = "UMAP", 
                          colour_by="cluster_type", text_by = "cluster_type", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
library(RColorBrewer)
plotNhoodGraphDA(scmilo, results, layout="UMAP",alpha = 1) +
  scale_color_continuous()
              

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
results <- annotateNhoods(scmilo, results, coldata_col = "cluster_type")
head(results)
results$cluster_type<-factor(results$cluster_type,levels = c("C0","C1","C2","C3"))
plotDAbeeswarm(results, group.by = "cluster_type",alpha =1)
#alpha值调高可以看所有的结果
