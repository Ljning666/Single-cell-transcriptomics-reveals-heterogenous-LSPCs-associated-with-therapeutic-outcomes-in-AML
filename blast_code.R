rm(list = ls())
library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(tidyverse)
library(Matrix)
aml_list<-list()
###GSE185991----
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\GSE185991\\")
exp_file<-list.files()
for(i in 1:length(exp_file)){
  data_dir<-paste0("F:\\study\\AML单细胞\\new_result\\result1\\data\\GSE185991\\",exp_file[i])
  aa<-Read10X(data.dir = data_dir, gene.column = 2)
  if(class(aa)=="list"){
    expression<-as.matrix(aa[[1]])
  }else{
    expression<-as.matrix(aa)
  }
  colnames(expression)<-unlist(lapply(strsplit(colnames(expression),"-"),function(x){return(x[1])}))
  colnames(expression)<-paste0(exp_file[i],"_",colnames(expression))
  exp_count<-CreateSeuratObject(counts = expression,min.cells = 3,min.features = 200,project=exp_file[i])
  dim(exp_count@assays$RNA)
  exp_count[["percent.mt"]]<-PercentageFeatureSet(exp_count,pattern="^MT-")
  ####过滤标准可根据样本测序情况调整
  exp_count <-subset(exp_count,subset=nFeature_RNA<6000 & nFeature_RNA>200 & nCount_RNA<40000 & nCount_RNA>400  & percent.mt<20)
  exp_count@meta.data$original <- exp_file[i]
  aml_list[[i]]<-exp_count
  print(i)
}
##OMIX002180----
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\OMIX002180\\data\\")
exp_file<-list.files()
for(i in 1:length(exp_file)){
  data_dir<-paste0("F:\\study\\AML单细胞\\new_result\\result1\\data\\OMIX002180\\data\\",exp_file[i])
  aa<-Read10X(data.dir = data_dir, gene.column = 2)
  expression<-as.matrix(aa)
  colnames(expression)<-unlist(lapply(strsplit(colnames(expression),"-"),function(x){return(x[1])}))
  colnames(expression)<-paste0(exp_file[i],"_",colnames(expression))
  exp_count<-CreateSeuratObject(counts = expression,min.cells = 3,min.features = 200,project=exp_file[i])
  dim(exp_count@assays$RNA)
  exp_count[["percent.mt"]]<-PercentageFeatureSet(exp_count,pattern="^MT-")
  ####过滤标准可根据样本测序情况调整
  exp_count <-subset(exp_count,subset=nFeature_RNA<6000 & nFeature_RNA>200 & nCount_RNA<40000 & nCount_RNA>400  & percent.mt<20)
  exp_count@meta.data$original<-exp_file[i]
  aml_list[[i+16]]<-exp_count
  print(i)
}
###health----
setwd("F:\\study\\AML单细胞\\AML_data\\GSE120221\\data\\")
exp_file<-list.files()
for(i in 1:length(exp_file)){
  data_dir<-paste0("F:\\study\\AML单细胞\\AML_data\\GSE120221\\data\\",exp_file[i])
  aa<-Read10X(data.dir = data_dir, gene.column = 2)
  expression<-as.matrix(aa)
  colnames(expression)<-unlist(lapply(strsplit(colnames(expression),"-"),function(x){return(x[1])}))
  colnames(expression)<-paste0(exp_file[i],"_",colnames(expression))
  exp_count<-CreateSeuratObject(counts = expression,min.cells = 3,min.features = 200,project=exp_file[i])
  exp_count[["percent.mt"]]<-PercentageFeatureSet(exp_count,pattern="^MT-")
  ####过滤标准可根据样本测序情况调整
  exp_count <-subset(exp_count,subset=nFeature_RNA<6000 & nFeature_RNA>200 & nCount_RNA<40000 & nCount_RNA>400  & percent.mt<20)
  exp_count@meta.data$original<-exp_file[i]
  aml_list[[i+19]]<-exp_count
  print(i)
}
###合并数据集----
names(aml_list) <- sapply(aml_list, function(x){unique(x@meta.data$original)})
exp_merge<-Reduce(merge,aml_list)
gc()
dim(exp_merge@assays$RNA)##27796 169001
saveRDS(exp_merge,file="F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\exp_merge_QC_Doublet.rds")#质控后的seurat对象

###剔除双核细胞----
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\")
DoubletFinder_res = read.table("DoubletFinder_res.txt",sep = "\t",header = T)
doublet<-DoubletFinder_res$cell[DoubletFinder_res$doubletFinder!="Singlet"]
cell_use<-setdiff(colnames(exp_merge),doublet)
exp_merge = subset(exp_merge,cells = cell_use)
dim(exp_merge@assays$RNA)##27796 162859
saveRDS(exp_merge,file="F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\exp_merge_QC.rds")#质控后的seurat对象
exp_merge<-readRDS("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\exp_merge_QC.rds")
exp_merge <- NormalizeData(exp_merge, verbose = FALSE)
###识别高可变基因----
gene<-list()
query<-as.data.frame(table(exp_merge@meta.data$original))
for (i in 1:nrow(query)) {
  scRNA_list_sample<-exp_merge[,which(exp_merge@meta.data$original==query[i,1])]
  scRNA_list_sample <- NormalizeData(scRNA_list_sample, verbose = FALSE)
  scRNA_list_sample <- FindVariableFeatures(scRNA_list_sample, selection.method = "vst", nfeatures =2000, verbose = FALSE)
  mito.genes <- grep(pattern = "^(TRAV|TRBV|TRGV|TRDV|IGHV|IGKV|IGLV|MT-|RP[LS])", rownames(scRNA_list_sample), value = TRUE)
  scRNA_list_sample@assays$RNA@meta.data$var.features[scRNA_list_sample@assays$RNA@meta.data$var.features%in%mito.genes]<-NA
  gene[[i]]<-scRNA_list_sample@assays$RNA@meta.data$var.features[!is.na(scRNA_list_sample@assays$RNA@meta.data$var.features)]
  print(i)
}
#所有样本top2000的高变基因取并集，然后排序，方差越大，排名越靠后，最后，所有数据集基因排序取平均，在取前1500
gene_sample_all<-as.data.frame(gene[[1]])
gene_sample_all[,2]<-rev(seq(1,nrow(gene_sample_all),1))
colnames(gene_sample_all)[1]<-'gene'
for (i in 2:nrow(query)) {
  gene_sample<-as.data.frame(gene[[i]])
  gene_sample[,2]<-rev(seq(1,nrow(gene_sample),1))
  colnames(gene_sample)[1]<-'gene'
  gene_sample_all<-merge(gene_sample_all,gene_sample,by='gene',all=T)
  print(i)
}
rownames(gene_sample_all)<-gene_sample_all[,1]
gene_sample_all<-gene_sample_all[,-1]
gene_sample_all[is.na(gene_sample_all)] <- 0
gene_mean<-as.data.frame(rownames(gene_sample_all))
gene_mean[,2]<-apply(gene_sample_all, 1, mean)
gene_mean<-gene_mean[order(gene_mean$V2),]
gene_mean[,3]<-rev(seq(1,nrow(gene_mean),1))
gene_mean_filter<-gene_mean[which(gene_mean$V3<=1500),]
colnames(gene_mean_filter)<-c('gene','mean','mean_order')
gene_mean_filter<-gene_mean_filter[order(gene_mean_filter$mean_order),]
var_gene<-as.character(gene_mean_filter$gene)
feature<-row.names(exp_merge)
feature[!feature%in%var_gene]<-NA
exp_merge <- FindVariableFeatures(exp_merge, selection.method = "vst", nfeatures =2000, verbose = FALSE)
exp_merge@assays$RNA@meta.data$var.features<-feature#替换高变基因
exp_merge <-ScaleData(exp_merge,features=exp_merge@assays$RNA@meta.data$var.features)
exp_merge <-RunPCA(exp_merge,features=exp_merge@assays$RNA@meta.data$var.features)
gc()
##harmony整合
system.time({exp_merge <- RunHarmony(exp_merge,max.iter.harmony = 50,group.by.vars = "original")})
#降维聚类
exp_merge <- RunUMAP(exp_merge, reduction = "harmony", dims = 1:30)
exp_merge <- RunTSNE(exp_merge, reduction = "harmony", dims = 1:30)
exp_merge <- FindNeighbors(exp_merge, reduction = "harmony", dims = 1:30)
saveRDS(exp_merge,file="F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\exp_merge_harmony.rds")#去批次后的seurat对象
###调参
library(Seurat)
library(dplyr)
library(data.table)
library(tidyverse)
library(Matrix)
exp_merge <- readRDS("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\exp_merge_harmony.rds")
exp_merge <- FindClusters(exp_merge, resolution = 0.5)
meta<-exp_merge@meta.data
meta[1:87753,8]<-"AML"
meta[87754:162859,8]<-"Health"
colnames(meta)[8]<-"source"
exp_merge<-AddMetaData(object = exp_merge,     #seurat对象
                   metadata =meta )   #需要添加的metadata
sum = exp_merge@meta.data %>% mutate(sourceBin = ifelse(exp_merge$source == 'AML',1,0)) %>% group_by(seurat_clusters) %>% 
  summarise(AML.pct = sum(sourceBin)/n()*100, HC.pct = 100-sum(sourceBin)/n()*100)

####marker注释----
DotPlot(exp_merge,group.by = "seurat_clusters",features=c("CD34","AVP","SPINK2","PRAME",
                                                      "MPO","AZU1","ELANE","CD38",
                                                      "PNASE2","CPA3","RUNX1","NREP",
                                                      "ARMH1","C1QBP","TRH",
                                                      "CLEC11A","HGF","STAB1",
                                                      "CCND2","LYZ","FCN1","CD14",
                                                      "HLA-DMB","HLA-DRA","HLA-DPB1",
                                                      "HBB","HBA1","CA1",
                                                      "CD3E","CD3D","CD3G","GNLY",
                                                      "CD19","CD79A","IGHM","IGHA1",
                                                      "IGHG1"))
m.cluster = sum$seurat_clusters[which(sum$HC.pct<45)]
exp_merge$Malignant_cluster = ifelse(exp_merge$seurat_clusters %in% m.cluster,'Leukemia-like','Normal-like')
#integrated$Malignant_cluster = ifelse(dat$source2=='HC' | dat$CellType %in% c('B','NaiveT','CTL','CTL_NK','Ery'),'NormalLike',dat$Malignant_cluster)
aa<-exp_merge@meta.data
aml<-aa[which(aa$source=="AML"),]
table(aml$Malignant_cluster)
aml$cell_name<-rownames(aml)
aaa<-aml[,c(10,9)]
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\")
write.table(aaa,file ="malignant_new.txt",row.names = F,sep = "\t",quote = F )

DimPlot(exp_merge,reduction = "umap",group.by = "seurat_clusters",
        label=T,raster=FALSE,label.size = 6)
DimPlot(exp_merge,reduction = "umap",group.by = "source",
        label=T,raster=FALSE,label.size = 6)
scRNA<-exp_merge
library(ggplot2)
library(dittoSeq)
library(ggrastr)
library(tidydr)
library(scatterpie)
library(ggrepel)
###画图umap
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(source = scRNA@meta.data$source)
colnames(df_umap)
label_umap <- df_umap %>%group_by(source) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$source
#ggplot作图----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = source),size =0.1,shape=16) +
  scale_color_manual(values=c("AML"="#BD2619","Health"="#EDEEEE"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label =source), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health")
ggsave("umap_source.pdf",width = 8,height = 6)
###
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\")
meta<-read.delim("meta.txt",header = T)
rownames(meta)<-meta$cell_name
malignant<-read.delim("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\malignant.txt",header = T)
meta<-merge(meta,malignant,by="cell_name")
meta[which(meta$cell_type%in%c('B cells','DCs','Erythroid','NK_T cells')),22]<-"Normal-like"
table(meta$cell_type,meta$Malignant_cluster)

meta1<-scRNA@meta.data
meta1[1:87753,9]<-meta[1:87753,22]
meta1[which(meta1$source=="Health"),9]<-"Normal"
table(meta1$Malignant_cluster)
table(meta$source)
scRNA<-AddMetaData(object = scRNA,     #seurat对象
                   metadata =meta1 )   #需要添加的metadata
###画图umap
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(malignant = scRNA@meta.data$Malignant_cluster)
colnames(df_umap)
label_umap <- df_umap %>%group_by(malignant) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$malignant
#ggplot作图----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = malignant),size =0.3,shape=16) +
  scale_color_manual(values=c("Leukemia-like"="#DF251F","Normal-like"="#92C7E9","Normal"="#EDEEEE"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label =malignant), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health")
ggsave("umap_malignant.pdf",width = 8,height = 6)
