##作图
rm(list = ls())
library(Seurat)
library(dplyr)
library(igraph)
library(monocle)
library(RColorBrewer)
library(grid)
library(ggrastr)
library(tidydr)
library(ggplot2)
library(ggrepel)
setwd("F:/study/AML单细胞/new_result/result1/malignant/lspc")
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
###cluster之间相关性热图----
av<-AverageExpression(lspc_scrna,group.by = "cluster_type",assays = "RNA")
av<-av[[1]]
head(av)
###选出标准差最大的1000个gene
cg<-names(tail(sort(apply(av,1,sd)),1000))
##查看细胞群相关性矩阵
View(cor(as.matrix(av[cg,]),method = "spearman"))
pheatmap::pheatmap(cor(as.matrix(av[cg,]),method = "spearman"))
####亚群注释分组marker----
#使用 AverageExpression 提取亚群基因的均值表达量:
marker<-c("CD14","CLEC7A","FCN1","LYZ","VCAN",
          "S100A8","S100A9","S100A10",
          "SOCS2","CD34","AVP","SPINK2","FAM30A","HOXA9",
          "CD99","CD44","VIM","MEF2C",
          "TOP2A","MKI67","CENPE","CENPF","HMMR","CCNB1","CENPA",
          "TK1","MYBL2","PCNA","E2F1","MCM2",
          "ELANE","AZU1","CTSG","PRTN3","MT1X","MZB1","MEG3","ARIH1")
##"BCL11A","XBP1","NFE2","MAX", "GLI3","TNFSF10","CXCL1","CXCL2", "CXCL3",
#"JUN","PTGS2","PTPRC","SPTBN1",
#"GATA2", "ELK4","REST","CREB1",'MYCN'
mean<-AverageExpression(lspc_scrna,group.by = "cluster_type",features=marker,slot = 'data')%>%data.frame() %>%  as.matrix()
htdf <- t(mean)
rownames(htdf)<-unlist(lapply(strsplit(rownames(htdf),"RNA."),function(x){return(x[2])}))
aa<-c("C12","C13","C2","C9","C1","C0","C4","C7",
      "C14","C15", "C10","C8","C11","C5","C3","C6")
htdf <- htdf[aa,marker]
library(pheatmap)
pheatmap(htdf,scale="column",cluster_cols = F,cluster_rows = F,
         column_title = "Marker genes",row_title = "",
         row_names_gp = gpar(fontsize = 10),row_names_side = 'right',
         border = T,color =colorRampPalette(colors = c("#225EA8","white","#DE221E"))(100), 
         column_names_side = 'bottom')
#lspc_scrna$cluster_type<-factor(lspc_scrna$cluster_type,
#levels = c("C12","C13","C2","C9","C1","C0","C4","C7",
#"C14","C15", "C10","C8","C11","C5","C3","C6"))
#DotPlot(lspc_scrna,group.by = "cluster_type",features=marker)+ 
#coord_flip()  #翻转

###umap展示----
umap = lspc_scrna@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cluster_type = paste0("C",lspc_scrna@meta.data$seurat_clusters)) # 注释后的label信息 ，改为cluster_type
head(umap)
cluster_type_mean <- umap %>%
  group_by(cluster_type) %>%
  summarise(
    UMAP_1 = mean(umap_1),
    UMAP_2 = mean(umap_2)
  )
colourCount <-length(unique(lspc_scrna@meta.data$cluster_type))
#color<-colorRampPalette(brewer.pal(8, "Paired"))(colourCount)
ggplot()+
  geom_point_rast(data=umap, aes(x= umap_1 , y = umap_2 ,color = cluster_type),size =0.1,shape=16,alpha=1) +
  scale_color_manual(values =c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = cluster_type), 
                  data = cluster_type_mean,size = 5,point.padding=unit(0.5, "lines"))
setwd("F:\\study\\AML单细胞\\new_result\\result1\\lspc\\lspc")
ggsave("umap_cluster.pdf",width = 7,height = 6)
###细胞类型注释----
new.cluster.ids <- c("LSPC-like","LSPC-like","S100hi GMP-like","GMP-like",
                     "LSPC-like","Cycling GMP-like","GMP-like",
                     "LSPC-like","Cycling GMP-like","S100hi GMP-like",
                     "Cycling GMP-like","Cycling GMP-like","Mono-like",
                     "S100hi GMP-like","LSPC-like","LSPC-like") 
names(new.cluster.ids) <- levels(lspc_scrna) 
lspc_scrna <- RenameIdents(lspc_scrna, new.cluster.ids)
###slice干性打分----
library(tidyverse)
###将所有数据集干性得分合并
###GSE185991
setwd("F:\\study\\AML单细胞\\AML_data\\GSE185991\\data\\")
exp_file<-list.files()
setwd("F:\\study\\AML单细胞\\AML_data\\GSE185991\\GSE185991_slice干性得分\\")
slice_file<-list.files()
#首先将第一个文件的数据赋值给它
slice<-read.table(slice_file[1],header = T)
slice$cell<-paste0(exp_file[1],"_",slice$cell)
#从第二个文件开始合并
for(i in 2:length(slice_file)){
  a1<-read.table(slice_file[i],header = T)
  a1$cell<-paste0(exp_file[i],"_",a1$cell)
  slice<-bind_rows(slice,a1)
  print(i)
}
###OMIX002180
setwd("F:\\study\\AML单细胞\\AML_data\\OMIX002180\\data\\")
exp_file<-list.files()
setwd("F:\\study\\AML单细胞\\AML_data\\OMIX002180\\OMIX002180_slice干性得分\\")
slice_file<-list.files()
#首先将第一个文件的数据赋值给它
slice1<-read.table(slice_file[1],header = T)
slice1$cell<-paste0(exp_file[1],"_",slice1$cell)
#从第二个文件开始合并
for(i in 2:length(slice_file)){
  a1<-read.table(slice_file[i],header = T)
  a1$cell<-paste0(exp_file[i],"_",a1$cell)
  slice1<-bind_rows(slice1,a1)
  print(i)
}
slice<-rbind(slice,slice1)
rownames(slice)<-slice$cell
scRNA<-AddMetaData(object = lspc_scrna,     #seurat对象
                   metadata =slice )   #需要添加的metadata
meta<-scRNA@meta.data
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
meta$cluster_type<-factor(meta$cluster_type,levels = c("C0","C1","C2","C3"))
ggboxplot(meta, x = "cluster_type", y = "slice",fill="cluster_type",
          palette =c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Slice",xlab = "",
          error.plot="pointrange",
          title = "")+
  stat_compare_means(method = "anova",label.y = 1) 
setwd("F:\\study\\AML单细胞\\new_result\\result1\\lspc\\lspc")
ggsave("slice.pdf",width = 5,height =4)
###cytotrace----
#定位python路径
library(reticulate)#使R能读懂python环境
Sys.setenv(RETICULATE_PYTHON="D:\\软件 (D)\\python.exe")
use_python("D:\\软件 (D)\\python.exe") 
py_available()
use_python("D:\\软件 (D)\\python.exe")
py_config()
py_available()
import("numpy")
import("sklearn")
import("scipy")
import("annoy")
import("scanoramaCT")
#导入包
library(CytoTRACE)
library(Seurat)
exp1 <- as.matrix(lspc_scrna@assays$RNA$counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- lspc_scrna$cluster_type
phenot <- as.character(phenot)
names(phenot) <- rownames(lspc_scrna@meta.data)
emb <- lspc_scrna@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = 'F:/study/AML单细胞/new_result/result1/lspc/lspc/cytotrace/')
plotCytoGenes(results, numOfGenes = 30, outputDir = 'F:/study/AML单细胞/new_result/result1/lspc/lspc/cytotrace/')

###monocle----
###随机抽样每个cluster一半细胞
a<-data.frame(colnames(lspc_scrna),lspc_scrna@meta.data$seurat_clusters)
colnames(a)<-c("cell","cluster")
set.seed(123)
cluster0<-a[which(a$cluster==0),]
c0<-sample(cluster0$cell,16971/3)
cluster1<-a[which(a$cluster==1),]
c1<-sample(cluster1$cell,15376/3)
cluster2<-a[which(a$cluster==2),]
c2<-sample(cluster2$cell,3543/3)
cluster3<-a[which(a$cluster==3),]
c3<-sample(cluster3$cell,3071/3)
unions <- function (...) {
  Reduce(union, list(...))
}
cell_use<-unions(c0,c1,c2,c3)
exp<-subset(lspc_scrna,cell=cell_use)
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
data <- as(as.matrix(exp@assays$RNA$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = exp@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
##估计size factor和离散度
HSMM <- estimateSizeFactors(HSMM)   
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3) ##cut 3
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。
##Step 1: choosing genes that define progress
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~cluster_type",cores = 10)
HSMM_ordering_genes <-row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
#ordering_genes <- row.names (subset(diff_test_res, qval < 1e-4)) ## 不要也写0.1 ，而是要写0.01。
HSMM <- setOrderingFilter(HSMM, HSMM_ordering_genes)
plot_ordering_genes(HSMM)
#Step 2: reducing the dimensionality of the data 
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
#Step 3: ordering the cells in pseudotime 
HSMM <- orderCells(HSMM)
setwd("F:\\study\\AML单细胞\\new_result\\result1\\lspc\\lspc\\monocle")
plot_cell_trajectory(HSMM,color_by="Pseudotime", size=0.1,show_backbone=TRUE)  
ggsave("trajectory_Pseudotime.pdf",width = 5,height = 4)
mycol<-c("#F1B0B0","#D81D20","#5C412E","#E6F4F8")
#plot_cell_trajectory(HSMM,color_by="State", size=0.1,show_backbone=TRUE)+
  scale_colour_manual(values = mycol)  
#ggsave("trajectory_State.pdf",width = 5,height = 4)
plot_cell_trajectory(HSMM,color_by="cluster_type", size=0.1,show_backbone=TRUE)+
  scale_color_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))
ggsave("trajectory_cluster_type.pdf",width = 5,height = 4)

plot_cell_trajectory(HSMM,color_by="cluster_type", size=0.1,show_backbone=TRUE) +
  facet_wrap("~cluster_type")+
  scale_color_manual(values =c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))
ggsave("trajectory_cluster_type_fen.pdf",width = 8,height = 6)

plot_cell_trajectory(HSMM,color_by="new_type", size=0.1,show_backbone=TRUE)+
  facet_wrap("~new_type")+
  scale_color_manual(values=c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"))
ggsave("trajectory_new_type.pdf",width = 8,height = 4)
##密度图
library(ggpubr) 
library(ggplot2)
df <- pData(HSMM)  ## pData(HSMM)取出的是HSMM对象中HSMM@phenoData@data的内容 View(df) 
ggplot(df, aes(Pseudotime, colour = cluster_type, fill=cluster_type)) +   
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ 
  scale_color_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))+
  scale_fill_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))
ggsave("trajectory密度分布图.pdf",width = 5,height = 4)
HSMM$time_point = ceiling(HSMM$Pseudotime)
table(HSMM$time_point)
df <- pData(HSMM)
p2 = ggplot(df, aes(x = time_point, fill = cluster_type)) +
  geom_bar( position = "fill") + # 百分比柱状图
  scale_fill_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))+#c("#efff0f","#ff0f0f","#0fc3ff","#ff42d6","#870fff","#7bff0f","#5f5d5f")
  #scale_x_discrete(limits = c("Tumor","Metastatic","Large","Small")) + #更改x轴顺序#sample_dsmt
  #scale_x_discrete(limits = c("Tumor","Metastatic","Non-metastatic")) + #更改x轴顺序#sample_type
  #scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Cluster")) +
  labs(title = "",
       x = "Pseudotime",
       y = "Fraction") +
  theme_classic()
#coord_flip() # 倒转x与y轴
p2
ggsave("Pseudotime_cluster_type柱状图.pdf",width = 8,height = 5)

###统计细胞占比----
Cellratio <- prop.table(table(Idents(lspc_scrna),lspc_scrna$Sample), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
library(tidyr)
freq <-spread(Cellratio, Var1, Freq)
colnames(freq)[1] <- 'sample'
setwd("F:\\study\\AML单细胞\\new_result\\result1\\lspc\\lspc\\")
write.csv(freq,file = "sample_cluster_type.csv",row.names = F)
####ggplot2展示----
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(tidydr)
###FAB----
lspc_scrna@meta.data$FAB[is.na(lspc_scrna@meta.data$FAB)]<-"NA"
umap = lspc_scrna@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(FAB = lspc_scrna@meta.data$FAB) # 注释后的label信息 ，改为cluster_type
time_mean <- umap %>%
  group_by(FAB) %>%
  summarise(
    UMAP_1 = mean(umap_1),
    UMAP_2 = mean(umap_2)
  )
ggplot()+
  geom_point_rast(data=umap, aes(x= umap_1 , y = umap_2 ,color = FAB),size =0.1,shape=16) +
  scale_color_manual(values=c("M1"="#F6AB58","M2"="#A6DA6B","M3"="#ABDCE0",
                              "M4"="#E66351","M5"="#3C8AD5","NA"="#FCE8B6"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = FAB), 
                  data = time_mean,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_FAB.pdf",width = 7,height = 6)

###type----
umap = lspc_scrna@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(type = lspc_scrna@meta.data$new_type) # 注释后的label信息 ，改为cluster_type
time_mean <- umap %>%
  group_by(type) %>%
  summarise(
    UMAP_1 = mean(umap_1),
    UMAP_2 = mean(umap_2)
  )
ggplot()+
  geom_point_rast(data=umap, aes(x= umap_1 , y = umap_2 ,color = type),size =0.3,shape=16) +
  scale_color_manual(values=c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = type), 
                  data = time_mean,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_new_type.pdf",width = 7,height = 6)
###细胞比例图----
plot_data = as.data.frame(lspc_scrna@meta.data)
library(ggplot2)
library(scales)
library(paletteer)
library(ggsci)
library(ggbreak)
###按照cluster_type画----
ggplot(plot_data, aes(x = new_type, fill =cluster_type)) +
  geom_bar( position = "fill",width = 0.5) + # 百分比柱状图
  scale_fill_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Type")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme_classic()
ggsave("type_cluster_type.pdf",width = 4,height = 3)

###按照type画----
ggplot(plot_data, aes(x = cluster_type, fill = new_type)) +
  geom_bar( position = "fill",width = 0.6) + # 百分比柱状图
  #scale_fill_brewer(palette = "Blues") +     # 调色板{RColorBrewer}
  scale_fill_manual(values = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Type")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme_classic()
ggsave("cluster_type_type.pdf",width = 6,height = 4)

###按照FAB画----
ggplot(plot_data, aes(x = cluster_type, fill = FAB)) +
  geom_bar( position = "fill",width = 0.6) + # 百分比柱状图
  #scale_fill_brewer(palette = "Blues") +     # 调色板{RColorBrewer}
  scale_fill_manual(values = c("M1"="#F6AB58","M2"="#A6DA6B","M3"="#ABDCE0",
                               "M4"="#E66351","M5"="#3C8AD5","NA"="#FCE8B6"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "FAB")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme_classic()
ggsave("cluster_type_FAB.pdf",width = 6,height = 4)
library(RColorBrewer)
###按照样本画----
colourCount <-  length(unique(plot_data$patient))
ggplot(plot_data, aes(x = cluster_type, fill = patient)) +
  geom_bar( position = "fill",width = 0.7) + # 百分比柱状图
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Paired"))(colourCount))+  
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Sample")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme_classic()
ggsave("cluster_type_patient.pdf",width = 6,height = 5)
###按照所有信息画----
aa<-c("PT06","PT07","PT12","PT01","PT02","PT11","PT13","PT17","Pt3","Pt9",
      "PT08","PT09", "PT10","PT15","PT18","Pt10")
plot_data$patient<-factor(plot_data$patient,levels = aa)
ggplot(plot_data,aes(x = patient, fill =cluster_type)) +
  geom_bar( position = "fill") + # 百分比柱状图
  scale_fill_manual(values = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40,vjust = 0.85,hjust = 0.75,size = 12))
ggsave("patient_cluster_type.pdf",width = 8,height =5)
###计算显著性----
library(tidyr)
count = plot_data %>% group_by(seurat_clusters,new_type) %>% summarise(count = n()) %>% 
  spread(key = new_type,value = count, fill = 0) %>% as.data.frame()
rownames(count) = count$seurat_clusters
count = count[,-1]
count$col_sum = apply(count,1,sum)
count[5,] = apply(count,2,sum)
roe = count[1:4,1:2]
for(i in 1:nrow(roe)){
  for(j in 1:ncol(roe)){
    roe[i,j] = count[i,j]/((count[i,3]*count[5,j])/count[5,3])
  }
}
###热图----
library(pheatmap)
library("RColorBrewer")
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
aa<-pheatmap(roe,cluster_rows = F,cluster_cols = F,display_numbers = TRUE,
             # 设置数值字体大小
             fontsize_number = 14,angle_col = 45,
             border_color =F, cellwidth = 60, cellheight =40)
save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(aa,"cluster_roe.pdf")
table(plot_data$cluster_type)
###多组一起展示----
library(reshape2)
a<-data.frame(plot_data$Sample,plot_data$new_type)
a<-unique(a)
for (i in 1:nrow(a)) {
  data1<-plot_data[which(plot_data$Sample==a[i,1]),]
  a[i,3]<-length(which(data1$seurat_clusters==0))/length(which(plot_data$Sample==a[i,1]))
  a[i,4]<-length(which(data1$seurat_clusters==1))/length(which(plot_data$Sample==a[i,1]))
  a[i,5]<-length(which(data1$seurat_clusters==2))/length(which(plot_data$Sample==a[i,1]))
  a[i,6]<-length(which(data1$seurat_clusters==3))/length(which(plot_data$Sample==a[i,1]))
}
colnames(a)<-c("Sample","type","C0","C1","C2","C3")
melted_data <- melt(a, id.vars = c("Sample", "type"), measure.vars = colnames(a)[3:ncol(a)])
colnames(melted_data)<-c("sample","type","cluster_type","fraction")
require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)
ggboxplot(melted_data,x ="cluster_type", y ="fraction",fill ="type",
          palette = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Fraction",xlab = "",
          error.plot="pointrange")+# Add global p-value
  stat_compare_means(aes(group = type),method = "t.test",label.y = 1)+
  theme(axis.text.x = element_text(angle = 40,vjust = 0.85,hjust = 0.75))
ggsave("cluster_type_type_箱式图.pdf",width = 8,height =5)

###识别marker----
marker <- FindAllMarkers(lspc_scrna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
c0<-marker[which(marker$cluster==0),]
c1<-marker[which(marker$cluster==1),]
c2<-marker[which(marker$cluster==2),]
c3<-marker[which(marker$cluster==3),]
marker_df = marker %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)
write.table(marker,file="marker.txt",col.names=T,sep="\t",row.names=T,quote=F)
###marker气泡图----
plot_data<-lspc_scrna@meta.data
library("scales")
library(ggsci)
bubble.df=as.matrix(lspc_scrna[["RNA"]]$data[marker_df$gene,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$cell_name=rownames(bubble.df)
bubble.df=merge(bubble.df,plot_data[,c("cell_name","cluster_type")],by = "cell_name")
bubble.df$cell_name=NULL
cluster_type_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$cluster_type)) {
  bubble.df_small=bubble.df%>%filter(cluster_type==i)
  for (j in marker_df$gene) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    cluster_type_v=append(cluster_type_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}

plotdf=data.frame(
  cluster_type=cluster_type_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
#plotdf$cluster_type=factor(plotdf$cluster_type,levels = sort(unique(plotdf$cluster_type)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(marker_df$gene)))
plotdf$exp=ifelse(plotdf$exp>1,1,plotdf$exp)
plotdf%>%ggplot(aes(x=cluster_type,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_colour_distiller(palette = "RdYlBu")+
  #scale_color_gradient(low="#FAFDFD",high = "#A3217B")+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
ggsave(filename = "bubble_marker.pdf",width = 10,height = 18,units = c("cm"))
###cyclone计算细胞周期----
a<-lspc_scrna@assays$RNA$counts
a[1:4,1:4]
dat<-lspc_scrna@assays$RNA$data
dat[1:4,1:4]
group_list=lspc_scrna$seurat_clusters
library(scran)
sce <- SingleCellExperiment(list(counts=dat)) 
sce
# scran包安装好后，会在exdata文件夹中找到附件文件
library(org.Hs.eg.db)
# syste,.file会列出文件所在的路径，下图就是exdata文件夹下的文件，看到除了小鼠还有人的相关的RDS数据。这个RDS其实和平常看到的Rdata差不多，只不过Rdata是针对多个对象，Rds是针对一个对象进行存储和读取
mm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                package="scran"))
# 将symbol转为ensembl基因
ensembl <- mapIds(org.Hs.eg.db, keys=rownames(sce), 
                  keytype="SYMBOL", column="ENSEMBL")
head(ensembl)
system.time(assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl))
# 这一过程会比较慢，用system.time计算一下时间看看
str(assigned) # 包含了phases、scores、normalized.scores三个元素
phase<-as.data.frame(assigned$phases)
rownames(phase)<-colnames(lspc_scrna)
colnames(phase)[1]<-"phase"
lspc_scrna<-AddMetaData(object = lspc_scrna,     #seurat对象
                       metadata =phase )   #需要添加的metadata
library(tidyr)
plot<-lspc_scrna@meta.data
plot$cluster_type=factor(plot$cluster_type,levels =c("C0","C1","C2","C3"))
ggplot(plot, aes(x = cluster_type, fill = phase)) +
  geom_bar( position = "fill",width = 0.5) + # 百分比柱状图
  #scale_fill_brewer(palette = "Blues") +     # 调色板{RColorBrewer}
  scale_fill_manual(values = c("G1"="#3F89CE","G2M"="#B1D7E9","S"="#FDEDD4"))+
  #scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Time")) +
  labs(title = "",
       x = "",
       y = "Fraction")+theme_classic()
#+coord_flip()
ggsave("cluster_phase_cyclone.pdf",width = 6,height = 4)
###seurat计算细胞周期----
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
lspc_scrna<- CellCycleScoring(lspc_scrna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
###每个cluster细胞周期占比
library(tidyr)
plot<-lspc_scrna@meta.data
plot$cluster_type=factor(plot$cluster_type,levels =c("C0","C1","C2","C3"))
ggplot(plot, aes(x = cluster_type, fill = Phase)) +
  geom_bar( position = "fill",width = 0.5) + # 百分比柱状图
  #scale_fill_brewer(palette = "Blues") +     # 调色板{RColorBrewer}
  scale_fill_manual(values = c("G1"="#3F89CE","G2M"="#B1D7E9","S"="#FDEDD4"))+
  #scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Time")) +
  labs(title = "",
       x = "",
       y = "Fraction")+theme_classic()
#+coord_flip()
ggsave("cluster_phase.pdf",width = 6,height = 4)
#加载R包
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(tidyverse)
###quiescence_geneset----
setwd("F:\\study\\AML单细胞\\AML_data\\GSE185991_OMIX002180单细胞数据提取LSC\\lspc\\malignant_slice_0.9\\new_result\\cellstate\\quiescence")
quiescence_geneset<-read.csv("quiescence_geneset.csv",header = T)
a<-list()
a[[1]]<-quiescence_geneset[1:40,1]
a[[2]]<-quiescence_geneset[1:23,2]
names(a)<-colnames(quiescence_geneset)
normalized_expr_gran<-as.matrix(GetAssayData(lspc_scrna,assay='RNA',slot='data'))
library(GSVA)
gsva_matrix<- gsva(normalized_expr_gran,a,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
meta<-lspc_scrna@meta.data
score<-as.data.frame(t(gsva_matrix))
data<-merge(meta,score,by="row.names")
data$cluster_type=factor(data$cluster_type,levels =c("C0","C1","C2","C3"))
ggboxplot(data, x = "cluster_type", y = colnames(data)[33], fill  = "cluster_type",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title ="Negative_G0_To_G1")+
  stat_compare_means(method = "anova",label.y = 1)+      # Add global p-value
  stat_compare_means(method = "wilcox.test", label = "p.signif", #设置星标
                     ref.group ="C1")  
setwd("F:\\study\\AML单细胞\\new_result\\result1\\malignant\\lspc")
ggsave(file="Negative_G0_To_G1.pdf",width = 5,height = 4)

ggboxplot(data, x = "cluster_type", y = colnames(data)[34], fill  = "cluster_type",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title ="G0Mhigh")+
  stat_compare_means(method = "anova",label.y = 1.5)+      # Add global p-value
  stat_compare_means(method = "wilcox.test", label = "p.signif", #设置星标
                     ref.group ="C1")  
ggsave(file="G0Mhigh.pdf",width = 5,height = 4)
###增殖毒性---
setwd("F:/study/AML单细胞/AML_data/GSE185991_OMIX002180单细胞数据提取LSC/lspc/malignant_slice_0.9/增殖-毒性")
geneset<-read.csv("geneset.csv",header = F)
rownames(geneset)<-geneset$V1
geneset<-geneset[,-1]
aa<-as.data.frame(t(geneset[1,1:39]))
bb<-as.data.frame(t(geneset[2,1:166]))
cc<-as.data.frame(t(geneset[3,1:22]))
features <- list(aa[,1],bb[,1],cc[,1])#示例，需要真是的基因symbol
names(features)<-rownames(geneset)
normalized_expr_gran<-as.matrix(GetAssayData(lspc_scrna,assay='RNA',slot='data'))
library(GSVA)
gsva_matrix<- gsva(normalized_expr_gran,features,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
meta<-lspc_scrna@meta.data
score<-as.data.frame(t(gsva_matrix))
data<-merge(meta,score,by="row.names")
data$cluster_type=factor(data$cluster_type,levels =c("C0","C1","C2","C3"))
ggboxplot(data, x = "cluster_type", y = colnames(data)[35], fill  = "cluster_type",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title ="Proliferation")+
  stat_compare_means(method = "anova",label.y = 1.2)+      # Add global p-value
  stat_compare_means(method = "wilcox.test", label = "p.signif", #设置星标
                     ref.group ="C1")  
ggsave(file="Proliferation.pdf",width = 5,height = 4)
###Hallmark
library(clusterProfiler)
library(org.Hs.eg.db)
library(qusage)
library(Seurat)
gene_set<-qusage::read.gmt('F://study//AML单细胞//AML_data//GSE185991_OMIX002180单细胞数据提取LSC//lspc//malignant_slice_0.9//new_result//hallmark//h.all.v7.4.symbols.gmt')
normalized_expr_gran<-as.matrix(GetAssayData(lspc_scrna,assay='RNA',slot='data'))
library(GSVA)
gsva_matrix<- gsva(normalized_expr_gran, gene_set,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
setwd("F:\\study\\AML单细胞\\new_result\\result1\\malignant\\lspc\\hallmark")
write.table(gsva_matrix,"hallmark_gsva.txt",col.names = T,row.names = T,sep = "\t",quote = F)
###读取数据
hallmark_process<-read.csv("F://study//AML单细胞//AML_data//GSE185991_OMIX002180单细胞数据提取LSC//lspc//malignant_slice_0.9//new_result//hallmark//hallmark_process.csv",header = T)
hallmark_process$hallmark_name<-gsub("HALLMARK_", "", hallmark_process$hallmark_name)
gsva_matrix<-read.delim("hallmark_gsva.txt",header = T)
a<-gsva_matrix
b<-a
gsva_matrix<-b
meta<-lspc_scrna@meta.data
meta$cluster_type=factor(meta$cluster_type,levels =c("C0","C1","C2","C3"))
meta<-meta[order(meta$cluster_type),]
gsva_matrix<-gsva_matrix[,meta$cell_name]
rownames(gsva_matrix)<-gsub("HALLMARK_", "",rownames(gsva_matrix))
gsva_matrix<-gsva_matrix[hallmark_process$hallmark_name,]
library(pheatmap)
annotation_col = data.frame(
  cluster = meta$cluster_type,
  type = meta$new_type)
row.names(annotation_col) <- rownames(meta)
annotation_row = data.frame(
  Biological_process = hallmark_process$Process.category)
row.names(annotation_row) <- hallmark_process$hallmark_name
cluster_color<-c("#F1B0B0","#D81D20","#5C412E","#E6F4F8")
names(cluster_color)<-c("C0","C1","C2","C3")
type_color<-c("#3082BD","#A2D39A")
names(type_color)<-c("CR","Refractory/Relapse")
Biological_process_color<-c("#E66351","#F6AB58","#FCE8B6","#A6DA6B",
                            "#1C9751","#ABDCE0","#7DC0F5","#528FAB")
names(Biological_process_color)<-c("cellular component","development","DNA damage",
                                   "immune","metabolic","pathway","proliferation","signaling")
ann_colors <- list(cluster=cluster_color, type= type_color,Biological_process=Biological_process_color) #颜色设置
gsva_matrix1<-t(apply(gsva_matrix, 1, function(x){(x-min(x))/(max(x)-min(x))})) 
#rownames(gsva_matrix1)<-gsub("HALLMARK_", "", rownames(gsva_matrix1))
pheatmap(gsva_matrix1,cluster_rows = F,cluster_cols = F,show_rownames =T,
         show_colnames = F,fontsize=12,
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors,annotation_names_row = F)

####均值----
###读取数据
hallmark_process<-read.csv("hallmark_process.csv",header = T)
hallmark_process$hallmark_name<-gsub("HALLMARK_", "", hallmark_process$hallmark_name)
gsva_matrix<-read.delim("hallmark_gsva.txt",header = T)
setwd("F:/study/AML单细胞/AML_data/GSE185991_OMIX002180单细胞数据提取LSC/lspc/malignant_slice_0.9/new_result")
meta<-lspc_scrna@meta.data
meta$seurat_clusters<-paste0("C",meta$seurat_clusters)
meta$seurat_clusters<-factor(meta$seurat_clusters,levels = c("C0","C1","C2","C3"))
meta<-meta[order(meta$seurat_clusters,decreasing = F),]
gsva_matrix<-gsva_matrix[hallmark_process$hallmark_name,]
gsva_matrix1<-t(apply(gsva_matrix, 1, function(x){(x-min(x))/(max(x)-min(x))})) 
library(pheatmap)
annotation_col = data.frame(
  cluster = c("C0","C1","C2","C3"))
row.names(annotation_col) <- c("C0","C1","C2","C3")
annotation_row = data.frame(
  Biological_process = hallmark_process$Process.category)
row.names(annotation_row) <- hallmark_process$hallmark_name
cluster_color<-c("#F1B0B0","#D81D20","#5C412E","#E6F4F8")
names(cluster_color)<-c("C0","C1","C2","C3")
Biological_process_color<-c("#E66351","#F6AB58","#FCE8B6","#A6DA6B",
                            "#1C9751","#ABDCE0","#7DC0F5","#528FAB")
names(Biological_process_color)<-c("cellular component","development","DNA damage",
                                   "immune","metabolic","pathway","proliferation","signaling")
ann_colors <- list(cluster=cluster_color,Biological_process=Biological_process_color) #颜色设置
C0_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='0'))
C0_gsva<-gsva_matrix1[,C0_cells]
C0_gsva_mean<-apply(C0_gsva, 1, mean)
C1_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='1'))
C1_gsva<-gsva_matrix1[,C1_cells]
C1_gsva_mean<-apply(C1_gsva, 1, mean)
C2_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='2'))
C2_gsva<-gsva_matrix1[,C2_cells]
C2_gsva_mean<-apply(C2_gsva, 1, mean)
C3_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='3'))
C3_gsva<-gsva_matrix1[,C3_cells]
C3_gsva_mean<-apply(C3_gsva, 1, mean)
gsva_mat_mean<-data.frame(C0=C0_gsva_mean,C1=C1_gsva_mean,C2=C2_gsva_mean,
                          C3=C3_gsva_mean)
pheatmap(gsva_mat_mean,cluster_rows = F,cluster_cols = F,show_rownames =T,
         show_colnames = F,fontsize=12,
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors,annotation_names_row = F)
###ZDHHC21
DotPlot(lspc_scrna,group.by = "seurat_clusters",features="ZDHHC21")
mat1=as.data.frame(lspc_scrna[["RNA"]]$data["ZDHHC21",])
colnames(mat1)="exp"
meta<-lspc_scrna@meta.data
mat3=merge(mat1,meta,by="row.names")
mat3$cluster_type=factor(mat3$cluster_type,levels =c("C0","C1","C2","C3"))
ggboxplot(mat3, x = "cluster_type", y = "exp", fill  = "cluster_type",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Expression",xlab = "",
          error.plot="pointrange",
          title ="ZDHHC21")+
  stat_compare_means(method = "anova",label.y = 1.5)+      # Add global p-value
  stat_compare_means(method = "wilcox.test", label = "p.signif", #设置星标
                     ref.group ="C1",label.y = 1.2)  
ggsave(file="ZDHHC21.pdf",width = 5,height = 4)
###氧化磷酸化、糖酵解----
gsva_matrix<-as.data.frame(t(gsva_matrix))
gsva_matrix$cell_name<-rownames(gsva_matrix)
meta<-lspc_scrna@meta.data
meta$seurat_clusters<-paste0("C",meta$seurat_clusters)
data_use<-merge(meta,gsva_matrix,by="cell_name")
data_use$seurat_clusters<-factor(data_use$seurat_clusters,levels = c("C0","C1","C2","C3"))
library(ggpubr)
ggboxplot(data_use, x = "seurat_clusters", y ="OXIDATIVE_PHOSPHORYLATION", fill = "seurat_clusters",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "OXIDATIVE_PHOSPHORYLATION")+
  stat_compare_means(method = "anova",label.y =0.85)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.8, label = "p.signif")  
ggsave(file="OXIDATIVE_PHOSPHORYLATION.pdf",width = 5,height = 4)

ggboxplot(data_use, x = "seurat_clusters", y ="OXIDATIVE_PHOSPHORYLATION", fill = "new_type",color = "seurat_clusters",
          palette = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "OXIDATIVE_PHOSPHORYLATION")+
  stat_compare_means(method = "anova",label.y =0.81)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.76, label = "p.signif")  
ggsave(file="OXIDATIVE_PHOSPHORYLATION_newtype.pdf",width = 6,height = 4)
####GLYCOLYSIS----
ggboxplot(data_use, x = "seurat_clusters", y ="GLYCOLYSIS", fill = "seurat_clusters",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "GLYCOLYSIS")+
  stat_compare_means(method = "anova",label.y =0.35)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.32, label = "p.signif")  
ggsave(file="GLYCOLYSIS.pdf",width = 5,height = 4)

ggboxplot(data_use, x = "seurat_clusters", y ="GLYCOLYSIS", fill = "new_type",color = "seurat_clusters",
          palette = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "GLYCOLYSIS")+
  stat_compare_means(method = "anova",label.y =0.33)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.32, label = "p.signif")  
ggsave(file="GLYCOLYSIS_new_type.pdf",width = 6,height = 4)
###缺氧
ggboxplot(data_use, x = "seurat_clusters", y ="HYPOXIA", fill = "seurat_clusters",
          palette = c("C0"="#F1B0B0","C1"="#D81D20","C2"="#5C412E","C3"="#E6F4F8"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "HYPOXIA")+
  stat_compare_means(method = "anova",label.y =0.35)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.32, label = "p.signif")  
ggsave(file="HYPOXIA.pdf",width = 5,height = 4)
ggboxplot(data_use, x = "seurat_clusters", y ="HYPOXIA", fill = "new_type",color = "seurat_clusters",
          palette = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Score",xlab = "",
          error.plot="pointrange",
          title = "HYPOXIA")+
  stat_compare_means(method = "anova",label.y =0.35)+
  stat_compare_means(method = "wilcox.test",ref.group = "C1",label.y = 0.32, label = "p.signif")  
ggsave(file="HYPOXIA_newtype.pdf",width = 6,height = 4)
###缺氧折线图
setwd("F:\\study\\AML单细胞\\AML_data\\GSE185991_OMIX002180单细胞数据提取LSC\\lspc\\malignant_slice_0.9\\new_result\\cellstate\\hypoxia")
geneset<-read.csv("hypoxia_signature.csv",header = T)
geneset <- geneset[,1:5]
colnames(geneset) <- c("Harris_Hypoxia","Semenza_HIF1_Targets","Winter_Hypoxia_Metagene","Elvidge_Hypoxia_up","Leonard_Hypoxia")
set<-list()
set[[1]]<-geneset[1:81,1]
set[[2]]<-geneset[1:36,2]
set[[3]]<-geneset[1:95,3]
set[[4]]<-geneset[1:172,4]
set[[5]]<-geneset[1:41,5]
names(set) <- c("Harris_Hypoxia","Semenza_HIF1_Targets","Winter_Hypoxia_Metagene","Elvidge_Hypoxia_up","Leonard_Hypoxia")
normalized_expr_gran<-as.matrix(GetAssayData(lspc_scrna,assay='RNA',slot='data'))
library(GSVA)
gsva_matrix<- gsva(normalized_expr_gran,set,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
C0_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='0'))
C0_gsva<-gsva_matrix[,C0_cells]
C0_gsva_mean<-apply(C0_gsva, 1, mean)
C1_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='1'))
C1_gsva<-gsva_matrix[,C1_cells]
C1_gsva_mean<-apply(C1_gsva, 1, mean)
C2_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='2'))
C2_gsva<-gsva_matrix[,C2_cells]
C2_gsva_mean<-apply(C2_gsva, 1, mean)
C3_cells<-rownames(subset(lspc_scrna@meta.data,seurat_clusters=='3'))
C3_gsva<-gsva_matrix[,C3_cells]
C3_gsva_mean<-apply(C3_gsva, 1, mean)
gsva_mat_mean<-data.frame(C0=C0_gsva_mean,C1=C1_gsva_mean,C2=C2_gsva_mean,
                          C3=C3_gsva_mean)
gsva_mat_mean <- as.data.frame(t(gsva_mat_mean))
gsva_mat_mean <-gsva_mat_mean[,c(1,3,4)]
gsva_mat_mean$cluster <- c("C0", "C1", "C2", "C3")
library(tidyr)
gsva_mat_mean <- gsva_mat_mean %>% pivot_longer(cols = c("Harris_Hypoxia", "Winter_Hypoxia_Metagene",
                                                         "Elvidge_Hypoxia_up"),
                                                names_to = 'pathway',
                                                values_to = 'score')

ggplot(gsva_mat_mean,aes(x=cluster,y=score,group = pathway,shape=pathway,bg = "white"))+
  geom_point(aes(color=pathway),show.legend = F)+
  geom_line(aes(color=pathway))+
  scale_color_manual(name='pathway',labels=c("Harris_Hypoxia","Semenza_HIF1_Targets",
                                             "Leonard_Hypoxia"),
                     values = c("#9DBBD9","#D6E4EE","#D075AF"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(color = "black"))
setwd("F:/study/AML单细胞/new_result/result1/malignant/lspc/hallmark")
ggsave(file="Hypoxia_折线图.pdf",width = 5,height = 3)
