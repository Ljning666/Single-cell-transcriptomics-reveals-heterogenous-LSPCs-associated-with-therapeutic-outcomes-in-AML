rm(list = ls())
#加载包#
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
#加载绘图数据#
hallmark_process<-read.csv("I:/课题3/AML单细胞/AML_data/GSE185991_OMIX002180单细胞数据提取LSC/hsc/malignant_slice_0.9/new_result/hallmark/hallmark_process.csv",header = T)
hallmark_process$hallmark_name<-gsub("HALLMARK_", "", hallmark_process$hallmark_name)
setwd("F:/AML单细胞/new_result/result2/malignant/lspc/hallmark")#设置工作路径
gsva_matrix<-read.delim("hallmark_gsva.txt",header = T)
rownames(gsva_matrix)<-gsub("HALLMARK_", "",rownames(gsva_matrix))
gsva_matrix[1:4,1:4]
library(Seurat)
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
lspc_scrna$cluster_type<-factor(lspc_scrna$cluster_type,levels = c("C0","C1","C2","C3"))
gsva_matrix<-gsva_matrix[hallmark_process$hallmark_name,]
gsva_matrix1<-t(apply(gsva_matrix, 1, function(x){(x-min(x))/(max(x)-min(x))})) 
gsva_matrix1[1:4,1:4]
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
         show_colnames = F,fontsize=12,number_format = "%.3f",display_numbers=T,
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors,annotation_names_row = F)
cir1<-gsva_mat_mean[rev(rownames(gsva_mat_mean)),]
cir1<-t(scale(t(cir1)))
a1<-cir1[1:3,]
a1<-a1[order(a1[,2],decreasing = T),]
a1<-as.data.frame(a1)
a2<-cir1[4:9,]
a2<-a2[order(a2[,2],decreasing = T),]
a2<-as.data.frame(a2)
a3<-cir1[10:12,]
a3<-a3[order(a3[,2],decreasing = T),]
a3<-as.data.frame(a3)
a4<-cir1[13:18,]
a4<-a4[order(a4[,2],decreasing = T),]
a4<-as.data.frame(a4)
a5<-cir1[19:25,]
a5<-a5[order(a5[,2],decreasing = T),]
a5<-as.data.frame(a5)
a6<-cir1[26:31,]
a6<-a6[order(a6[,2],decreasing = T),]
a6<-as.data.frame(a6)
a7<-cir1[32:36,]
a7<-a7[order(a7[,2],decreasing = T),]
a7<-as.data.frame(a7)
a8<-cir1[37:50,]
a8<-a8[order(a8[,2],decreasing = T),]
a8<-as.data.frame(a8)
library(dplyr)
cir2<-bind_rows(a1,a2,a3,a4,a5,a6,a7,a8)
range(cir2)
#计算数据大小范围
range(cir1)
#绘制基础环形热图：
mycol<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7)

mycol=colorRamp2(c(-1,0,1),c("#4575B4","#FFFFBF","#D73027"))
circos.heatmap(cir2,col=mycol)
circos.clear()#绘制完成后需要使用此函数完全清除布局
#在circos.heatmap()中添加参数进行环形热图的调整和美化：

circos.par(gap.after=c(40)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽

circos.heatmap(cir1,col=mycol,
               
               #dend.side="inside", #聚类放在环形内测,控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               
               rownames.side="outside", #基因名放在环形外侧,控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               
               rownames.col="black",
               
               rownames.cex=0.5, #字体大小
               
               rownames.font=0.5, #字体粗细
               
               bg.border = "lightgrey",#背景边缘颜色
               
               #show.sector.labels = T，#显示分的区域的标签
               
               cluster=FALSE) #cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类

#circos.clear()#画完图结束。一定要运行这个，不然后续画图会叠加
#添加图例标签等

#library(ComplexHeatmap)

#install.packages("gridBase")

library(gridBase)

lg=Legend(title="Score",col_fun=mycol,
          
          direction = c("vertical"),
          
          #title_position= c('topcenter')，
          
)

#grid.draw(lg)

draw(lg, x = unit(0.95,"npc"), y = unit(0.5,"npc"), just = c("right","center"))#画在右边

#添加列名：

circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  
  if(CELL_META$sector.numeric.index==1){   #if(CELL_META$sector.numeric.index == 3) { # the last sector
    
    cn=colnames(cir1)
    
    n=length(cn)
    
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.5,"mm"),#x坐标
                
                1:n+5,#调整y坐标
                
                cn,cex=0.6,adj=c(0,0.5),facing="inside")}
  
},bg.border=NA)

circos.clear()
#分组热图绘制#

library(dendextend)

#但如果矩阵数据分组，可用split参数来指定分类变量

ann_row = data.frame(pathway=c(rep("cellular component",3),rep("developmen",6),
                               rep("DNA damage",3),rep("immune",6),
                               rep("metabolic",7),rep("pathway",6),
                               rep("proliferation",5),rep("signaling",14)))#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir2)

ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix

#分组绘图

circos.par(gap.after=c(5,5,5,5,5,5,5,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir2,col=mycol,
               
               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               
               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               
               track.height = 0.2, #轨道的高度，数值越大圆环越粗
               
               rownames.col="black",
               
               bg.border="lightgrey", #背景边缘颜色
               
               split = ann_row,#用行注释分裂热图
               
               show.sector.labels = T,
               
               rownames.cex=0.5,#字体大小
               
               rownames.font=0.5,#字体粗细
               
               cluster=FALSE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               
)
#图例与列名设置

library(gridBase)

lg=Legend(title="Score",col_fun=mycol,
          
          direction = c("vertical"),
          
          #title_position= c('topcenter')，
          
)

#grid.draw(lg)

draw(lg, x = unit(0.95,"npc"), y = unit(0.5,"npc"), just = c("right","center"))#画在右边

circos.clear()
