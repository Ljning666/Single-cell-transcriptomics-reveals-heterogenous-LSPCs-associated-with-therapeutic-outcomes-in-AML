rm(list=ls())
gc()
###读取数据
library(Seurat)
library(dplyr)
library(data.table)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(dittoSeq)
library(ggrastr)
library(tidydr)
library(scatterpie)
library(ggrepel)
setwd("F:\\study\\AML单细胞\\new_result\\result1\\data\\AML_health\\")
scRNA<-readRDS("scrna_aml.rds")
meta<-scRNA@meta.data
table(meta$cell_type,meta$Malignant_cluster)
#meta$cell_type<-as.character(meta$cell_type)
meta[which(meta$Malignant_cluster=="Leukemia-like"),13]<-"AML_Blasts"
meta[which(meta$Malignant_cluster=="Leukemia-like"),10]<-"AML_Blasts"
meta[which(meta$type=="refractory"),23]<-"Refractory/Relapse"
meta[which(meta$type=="relapse"),23]<-"Refractory/Relapse"
table(meta$type)
scRNA<-AddMetaData(object = scRNA,     #seurat对象
                   metadata =meta)   #需要添加的metadata
###画图umap_cell_type
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = scRNA@meta.data$cell_type)
colnames(df_umap)
label_umap <- df_umap %>%group_by(cell_type) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$cell_type
###umap----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = cell_type),size =0.1,shape=16) +
  scale_color_manual(values=c("AML_Blasts"="#BD0026","HSPCs"= "#4292C4","Monocytes"= "#E3C496","GMPs"="#FD8D3C", 
                              "NK_T cells"= "#FFFFB2" ,"B cells"= "#751B68", 
                              "DCs"= "#F7C77C","Erythroid"="#EA4B30"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = cell_type), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_celltypeold.pdf",width = 8,height = 6)
###画图umap_newcelltype
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(newcelltype = scRNA@meta.data$newcelltype)
colnames(df_umap)
label_umap <- df_umap %>%group_by(newcelltype) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$newcelltype
###umap----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = newcelltype),size =0.1,shape=16) +
  scale_color_manual(values=c("AML_Blasts"="#BD0026","HSPCs"= "#4292C4","Monocytes"= "#E3C496","GMPs"="#FD8D3C", 
                              "NK_T cells"= "#FFFFB2" ,"B cells"= "#751B68", 
                              "DCs"= "#F7C77C","Erythroid"="#EA4B30"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = newcelltype), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_newcelltype.pdf",width = 8,height = 6)

###画图umap_patient
  color<-colorRampPalette(brewer.pal(8, "Spectral"))(16)
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(patient = scRNA@meta.data$patient)
colnames(df_umap)
label_umap <- df_umap %>%group_by(patient) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$patient
###umap----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = patient),size =0.2,shape=16) +
  scale_color_manual(values =c("PT01"="#B2B2A5","PT02"="#F05232","PT06"="#80132A",
                               "PT07"="#A53923","PT08"="#FABBA1","PT09"="#F59EB5",
                               "PT10"="#F069A0","PT11"="#772775","PT12"="#DFD4B3",
                               "PT13"="#FF9933","PT15"="#E4C494","PT17"="#4D134B",
                               "PT18"="#C8E3B4","Pt3"="#DADAEA","Pt9"="#6B52A2",
                               "Pt10"="#FED978"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = patient), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_patient.pdf",width = 8,height = 6)

###画图umap_type
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(type = scRNA@meta.data$type)
colnames(df_umap)
label_umap <- df_umap %>%group_by(type) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
rownames(label_umap ) <- label_umap$type
###umap----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = type),size =0.1,shape=16) +
  scale_color_manual(values =c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = type), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_type.pdf",width = 8,height = 6)

###画图umap_FAB
df_umap <- scRNA@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(FAB = scRNA@meta.data$FAB)
colnames(df_umap)
label_umap <- df_umap %>%group_by(FAB) %>%
  summarise(UMAP_1 = median(umap_1),
            UMAP_2 = median(umap_2))%>%
  as.data.frame()
label_umap[4,1]<-"Na"
rownames(label_umap ) <- label_umap$FAB
###umap----
ggplot()+
  geom_point_rast(data=df_umap, aes(x= umap_1 , y = umap_2 ,color = FAB),size =0.1,shape=16) +
  scale_color_manual(values =c("M2"="#f4d19b",
                               "M4"="#b17aa2","M5"="#fdfce5","Na"="#f9f3df"))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=16), #设置legend标签的大小
    legend.key.size=unit(0.8,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_text_repel(aes(x= UMAP_1 , y = UMAP_2 ,label = FAB), 
                  data = label_umap,size = 5,point.padding=unit(0.5, "lines"))
ggsave("umap_FAB.pdf",width = 8,height = 6)
###
plot_data = as.data.frame(scRNA@meta.data)
library(ggplot2)
library(scales)
library(paletteer)
library(ggsci)
library(ggbreak)
###按照所有16个患者画
aa<-c("PT01","PT02","PT06","PT07","PT08","PT09", "PT10",
      "PT11","PT12","PT13","PT15","PT17","PT18","Pt3","Pt9","Pt10")
plot<-plot_data[plot_data$patient%in%aa,]
plot$patient<-factor(plot$patient,levels = aa)
plot <- plot[order(plot$patient), ]
plot$cell_type<-factor(plot$cell_type,levels = c("AML_Blasts","HSPCs","NK_T cells","GMPs",
                                                 "Monocytes","DCs",
                                                 "Erythroid","B cells"))
ggplot(plot, aes(x = patient, fill = cell_type)) +
  geom_bar( position = "fill") + # 百分比柱状图
  scale_fill_manual(values = c("AML_Blasts"="#BD0026","HSPCs"= "#4292C4","Monocytes"= "#E3C496","GMPs"="#FD8D3C", 
                               "NK_T cells"= "#FFFFB2" ,"B cells"= "#751B68", 
                               "DCs"= "#F7C77C","Erythroid"="#EA4B30"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Cell type")) +
  labs(title = "",
       x = "Patient",
       y = "Fraction") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "#000000",size = 1),
        axis.text = element_text(colour = "#000000" ,size = 8),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75,size = 11), ##就是这里
        axis.ticks = element_line(colour = "#000000" ,size = 1) ,
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0,0,0),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("bili_patient_dx.pdf",width = 10,height = 6)
plot_data$type<-factor(plot_data$type,levels = c("CR","Refractory/Relapse"))
plot_data$cell_type<-factor(plot_data$cell_type,levels = c("AML_Blasts","HSPCs","NK_T cells","GMPs",
                                                           "Monocytes","DCs",
                                                           "Erythroid","B cells"))
ggplot(plot_data, aes(x = type, fill = cell_type)) +
  geom_bar( position = "fill",width = 0.5) + # 百分比柱状图
  scale_fill_manual(values = c("AML_Blasts"="#BD0026","HSPCs"= "#4292C4","Monocytes"= "#E3C496","GMPs"="#FD8D3C", 
                               "NK_T cells"= "#FFFFB2" ,"B cells"= "#751B68", 
                               "DCs"= "#F7C77C","Erythroid"="#EA4B30"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Cell type")) +
  labs(title = "",
       x = "",
       y = "Fraction") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "#000000",size = 1),
        axis.text = element_text(colour = "#000000" ,size = 12),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75,size = 12), ##就是这里
        axis.ticks = element_line(colour = "#000000" ,size = 1) ,
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0,0,0),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("bili_type_celltype.pdf",width = 5,height = 4)
####按照type排序
aa<-c("PT06","PT07","PT12","PT01","PT02","PT11","PT13","PT17","Pt3","Pt9",
      "PT08","PT09", "PT10","PT15","PT18","Pt10")
plot<-plot_data[plot_data$patient%in%aa,]
plot$patient<-factor(plot$patient,levels = aa)
plot <- plot[order(plot$patient), ]
plot$cell_type<-factor(plot$cell_type,levels = c("AML_Blasts","HSPCs","NK_T cells","GMPs",
                                                 "Monocytes","DCs",
                                                 "Erythroid","B cells"))
ggplot(plot, aes(x = patient, fill = cell_type)) +
  geom_bar( position = "fill") + # 百分比柱状图
  scale_fill_manual(values = c("AML_Blasts"="#BD0026","HSPCs"= "#4292C4","Monocytes"= "#E3C496","GMPs"="#FD8D3C", 
                               "NK_T cells"= "#FFFFB2" ,"B cells"= "#751B68", 
                               "DCs"= "#F7C77C","Erythroid"="#EA4B30"))+
  scale_y_continuous(labels = percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Cell type")) +
  labs(title = "",
       x = "Patient",
       y = "Fraction") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "#000000",size = 1),
        axis.text = element_text(colour = "#000000" ,size = 12),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75,size = 12), ##就是这里
        axis.ticks = element_line(colour = "#000000" ,size = 1) ,
        axis.ticks.length = unit(1,'mm'),
        plot.margin = unit(c(0.5,0,0,0),"cm"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("bili_patient_dx_排序.pdf",width = 10,height = 6)
###细胞类型箱式图----
a<-data.frame(plot_data$patient,plot_data$type)
a<-unique(a)
for (i in 1:nrow(a)) {
  data1<-plot_data[which(plot_data$patient==a[i,1]),]
  a[i,3]<-length(which(data1$cell_type=="B cells"))/length(which(plot_data$patient==a[i,1]))
  a[i,4]<-length(which(data1$cell_type=="Monocytes"))/length(which(plot_data$patient==a[i,1]))
  a[i,5]<-length(which(data1$cell_type=="NK_T cells"))/length(which(plot_data$patient==a[i,1]))
  a[i,6]<-length(which(data1$cell_type=="HSPCs"))/length(which(plot_data$patient==a[i,1]))
  a[i,7]<-length(which(data1$cell_type=="Erythroid"))/length(which(plot_data$patient==a[i,1]))
  a[i,8]<-length(which(data1$cell_type=="DCs"))/length(which(plot_data$patient==a[i,1]))
  a[i,9]<-length(which(data1$cell_type=="GMPs"))/length(which(plot_data$patient==a[i,1]))
  a[i,10]<-length(which(data1$cell_type=="AML_Blasts"))/length(which(plot_data$patient==a[i,1]))
}
colnames(a)<-c("patient","type","B cells","Monocytes",
               "NK_T cells","HSPCs","Erythroid","DCs","GMPs","AML_Blasts")
melted_data <- melt(a,id.vars = c("patient", "type"), measure.vars = colnames(a)[3:ncol(a)])
colnames(melted_data)<-c("patient","type","cell_type","fraction")
require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)
library(ggprism)#提供了GraphPad prism风格的主题和颜色，主要用于美化我们的图形
library(rstatix)
melted_data$type<-factor(melted_data$type,
                         levels = x_level)
melted_data$cell_type<-factor(melted_data$cell_type,levels = c("AML_Blasts","HSPCs","NK_T cells","GMPs",
                                                               "Monocytes","DCs",
                                                               "Erythroid","B cells"))
my_comparisons<-list(c("CR","Refractory/Relapse"))
###画在一起
ggboxplot(melted_data,x ="cell_type", y ="fraction",fill ="type",
          palette = c("CR"="#3082BD","Refractory/Relapse"="#A2D39A"),
          bxp.errorbar = TRUE,outlier.shape = NA,
          bxp.errorbar.width =0.3,
          ylab="Fraction",xlab = "",
          error.plot="pointrange")+      # Add global p-value
  stat_compare_means(method = "anova",label.y = 1)     # Add global p-value
#stat_compare_means(aes(group = type),method = "wilcox.test",label = "p.signif") 
ggsave('cell_type_all.pdf',width = 7,height = 4)
###识别marker----
scRNA<- JoinLayers(scRNA)
scRNA$cell_type<-as.factor(scRNA$cell_type)
scRNA@active.ident<-scRNA$cell_type
marker <- FindAllMarkers(scRNA,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

