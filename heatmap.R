rm(list = ls())
options(stringsAsFactors = F)
library(pheatmap)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
exp <-  read.csv("E:Data/meth+400.csv",row.names = 1)

A <- as.matrix(exp)
for(i in 1:nrow(A)) A[i, ] <- scale(log(unlist(A[i, ] + 1), 1.5))#标准化处理
scale(A[,1:363])

annotation_col = data.frame(
  group = c( rep("Differentiated",98),rep("Immunoreactive",79),
             rep("Proliferative",100),rep("Mesenchymal",86)))
row.names(annotation_col) <- colnames(A)
groupcolor <- c("#85B22E","#5F80B4","#E29827","#57C3F3")
names(groupcolor) <- c("Differentiated","Immunoreactive","Proliferative","Mesenchymal")
ann_color <- list(group=groupcolor)#颜色设置

pheatmap(A,name = "rnorm",cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         legend = T,
         
         show_colnames = F,border_color = NA,show_rownames =T,
         annotation_col = annotation_col,annotation_colors = ann_color)


