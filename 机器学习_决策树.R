# 加载数据集
rm(list = ls())
library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
load('key_train_exprSet.Rdata')
set.seed(123456)
library(caret)

#######
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]
######
######自定义
#gene<-read.table("PANoptosis.txt")$V1
rt=rt[rownames(rt) %in% gene,]
rt=as.data.frame(t(rt))
group<-anno$group
data<-rt
dim(data)
colnames(data)
aggr(data)
set.seed(123456)
data2<-as.data.frame(data)
#建模
#mod1<-rpart(as.factor(group)~.,data = data2,method = "class",cp=0.000001)
mod1<-rpart(as.factor(group)~.,data = data2,method = "class")
#显示重要性
importances <- varImp(mod1)
importances %>%
  arrange(desc(Overall))

mod1$cp



write.csv(importances,"机器学习_决策树.csv")
importances$Gene<-rownames(importances)
importances=importances[,c(2,1)]
names(importances)=c("Gene","importance")
library(ggplot2)
library(ggpubr)
#展示前20个基因的重要性
importances<-dplyr::arrange(importances, desc(importance))
#af=importance[1:20,]
af=importances[1:20,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             dot.size = 3,                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_bw()         ,               # ggplot2 theme
             rotate=TRUE                                       )#翻转坐标轴 
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +#颜色
  grids()   
#保存图片
pdf(file="importance_决策树.pdf", width=3.5, height=4.5)
print(p1)
dev.off()
save(af,file = "机器学习_决策树.Rdata")
plotcp(mod1)

#查看模型最低CP值
#Complexity parameter是决策树每一次分裂时候最小的提升量，用于平衡模型精确度于复杂度
plotcp(mod1)
#模型优化（取最低CP值）
mod1<-rpart(as.factor(group)~.,data = data2,method = "class",cp=0.01000000)
visTree(mod1,main = "Decision Tree",height = "600px",
        colorY = c("greenYellow","hotPink","yellow"),legendWidth=0.2,legendNcol=2,  nodesFontSize = 16,edgesFontSize = 10,)
