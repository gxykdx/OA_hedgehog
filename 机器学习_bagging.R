rm(list = ls())
# 依赖包
library(tidyverse)       # for data wrangling
library(ggplot2)     # for awesome plotting
library(doParallel)  # for parallel backend to foreach
library(foreach)     # for parallel processing with for loops

# 建模包
library(caret)       # for general model fitting
library(rpart)       # for fitting decision trees
library(ipred)
set.seed(123456)
# for fitting bagged decision trees
#数据准备
load('key_train_exprSet.Rdata')
#set.seed(123456)
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
#rt=read.csv("diabetes.csv",header = T)
#model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)


####################################################
data$sample<-rownames(data)
anno$sample<-rownames(anno)
rt<-merge(data,anno,by="sample")
rownames(rt)<-rt$sample
rt<-rt[,2:45]
rownames(rt)<-NULL
rt$group=factor(rt$group,levels = c("control","OA"),labels = c("control","OA"))
####################################################

rt_bag1 <- bagging(
  formula = group ~ .,
  data = rt,
  nbagg = 100,  
  coob = TRUE,
  control = rpart.control(minsplit = 2, cp = 0)
)

rt_bag1
#vip(rt_bag1,aesthetics = list(fill=mycol))+theme_bw()

rt_bag2 <- train(
  group ~ .,
  data = rt,
  method = "treebag",
  trControl = trainControl(method = "cv", number = 10),
  nbagg = 200,  
  control = rpart.control(minsplit = 2, cp = 0)
)
rt_bag2
a<-rt_bag2$results
#vip(rt_bag2,aesthetics = list(fill=mycol))+
  theme_bw()

rt_bag2
importance <- varImp(rt_bag2)
important <- as.data.frame(importance$importance) 

importances<-important
write.csv(importances,"机器学习_bagging.csv")
importances$Gene<-rownames(importances)
importances=importances[,c(2,1)]
names(importances)=c("Gene","importance")
library(ggplot2)
library(ggpubr)
#展示前20个基因的重要性
importances<-dplyr::arrange(importances, desc(importance))
#order(importances$importance,decreasing = T)
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
pdf(file="importance_bagging.pdf", width=3.5, height=4.5)
print(p1)
dev.off()
save(af,file = "机器学习_bagging.Rdata")
