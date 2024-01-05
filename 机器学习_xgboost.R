library(xgboost)
library(caret)
library(tidyverse)
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
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)

rm(list = ls())
#set.seed(123)
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
# Fitting model(用caret实现)
TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)
#model<- train(x=data,y=anno$group,  method = "xgbTree", trControl = TrainControl,verbose = FALSE)

plot(varImp(model))
importance <- varImp(model)
head(importance)
important <- as.data.frame(importance$importance) 
a<-important
varimpdf <- data.frame(var = row.names(a),
                       impor = a[,1])


ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (XGBoost)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 3))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))

importances<-important
write.csv(importances,"机器学习_XGBoost.csv")
importances$Gene<-rownames(importances)
importances=importances[,c(2,1)]
names(importances)=c("Gene","importance")
library(ggplot2)
library(ggpubr)
#展示前20个基因的重要性
importances<-dplyr::arrange(importances, desc(importance))
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
pdf(file="importance_xgboost.pdf", width=3.5, height=4.5)
print(p1)
dev.off()
save(af,file = "机器学习_xgboost.Rdata")
