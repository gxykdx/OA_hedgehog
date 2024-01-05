rm(list = ls())
### 富集（运气做）
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
# 牢记样本顺序
load('GSE55235.Rdata')
pdata1=pdata
group_list1=c(rep('control',10),rep('OA',10))
load('GSE55457.Rdata')
pdata2=pdata
group_list2=c(rep('control',10),rep('OA',10))
load('GSE82107.Rdata')
pdata3=pdata
group_list3=c(rep('control',7),rep('OA',10))
load('GSE51588.Rdata')
pdata4=pdata
group_list4=c(rep('control',10),rep('OA',40))
load('GSE57218.Rdata')
pdata5=pdata
group_list5=c(rep('control',2),rep('OA',14),rep('control',1),rep('OA',19),rep('control',4))
load('GSE98918.Rdata')
pdata6=pdata
group_list6=c(rep('control',12),rep('OA',12))
load('GSE117999.Rdata')
pdata7=pdata
group_list7=c(rep('control',10),rep('OA',10))
load('GSE12021.Rdata')
pdata8=pdata
group_list8=c(rep('control',4),rep('OA',10),rep('control',5))



#group_list=c(group_list1,group_list2)
#group_list=c(group_list1,group_list2,group_list3)
group_list=c(group_list1,group_list2,group_list3,group_list4,group_list5,group_list6,group_list7,group_list8)
group_list=factor(group_list,levels=c('control','OA'))
#3的数据集，提前2个
a=grep('GSE55235',colnames(rt))
b=grep('GSE55457',colnames(rt))
c=grep('GSE82107',colnames(rt))
d=grep('GSE51588',colnames(rt))
e=grep('GSE57218',colnames(rt))
f=grep('GSE98918',colnames(rt))
g=grep('GSE117999',colnames(rt))
h=grep('GSE12021',colnames(rt))
exprSet=rt[,c(a,b,c,d,e,f,g,h)]
#exprSet=rt[,c(a,b,c)]
#exprSet=rt[,c(a,b)]

library(limma)
design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#tumor处需修改case
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 
write.csv(allDiff,"allDiff.csv")

##########
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


######自定义
#gene<-read.table("pan.txt")$V1




alldiff_mito=allDiff[rownames(allDiff) %in% gene,]
write.csv(alldiff_mito,file ='alldiff_mito.csv',quote = F)
alldiff_sig=alldiff_mito[alldiff_mito$adj.P.Val<0.05,]
write.csv(alldiff_sig,file ='alldiff_sig.csv',quote = F)

### 作热图
rt1=exprSet[,group_list=='control']
rt2=exprSet[,group_list=='OA']

rt=cbind(rt1,rt2)

anno=data.frame(row.names = colnames(rt),group=c(rep('control',75),rep('OA',135)))
save(rt,anno,file ='key_train_exprSet.Rdata')


rt=rt[rownames(rt) %in% rownames(alldiff_mito),]


pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100),show_colnames = F,annotation_col = anno)

## 保存整理好的这些文件
### 火山图

library(ggplot2)

alldiff_mito$label=ifelse(alldiff_mito$adj.P.Val<0.05,rownames(alldiff_mito),'')


ggplot(alldiff_mito,aes(logFC, -log10(adj.P.Val)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-0.2,0.2,0.5,-0.5), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(adj.P.Val), color= -log10(adj.P.Val)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  #scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+# 调整主题和图例位置：
  theme(panel.grid = element_blank()
  )+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10(adj.P.Val)"),
         size = "none")+
  # 添加标签：
  geom_text(aes(label=label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust=1)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)")

### 箱线图
rt=rt[rownames(alldiff_mito),]
rt=as.data.frame(t(rt))
rt$group=anno$group
ss=intersect(gene,colnames(rt))
rt=rt[,c(ss,'group')]

library(tidyverse)
rt2=tidyr::pivot_longer(rt,cols = -c('group'),names_to = "gene",values_to = 'expression',
)


library(ggpubr)

ggboxplot(
  rt2,
  x = "gene",
  y = "expression",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "expression", palette=c("lightblue","#fe8180")
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))

### 相关性
load('key_train_exprSet.Rdata')


#####
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


######自定义
#gene<-read.table("pan.txt")$V1


rt=rt[rownames(rt) %in% gene,]
rt=as.data.frame(t(rt))

library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

rt=rt[anno$group=='OA',]
library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

sig_mito=read.csv('alldiff_sig.csv')
sig_mito=sig_mito$X

rt=rt[,sig_mito]
library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "#fe8180")
)

library(ggstatsplot)
ggstatsplot::ggscatterstats(rt,x='PRKX','GLI3')
ggstatsplot::ggscatterstats(rt,x='RPS27A','TOMM20')


###


# 加载数据集
rm(list = ls())
load('key_train_exprSet.Rdata')
set.seed(123456)


library(caret)


#######
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


######
######自定义
#gene<-read.table("pan.txt")$V1


rt=rt[rownames(rt) %in% gene,]
rt=as.data.frame(t(rt))

rt$group=anno$group
rt$group=factor(rt$group)


control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
#set.seed(2450836)
# 执行SVM-RFE算法
results <- rfe(rt[,1:43], 
               rt[,44],
               rfeControl = control,
               method = "svmRadial")


# 结果分析
print(results)
# 列出选择的变量集
rfe_gene=predictors(results)[1:20]
save(rfe_gene,file = "机器学习_SVM.Rdata")
write.csv(rfe_gene,"机器学习_svm-rfe.csv")
# 绘图
#plot(results,type=c("OA","control"))
plot(results,type=c("g","o"))
######################################
#############随机森林（减少变量）#################
####################################
#install.packages("randomForest")
library(randomForest)
#set.seed(2450836)
#!!!数字代表基因数量
x=rt[,1:43]
rf=randomForest(x = x,y = rt$group,ntree = 1000)

#rf=randomForest(group~., data=data, ntree=500)

plot(rf, main="Random forest", lwd=2)
#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
#set.seed(2450836)
rf2=randomForest(x = x,y = rt$group, ntree=optionTrees)
#查看基因的重要性
importance=importance(x=rf2)
dev.off()
#绘制基因的重要性图
varImpPlot(rf2, main="")

###########################图片美化#########################
#查看基因的重要性
#绘制基因的重要性图
importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")
#展示前20个基因的重要性
importance<-dplyr::arrange(importance, desc(importance))
#importance[order(importance$importance),]
af=importance[1:20,]
library(ggpubr)
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
pdf(file="importance_随机森林.pdf", width=3.5, height=4.5)
print(p1)
dev.off()
save(af,file = "机器学习_随机森林.Rdata")
write.csv(importance,"随机森林.csv")



#挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
# 可以调节！
rfGenes=names(rfGenes[1:20])     #挑选重要性评分大于0的基因

gene=intersect(rfe_gene,rfGenes)
gene

save(gene,file ='gene_rfe_rf.Rdata')



#rm(list = ls())
### 多因素逻辑回归
load('gene_rfe_rf.Rdata')
gene<-c("GLI3","CSNK1G2","WNT5A","WNT5B","PRKX","RAB23","BMP2","FBXW11")
gene<-c("GLI3","CSNK1G2","WNT5A","PRKX","RAB23","PRKACA")
gene<-c("GLI3","CSNK1G2","WNT5A","WNT5B","PRKX","RAB23","BMP2","FBXW11","PRKACA","BMP8B")
gene<-c("GLI3","CSNK1G2","WNT5A","WNT5B","PRKX","RAB23","FBXW11","PRKACA","BMP8B")


rt=rt[,c(gene,'group')]
rt$group=ifelse(rt$group=='control',0,1)


glm=glm(group ~.,data = rt,family= binomial(link='logit'))

OR=exp(glm$coefficients)[2:(length(gene)+1)]
OR_CI=exp(confint(glm,level = 0.95))
OR.95L=OR_CI[2:(length(gene)+1),1]
OR.95H=OR_CI[2:(length(gene)+1),2]


df=data.frame(gene=gene,OR=OR,OR.95L=OR.95L,OR.95H=OR.95H)


#install.packages("forestploter")
library(forestploter)

df$`OR (95% CI)` <- ifelse(is.na(df$OR), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   df$OR,df$OR.95L,df$OR.95H))


df$"" = paste(rep(" ",57), collapse = " ")
plot <- forest(df[, c(1, 6, 5)],
               est = df$OR,
               lower =df$OR.95L,
               upper =df$OR.95H,ref_line = 1,
               ci_column = 2)
plot

pred=predict(glm,rt,type = 'response')
pred_df=as.data.frame(pred)
save(pred_df,file ='pred_logistic.Rdata')



library(pROC)
plot.roc(rt$group,pred,print.auc=T,print.thres=T)

###bootstrap验证
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data) 
  # 得到预测的ROC
  glmROC=roc(rt$group,as.numeric(glm.pred))
  plot.roc(rt$group,glm.pred,add=T,col='#0072b5')
  # 返回ROC的AUC
  return(glmROC$auc)
} 

library(boot) 
set.seed(123456) 


# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))

rm(outtab)

###boot 灵敏度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(sen)
} 

library(boot) 
set.seed(123456789) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))



###boot 特异度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(spe)
} 

library(boot) 
set.seed(123456789) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))




#install.packages('regplot')
#BiocManager::install("vioplot")
#library(vioplot)
library(regplot)
regplot(glm,observation = rt[2,],interval = 'confidence',clickable = T,title = 'Nomogram')


library(rms)
glm2=lrm(group~., data = rt,x=T,y=T)

cal<-calibrate(glm2, method = 'boot', B=1000, data = rt)
par(oma=c(3,3,3,3)) 
par(mar=c(6,5,4,3) + 0.1) 
plot(cal,
     xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",main = "Calibration Curve")

#install.packages('set')
#devtools::install_github("yikeshu0611/do")
#if (!require(ggDCA)) {
#  devtools::install_github("yikeshu0611/ggDCA")
#}
#install.packages("rmda")
library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)

#注：若是devtools::install_github('yikeshu0611/ggDCA')也报错，可先运行：
#options(unzip ='internal')




set.seed(123456)
colnames(rt)

model1 <- decision_curve(group ~ GLI3, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
model2 <- decision_curve(group ~ WNT5B, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
colnames(rt)
#!!!!!!!!!!!!!!!!!!
model0 <- decision_curve(group ~ GLI3 +CSNK1G2+ WNT5A +WNT5B+ PRKX+RAB23+FBXW11+PRKACA+BMP8B,
                         data = rt,thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)

par(oma=c(1,1,1,1)) 
par(mar=c(6,5,4,3) + 0.1) 
plot_decision_curve( list(model1,model2,model0),
                     # curve.names = c("Baseline model", "Full model"),
                     # col = c("blue", "red"),
                     confidence.intervals = F,  #remove confidence intervals
                     cost.benefit.axis = FALSE, #remove cost benefit axis
                     legend.position = "bottomleft") #add the legend

