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

library(limma)
design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#tumor处需修改case
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 


##########
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


######自定义
#gene<-read.table("pan.txt")$V1

alldiff_mito=allDiff[rownames(allDiff) %in% gene,]
#write.csv(alldiff_mito,file ='alldiff_mito.csv',quote = F)
alldiff_sig=alldiff_mito[alldiff_mito$adj.P.Val<0.05,]
#write.csv(alldiff_sig,file ='alldiff_sig.csv',quote = F)

### 作热图
rt1=exprSet[,group_list=='control']
rt2=exprSet[,group_list=='OA']

rt=cbind(rt1,rt2)

anno=data.frame(row.names = colnames(rt),group=c(rep('control',75),rep('OA',135)))
#save(rt,anno,file ='key_train_exprSet.Rdata')


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

library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)
library(patchwork)
library(limma)
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
library(limma)
#rm(list = ls())
getwd()
#setwd("E:/麻醉/WGCNA_Pro_top5000")
GSE="ND"
#rt=read.table(paste0(GSE,".txt"),sep="\t",header=T,check.names=F)
#load("F:/顾客/肺动脉高压/GSE117261/03_差异分析/差异分析_input.Rdata")
ND<-gene
exp<-exprSet
exp<-exp[ND,]
rt<-exp
rt$geneNames<-rownames(rt)
rt<-na.omit(rt)
rt<-rt[,c(211,1:210)]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
data<-rt
data=t(data)
#data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
#控制组放置最前
#分组
sample<-anno
data=data[rownames(sample),]

afcon=as.matrix(table(sample[,1]))[1,1]
afcon=as.vector(afcon)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(as.numeric(group))
rownames(group)=rownames(data)
colnames(group)="Type"
input <- as.data.frame(cbind(group,data))
input$Type=as.factor(input$Type)
#采用十折交叉验证
svmRFE(input, k = 10, halve.above = 100) #分割数据，分配随机数
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件，AvgRank为按 10 次折叠的平均排名排序
write.csv(top.features,"机器学习_feature_svm.csv")

# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前40个变量进行SVM模型构建，然后导入已经运行好的结果
featsweep = lapply(1:44, FeatSweep.wrap, results, input) 

# 画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("svm-error.pdf",width = 4,height = 4)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()

pdf("svm-accuracy.pdf",width = 4,height = 4)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()

# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 
