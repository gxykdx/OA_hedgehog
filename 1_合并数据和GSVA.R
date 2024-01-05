#GSE55235
##——————————————————————————数据下载和整理——————————————————————————————————
if(! require("GEOquery")) BiocManager::install("GEOquery",update=F,ask=F)

rm(list = ls())
library(GEOquery)
##下载并且了解你的数据
gset=getGEO('GSE55235',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,1:20]
pdata<-pdata[1:20,]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',10),rep('OA',10))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL96-57554.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
anno = anno[!anno$`Gene Symbol`== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE55235.Rdata')
write.table(exprSet,file ='GSE55235.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))



#GSE55457--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE55457',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,-(11:23)]
pdata<-pdata[-(11:23),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',10),rep('OA',10))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL96-57554.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
anno = anno[!anno$`Gene Symbol`== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE55457.Rdata')
write.table(exprSet,file ='GSE55457.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))






#GSE82107--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE82107',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-(11:23)]
#pdata<-pdata[-(11:23),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',7),rep('OA',10))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL570-55999.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
anno = anno[!anno$`Gene Symbol`== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE82107.Rdata')
write.table(exprSet,file ='GSE82107.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))





#GSE51588--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE51588',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-(11:23)]
#pdata<-pdata[-(11:23),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',10),rep('OA',40))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL13497-9755.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,7)] 
#!!!!
anno = anno[!anno$GENE_SYMBOL== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE51588.Rdata')
write.table(exprSet,file ='GSE51588.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))





#GSE1919--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE1919',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,(1:10)]
pdata<-pdata[(1:10),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',5),rep('OA',5))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL91-30375.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
#anno = anno[!anno$GENE_SYMBOL== "",] 
anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE1919.Rdata')
write.table(exprSet,file ='GSE1919.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))






#GSE12021-GPL96--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE12021',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,-c(5,8,9,10,11,12,14,15,22,23,24,26)]
pdata<-pdata[-c(5,8,9,10,11,12,14,15,22,23,24,26),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',4),rep('OA',10),rep('control',5))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL96-57554.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
#anno = anno[!anno$GENE_SYMBOL== "",] 
anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE12021-GPL96.Rdata')
write.table(exprSet,file ='GSE12021-GPL96.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))




#GSE12021-GPL97--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE12021',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[2]])
pdata=pData(gset[[2]])

library(data.table)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('OA',10),rep('control',4))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL97-17394.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#!!!!
#anno = anno[!anno$GENE_SYMBOL== "",] 
anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE12021-GPL97.Rdata')
write.table(exprSet,file ='GSE12021-GPL97.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))




#GSE48556--------------------------
############
##下载并且了解你的数据
gset=getGEO('GSE48556',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
table(pdata$source_name_ch1)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
#pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('OA',106),rep('control',33))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL6947-13512.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,14)] 
#!!!!
anno = anno[!anno$Symbol== "",] 
#anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE48556.Rdata')
write.table(exprSet,file ='GSE48556.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))





#GSE89408--------------------------
############
##下载并且了解你的数据
rm(list = ls())
gset=getGEO('GSE89408',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
#table(pdata$source_name_ch1)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
#pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',106),rep('control',33))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL6947-13512.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,14)] 
#!!!!
anno = anno[!anno$Symbol== "",] 
#anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE89408.Rdata')
write.table(exprSet,file ='GSE89408.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))





#GSE117999--------------------------
############
##下载并且了解你的数据
rm(list = ls())
library(readxl)
library(tidyverse)
library(GEOquery)
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#rm(list = ls())
#gset=getGEO('GSE117999',destdir = '.',getGPL = F)
gset=getGEO(filename="GSE117999_series_matrix.txt.gz",getGPL = F)
#
#exprSet=exprs(gset[[1]])
#exprSet=gset[[1]]

exprSet <- exprs(gset) %>% 
  as.data.frame()  # 提取表达矩阵

##有空值，需剔除
exprSet <- exprSet[,-c(9,11,15,19)]
##有空值，需剔除
#pdata=pData(gset[[1]])
pdata <- pData(gset)
pdata <- pdata[-c(9,11,15,19),]
library(data.table)
#table(pdata$source_name_ch1)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
#pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',10),rep('OA',10))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
#exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL20844-88004.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,10)] 
#!!!!
anno = anno[!anno$GENE_SYMBOL== "",] 
#anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE117999.Rdata')
write.table(exprSet,file ='GSE117999.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)


######自定义
gene<-read.table("pan.txt")$V1



alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))





#GSE98918--------------------------
############
##下载并且了解你的数据
rm(list = ls())
library(readxl)
library(tidyverse)
library(GEOquery)
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#rm(list = ls())
#gset=getGEO('GSE117999',destdir = '.',getGPL = F)
gset=getGEO(filename="GSE98918_series_matrix.txt.gz",getGPL = F)
#
#exprSet=exprs(gset[[1]])
#exprSet=gset[[1]]

exprSet <- exprs(gset) %>% 
  as.data.frame()  # 提取表达矩阵

##有空值，需剔除
#exprSet <- exprSet[,-c(9,11,15,19)]
##有空值，需剔除
#pdata=pData(gset[[1]])
pdata <- pData(gset)
#pdata <- pdata[-c(9,11,15,19),]
library(data.table)
#table(pdata$source_name_ch1)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
#pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',12),rep('OA',12))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
#exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL20844-88004.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,10)] 
#!!!!
anno = anno[!anno$GENE_SYMBOL== "",] 
#anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE98918.Rdata')
write.table(exprSet,file ='GSE98918.txt',sep = '\t',col.names = NA,quote = F)



#GSE57218--------------------------
############
##下载并且了解你的数据
rm(list = ls())
library(readxl)
library(tidyverse)
library(GEOquery)
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#rm(list = ls())
#gset=getGEO('GSE117999',destdir = '.',getGPL = F)
gset=getGEO(filename="GSE57218_series_matrix.txt.gz",getGPL = F)
#
#exprSet=exprs(gset[[1]])
#exprSet=gset[[1]]

exprSet <- exprs(gset) %>% 
  as.data.frame()  # 提取表达矩阵

##有空值，需剔除
#exprSet <- exprSet[,-c(9,11,15,19)]
##有空值，需剔除
#pdata=pData(gset[[1]])
pdata <- pData(gset)
pdata = filter(pdata, pdata$`disease state:ch1` !='Preserved')
exprSet <-exprSet[,-c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69)]
#pdata <- pdata[-c(9,11,15,19),]
library(data.table)
#table(pdata$source_name_ch1)
#exprSet=fread('./GSE55235_RAW.txt',data.table = F)
#rownames(exprSet)=exprSet$ID_REF
#exprSet=exprSet[,-c(1,4,5,6,7,8,10,11,18,19,20,22)]
#pdata<-pdata[-c(1,4,5,6,7,8,10,11,18,19,20,22),]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('control',2),rep('OA',14),rep('control',1),rep('OA',19),rep('control',4))
group_list=factor(group_list,levels = c('control','OA'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
#exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL6947-13512.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,14)] 
#!!!!
anno = anno[!anno$Symbol== "",] 
#anno = anno[!anno$`Gene Symbol`== "",]
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE57218.Rdata')
write.table(exprSet,file ='GSE57218.txt',sep = '\t',col.names = NA,quote = F)
















#######去批次
rm(list = ls())
library(limma)
# 没有先安装
#BiocManager::install('sva',ask = F,update = F)
library(sva)

# 运行以下代码去批次
mergeFile="merge.preNorm.txt"            #合并后的文件名称
normalizeFile="merge.normalzie.txt"      #矫正后的文件名称
# 一定是在昨天输出这两个txt的文件夹下
files=c('GSE55235.txt','GSE55457.txt','GSE82107.txt','GSE51588.txt','GSE57218.txt','GSE98918.txt',
        'GSE117999.txt','GSE12021.txt')        #输入文件名称
#获取交集基因
length(files)
library(data.table)

geneList=list()
for(i in 1:length(files)){
  fileName=files[i]
  rt=fread(fileName,data.table = F)
  header=unlist(strsplit(fileName, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
# 共同基因
intersectGenes=Reduce(intersect, geneList)
library(limma)
#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  fileName=files[i]
  header=unlist(strsplit(fileName, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=fread(fileName,data.table = F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对数值大的数据取log2
  #qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  #LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  #if(LogC){
  #rt[rt<0]=0
  #rt=log2(rt+1)}
  #rt=normalizeBetweenArrays(rt)
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(header[1],ncol(rt)))
}

allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)

#对数据进行批次矫正，输出矫正后的结果
library(sva)
normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)

#合并前的PCA###################
#install.packages('ggplot2')
library(ggplot2)        #引用包
#读取输入文件,提取数据
rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

ggplot(data=PCA,aes(x=PC1,y=PC2))+
  geom_point(aes(color=Type,shape=Type),
             size=5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw()+
  scale_shape_manual(values=15:22)


#合并后PCA###################
#读取输入文件,提取数据
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图


ggplot(data=PCA,aes(x=PC1,y=PC2))+
  geom_point(aes(color=Type,shape=Type),
             size=5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw()+
  scale_colour_manual(values=bioCol)+
  scale_shape_manual(values=15:22)



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


########
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


mito=list(gene)
names(mito)='Hedgehog'

library(GSVA)
exprSet=as.matrix(exprSet)
score=gsva(exprSet,mito,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
score
max(score)
min(score)


score=as.data.frame(t(score))

score$group=group_list
library(ggpubr)

ggboxplot(score,x='group',y='Hedgehog',
          color = "black",
          fill = "group",
          palette=c("lightblue","#fe8180"))+stat_compare_means()
