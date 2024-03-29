#引用包
library(randomForest)
rm(list = ls())
set.seed(123)
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
#随机森林树

rf=randomForest(as.factor(group)~., data=data, ntree=1000)
pdf(file="森林.pdf", width=5, height=5)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要性
importance=importance(x=rf2)

#绘制基因的重要性图
pdf(file="GeneIm.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

#挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>2])     #挑选重要性评分大于2的基因
#rfGenes=names(rfGenes[1:30])         #挑选重要性评分最高的30个基因
write.table(rfGenes, file="随机森林Genes.txt", sep="\t", quote=F, col.names=F, row.names=F)

#输出重要基因的表达量
sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="imGeneExp.txt", sep="\t", quote=F, col.names=F)