library(glmnet)
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

x=as.matrix(rt)
#y=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(rt))
y<-group

#2.岭回归
#在glmnet函数中设置参数alpha=0即执行岭回归。

#岭回归
ridge <- glmnet(x, y, family="binomial", nlambda=1000, alpha=0) 
#10折交叉验证
ridgeCV <- cv.glmnet(x, y, family = "binomial",
                     type.measure = "deviance",  ##type.measure = c("default", 连续型变量选"mse",二分类变量选 "deviance", "class", "auc", "mae", "C")
                     alpha=0,
                     nfolds = 10)
#绘制CV曲线图，选择最佳lambda值
par(mfrow=c(1,2))
plot(ridge, xvar="lambda")
plot(ridgeCV)
dev.off()
#查看交叉验证结果
ridgeCV

#1se模型--提取简洁模型的参数系数
#se_lambda<-lassoCV$lambda.1se #求出最小值一个标准误的λ值
se_lambda<-ridgeCV$lambda.min  #求出最优值一个标准误的λ值
#se_coef<-coef(lassoCV, s = "lambda.1se")##λ=最小值一个标准误时各变量的系数值
se_coef<-coef(ridgeCV, s = "lambda.min")##λ=最小值一个标准误时各变量的系数值
se_coef
index<-which(se_coef!=0)#非零系数
coef<-se_coef[index][-1]#对应回归系数
diffvariables=row.names(se_coef)[index][-1]#非零变量
#lasso.result.se<-cbind(diffvariables,coef)#输出结果
write.table(diffvariables, file="ridge.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
save(diffvariables,file = "ridge.Rdata")
