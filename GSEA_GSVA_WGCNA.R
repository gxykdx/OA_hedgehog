rm(list = ls())
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(GSVA)

set.seed(123456)
load('allDiff.Rdata')
deg<-allDiff
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

###3.GSEA分析
GSinput="hedgehog"
geneset=read.table(paste0(GSinput,".txt"),sep="\t",header=F,check.names=F)
geneset$name=rep(GSinput,nrow(geneset))
geneset=geneset[,c(2,1)]
colnames(geneset)=c("term","gene")
rownames(geneset)=geneset$gene
Ensembl_ID <- bitr(geneset$gene, fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
kegmt=data.frame(Ensembl_ID ,geneset[match(Ensembl_ID$SYMBOL,geneset$gene),])
kegmt=kegmt[,c(3,2)]
colnames(kegmt)=c("term","gene")


#开始GSEA富集分析
kk2<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 1,pAdjustMethod = "none") #GSEA分析

#保存GSEA结果
GSEAOUT=as.data.frame(kk2@result)
write.table(GSEAOUT,file="4.GSEAOUT.xls",sep="\t",quote=F,col.names=T)

#保存

library(GseaVis)
kk2<-kk2@result
#自定义配色1：
#which(grepl("REACTOME_SIGNALING_BY_HEDGEHOG",gg$ID))
gseaNb(object = kk2,
       geneSetID = kk2@result$ID[1],
       subPlot = 3,
       addPval = T,
       pvalX = 0.95,
       pvalY = 0.8,
       curveCol = c("#FF9999","#99CC00"),
       htCol = c("#FF9999","#99CC00"),
       rankCol = c("#FF9999", "white", "#99CC00"))

#pdf(file=paste0("4.", "GSEA.pdf"),height=8,width=8)
#gseaplot2(kk2,geneSetID = rownames(kk2@result),
          title = "Ferroptosis-associated",  #设置title
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, #展示上3部分
          pvalue_table = T)# 显示p值
dev.off()




#############################################################################################################5.GSVA-WGCNA

###准备
getwd()
afdir <- paste0(getwd(),"/5.WGCNA")           #路径必须中文
dir.create(afdir)

afgmt=c(GSinput,"NA",geneset$gene)
afgmt=as.matrix(afgmt)
write.table(t(afgmt),file="afGMT.gmt",sep="\t",quote=F,col.names=F,row.names = F)
afgmt=getGmt("afGMT.gmt")

load("F:/顾客/OA_pan/key_train_exprSet.Rdata")
data<-rt

GSVA <- gsva(expr=as.matrix(data), afgmt, kcdf="Gaussian",method = "gsva",parallel.sz=1)
GSVA=t(GSVA)

#输出GSVA得到的ES值
GSVAOUT=rbind(ID=colnames(GSVA),GSVA)
write.table(GSVAOUT,file=paste0("5.","GSVA.ES", ".xls"),sep="\t",quote=F,col.names = F)

sample<-anno
colnames(sample)<-c("V2")

traitData=cbind(sample,GSVA)
#Hedgehog
colnames(traitData)<-c("V2","Hedgehog")

traitData[,3]=traitData[,1]
traitData[,1]=ifelse(traitData[,1]=="control",1,0)
traitData[,3]=ifelse(traitData[,3]=="OA",1,0)
traitData=traitData[,c(1,3,2)]
#修改性状名称
colnames(traitData)=c("control","OA","Hedgehogscore")
cluster<-traitData

gc()
library(WGCNA)
# 可换5000
WGCNA_matrix = t(rt[order(apply(rt,1,mad), decreasing = T)[1:5000],])
datExpr0 <- WGCNA_matrix  ## top mad genes
datExpr0 <- as.data.frame(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # 返回TRUE则继续
# （可选）如果存在太多的缺失值
if (!gsg$allOK){
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## 样本过滤前
sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

## 根据图片挑选cutheight,我们不丢样本了!!!!!
clust = cutreeStatic(sampleTree, cutHeight = 90, minSize = 10)
table(clust) # 0代表切除的，1代表保留的
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

## 更新anno
cluster=cluster[rownames(datExpr),,drop=F]

## 样本过滤后
sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers (after)", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


datTraits =cluster
datExpr=datExpr[rownames(datTraits),]
sampleNames = rownames(datExpr)
# 能全部对上
traitRows = match(sampleNames, rownames(datTraits))  

###power值散点图
#allowWGCNAThreads()
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize = 100000)

dev.off()
par(mfrow = c(1,2))
cex1 = 0.85
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()
###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
# 发现不合适就自定义，此处我自定义!!!!!
#softPower=4
adjacency = adjacency(datExpr, power = softPower)

net = blockwiseModules(datExpr, power = softPower,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PRTOM",
                       verbose = 3)
# 显示模块数量以及各自包含的基因数目
# 0表示未分入任何模块的基因
# 1是最大的模块，往后依次降序排列，分别对应各自模块的基因
table(net$colors)


mergedColors = labels2colors(net$colors)
mergedColors
##手动保存
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用color labels重新计算MEs（Module Eigengenes:模块的第一主成分）
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #（这是重点）计算ME和表型相关性
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# 设置热图上的文字（两行数字：第一行是模块与各种表型的相关系数；
# 第二行是p值）
# signif 取有效数字
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 然后对moduleTraitCor画热图
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


##把turquoise挑出来
# 选择导出模块
module = c('turquoise')
# 选择模块中基因/探针
probes = names(datExpr)
inModule = (moduleColors %in% module)
modProbes_POS = probes[inModule]
#modprobes可以后续分析
modProbes_POS
write.table(modProbes_POS,file ='turquoise.txt',row.names = F,col.names = F,quote=F)


########FGSEA
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'BP')
# 变list
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(BP)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  BP_x <- fgsea(pathways = BP, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

# 获得各亚型特征的通路
## cluster1的marker基因和fgsea
# 再次确认
identical(colnames(rt),rownames(cluster))

library(limma)
# 先分组信息使用cluster1列

group_list=cluster$Cluster1
group_list=ifelse(group_list==0,'other','Cluster1')
group_list=factor(group_list,levels=c('other','Cluster1'))
design=model.matrix(~ group_list)

fit=lmFit(rt,design)
fit=eBayes(fit) 
c1_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

# 使用预置函数直接进行fgsea
library(dplyr)
BP_c1 <- preranked_BP(c1_allDiff)
sig_BP_c1 <- BP_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_BP_c1 <- sig_BP_c1[order(sig_BP_c1$NES,decreasing = T),]

## cluster2类的marker基因和fgsea
# 先分组信息使用cluster2列

group_list=cluster$Cluster2
group_list=ifelse(group_list==0,'other','Cluster2')
group_list=factor(group_list,levels=c('other','Cluster2'))
design=model.matrix(~ group_list)


fit=lmFit(rt,design)
fit=eBayes(fit) 
c2_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1)



# 使用预置函数直接进行fgsea
library(dplyr)
BP_c2 <- preranked_BP(c2_allDiff)
sig_BP_c2 <- BP_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_BP_c2 <- sig_BP_c2[order(sig_BP_c2$NES,decreasing = T),]



## cluster3类的marker基因和fgsea
# 先分组信息使用cluster2列
#group_list=cluster$Cluster3
#group_list=ifelse(group_list==0,'other','Cluster1')
#group_list=factor(group_list,levels=c('other','Cluster1'))
#design=model.matrix(~ group_list)


#fit=lmFit(rt,design)
#fit=eBayes(fit) 
#c3_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1)



# 使用预置函数直接进行fgsea,不用运行
#library(dplyr)
#BP_c3 <- preranked_BP(c3_allDiff)
#sig_BP_c3 <- BP_c3 %>% filter(abs(NES)>1 & pval<0.01)
#sig_BP_c3 <- sig_BP_c3[order(sig_BP_c3$NES,decreasing = T),]




# 合并作图
##A和B合并 suffixes后缀，一个是A，一个是B

merged <- merge(sig_BP_c1[sample(1:nrow(sig_BP_c1),20),c(1,5)], sig_BP_c2[sample(1:nrow(sig_BP_c2),20),c(1,5)], by = "pathway" , all = T,
                suffixes = c(".c1",".c2"))

## 重复合并，合并C
#merged2 <- merge(merged,sig_BP_c3[sample(1:nrow(sig_BP_c3),20),c(1,5)], by = "pathway" , all = T,
suffixes = c(".c1",".mixed"))

#merged3=merged2
merged3=merged
#colnames(merged3) <- c("pathway","Cluster1", "Cluster2",'Cluster3')
colnames(merged3) <- c("pathway","Cluster1", "Cluster2")
merged3[is.na(merged3)] <- 0
merged3=as.data.frame(merged3)
rownames(merged3) <- merged3$pathway
rownames(merged3)<- tolower(rownames(merged3))
merged3[,1] <- NULL

library(pheatmap)
pheatmap(merged3, cluster_cols = F, cluster_rows = T, border_color=NA, 
         cellwidth =30,  color = colorRampPalette(c(rep("Darkblue",1), "white", rep("red",1)))(1000) ,
         breaks = seq(-2,2,length.out = 1000), main = "Biological process enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
#save(c1_allDiff,c2_allDiff,c3_allDiff,file ='alldiff_cluster.Rdata')
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster.Rdata')

####### 输出protemap
load('alldiff_cluster.Rdata')
c2_sig=c2_allDiff[c2_allDiff$logFC>0.5 & c2_allDiff$adj.P.Val<0.05,]
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
gene = bitr(rownames(c2_sig), fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(SYMBOL = rownames(c2_sig),logFC=c2_sig$logFC)

gene_df <- merge(gene,gene_df,by="SYMBOL")
gene_df <- dplyr::distinct(gene_df,UNIPROT,.keep_all=TRUE)

prote=gene_df[,-1]
write.table(prote,file ='prote_c2.tsv',sep = '\t',row.names = F,col.names = F,quote = F)


c3_sig=c3_allDiff[c3_allDiff$logFC>0.5 & c3_allDiff$adj.P.Val<0.05,]
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
gene = bitr(rownames(c3_sig), fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(SYMBOL = rownames(c3_sig),logFC=c3_sig$logFC)

gene_df <- merge(gene,gene_df,by="SYMBOL")
gene_df <- dplyr::distinct(gene_df,UNIPROT,.keep_all=TRUE)

prote=gene_df[,-1]
write.table(prote,file ='prote_c3.tsv',sep = '\t',row.names = F,col.names = F,quote = F)
