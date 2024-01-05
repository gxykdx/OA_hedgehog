rm(list = ls())
load('key_train_exprSet.Rdata')
cluster=read.table('./cluster.txt',row.names = 1)
colnames(cluster)='Cluster'

cluster$Cluster=paste0('Cluster',cluster$Cluster)

rt=rt[,rownames(cluster)]

cluster$Cluster1=ifelse(cluster$Cluster=='Cluster1',1,0)
cluster$Cluster2=ifelse(cluster$Cluster=='Cluster2',1,0)
#cluster$Cluster3=ifelse(cluster$Cluster=='Cluster3',1,0)
cluster$Cluster=NULL


########FGSEA
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "C2")
# 变list
KEGG <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(KEGG)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_KEGG <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  KEGG_x <- fgsea(pathways = KEGG, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000)
  
  KEGG_x$pathway<-gsub("KEGG_","",KEGG_x$pathway)
  KEGG_x$pathway<-gsub("_"," ",KEGG_x$pathway)
  return(KEGG_x)
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
KEGG_c1 <- preranked_KEGG(c1_allDiff)
sig_KEGG_c1 <- KEGG_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_KEGG_c1 <- sig_KEGG_c1[order(sig_KEGG_c1$NES,decreasing = T),]

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
KEGG_c2 <- preranked_KEGG(c2_allDiff)
sig_KEGG_c2 <- KEGG_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_KEGG_c2 <- sig_KEGG_c2[order(sig_KEGG_c2$NES,decreasing = T),]



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

merged <- merge(sig_KEGG_c1[sample(1:nrow(sig_KEGG_c1),20),c(1,5)], sig_KEGG_c2[sample(1:nrow(sig_KEGG_c2),20),c(1,5)], by = "pathway" , all = T,
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
         breaks = seq(-2,2,length.out = 1000), main = "KEGG enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
#save(c1_allDiff,c2_allDiff,c3_allDiff,file ='alldiff_cluster.Rdata')
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster_kegg.Rdata')


########FGSEA
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'CC')
# 变list
CC <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(CC)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_CC <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  CC_x <- fgsea(pathways = CC, 
                  stats = ranks,
                  minSize=10,
                  maxSize=500,
                  nperm=1000)
  
  CC_x$pathway<-gsub("GOCC_","",CC_x$pathway)
  CC_x$pathway<-gsub("_"," ",CC_x$pathway)
  return(CC_x)
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
CC_c1 <- preranked_CC(c1_allDiff)
sig_CC_c1 <- CC_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_CC_c1 <- sig_CC_c1[order(sig_CC_c1$NES,decreasing = T),]

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
CC_c2 <- preranked_CC(c2_allDiff)
sig_CC_c2 <- KEGG_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_CC_c2 <- sig_CC_c2[order(sig_CC_c2$NES,decreasing = T),]



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

merged <- merge(sig_CC_c1[sample(1:nrow(sig_CC_c1),20),c(1,5)], sig_CC_c2[sample(1:nrow(sig_CC_c2),20),c(1,5)], by = "pathway" , all = T,
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
         breaks = seq(-2,2,length.out = 1000), main = "Cellular component enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
#save(c1_allDiff,c2_allDiff,c3_allDiff,file ='alldiff_cluster.Rdata')
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster_Cellular_Component.Rdata')




########FGSEA
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'MF')
# 变list
MF <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(MF)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_MF <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  MF_x <- fgsea(pathways = MF, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000)
  
  MF_x$pathway<-gsub("GOMF_","",MF_x$pathway)
  MF_x$pathway<-gsub("_"," ",MF_x$pathway)
  return(MF_x)
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
MF_c1 <- preranked_MF(c1_allDiff)
sig_MF_c1 <- MF_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_MF_c1 <- sig_MF_c1[order(sig_MF_c1$NES,decreasing = T),]

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
MF_c2 <- preranked_MF(c2_allDiff)
sig_MF_c2 <- MF_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_MF_c2 <- sig_MF_c2[order(sig_MF_c2$NES,decreasing = T),]



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

merged <- merge(sig_MF_c1[sample(1:nrow(sig_MF_c1),20),c(1,5)], sig_MF_c2[sample(1:nrow(sig_MF_c2),20),c(1,5)], by = "pathway" , all = T,
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
merged4<-merged3
merged3<-merged3[-24,]
pheatmap(merged3, cluster_cols = F, cluster_rows = T, border_color=NA, 
         cellwidth =30,  color = colorRampPalette(c(rep("Darkblue",1), "white", rep("red",1)))(1000) ,
         breaks = seq(-2,2,length.out = 1000), main = "Molecular function enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
#save(c1_allDiff,c2_allDiff,c3_allDiff,file ='alldiff_cluster.Rdata')
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster_Molecular_function.Rdata')


########FGSEA
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "H")
# 变list
hallmak <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(hallmak)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_hallmak <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  hallmak_x <- fgsea(pathways = hallmak, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000)
  
  hallmak_x$pathway<-gsub("HALLMARK_","",hallmak_x$pathway)
  hallmak_x$pathway<-gsub("_"," ",hallmak_x$pathway)
  return(hallmak_x)
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
hallmak_c1 <- preranked_hallmak(c1_allDiff)
sig_hallmak_c1 <- hallmak_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_hallmak_c1 <- sig_hallmak_c1[order(sig_hallmak_c1$NES,decreasing = T),]

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
hallmak_c2 <- preranked_hallmak(c2_allDiff)
sig_hallmak_c2 <- hallmak_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_hallmak_c2 <- sig_hallmak_c2[order(sig_hallmak_c2$NES,decreasing = T),]



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

merged <- merge(sig_hallmak_c1[sample(1:nrow(sig_hallmak_c1),20),c(1,5)], sig_hallmak_c2[sample(1:nrow(sig_hallmak_c2),20),c(1,5)], by = "pathway" , all = T,
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
         breaks = seq(-2,2,length.out = 1000), main = "HALLMARK enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
#save(c1_allDiff,c2_allDiff,c3_allDiff,file ='alldiff_cluster.Rdata')
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster_hallmark.Rdata')
