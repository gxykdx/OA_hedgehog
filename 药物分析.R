rm(list = ls())
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(stringr)
library(stringi)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(patchwork)
library(ggplotify)
load('key_train_exprSet.Rdata')
cluster=read.table('./cluster.txt',row.names = 1)
colnames(cluster)='Cluster'

cluster$Cluster=paste0('Cluster',cluster$Cluster)

rt=rt[,rownames(cluster)]
exprSet<-rt
group_list<-cluster$Cluster
group_list=factor(group_list,levels=c('Cluster2','Cluster1'))
#cluster$Cluster1=ifelse(cluster$Cluster=='Cluster1',1,0)
#cluster$Cluster2=ifelse(cluster$Cluster=='Cluster2',1,0)

library(limma)
design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#tumor处需修改case
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 
write.csv(allDiff,"cluster_Diff.csv")


library(ggplot2)

#alldiff_mito$label=ifelse(alldiff_mito$adj.P.Val<0.05,rownames(alldiff_mito),'')


ggplot(allDiff,aes(logFC, -log10(adj.P.Val)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  #geom_vline(xintercept = c(-0.2,0.2,0.5,-0.5), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(0), linetype = "dashed", color = "#999999")+
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
  #geom_text(aes(label=label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust=1)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)")

library(PharmacoGx)
library(parallel)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 提取药物处理矩阵
load("CMAP_gene_signatures.RData")
camp_sig <- CMAP.genePerturbations[,,c("tstat")] %>% data.frame()

# 基因名转换
camp_sig$ENTREZID <- do.call(rbind, strsplit(rownames(camp_sig),'\\.'))[,2]

SYMBOL <- bitr(camp_sig$ENTREZID, fromType = "ENTREZID",
               toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

camp_sig <- merge(SYMBOL, camp_sig, by = "ENTREZID")
camp_sig <- column_to_rownames(camp_sig, var = "SYMBOL"); camp_sig <- camp_sig[,-1]

# 保存数据
saveRDS(camp_sig, "camp_sig.rds")


dis_sig<-allDiff
# 提取变化倍数最大的基因
# 例文（Yang et al）提到疾病分子特征数量选择100时可以获得较好的预测性能。但这篇研究是基于LINCS数据，cmap数据的维度与LINCS差别明显。这里我们建议疾病分子特征数量可以稍微多一点，以下演示使用top300基因进行XSum分析
dis_sig <- rbind(top_n(dis_sig, 150, logFC), top_n(dis_sig, -150, logFC))

# 将logfc转成二分类变量（logfc>0的基因用1表示，logfc小于0的基因用-1表示）
# 使用XSum时不需要考虑差异基因的差异倍数，这步分析是为了让大家更好的理解，并不是并要的

dis_sig$logFC[dis_sig$logFC>0] <- 1; dis_sig$logFC[dis_sig$logFC<0] <- -1
dis_sig$id<-rownames(dis_sig)
rownames(dis_sig) <- NULL

# 保存结果
write.table(dis_sig, "dis_sig.csv", sep=",", quote=F, row.names=F, col.names=T)





  
rm(list = ls())
# 读入drug signature
camp_sig <- readRDS("camp_sig.rds")

# 读入disease signature
dis_sig <- read.csv('dis_sig.csv', sep=',', header=TRUE)

# 读入XSum函数
source("Core_function.R")

# 选择XSum的topN（Yang et al的研究提到topN选择200效果可能比较好，但这个结论可能不适用与cmap的数据，这里我们选择topN = 500）
XLogFC <- eXtremeLogFC(camp_sig, N = 500)

up_gene <- dis_sig$id[dis_sig$logFC == 1]
dn_gene <- dis_sig$id[dis_sig$logFC == -1]

xsum <- data.frame(score=XSum(XLogFC, up_gene, dn_gene))
xsum <- rownames_to_column(xsum, var = "id")

# 把结果标准化至-1到1（这步也可不做）
#############################################################################
xsum_pos <- xsum[xsum$score>0,]
xsum_pos$score <- xsum_pos$score/max(xsum_pos$score)

xsum_neg <- xsum[xsum$score<0,]
xsum_neg$score <- xsum_neg$score/min(xsum_neg$score) * -1

xsum <- rbind(xsum_pos, xsum_neg)
#############################################################################


# 将结果从低到高排序
xsum <- xsum[order(xsum$score),]
head(xsum)

xsum$number <- 1:nrow(xsum)

# 突出显示top5的药物，标出药物名
select <- xsum[1:10,]

###基于CMap的理论（已经被大量研究验证），这里我们得到的分数越低，这个药物越有可能逆转疾病的分子特征，
#理论上更有可能具有治疗该疾病的能力。该演示数据结果提示 X4.5.dianilinophthalimide 是最有可能治疗肝癌的药物

# 开始画图
ggplot(xsum, aes(number,score))+
  geom_point(size=0.5, color="grey50") + 
  geom_point(data = select, alpha = 1, 
             size = 2, color = "#5ec7dd") + 
  geom_label_repel(data = select, aes(label=id), 
                   color = "white",
                   alpha = 1, point.padding = 1, 
                   size = 5, fill = "#009bc7",
                   segment.size = 1, nudge_x=-0.5, 
                   segment.color = "grey50",
                   direction = "x",
                   hjust = 1) + 
  theme_classic()

ggsave("CMAP_XSum.pdf", width = 9, height = 5)
