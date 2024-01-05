rm(list = ls())
load('key_train_exprSet.Rdata')
rt=rt[,anno$group=='OA']
exprSet<-rt
cluster<-read.table("cluster.txt")
rownames(cluster)<-cluster$V1
#load("F:/顾客/OA_pan/key_train_exprSet.Rdata")
group_list<-as.character(cluster$V2)
group_list=factor(group_list,levels=c('1','2'))
########
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]


######自定义
######自定义
#gene<-read.table("pan.txt")$V1
#gene<-read.table("NR.txt")$V1


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


library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
rm(list = ls())

# 读入前面生成的通路中的基因列表
(load("hallmark.gs.RData")) #保存在当前文件夹
load('key_train_exprSet.Rdata')

rt=rt[,anno$group=='OA']
#exprSet<-rt
cluster<-read.table("cluster.txt")
rownames(cluster)<-cluster$V1

exp<-rt
gsym.expr <- exp
head(gsym.expr)

# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(gsym.expr), gs)

# 预览GSVA分析返回的矩阵
head(gsva_es)

# 把通路的表达量保存到文件
write.csv(gsva_es, "gsva_output_cluster.csv", quote = F)
#```

#到这里，GSVA分析就结束了。

#￥接下来，你可能要对这些通路进行进一步的分析，此处以差异表达分析为例。

## 通路的差异表达分析

#```{r}
# 分组
group_list <- cluster
head(group_list)

group_list$V2<-as.character(group_list$V2)
group_list<-as.data.frame(group_list)
group_list$V2[group_list$V2 == '1']<- c("cluster1") #将Group列所有case换成0
group_list$V2[group_list$V2 == '2']<- c("cluster2") #将Group列所有case换成0
group_list<-group_list[,-1]
# 设置对比
design <- model.matrix(~ 0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(gsva_es)
design

# 构建差异比较矩阵
#contrast.matrix <- makeContrasts(cluster1-cluster2, levels = design)
contrast.matrix <- makeContrasts(cluster2-cluster1, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#把通路的limma分析结果保存到文件
write.csv(x, "gsva_limma_cluster.csv", quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, "easy_input2_for39bar_cluster.csv", quote = F, row.names = F)
#```

## 开始画图

#t为柱子的长度（图中横轴），行名为图中pathway的名称。

#如果想做更多细节上的调整，可参考FigureYa39bar。

#```{r, fig.width=6, fig.height=8}
df <- read.csv("easy_input2_for39bar_cluster.csv")
head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('lightblue', 'snow3', '#fe8180'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= -0.01, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","black","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, cluster2 \n versus cluster1")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴


ggsave("gsva_cluster.pdf", width = 6, height = 8)
