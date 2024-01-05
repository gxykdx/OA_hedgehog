rm(list = ls())
af <- readRDS("F:/gk/OA_hedgehog/af.rds")
scRNA<-af

# 先装包
library(RColorBrewer) 
library(viridis)
#install.packages("wesanderson")
library(wesanderson)
library(Seurat)
####COLOR PALETTE FOR PLOTS ####
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector

col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]



##### QUALITY FILTERING #####
View(scRNA@meta.data)
scRNA$tissue_type=stringr::str_remove(scRNA$orig.ident,pattern = '[0-9]')####将数字删掉
View(scRNA@meta.data)
scRNA$tissue_type=factor(scRNA$tissue_type,levels = c('N','O'))
length(scRNA@active.ident)
#!!!
# 人：^MT-
# 人用下面的！！！！！，不要用^mt-！！！！！！！！！！！：
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt.")
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
RidgePlot(scRNA, features="nFeature_RNA")
RidgePlot(scRNA, features="nCount_RNA")
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
scRNA<- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
length(scRNA@active.ident)


##### PRE-PROCESSING #####
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)


#################################################来自挑圈联靠#############################################
# 这里我们展示一下前10个改变最大的基因
top10 <- head(VariableFeatures(scRNA), 10)
# 对高可变基因进行可视化
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#下面可视化的解读则是，左边和右边的图表达的是一个东西，只不过右边的图把表达差异在各个样本中最显著的基因表示了出来，然后高可变基因用红色表示
#################################################来自挑圈联靠#############################################


####对高变基因进行归一化处理
length(scRNA@active.ident)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")


#################################################来自挑圈联靠#############################################
#当然一个PC内包含的基因，我们可以用热图的形式来展现出来
DimHeatmap(scRNA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 3, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 4, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 5, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 1:12, cells = 500, balanced = TRUE)
#DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)#1：4是热图分成几个，ncol是分为几列。
#################################################来自挑圈联靠#############################################


ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:6)

#################################################来自挑圈联靠#############################################
#接下来我们既然已经对庞大的数据进行了降维（也就是聚堆）的形式，那么我们究竟要选择几个PC来代表这么庞大的数据呢，肯定不能都选，否则我的数据量还是这些，就没有我们前面一直渗透的缩小缩小再缩小的含义了
pbmc <- JackStraw(scRNA, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#上面两个代码是通过不同的方式来帮助我们选择PC的数目，并且分别都对应不同的可视化
JackStrawPlot(pbmc, dims = 1:20)
#################################################来自挑圈联靠#############################################


#scRNA_1 <- FindNeighbors(scRNA, dims = 1:5)
#scRNA_2 <- FindNeighbors(scRNA, dims = 1:6)
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:6)
DimPlot(scRNA, order=T, group.by = "orig.ident", label=F, cols = col_vector)
DimPlot(scRNA, order=T, group.by = "seurat_clusters", label=T, cols = col_vector)
DimPlot(scRNA,split.by = 'tissue_type', cols = col_vector)

# 人种注意都大写不要小写！！！！！！！！！！！！！！！！！！！！！
FeaturePlot(scRNA,features = 'GAS1',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
#FeaturePlot(scRNA,features = 'GABARAPL1',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
VlnPlot(scRNA,features = 'WNT5A')
VlnPlot(scRNA,features = 'WNT5B')
Idents(scRNA)=scRNA$seurat_clusters
VlnPlot(scRNA,features = 'WNT5A')

library(patchwork)
p1<-FeaturePlot(scRNA,features = 'GLI3',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p2<-FeaturePlot(scRNA,features = 'CSNK1G2',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p3<-FeaturePlot(scRNA,features = 'WNT5A',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p4<-FeaturePlot(scRNA,features = 'WNT5B',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p5<-FeaturePlot(scRNA,features = 'PRKX',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p6<-FeaturePlot(scRNA,features = 'RAB23',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p7<-FeaturePlot(scRNA,features = 'FBXW11',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p8<-FeaturePlot(scRNA,features = 'PRKACA',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
p9<-FeaturePlot(scRNA,features = 'BMP8B',order = T,label = T,pt.size = 1,split.by = 'tissue_type')
#p1+p2+p3+p4+p5+p6+p7+p8+p9
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9)


save(scRNA,file ='scRNA.Rdata')


#rm(list = ls())
load('scRNA.Rdata')

#BiocManager::install("SingleR")
#BiocManager::install("celldex")
library(celldex)
library(SingleR)
refdata <- SingleR::HumanPrimaryCellAtlasData()
#refdata <- SingleR::MouseRNAseqData()
library(Seurat)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

Idents(scRNA)=scRNA$celltype
DimPlot(scRNA,split.by = 'tissue_type', cols = col_vector)
saveRDS(scRNA,file ='scRNA_anno.RDS')


rm(list = ls())
scRNA=readRDS('./scRNA_anno.RDS')
DimPlot(scRNA,split.by = 'tissue_type', cols = col_vector)
library(RColorBrewer) 
library(viridis)
library(wesanderson)
####COLOR PALETTE FOR PLOTS ####
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector

col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]


###
##gene=read.table('../REACTOME_MITOPHAGY.v7.5.1.gmt')
##gene=as.data.frame(t(gene))
##gene=gene[3:31,]

gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]



# 以下是为了人基因转鼠基因，人单细胞不用跑下面三行
# ！！！！！！人种不要多跑！！！！！！！！！！！
#library(Hmisc)
#gene=tolower(gene)
#gene=capitalize(gene)


Mitoscore=list(gene)
names(Mitoscore)='Mitoscore'

library(Seurat)
scRNA <- AddModuleScore(object = scRNA, features = Mitoscore, name ='Mitoscore')

library(ggplot2)

VlnPlot(scRNA,features ='Mitoscore1', pt.size = 0,split.by = 'tissue_type') 
FeaturePlot(scRNA,features = 'Mitoscore1',order = F,split.by = 'tissue_type',label = T,cols =viridis(10))


## 显著性
epi<- subset(scRNA, celltype %in% c("Chondrocytes"))
epi_sham<- subset(epi, tissue_type %in%  c("N"))
epi_MCAO<- subset(epi, tissue_type %in%  c("O"))
wilcox.test(epi_sham$Mitoscore1, epi_MCAO$Mitoscore1, alternative = "two.sided") #W = 210385, p-value = 0.0114

## 显著性
epi<- subset(scRNA, celltype %in% c("Tissue_stem_cells"))
epi_sham<- subset(epi, tissue_type %in%  c("N"))
epi_MCAO<- subset(epi, tissue_type %in%  c("O"))
wilcox.test(epi_sham$Mitoscore1, epi_MCAO$Mitoscore1, alternative = "two.sided") #W = 353, p-value = 0.5471


## 显著性
epi<- subset(scRNA, celltype %in% c("Macrophage"))
epi_sham<- subset(epi, tissue_type %in%  c("N"))
epi_MCAO<- subset(epi, tissue_type %in%  c("O"))
wilcox.test(epi_sham$Mitoscore1, epi_MCAO$Mitoscore1, alternative = "two.sided") #W = 107, p-value = 0.7825


## 显著性
#endo<- subset(scRNA, celltype %in% c("Endothelial cells"))
#endo_sham<- subset(endo, tissue_type %in%  c("sham"))
#endo_MCAO<- subset(endo, tissue_type %in%  c("MCAO"))
#wilcox.test(endo_sham$Mitoscore1, endo_MCAO$Mitoscore1, alternative = "two.sided") #p-value < 2.2e-16


DimPlot(scRNA, order=T, group.by = "celltype", split.by = 'tissue_type',label=T, cols = col_vector)


##4. 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)

#挑出Endothelial cells，这里你可以换

dir.create("pseudotime")
##提取细胞子集，Endothelial cells可以换！！！
scRNAsub <- subset(scRNA, celltype %in% c('Chondrocytes'))
dim(scRNAsub)
set.seed(123456)
#a=sample(1:ncol(scRNAsub),2000,replace = F)
a=sample(1:ncol(scRNAsub),1357,replace = F)
scRNAsub=scRNAsub[,a]
dim(scRNAsub)
##重新降维聚类,提取子集后必须重新降维
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")


#################################################来自挑圈联靠#############################################
#当然一个PC内包含的基因，我们可以用热图的形式来展现出来
DimHeatmap(scRNA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 3, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 4, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 1:12, cells = 500, balanced = TRUE)
#################################################来自挑圈联靠#############################################


pc.num=1:8
##细胞聚类
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub )

##都不用更改代码
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
DimPlot(scRNAsub, reduction = "umap",split.by = 'tissue_type',label = T,cols = col_vector) 
VlnPlot(scRNAsub,features = 'Mitoscore1',cols = col_vector)

VlnPlot(scRNAsub,features = 'Mitoscore1',group.by = 'tissue_type')
#VlnPlot(scRNAsub,features =gene[11:20],group.by = 'tissue_type')

FeaturePlot(scRNAsub,features = 'Mitoscore1',split.by = 'tissue_type',label = T,cols = c('#dadada','#bc3c29'))

## 显著性
#epi<- subset(scRNAsub, celltype %in% c("Chondrocytes"))
#epi_sham<- subset(epi, tissue_type %in%  c("N"))
#epi_MCAO<- subset(epi, tissue_type %in%  c("O"))
#wilcox.test(epi_sham$Mitoscore1, epi_MCAO$Mitoscore1, alternative = "two.sided") #W = 107, p-value = 0.7825


saveRDS(scRNAsub,file ='./pseudotime/scRNA_endo.RDS')

#setwd('./scRNA/')

rm(list = ls())
scRNAsub=readRDS('./pseudotime/scRNA_endo.RDS')
dim(scRNAsub)
## 拟时序分析 先关掉再开
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNAsub@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)


#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1

##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot2

plot3 <-  plot_cell_trajectory(mycds, color_by = "tissue_type")
plot3

plot4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot4

##合并出图
plotc <- plot1|plot2|plot3|plot4
plotc

########
gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]

library(Hmisc)
gene=tolower(gene)
gene=capitalize(gene)

# 展现趋势统一的基因
p1 <- plot_genes_in_pseudotime(mycds['GLI3',], color_by = "seurat_clusters")
p2 <- plot_genes_in_pseudotime(mycds['CSNK1G2',], color_by = "seurat_clusters")
p3 <- plot_genes_in_pseudotime(mycds['WNT5A',], color_by = "seurat_clusters")
p4 <- plot_genes_in_pseudotime(mycds['WNT5B',], color_by = "seurat_clusters")
p5 <- plot_genes_in_pseudotime(mycds['PRKX',], color_by = "seurat_clusters")
p6 <- plot_genes_in_pseudotime(mycds['RAB23',], color_by = "seurat_clusters")
p7 <- plot_genes_in_pseudotime(mycds['FBXW11',], color_by = "seurat_clusters")
p8 <- plot_genes_in_pseudotime(mycds['PRKACA',], color_by = "seurat_clusters")
p9 <- plot_genes_in_pseudotime(mycds['BMP8B',], color_by = "seurat_clusters")
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9)


s.genes <- c("GLI3","CSNK1G2","WNT5A","WNT5B","PRKX","RAB23","FBXW11","PRKACA","BMP8B")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "seurat_clusters", color_by = "seurat_clusters")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "seurat_clusters", color_by = "seurat_clusters")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "seurat_clusters")
plotc <- p1|p2|p3
plotc
#ggsave("pseudotime/genes_visual.png", plot = plotc, width = 8, height = 4.5)







#######################
## 细胞通讯分析

rm(list = ls())
scRNA=readRDS('./scRNA_anno.RDS')

###
##gene=read.table('../REACTOME_MITOPHAGY.v7.5.1.gmt')
##gene=as.data.frame(t(gene))
##gene=gene[3:31,]

gene=read.table('KEGG_HEDGEHOG_SIGNALING_PATHWAY.v2022.1.Hs.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:58,]



# 以下是为了人基因转鼠基因，人单细胞不用跑下面三行
# ！！！！！！人种不要多跑！！！！！！！！！！！
#library(Hmisc)
#gene=tolower(gene)
#gene=capitalize(gene)


Mitoscore=list(gene)
names(Mitoscore)='Mitoscore'

library(Seurat)
scRNA <- AddModuleScore(object = scRNA, features = Mitoscore, name ='Mitoscore')

# 先定义细胞
scRNAendo=subset(scRNA,celltype %in% 'Chondrocytes')
#scRNAendo=subset(scRNA,celltype %in% Chondrocytes)
summary(scRNAendo$Mitoscore1)
scRNAendo$Mitogroup='Mitophagy_median'

for (i in 1:nrow(scRNAendo@meta.data)) {
  if (scRNAendo$Mitoscore1[i]>0.01876) {
    scRNAendo$Mitogroup[i]='Mitophagy_high'
  }
  if (scRNAendo$Mitoscore1[i]<-0.07623) {
    scRNAendo$Mitogroup[i]='Mitophagy_low'
  }
}

scRNA$Mitogroup=''
# 排除内皮的其他细胞
scRNAother=subset(scRNA, celltype != 'Chondrocytes')
#scRNAother<-scRNA
scRNA_chat=merge(scRNAendo,c(scRNAother))

scRNA_chat$Mitocell=paste0(scRNA_chat$Mitogroup,scRNA_chat$celltype)

## 用最新的R
# BiocManager::install("sqjin/CellChat")
library(CellChat)

# 选取病鼠，如果内存不够，再进一步减少细胞数，例如随机抽2000个
scRNA_chat <- subset(scRNA_chat, orig.ident=='O')

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "Mitocell")

#CellChatDB <- CellChatDB.mouse
CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

unique(cellchat@idents)

# 等待
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

write.csv(df.net,file ='cellchat.csv',quote=F)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat <- aggregateNet(cellchat)
dev.off()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
table(scRNA_chat$Mitocell)
p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('Mitophagy_highChondrocytes','Mitophagy_lowChondrocytes'),
                           targets.use = c('Macrophage','Tissue_stem_cells'),
                           
                           remove.isolate = FALSE)+coord_flip()


p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('Mitophagy_lowChondrocytes'),
                           targets.use = c('Macrophage'),
                           
                           remove.isolate = FALSE)+coord_flip()



p_bubble


#必须把上一个图关掉
dev.off()
##netVisual_aggregate(cellchat, signaling = 'GAS',  
                    sources.use = c('Mitophagy_lowChondrocytes'),
                    targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'))

#####去df.net查看通路名称
netVisual_aggregate(cellchat, signaling = 'ANGPTL',  
                    sources.use = c('Mitophagy_lowEndothelial cells'),
                    targets.use = c('Macrophage','Tissue_stem_cells'))


netVisual_aggregate(cellchat, signaling = 'SPP1',  
                    sources.use = c('Mitophagy_lowEndothelial cells'),
                    targets.use = c('Macrophage'))



#必须把上一个图关掉
dev.off()
netVisual_aggregate(cellchat, signaling = 'PTN',  
                    sources.use = c('Mitophagy_highEndothelial cells','Mitophagy_lowEndothelial cells'),
                    targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'))


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = 'ncWNT', width = 8, height = 2.5, font.size = 10)


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

