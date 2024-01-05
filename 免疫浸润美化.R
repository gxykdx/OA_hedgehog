library(IOBR)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)

###细胞比例
cibersort=as.data.frame(ciber)
#rownames(cibersort)=cibersort[,1]
#immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
immune<-cibersort
data=as.matrix(immune)
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=anno
#colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$group),]
gaps=c(1, as.vector(cumsum(table(data$group))))
xlabels=levels(factor(data$group))
data$group=factor(data$group, levels=c(cluster[1,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data,id.vars=c("group","GSM"))
colnames(data)=c("group","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
#定义细胞浸润的颜色
colaa=colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell cycle phase',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # 图例："left" 左, "right" 右,  "bottom" 下, "top" 上
  scale_fill_manual(values=colaa)+
  theme_bw()+
  xlab(NULL)+
  theme(axis.text.x  = element_blank())+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))+
  facet_grid(. ~ group,scales="free")+
  theme(strip.text.x = element_text(size = 20,colour = "white"))+       #分面字体颜色
  theme(strip.background.x = element_rect(fill = c("grey"), colour = "black")) #分面颜色
ggsave("immune.ration.pdf",width = 15,height = 8)        #输出图片



###细胞比列差异箱式

###细胞比列差异箱式
cibersort=as.data.frame(ciber)
#rownames(cibersort)=cibersort[,1]
#immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
immune<-cibersort
data=as.matrix(immune)
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=anno
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data=melt(data,id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
data$Type=factor(data$Type, c(cluster[1,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression",  fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme_bw()+
  rotate_x_text(50)

pdf(file="immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()


########################################################
library(IOBR)
library(limma)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)

###细胞比例
cibersort=as.data.frame(mcp)
#rownames(cibersort)=cibersort[,1]
#immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
immune<-cibersort
data=as.matrix(immune)
colnames(immune)=gsub("_MCPcounter"," ",colnames(immune))
colnames(data)=gsub("_MCPcounter"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=anno
#colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$group),]
gaps=c(1, as.vector(cumsum(table(data$group))))
xlabels=levels(factor(data$group))
data$group=factor(data$group, levels=c(cluster[1,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data,id.vars=c("group","GSM"))
colnames(data)=c("group","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
#定义细胞浸润的颜色
colaa=colaa=col=c(pal_aaas()(10),pal_jco()(10),pal_nejm()(8),pal_d3()(10),pal_jama()(6))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell cycle phase',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # 图例："left" 左, "right" 右,  "bottom" 下, "top" 上
  scale_fill_manual(values=colaa)+
  theme_bw()+
  xlab(NULL)+
  theme(axis.text.x  = element_blank())+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))+
  facet_grid(. ~ group,scales="free")+
  theme(strip.text.x = element_text(size = 20,colour = "white"))+       #分面字体颜色
  theme(strip.background.x = element_rect(fill = c("grey"), colour = "black")) #分面颜色
ggsave("immune.ration-mcp.pdf",width = 15,height = 4)        #输出图片



###细胞比列差异箱式

###细胞比列差异箱式
cibersort=as.data.frame(mcp)
cibersort=cibersort[,! colnames(cibersort) %in% c("Endothelial_cells_MCPcounter","Fibroblasts_MCPcounter" )]
#rownames(cibersort)=cibersort[,1]
#immune=cibersort[,-1]
#immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
immune<-cibersort
data=as.matrix(immune)
colnames(immune)=gsub("_MCPcounter"," ",colnames(immune))
colnames(data)=gsub("_MCPcounter"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=anno
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data=melt(data,id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
group=levels(factor(data$Type))
data$Type=factor(data$Type, c(cluster[1,1],cluster[nrow(cluster),1]))
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression",  fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Type", 
                  width=0.8,
                  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+
  stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme_bw()+
  rotate_x_text(50)

pdf(file="immune.diff-mcp.pdf", width=5, height=5)
print(boxplot)
dev.off()

