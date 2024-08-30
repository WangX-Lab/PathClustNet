##########################################1.algorithm#############################################################
setwd("/Users/moonly/Desktop/Pan-cancer/pan3")
rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#Step1:Top5000
#输入pancan表达数据
expr <- read.csv("/Users/moonly/Desktop/Pan-cancer/pan_tumor_exp.csv",header = T, stringsAsFactors = F)
expr[1:3,1:3]
rownames(expr) <- expr[,1]
expr <- expr[,-1]

na_counts <- rowSums(is.na(expr))
exp6 <- expr[na_counts <= 100, ]
dim(exp6)
exp7 <- expr[na_counts <= 10, ]
dim(exp7)
exp8 <- expr[na_counts == 0, ]
dim(exp8)
table(is.na(exp8))

exp7 <- na.omit(exp7)

write.csv(exp8,"exp(rowSums(is.na)>0).csv")
rm(exp7)
#筛选top5000
exp_ratio = exp8
SD = apply(exp_ratio,1 ,function(x) sd(x,na.rm=TRUE))#计算每一条通路的SD，用于筛选top通路
write.csv(SD ,"TCGA_ratio_SD_exp.csv")

exp_SD = SD
class(exp_SD)
class(SD)
head(SD)
exp_SD <- as.data.frame(exp_SD)
exp_SD[,2] <- rownames(exp_SD)
names(exp_SD) = c("samp_SD","exp")
SD_exp_order1 = exp_SD[order(-exp_SD$samp_SD),]#降序
#取---------------top5000的通路---------------------------------------
exp_5000 = SD_exp_order1[1:5000,]#保证匹配为5000
exp_5000_matrix = exp_ratio[(rownames(exp_ratio) %in% exp_5000[,2]),]
write.csv(exp_5000_matrix,"exp_5000_matrix.csv")

rm(na_counts)

#Step2:spearman
result1 <- cor(t(exp_5000_matrix), method = "pearson", use = "complete.obs")
result1[1:3,1:3]
write.csv(result1,"spearman results1.csv")

#Step3:Consensus聚类


library(ConsensusClusterPlus)

d = as.matrix(result1)
class(d)

result2 <- ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title="pan3(5000)",clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

#Step4:根据聚类结果构建基因集
####ID 转换
library(org.Hs.eg.db)
k=keys(org.Hs.eg.db,keytype = "ENSEMBL")
head(k,5)
list=select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
dim(list)

head(list,5)

ID <- read.csv("ID .csv",header = T, stringsAsFactors = F)
ID
ID_list=list[match(ID$ENSEMBL,list$ENSEMBL), ]
ID_list

write.csv(ID_list,"ID_trans.csv")

# Step5:ssgsea打分 ----------------------------------------------------------------

rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/")
library("BiocManager")
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")


library(GSVA)

genesets = read.csv("genesets from GOBP(ID)3.csv", stringsAsFactors = FALSE,header = FALSE)
colnames(genesets) = genesets[1,]
genesets = as.data.frame(t(genesets))
results=c()  

for(i in 1:dim(genesets)[1]){
  if (genesets[i,dim(genesets)[2]]=="") {geneset=list(as.character(genesets[i,2:length(genesets[i,])])[-which(genesets[i,2:length(genesets[i,])]=="")])}
  else geneset=list(as.character(genesets[i,2:length(genesets[i,2:length(genesets[i,])])]))
  names(geneset)=genesets[i,1]
  results=c(results,geneset)
}

geneSet = results

rna_seq = read.csv("exp(rowSums(is.na)>0).csv",stringsAsFactors = F,header = T,check.names = F,row.names = 1)

rna_seq = as.matrix(rna_seq)
rna_seq[1:3,1:3]
res = gsva(rna_seq,geneSet,method="ssgsea",ssgsea.norm = TRUE,verbose = TRUE)

colnames(res) = gsub(colnames(res),pattern=".",replacement="-",fixed = TRUE)

#ssgsea_score = t(res)

res = as.data.frame(res)
res[1:3,1:3]
write.csv(res, "ssgsea_pathway.csv")

# Step6:层次聚类

rm(list = ls())

ssgsea = read.csv("ssgsea_pathway.csv",header=F,encoding="UTF-8")
ssgsea[1:3,1:2]
rownames(ssgsea) = as.character(ssgsea[,1])

ssgsea = ssgsea[,-1]

colnames(ssgsea) = t(ssgsea[1,])

ssgsea = ssgsea[-1,]
#ssgsea = res
ssgsea = as.matrix(ssgsea)

#将表达值转换为数值型

ssgsea = matrix(as.numeric(ssgsea),nrow=nrow(ssgsea),ncol=ncol(ssgsea),dimnames=list(rownames(ssgsea),colnames(ssgsea)))

#scale标准化

df=scale(t(ssgsea))  ##由于scale是对每一列标准化，要对变量标准化，所以这里要转置

d=dist(df,method = "euclidean")

sample.hc = hclust(d,method="ward.D2")


sample.id <- cutree(sample.hc,4)   #k

write.csv(sample.id,"cluster_2list.csv")

b = read.csv("cluster_2list.csv",stringsAsFactor=FALSE)

library(pheatmap)

b[,2] = paste("cluster",b[,2],sep = "")

anno_col = data.frame(cluster=factor(b[,2]))

rownames(anno_col)=as.character(b[,1])
table(anno_col)
ann_colors = list(cluster = c(cluster1="#d9a62e", cluster2="#824880", cluster3="#cd6234", cluster4="#F8C9C1"))#cluster5="#C0BBBE"

pdf("heatmap.pdf",width=15,height=10)

heatmap = pheatmap(ssgsea,scale = 'row',cellheight = 12,show_colnames = FALSE,color=colorRampPalette(c("blue2", "white", "red"))(20),legend=F,
                   clustering_distance_cols = "euclidean",cluster_rows=FALSE,annotation_col=anno_col, annotation_colors = ann_colors,
                   clustering_method = "ward.D2",cutree_cols=2)

dev.off()