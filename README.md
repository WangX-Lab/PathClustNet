# PathClustNet
## Step 0: 设置路径 
```R
setwd("/Users/moonly/Desktop/Pan-cancer/pan3")

rm(list = ls())

gc()

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
```

## Step1:取泛癌TOP5000的基因表达数据

输入pancan表达数据
![image](https://github.com/user-attachments/assets/9ab10089-50cd-4635-bd8b-cabecf3b62aa)

```R

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

```

![image](https://github.com/user-attachments/assets/e4d3fb26-a4ad-4c66-85c9-487a35c420b1)

筛选top5000

```R
rm(exp7)
exp_ratio = exp8
SD = apply(exp_ratio,1 ,function(x) sd(x,na.rm=TRUE))
write.csv(SD ,"TCGA_ratio_SD_exp.csv")
exp_SD = SD
class(exp_SD)
class(SD)
head(SD)
exp_SD <- as.data.frame(exp_SD)
exp_SD[,2] <- rownames(exp_SD)
names(exp_SD) = c("samp_SD","exp")
SD_exp_order1 = exp_SD[order(-exp_SD$samp_SD),]
exp_5000 = SD_exp_order1[1:5000,]
exp_5000_matrix = exp_ratio[(rownames(exp_ratio) %in% exp_5000[,2]),]
write.csv(exp_5000_matrix,"exp_5000_matrix.csv")
```

![image](https://github.com/user-attachments/assets/37b37193-68ad-4f93-856f-88792d845897)

## Step2:spearman

```R

result1 <- cor(t(exp_5000_matrix), method = "pearson", use = "complete.obs")
result1[1:3,1:3]
write.csv(result1,"spearman results1.csv")

```
![image](https://github.com/user-attachments/assets/7b0bccfb-2b4b-4e4e-9df7-9733eef0e1e0)

## Step3:Consensus聚类

```R

library(ConsensusClusterPlus)
d = as.matrix(result1)
class(d)
result2 <- ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title="pan3(5000)",clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

```

Consensus 聚类结果会储存在一个文件夹中
![image](https://github.com/user-attachments/assets/f246417f-b378-4fd5-a4c1-dfe003c41ce3)



根据聚类结果确认最优类为8类
![image](https://github.com/user-attachments/assets/6e7e91fc-6c69-417b-a7a8-37505bd998e4)


在Consensus聚类结果文件夹中找到k=8的基因聚类结果
![image](https://github.com/user-attachments/assets/3c10fc6d-6c30-47ad-a03a-0aec9aef4319)


将上述聚类结果整理成genesets表格，格式如下：
![image](https://github.com/user-attachments/assets/3bba83b1-fc29-4d46-9945-613d3652712b)


## Step4:根据聚类结果构建基因集
ID 转换

```R

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

```

将以上8类基因集中的每一类基因集进行基因集富集分析（GSEA），同类合并，整理结果如下，显示为6类（基质、发育、免疫、细胞循环、神经、代谢）通路基因集：
![image](https://github.com/user-attachments/assets/90a74147-e49b-4467-9d9f-2f9b0ce00b8a)

随机挑选每类基因集里面的5类通路，从GOBP中找到每类通路所对应的基因，构建基因集，得到30个通路的基因集矩阵，示例如下：
![image](https://github.com/user-attachments/assets/b36da6cc-7845-4ef0-930a-3fe81bad4b17)



## Step5:ssgsea打分 
接下来，我们对通路富集矩阵与exp表达数据进行ssgsea打分

```R

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
res = as.data.frame(res)
res[1:3,1:3]
write.csv(res, "ssgsea_pathway.csv")

```
![image](https://github.com/user-attachments/assets/511bafea-e652-4c13-9cd9-e36a0b99e618)



## Step6:层次聚类

```R

rm(list = ls())
ssgsea = read.csv("ssgsea_pathway.csv",header=F,encoding="UTF-8")
ssgsea[1:3,1:2]
rownames(ssgsea) = as.character(ssgsea[,1])
ssgsea = ssgsea[,-1]
colnames(ssgsea) = t(ssgsea[1,])
ssgsea = ssgsea[-1,]
ssgsea = as.matrix(ssgsea)
ssgsea = matrix(as.numeric(ssgsea),nrow=nrow(ssgsea),ncol=ncol(ssgsea),dimnames=list(rownames(ssgsea),colnames(ssgsea)))
df=scale(t(ssgsea)) 
d=dist(df,method = "euclidean")
sample.hc = hclust(d,method="ward.D2")
sample.id <- cutree(sample.hc,4)
write.csv(sample.id,"cluster_2list.csv")
b = read.csv("cluster_2list.csv",stringsAsFactor=FALSE)
library(pheatmap)
b[,2] = paste("cluster",b[,2],sep = "")
anno_col = data.frame(cluster=factor(b[,2]))
rownames(anno_col)=as.character(b[,1])
table(anno_col)
ann_colors = list(cluster = c(cluster1="#d9a62e", cluster2="#824880", cluster3="#cd6234", cluster4="#F8C9C1"))
pdf("heatmap.pdf",width=15,height=10)
heatmap = pheatmap(ssgsea,scale = 'row',cellheight = 12,show_colnames = FALSE,color=colorRampPalette(c("blue2", "white", "red"))(20),legend=F,
                   clustering_distance_cols = "euclidean",cluster_rows=FALSE,annotation_col=anno_col, annotation_colors = ann_colors,
                   clustering_method = "ward.D2",cutree_cols=2)
dev.off()

```

![image](https://github.com/user-attachments/assets/26ea2adf-26b0-4560-b51a-f6b50be608bf)


