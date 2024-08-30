# PATHWAY GC
## Step 0: 设置路径 
setwd("/Users/moonly/Desktop/Pan-cancer/pan3")

rm(list = ls())

gc()

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

## Step1:取泛癌TOP5000的基因表达数据

#输入pancan表达数据
![image](https://github.com/user-attachments/assets/312b1092-4260-491b-9b4f-3f70946affd0)

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

![2 exp(rowSums(is na) 0)](https://github.com/user-attachments/assets/4ea84f9d-3d69-4e1f-9eec-ed99489d49bc)

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

SD_exp_order1 = exp_SD[order(-exp_SD$samp_SD),]

exp_5000 = SD_exp_order1[1:5000,]#保证匹配为5000

exp_5000_matrix = exp_ratio[(rownames(exp_ratio) %in% exp_5000[,2]),]

write.csv(exp_5000_matrix,"exp_5000_matrix.csv")
![image](https://github.com/user-attachments/assets/fbacbea2-f8d4-495c-9ce2-b3df624cfcdc)

rm(na_counts)

## Step2:spearman
result1 <- cor(t(exp_5000_matrix), method = "pearson", use = "complete.obs")

result1[1:3,1:3]

write.csv(result1,"spearman results1.csv")
![4 spearman_result1](https://github.com/user-attachments/assets/1d0c49b7-8c7f-4494-9b92-99767f75cdc3)

## Step3:Consensus聚类


library(ConsensusClusterPlus)

d = as.matrix(result1)

class(d)

result2 <- ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title="pan3(5000)",clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="pdf",writeTable=TRUE)

#Consensus 聚类结果会储存在一个文件夹中
![image](https://github.com/user-attachments/assets/f5916721-1ca3-4ae6-824b-98667609aa8e)

#根据聚类结果确认最优类为8类
![image](https://github.com/user-attachments/assets/2915df03-4d8a-43e4-883f-071775c6cd0f)

#在Consensus聚类结果文件夹中找到k=8的基因聚类结果
![image](https://github.com/user-attachments/assets/d3bbd222-893a-464a-a796-eea481bd75f9)

#将上述聚类结果整理成genesets表格，格式如下：
![image](https://github.com/user-attachments/assets/ab1aa94d-8891-4681-bced-fd218009870d)


## Step4:根据聚类结果构建基因集
#ID 转换

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
#将以上8类基因集中的每一类基因集进行基因集富集分析（GSEA），同类合并，整理结果如下，显示为6类（基质、发育、免疫、细胞循环、神经、代谢）通路基因集：
![image](https://github.com/user-attachments/assets/114f5250-0a91-4e51-859f-1f44a6d34f9a)

#随机挑选每类基因集里面的5类通路，从GOBP中找到每类通路所对应的基因，构建基因集，得到30个通路的基因集矩阵，示例如下：
![image](https://github.com/user-attachments/assets/a991dd1f-1de4-4500-8c7d-6b59ed8e6d3d)


## Step5:ssgsea打分 
#接下来，我们对通路富集矩阵与exp表达数据进行ssgsea打分
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

```R

for(i in 1:dim(genesets)[1]){
  if (genesets[i,dim(genesets)[2]]=="") {geneset=list(as.character(genesets[i,2:length(genesets[i,])])[-which(genesets[i,2:length(genesets[i,])]=="")])}
  else geneset=list(as.character(genesets[i,2:length(genesets[i,2:length(genesets[i,])])]))
  names(geneset)=genesets[i,1]
  results=c(results,geneset)
}
```

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
![image](https://github.com/user-attachments/assets/558a7243-f85e-4675-a5bc-10d3f32845cf)

## Step6:层次聚类

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

```R
heatmap = pheatmap(ssgsea,scale = 'row',cellheight = 12,show_colnames = FALSE,color=colorRampPalette(c("blue2", "white", "red"))(20),legend=F,
                   clustering_distance_cols = "euclidean",cluster_rows=FALSE,annotation_col=anno_col, annotation_colors = ann_colors,
                   clustering_method = "ward.D2",cutree_cols=2)

dev.off()
```
![image](https://github.com/user-attachments/assets/9c73f14c-b21d-48ac-b38f-54febf12e5a1)
