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

#######################################1.2 survival分析##################################

rm(list = ls())
gc()
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

#setwd("/Users/moonly/Desktop/Pan-cancer/all sample/survival")
library(readxl)
cluster <- read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/cluster.csv",header=TRUE, stringsAsFactors=FALSE,check.names=FALSE)
head(cluster)
cluster[,1] = gsub(cluster[,1],pattern=".",replacement="-",fixed = TRUE)
cluster[,2] = gsub(cluster[,2],pattern=".",replacement="-",fixed = TRUE)
cluster[,3] = gsub(cluster[,3],pattern=".",replacement="-",fixed = TRUE)
cluster[,4] = gsub(cluster[,4],pattern=".",replacement="-",fixed = TRUE)
#cluster[,5] = gsub(cluster[,5],pattern=".",replacement="-",fixed = TRUE)

#cluster1 = b[which(b$cluster= 'cluster1'),]

c1 = cluster[,1]
c2 = cluster[,2]
c3 = cluster[,3]
c4 = cluster[,4]
#c5 = cluster[,5]

class(c1)

c1[c1==""] <- NA
c2[c2==""] <- NA
c3[c3==""] <- NA
c4[c4==""] <- NA
#c5[c5==""] <- NA

c1 <- na.omit(c1)
c2 <- na.omit(c2)
c3 <- na.omit(c3)
c4 <- na.omit(c4)
#c5 <- na.omit(c5)



##临床数据
clinic<-read.csv("/Users/moonly/Desktop/Pan-cancer/pan_clin.csv", header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)

#-----------------------OS----------------------------------------
clinic[which(as.numeric(clinic$OS.time)>3650),"OS"] = 0
clinic[which(as.numeric(clinic$OS.time)>3650),"OS.time"] = 3650
clinic[which(as.numeric(clinic$DSS.time)>3650),"DSS"] = 0
clinic[which(as.numeric(clinic$DSS.time)>3650),"DSS.time"] = 3650
clinic[which(as.numeric(clinic$DFI.time)>3650),"DFI"] = 0
clinic[which(as.numeric(clinic$DFI.time)>3650),"DFI.time"] = 3650
clinic[which(as.numeric(clinic$PFI.time)>3650),"PFI"] = 0
clinic[which(as.numeric(clinic$PFI.time)>3650),"PFI.time"] = 3650


clinic$type = ""      # 添加列名为type 的一列
clinic[clinic$sample %in% c1, "type"] = "cluster1"
clinic[clinic$sample %in% c2, "type"] = "cluster2"
clinic[clinic$sample %in% c3, "type"] = "cluster3"
clinic[clinic$sample %in% c4, "type"] = "cluster4"
#clinic[clinic$sample %in% c5, "type"] = "cluster5"


#2.1顺便跑个分布
colnames(clinic)
head(clinic)
a <- table(clinic$`cancer type abbreviation`,clinic$type)
write.csv(a,"c4_distribution.csv")


clinic = clinic[clinic$type != "",]
clinic$type = factor(clinic$type, level = c("cluster1", "cluster2", "cluster3","cluster4")) #"cluster5"


clinic1 = subset(clinic, !is.na(OS.time) & !is.na(OS))
clinic1$OS.time = clinic1$OS.time/30
cluster_A1 = length(which(clinic1$type == "cluster1"))
cluster_A2 = length(which(clinic1$type == "cluster2"))
cluster_A3 = length(which(clinic1$type == "cluster3"))
cluster_A4 = length(which(clinic1$type == "cluster4"))
#cluster_A5 = length(which(clinic1$type == "cluster5"))

clinic2 = subset(clinic, !is.na(DSS.time) & !is.na(DSS))
clinic2$DSS.time = clinic2$DSS.time/30
cluster_B1 = length(which(clinic2$type == "cluster1"))
cluster_B2 = length(which(clinic2$type == "cluster2"))
cluster_B3 = length(which(clinic2$type == "cluster3"))
cluster_B4 = length(which(clinic2$type == "cluster4"))
#cluster_B5 = length(which(clinic2$type == "cluster5"))

clinic3 = subset(clinic, !is.na(DFI.time) & !is.na(DFI))
clinic3$DFI.time = clinic3$DFI.time/30
cluster_C1 = length(which(clinic3$type == "cluster1"))
cluster_C2 = length(which(clinic3$type == "cluster2"))
cluster_C3 = length(which(clinic3$type == "cluster3"))
cluster_C4 = length(which(clinic3$type == "cluster4"))
#cluster_C5 = length(which(clinic3$type == "cluster5"))

clinic4 = subset(clinic, !is.na(PFI.time) & !is.na(PFI))
clinic4$PFI.time = clinic4$PFI.time/30
cluster_D1 = length(which(clinic4$type == "cluster1"))
cluster_D2 = length(which(clinic4$type == "cluster2"))
cluster_D3 = length(which(clinic4$type == "cluster3"))
cluster_D4 = length(which(clinic4$type == "cluster4"))
#cluster_D5 = length(which(clinic4$type == "cluster5"))

library(survival)

sA = survfit(Surv(OS.time, OS) ~ type, data = clinic1)
s1A = survdiff(Surv(OS.time, OS) ~ type, data = clinic1)
pvalueA = format(1 - pchisq(s1A$chisq, length(s1A$n) - 1),scientific=TRUE,digit=3)[[1]]

sB = survfit(Surv(DSS.time, DSS) ~ type, data = clinic2)
s1B = survdiff(Surv(DSS.time, DSS) ~ type, data = clinic2)
pvalueB = format(1 - pchisq(s1B$chisq, length(s1B$n) - 1),scientific=TRUE,digit=3)[[1]]

sC = survfit(Surv(DFI.time, DFI) ~ type, data = clinic3)
s1C = survdiff(Surv(DFI.time, DFI) ~ type, data = clinic3)
pvalueC = format(1 - pchisq(s1C$chisq, length(s1C$n) - 1),scientific=TRUE,digit=3)[[1]]


sD = survfit(Surv(PFI.time, PFI) ~ type, data = clinic4)
s1D = survdiff(Surv(PFI.time, PFI) ~ type, data = clinic4)
pvalueD = format(1 - pchisq(s1D$chisq, length(s1D$n) - 1),scientific=TRUE,digit=3)[[1]]

library(survminer)
library(ggplot2)
#install.packages("ggthemes")
library(ggthemes)
library(survival)
library(survminer)
library(patchwork)
library(ggthemes)
p1 = ggsurvplot(sA,title = paste0("OS"),legend = c(0.2,0.15),legend.title = "", legend.labs = c(paste0("c1 (n=", cluster_A1, ")"),paste0("c2 (n=", cluster_A2, ")"),paste0("c3 (n=", cluster_A3, ")"),paste0("c4 (n=", cluster_A4, ")")), xlab = "OS (months)",ylab = "",font.tickslab = c(13, "plain", "black"),font.xlab = c(12, "plain", "black"),font.ylab = c(12, "plain", "black"),font.legend = c(7, "plain", "black"),font.main = c(15, "bold", "black"), pval = paste("P = ",pvalueA), pval.size = 5, pval.coord=c(30,1.0),palette = c("#314a7b","#a0403b","#e6b422","#006e54"), ggtheme = theme_base())

p2 = ggsurvplot(sB, title = paste0("DSS"),legend = c(0.2,0.15),legend.title = "", legend.labs = c(paste0("c1 (n=", cluster_B1, ")"),paste0("c2 (n=", cluster_B2, ")"),paste0("c3 (n=", cluster_B3, ")"),paste0("c4 (n=", cluster_B4, ")")), xlab = "DSS (months)",ylab = "",font.tickslab = c(13, "plain", "black"),font.xlab = c(12, "plain", "black"),font.ylab = c(12, "plain", "black"),font.legend = c(7, "plain", "black"),font.main = c(15, "bold", "black"), pval = paste("P = ",pvalueB), pval.size = 5, pval.coord=c(30,1.0),palette = c("#314a7b","#a0403b","#e6b422","#006e54"), ggtheme = theme_base())

p3 = ggsurvplot(sC, title = paste0("DFI"),legend = c(0.2,0.15),legend.title = "", legend.labs = c(paste0("c1 (n=", cluster_C1, ")"),paste0("c2 (n=", cluster_C2, ")"),paste0("c3 (n=", cluster_C3, ")"),paste0("c4 (n=", cluster_C4, ")")), xlab = "DFI (months)",ylab = "",font.tickslab = c(13, "plain", "black"),font.xlab = c(12, "plain", "black"),font.ylab = c(12, "plain", "black"),font.legend = c(7, "plain", "black"),font.main = c(15, "bold", "black"), pval = paste("P = ",pvalueC), pval.size = 5, pval.coord=c(30,1.0),palette = c("#314a7b","#a0403b","#e6b422","#006e54"),ggtheme = theme_base())

p4 = ggsurvplot(sD, title = paste0("PFI"),legend = c(0.2,0.15),legend.title = "", legend.labs = c(paste0("c1 (n=", cluster_D1, ")"),paste0("c2 (n=", cluster_D2, ")"),paste0("c3 (n=", cluster_D3, ")"),paste0("c4 (n=", cluster_D4, ")")), xlab = "PFI (months)",ylab = "",font.tickslab = c(13, "plain", "black"),font.xlab = c(12, "plain", "black"),font.ylab = c(12, "plain", "black"),font.legend = c(7, "plain", "black"),font.main = c(15, "bold", "black"), pval = paste("P = ",pvalueD), pval.size = 5, pval.coord=c(30,1.0),palette = c("#314a7b","#a0403b","#e6b422","#006e54"), ggtheme = theme_base())

p <- p1$plot+p2$plot+p3$plot+p4$plot+plot_layout(nrow = 1)

pdf(file = paste0("Survival(4类-10years).pdf"), width = 18, height = 4.5)
print(p)
dev.off()



######################################1.3 therapy(stage/grade/response)###################################################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/lab/Pan-cancer/pan-results/VS/clinical/therapy/")
drug = read.csv("pancan_clinical_use.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
drug = drug[,-5]
drug = drug[,-5]
dim(drug)
# [1] 12591    7

sum(duplicated(drug[,1]))
# [1] 0

score = read.csv("/Users/moonly/Desktop/lab/Pan-cancer/pan-results/打分聚类生存/cluster_2list.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
colnames(score)[1] <- "sample"
cordata = merge(score, drug, by.x = "sample", by.y = "sample")
dim(cordata)
# [1] 10295     8
View(head(cordata))

result <- xtabs(~cluster + treatment_outcome, data = cordata)
rowsums <- margin.table(result, 1) #计算每一行的和
# perc <- prop.table(result,1)*100 #计算每一行的比例/百分比
res <- cbind(result,rowsums)
write.csv(res, "treatment_num.csv")

result1 <- xtabs(~subtype + tumor_stage, data = cordata)
# rowsums <- margin.table(result, 1) #计算每一行的和
# perc <- prop.table(result,1)*100 #计算每一行的比例/百分比
# res <- cbind(result,rowsums)
write.csv(result1, "tumor_stage_num.csv")

result2 <- xtabs(~subtype + clinical_stage, data = cordata)
write.csv(result2, "clinical_stage_num.csv")

result3 <- xtabs(~subtype + histological_grade, data = cordata)
write.csv(result3, "grade_num.csv")

####chisq test
#treatment
a = matrix(c(1356, 317, 764, 186, 1542, 466, 440, 26), nrow = 2)
View(a)
sum(a)
chisq.test(a)
chisq.test(a)$p.value
# [1]  3.016107e-16

#grade
a1 = matrix(c(368, 665, 640, 685, 725, 719, 198, 197), nrow = 2)
View(a1)
sum(a1)
chisq.test(a1)
chisq.test(a1)$p.value
# 4.389605e-13

#fisher
# fisher.test(a1)
#tumor_stage
a2 = matrix(c(1461, 960, 776, 648, 1827, 823, 148, 110), nrow = 2)
View(a2)
sum(a2)
chisq.test(a2)
chisq.test(a2)$p.value
# [1] 4.075645e-20

#clinical_stage
a3 = matrix(c(451, 448, 210, 510, 233, 83, 260, 96), nrow = 2)
View(a3)
sum(a3)
chisq.test(a3)
chisq.test(a3)$p.value
# [1] 1.355415e-58



#########################################3 Estimate(immune scores\stromal scores\tumor purity)####################################
# R package: estimate
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE) #lib.loc ="C:/Users/heyin/Documents/R/win-library/4.0", fields = c("estimate", "Version")
#若安装出现问题
# options(download.file.method = 'libcurl') 
# options(url.method='libcurl') 
library(estimate)
# TCGA -----------------------------------------------------------------------------------------------------------------------
library(estimate)

# dir.create("D:/pan_analysis_files/estimate/")
# files = list.files("C:/Users/Administrator/Desktop/skcm/rna_seq_csv/")
# cancer_type <-c ("GSE22155")
# for (i in 1:1){
# dir.create(paste0("C:/Users/Administrator/Desktop/skcm/estimate/", cancer_type[[i]], "/"))
# setwd(paste0("C:/Users/Administrator/Desktop/skcm/estimate/", cancer_type[[i]], "/"))

setwd("/Users/moonly/Desktop/lab/Pan-cancer/pan-results/Estimate")
a = read.csv("./exp/GSE105261_expr.csv", header=TRUE, stringsAsFactors=FALSE,check.names=FALSE)
dim(a)
[1] 16599 10324

write.table(a, "exp_pancancer.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

filterCommonGenes(input.f="exp_pancancer.txt", output.f="estimate_pancancer.gct", id="GeneSymbol") #注意id是symbol还是EntrezID
# filterCommonGenes(input.f=paste0("estimate_", cancer_type[[i]], ".txt"), output.f=paste0("estimate_", cancer_type[[i]], ".gct"), id="GeneSymbol")
estimateScore(input.ds="estimate_pancancer.gct", output.ds= "estimate_pancancer_score.gct", platform="affymetrix")

scores = read.table("estimate_pancancer_score.gct", skip=2, header=TRUE)
rownames(scores) <- scores[,1]
scores=scores[,-c(1,2)]
colnames(scores) = gsub(pattern=".", replacement="-", colnames(scores), fixed=TRUE)
write.csv(scores, "estimatescore_pancancer.csv")
scores_t = t(scores)
dim(scores_t)
[1] 10323     4

write.csv(scores_t, "T_estimatescore_pancancer.csv")

# }


#################Estimate_wilcox####################################
rm(list=ls())
setwd("/Users/moonly/Desktop/lab/Pan-cancer/pan-results/Estimate")
ki<-read.csv("T_estimatescore_pancancer.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(ki)[1]=c("sample")
group<-read.csv("cluster.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(group)[1]=c("sample")


c1 = group[,1]
c2 = group[,2]
c3 = group[,3]
c4 = group[,4]

c1[c1==""] <- NA
c2[c2==""] <- NA
c3[c3==""] <- NA
c4[c4==""] <- NA

c1 <- na.omit(c1)
c2 <- na.omit(c2)
c3 <- na.omit(c3)
c4 <- na.omit(c4)
c1 <- as.data.frame(c1)
c1[,2] = paste("cluster",1, sep = " ")
colnames(c1)[1] <- c("sample")
colnames(c1)[2] <- c("cluster")
head(c1)

c2 <- as.data.frame(c2)
c2[,2] = paste("cluster",2, sep = " ")
colnames(c2)[1] <- c("sample")
colnames(c2)[2] <- c("cluster")
head(c2)

c3 <- as.data.frame(c3)
c3[,2] = paste("cluster",3, sep = " ")
colnames(c3)[1] <- c("sample")
colnames(c3)[2] <- c("cluster")
head(c3)

c4 <- as.data.frame(c4)
c4[,2] = paste("cluster",4, sep = " ")
colnames(c4)[1] <- c("sample")
colnames(c4)[2] <- c("cluster")
head(c4)

c5 <- rbind(c1,c2)
c6 <- rbind(c5,c3)
c7 <- rbind(c6,c4)
dim(c7)
########合并
cordata<-merge(group,ki,by.x="sample",by.y="sample")
View(head(cordata))
dim(cordata)
[1] 10323     6

write.csv(cordata,"estimate_combine_new.csv",row.names=FALSE)

##########KW test######################
a = read.csv("estimate_combine_new.csv",check.names = FALSE, stringsAsFactors=FALSE,header=TRUE)
dim(a)
[1] 10323    6

pv=c()
meanresults = matrix()
for (i in 3:dim(a)[2]){
  c1=subset(a[,i], a[,2]==1) 
  c2=subset(a[,i], a[,2]==2) 
  #c3=subset(a[,i], a[,2]=="St3") 
  #c4=subset(a[,i], a[,2]=="St4") 
  # c5=subset(a[,i], a[,2]=="St5")
  pv1<- kruskal.test(list(c1, c2))$p.value  #KW
  pv2<-cbind(colnames(a)[i],pv1)
  pv<-rbind(pv,pv2)
  #mean value
  meanresult <- aggregate(a[,i],by=list(a$cluster),FUN=mean,na.rm=T)
  colnames(meanresult) = c("",colnames(a)[i])
  meanresults=cbind(meanresults,meanresult[2])
}
result=t(meanresults)
#meanresults1=as.matrix(meanresults)
#colnames(meanresults1)=NULL
res=cbind(pv, result[-1,])
colnames(res) = c("esti_score","pv_kw","cluster1","cluster2","cluster3","cluster4")
reS = res[sapply(c("ESTIMATEScore","TumorPurity","ImmuneScore", "StromalScore"), function(x) {which(rownames(res) == x)}), ] 
View(reS)

write.csv(reS,"Estimate_kw_new.csv",row.names = FALSE)


#########wilcox###################################
cordata <- read.csv("estimate_combine_new.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

St = c("cluster1","cluster2","cluster3","cluster4")
resj = c()
for(j in 1){
  resk = c()
  for(k in (j+1):2){
    pv=c() 
    result=c()
    for(f in 3:dim(cordata)[2]){ 
      high_pos = which(cordata[,2]==1)##高的组所对应样本位置
      low_pos = which(cordata[,2]==2)#低的组所对应样本位置
      
      pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
      pv2 = cbind(colnames(cordata)[f],pvG,pvL)
      pv<-rbind(pv,pv2)   
      rownames(pv)=NULL  
    }
    ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
    pvalue1<-pv[ord,]
    fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
    pv2<-cbind(pvalue1,fdr)
    colnames(pv2)<-c("esti_score","pvG","pvL","fdr")
    result=rbind(result,pv2)
    resu1 = result[sapply(c("ESTIMATEScore","TumorPurity","ImmuneScore", "StromalScore"), function(x) {which(result[,1] == x)}), ] #reorder as the names
    rownames(resu1) = resu1[,1]
    resu1 = resu1[,-1]
    colnames(resu1)[1] = paste0("pvG_",St[j]) 
    colnames(resu1)[2] = paste0("pvL_",St[k])   
    resk = cbind(resk, resu1)   
  }
  resj = cbind(resj, resk) 
}
View(resj)
dim(resj)
[1]  4 18

write.csv(resj,"Estimate_wilcox_new.csv")









############################################4 genome-instability############################################
#wilcox
setwd("/data/bioinfo2021/XML/pan-cancer/genome-instability")
rm(list=ls())
ki<-read.csv("pancan_CNA.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
head(ki)
dim(ki)
#[1] 8782    3 #TMB
#[1] 8948    2 #HRD
#[1] 10323    3 #DEPTH
#[1] 5560    4 #CNA
head(ki)
ki[,1] = substr(ki[,1],1,15)
colnames(ki)[1]=c("sample")
head(ki)
write.csv(ki,"pancan_CNA.csv",row.names=FALSE)

group<-read.csv("/data/bioinfo2021/XML/pan-cancer/cluster_2list.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
head(ki)
head(group)

#合并
cordata<-merge(group,ki,by.x="sample",by.y="sample")
dim(cordata)
head(cordata)
#[1] 8319    3 #TMB
#[1] 8942    3 #HRD
#[1] 10323    3 #DEPTH
#[1] 5560    4 #CNA
head(cordata)

St <- c("cluster 1", "cluster 2", "cluster 3", "cluster 4")
resj = c()

for(j in 1:3){
  pv = c()
  for(k in (j+1):4){
    high_pos <- which(cordata[,2] == St[j])##高的组所对应样本位置
    low_pos <- which(cordata[,2] == St[k])#低的组所对应样本位置
    
    pvG <- wilcox.test(as.numeric(cordata[,3][high_pos]), as.numeric(cordata[,3][low_pos]), alternative = "greater")$p.value
    pvL <- wilcox.test(as.numeric(cordata[,3][high_pos]), as.numeric(cordata[,3][low_pos]), alternative = "less")$p.value
    pv2 <- cbind(colnames(cordata)[3], pvG, pvL) 
    colnames(pv2)[2] <- paste0("pvG_", St[j]) 
    colnames(pv2)[3] <- paste0("pvL_", St[k])    
    pv = cbind(pv, pv2) 
    
  }
  resj = cbind(resj, pv) 
}

dim(resj)
#[1]  1 18
head(resj)
write.csv(resj,"HRD_wilcox_new.csv")

##########KW test######################
rm(list = ls())
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/genome-instability")
GI = c("TMB","HRD","depth","depth2","CNA")
for(j in 1:4){
  ki<-read.csv(paste0("pancan_",GI[j],".csv"),stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
  #colnames(ki)=c("sample",GI[j])
  group<-read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/cluster_2list.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
  a <-merge(group,ki,by.x="sample",by.y="sample")
  i=3
  c1=subset(a[,i], a[,2]=="cluster 1") 
  c2=subset(a[,i], a[,2]=="cluster 2") 
  c3=subset(a[,i], a[,2]=="cluster 3") 
  c4=subset(a[,i], a[,2]=="cluster 4")
  # c5=subset(a[,i], a[,2]=="St5")
  pv1<- kruskal.test(list(c1, c2, c3, c4))$p.value  #KW
  pv2<-cbind(colnames(a)[i],pv1)
  colnames(pv2) = c(GI[j],"value")
  
  #mean value
  meanresult <- aggregate(a[,i],by=list(a$cluster),FUN=mean,na.rm=T)
  colnames(meanresult) = c(GI[j],"value")
  result = rbind(pv2, meanresult)
  
  write.csv(result,paste0("result/pancan_", GI[j], "_kw.csv"),row.names = FALSE)
  write.csv(a,paste0("pancan_", GI[j], "_combine.csv"),row.names = FALSE)
}

#########wilcox###################################
GI = c("TMB","HRD","depth","depth2","CNA")
for(i in 1:4){
  cordata <- read.csv(paste0("pancan_",GI[i],"_combine.csv"),stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
  St = c("cluster 1","cluster 2","cluster 3","cluster 4")
  resj = c()
  for(j in 1:3){
    resk = c()
    for(k in (j+1):4){
      high_pos = which(cordata[,2]==St[j])##高的组所对应样本位置
      low_pos = which(cordata[,2]==St[k])#低的组所对应样本位置
      f = 3
      pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
      resu1 = cbind(colnames(cordata)[f],pvG,pvL)
      rownames(resu1) = resu1[,1]
      resu1 = resu1[,-1]
      resu = t(resu1)
      colnames(resu)[1] = paste0("pvG_",St[j]) 
      colnames(resu)[2] = paste0("pvL_",St[k])   
      resk = cbind(resk, resu)   
    }
    resj = cbind(resj, resk) 
  }
  rownames(resj) =  GI[i]
  write.csv(resj,paste0("result2/pancan_", GI[i], "_wilcox.csv"))
}

#################IFN#########

geneset = read.csv("IFN.csv", stringsAsFactors=FALSE,check.names=FALSE,header = TRUE)
rna_seq = read.csv("pancan_exp_onlytumor(ID).csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, row.names=1)
dim(rna_seq)
[1] 16599 10323
rna_seq = as.matrix(rna_seq)

#gsea score
library(GSVA)
memory.limit(10000000)
result = gsva(rna_seq, list(geneset[,1]), method="ssgsea", ssgsea.norm=TRUE, verbose=TRUE)
n <- t(result)
write.csv(n, "ssGSEA(IFN).csv")

###合并
ki<-read.csv("ssGSEA(IFN).csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(ki)=c("sample","IFN")
group<-read.csv("pancan_6_cluster4_st.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
#colnames(group)[1]=c("sample")

cordata<-merge(group,ki,by.x="sample",by.y="sample")
View(head(cordata))
dim(cordata)

write.csv(cordata,"IFN_ssGSEA_combine_new.csv",row.names=FALSE)

##########KW test######################
a = read.csv("IFN_ssGSEA_combine_new.csv",check.names = FALSE, stringsAsFactors=FALSE,header=TRUE)
dim(a)

# pv=c()
# meanresults = matrix()
# for (i in 3:dim(a)[2]){
i=3
c1=subset(a[,i], a[,2]=="cluster1") 
c2=subset(a[,i], a[,2]=="cluster2") 
c3=subset(a[,i], a[,2]=="cluster3") 
c4=subset(a[,i], a[,2]=="cluster4") 
# c5=subset(a[,i], a[,2]=="St5")
pv1<- kruskal.test(list(c1, c2, c3, c4))$p.value  #KW
pv2<-cbind(colnames(a)[i],pv1)
colnames(pv2) = c("IFN","value")

#mean value
meanresult <- aggregate(a[,i],by=list(a$subtype),FUN=mean,na.rm=T)
colnames(meanresult) = c("IFN","value")
result = rbind(pv2, meanresult)
View(result)

write.csv(result,"IFN_kw_new.csv",row.names = FALSE)

#########wilcox###################################
cordata <- read.csv("IFN_ssGSEA_combine_new.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
St = c("cluster1","cluster2","cluster3","cluster4")
resj = c()
for(j in 1:3){
  resk = c()
  for(k in (j+1):4){
    # pv=c() 
    # result=c()
    # for(f in 3:dim(cordata)[2]){ 
    high_pos = which(cordata[,2]==St[j])##高的组所对应样本位置
    low_pos = which(cordata[,2]==St[k])#低的组所对应样本位置
    f = 3
    pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
    pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
    resu1 = cbind(colnames(cordata)[f],pvG,pvL)
    rownames(resu1) = resu1[,1]
    resu1 = resu1[,-1]
    resu = t(resu1)
    colnames(resu)[1] = paste0("pvG_",St[j]) 
    colnames(resu)[2] = paste0("pvL_",St[k])   
    resk = cbind(resk, resu)   
  }
  resj = cbind(resj, resk) 
}
View(resj)
dim(resj)
# [1]  1 12

rownames(resj) = "IFN"
write.csv(resj,"IFN_wilcox_new.csv")

##############4 E##########

library(maftools)

setwd("E:/GISTIC2/cluster4/results")

g <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt", gisticAmpGenesFile="amp_genes.conf_90.txt", gisticDelGenesFile="del_genes.conf_90.txt", gisticScoresFile="scores.gistic", isTCGA=TRUE)

gisticChromPlot(gistic=g, markBands="all")


pdf("ChromPlot.pdf", width = 8, height = 4, onefile = TRUE)
gisticChromPlot(
  gistic = g,
  fdrCutOff = 0.1,
  txtSize = 0.8,
  cytobandTxtSize = 0.5,
  color = c("#D95F02", "#1B9E77"),
  markBands = "all",
  ref.build = "hg38",
  y_lims = c(-0.5, 0.5)
)
dev.off()


#读入GISTIC文件

setwd("E:/GISTIC2/cluster4/results")

cancer.gistic = readGistic(gisticAllLesionsFile = "all_lesions.conf_90.txt",
                           gisticAmpGenesFile = "amp_genes.conf_90.txt", 
                           gisticDelGenesFile = "del_genes.conf_90.txt", 
                           gisticScoresFile = "scores.gistic", 
                           isTCGA = T)

## -Processing Gistic files..
## --Processing amp_genes.conf_99.txt
## --Processing del_genes.conf_99.txt
## --Processing scores.gistic
## --Summarizing by samples

cancer.gistic #查看信息

#进行统计
getSampleSummary(cancer.gistic)

getGeneSummary(cancer.gistic)

getCytobandSummary(cancer.gistic)

write.GisticSummary(gistic=cancer.gistic, basename="cancer_gistic2")

##绘图
#genome plot
a <- getCytobandSummary(cancer.gistic)$Cytoband #提取Cytoband(已经按照q-value排序)
pdf(file= "E:/GISTIC2/gistic_c4.pdf",width=7,height=4)
gisticChromPlot(gistic = cancer.gistic, ref.build = "hg38", fdrCutOff = 0.05, markBands = a[1:3], y_lims = c(-1,1)) #markBands：标注染色体区带 , markBands = a[1:12],cytobandTxtSize = 0.6(default), y_lims = NULL, y_lims = c(-1,1)
dev.off()


###########################Oncogenic##########################################
##############GSVA打分##################
setwd("/Users/moonly/Desktop/pan_cancer/pancan/DDR&oncogenic pathway/oncogenic/")
library(GSVA)
gs = read.csv("oncogenic pathway(id).csv", stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)# immune signature --- Gene ID
a <- read.csv("/Users/moonly/Desktop/pan_cancer/pancan/pancan_exp_onlytumor(ID).csv",check.names=F,row.names=1,header=T,stringsAsFactors = FALSE)
a = as.matrix(a)
dim(a)

memory.limit(10000000)#扩容
result=c()
for(i in 1:dim(gs)[2]){ 
  gs1 = as.list(gs[i])
  gs1 = lapply(gs1, function(x) x[!is.na(x)])   # Gene symbol #!is.na(x)
  res = gsva(a, gs1, method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
  result<-rbind(result,res)
}
n <- t(result)
write.csv(n ,"oncogenic_gsva.csv")

###################合并##############
setwd("/Users/moonly/Desktop/pan_cancer/pancan/DDR&oncogenic pathway/oncogenic/")
rm(list=ls())
ki<-read.csv("oncogenic_gsva.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(ki)[1]=c("sample")
group<-read.csv("/Users/moonly/Desktop/pan_cancer/pancan/pancan_6_cluster4_st.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

cordata<-merge(group,ki,by.x="sample",by.y="sample")
View(head(cordata))
dim(cordata)
[1] 10323    15

write.csv(cordata,"oncogenic_ssGSEA_combine_new.csv",row.names=FALSE)

##########KW test######################
rm(list=ls())
a = read.csv("oncogenic_ssGSEA_combine_new.csv",check.names = FALSE, stringsAsFactors=FALSE,header=TRUE)

pv=c()
meanresults = matrix()
for (i in 3:dim(a)[2]){
  c1=subset(a[,i], a[,2]=="cluster 1") 
  c2=subset(a[,i], a[,2]=="cluster 2") 
  c3=subset(a[,i], a[,2]=="cluster 3") 
  c4=subset(a[,i], a[,2]=="cluster 4") 
  # c5=subset(a[,i], a[,2]=="St5")
  pv1<- kruskal.test(list(c1, c2, c3, c4))$p.value  #KW
  pv2<-cbind(colnames(a)[i],pv1)
  pv<-rbind(pv,pv2)
  #mean value
  meanresult <- aggregate(a[,i],by=list(a$subtype),FUN=mean,na.rm=T)
  colnames(meanresult) = c("",colnames(a)[i])
  meanresults=cbind(meanresults,meanresult[2])
}

result=t(meanresults)
res=cbind(pv, result[-1,])
colnames(res) = c("pathway","pv_kw","cluster1","cluster2","cluster3","cluster4")
View(res)

write.csv(res,"oncogenic_kw_new.csv",row.names = FALSE)


#########wilcox#####
rm(list=ls())
cordata <- read.csv("oncogenic_ssGSEA_combine_new.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

St = c("cluster 1","cluster 2","cluster 3","cluster 4")
resj = c()
for(j in 1:3){
  resk = c()
  for(k in (j+1):4){
    pv=c() 
    result=c()
    for(f in 3:dim(cordata)[2]){ 
      high_pos = which(cordata[,2]==St[j])##高的组所对应样本位置
      low_pos = which(cordata[,2]==St[k])#低的组所对应样本位置
      
      pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
      pv2 = cbind(colnames(cordata)[f],pvG,pvL)
      pv<-rbind(pv,pv2)   
      rownames(pv)=NULL  
    }
    ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
    pvalue1<-pv[ord,]
    fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
    pv2<-cbind(pvalue1,fdr)
    colnames(pv2)<-c("pathway","pvG","pvL","fdr")
    result=rbind(result,pv2)
    resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
    rownames(resu1) = resu1[,1]
    resu1 = resu1[,-1]
    colnames(resu1)[1] = paste0("pvG_",St[j]) 
    colnames(resu1)[2] = paste0("pvL_",St[k])   
    resk = cbind(resk, resu1)   
  }
  resj = cbind(resj, resk) 
}
View(resj)
dim(resj)
#[1] 13 18
write.csv(resj,"oncogenic_wilcox_new.csv")

##################################5 Biological####################################

##############################################################bio_process######################################
##############GSVA打分######
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/bio_process")
library(GSVA)
gc()
gs = read.csv("genesets.csv", stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)# immune signature --- Gene ID
a <- read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/exp/exp(rowSums(is.na)>0).csv",check.names=F,row.names=1,header=T,stringsAsFactors = FALSE)
a = as.matrix(a)
dim(a)
#[1] 16335 10323


result=c()
for(i in 1:dim(gs)[2]){ 
  gs1 = as.list(gs[i])
  gs1 = lapply(gs1, function(x) x[!is.na(x)])   # Gene symbol #!is.na(x)
  res = gsva(a, gs1, method="ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
  result<-rbind(result,res)
}
n <- t(result)
dim(n)
#[1] 10323    36
n <- as.data.frame(n)
colnames(n)
colnames(n)[20] <- "TGF-β signaling"
write.csv(n,"ssgsea.csv")


#####合并
#rm(list=ls())
#setwd("/data/bioinfo2021/XML/pan-cancer/bio_process/")

ki<-read.csv("ssgsea.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

colnames(ki)[1]=c("sample")
ki[,1] = gsub(pattern=".", replacement="-", ki[,1] , fixed=TRUE)


group<-read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/打分聚类生存/cluster_2list.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
#colnames(group)[1]=c("sample")

cordata<-merge(group,ki,by.x="sample",by.y="sample")
cordata[1:4,1:4]
dim(cordata)
#[1] 10323    38

write.csv(cordata,"ssGSEA_combine_new.csv",row.names=FALSE)

##########KW test######################
a = read.csv("ssGSEA_combine_new.csv",check.names = FALSE, stringsAsFactors=FALSE,header=TRUE)
dim(a)
#[1] 10323    28

pv=c()
meanresults = matrix()
for (i in 3:dim(a)[2]){
  c1=subset(a[,i], a[,2]=="cluster 1") 
  c2=subset(a[,i], a[,2]=="cluster 2") 
  c3=subset(a[,i], a[,2]=="cluster 3") 
  c4=subset(a[,i], a[,2]=="cluster 4")
  
  pv1<- kruskal.test(list(c1, c2, c3, c4))$p.value  #KW
  pv2<-cbind(colnames(a)[i],pv1)
  pv<-rbind(pv,pv2)
  #mean value
  meanresult <- aggregate(a[,i],by=list(a$cluster),FUN=mean,na.rm=T)
  colnames(meanresult) = c("",colnames(a)[i])
  meanresults=cbind(meanresults,meanresult[2])
}

result=t(meanresults)
#meanresults1=as.matrix(meanresults)
#colnames(meanresults1)=NULL
res=cbind(pv, result[-1,])
colnames(res) = c("geneset","pv_kw","cluster 1","cluster 2","cluster 3","cluster 4")
#reS = res[sapply(c("Antigen processing and presentation","Apoptosis","JAK-STAT signaling", "Invasion", "Metastasis","Differentiation","Inflammation","DNA damage", "DNA repair", "Mismatch repair","Homologous recombination","Stemness","stenmess_paper", "Proliferation", "Proliferation_paper","Cell cycle","EMT","EMT_paper", "p53 signaling", "TGF-β signaling","Wnt signaling","Notch signaling","PI3K-Akt signaling", "Hedgehog signaling", "Angiogenesis","Quiescence"), function(x) {which(rownames(res) == x)}), ] 
res[1:3,1:3]
dim(res)
#write.csv(reS,"geneset_kw_new.csv",row.names = FALSE)
write.csv(res,"geneset_kw_new.csv",row.names = FALSE)

#########wilcox###################################
cordata <- read.csv("ssGSEA_combine_new.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
St = c("cluster 1","cluster 2","cluster 3","cluster 4")
resj = c()
for(j in 1:3){
  resk = c()
  for(k in (j+1):4){
    pv=c() 
    result=c()
    for(f in 3:dim(cordata)[2]){ 
      high_pos = which(cordata[,2]==St[j])##高的组所对应样本位置
      low_pos = which(cordata[,2]==St[k])#低的组所对应样本位置
      
      pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
      pv2 = cbind(colnames(cordata)[f],pvG,pvL)
      pv<-rbind(pv,pv2)   
      rownames(pv)=NULL  
    }
    ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
    pvalue1<-pv[ord,]
    fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
    pv2<-cbind(pvalue1,fdr)
    colnames(pv2)<-c("pathway","pvG","pvL","fdr")
    result=rbind(result,pv2)
    #resu1 = result[sapply(c("Antigen processing and presentation","Apoptosis","JAK-STAT signaling", "Invasion", "Metastasis","Differentiation","Inflammation","DNA damage", "DNA repair", "Mismatch repair","Homologous recombination","Stemness","stenmess_paper", "Proliferation", "Proliferation_paper","Cell cycle","EMT","EMT_paper", "p53 signaling", "TGF-β signaling","Wnt signaling","Notch signaling","PI3K-Akt signaling", "Hedgehog signaling", "Angiogenesis","Quiescence"), function(x) {which(result[,1] == x)}), ] #reorder as the names
    resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
    rownames(resu1) = resu1[,1]
    resu1 = resu1[,-1]
    colnames(resu1)[1] = paste0("pvG_",St[j]) 
    colnames(resu1)[2] = paste0("pvL_",St[k])   
    resk = cbind(resk, resu1)   
  }
  resj = cbind(resj, resk) 
}
resj[1:3,1:3]
dim(resj)
#[1] 26 18

write.csv(resj,"pathway_wilcox_new.csv")
#########################################6 DEG gene####################################
#marker gene
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/DEG")


rm(list = ls())

c1 <- c(1,1,1,2,2,3)
c2 <- c(2,3,4,3,4,4)

id_df <- data.frame(id1 = c1,id2 = c2)

cluster <- read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/打分聚类生存/cluster.csv",header=TRUE, stringsAsFactors=FALSE,check.names=FALSE)

expr <- read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/exp/exp(rowSums(is.na)>0).csv",header=TRUE, row.names = 1, stringsAsFactors=FALSE,check.names=FALSE)
expr[1:3,1:3]
colnames(expr) = gsub(colnames(expr),pattern=".",replacement="-",fixed = TRUE)
cluster[,1] = gsub(cluster[,1],pattern=".",replacement="-",fixed = TRUE)
cluster[,2] = gsub(cluster[,2],pattern=".",replacement="-",fixed = TRUE)
cluster[,3] = gsub(cluster[,3],pattern=".",replacement="-",fixed = TRUE)
cluster[,4] = gsub(cluster[,4],pattern=".",replacement="-",fixed = TRUE)
#cluster[,5] = gsub(cluster[,5],pattern=".",replacement="-",fixed = TRUE)


expr[1:3,1:3]
cluster[1:3,1:3]


for(j in 1:length(c1)) {
  tag_a <- id_df[j,1]
  tag_b <- id_df[j,2]
  sample_a <- as.character(cluster[,tag_a])
  sample_b <- as.character(cluster[,tag_b])
  
  sample_a[sample_a == ""] <- NA
  sample_a <- na.omit(sample_a)
  
  sample_b[sample_b == ""] <- NA
  sample_b <- na.omit(sample_b)
  
  expr_a <- expr[,sample_a]
  expr_b <- expr[,sample_b]
  
  resC = c()
  for (i in 1:dim(expr_a)[1]) {
    x = as.numeric(expr_a[i,]); y = as.numeric(expr_b[i,])
    p.value = t.test(x, y, var.equal = TRUE, na.rm=TRUE)$p.value
    log2FC = mean(x, na.rm=TRUE)-mean(y, na.rm=TRUE)
    res = cbind(rownames(expr_a)[i], p.value, log2FC)
    resC = rbind(resC, res)
  }
  
  resC = resC[order(as.numeric(resC[,2])),]
  fdr = rep(1,dim(resC)[1]); for(i in 1:dim(resC)[1]) fdr[i] = as.numeric(resC[i,2])*dim(resC)[1]/i  ###多重检验校正
  resC = as.data.frame(cbind(resC, fdr))
  colnames(resC)[1] <- 'gene_symbol'
  
  write.csv(resC, paste0("DEGs_c",tag_a,"_c",tag_b,".csv"), row.names=FALSE)
  
  resC$log2FC <- as.numeric(resC$log2FC)
  resC$fdr <- as.numeric(resC$fdr)
  
  up <- subset(resC, log2FC > 0.585 & fdr < 0.05)
  down <- subset(resC, log2FC < -0.585 & fdr < 0.05)
  
  write.csv(up, paste0("up_DEGs_c",tag_a,"_c",tag_b,".csv"), row.names=FALSE)
  write.csv(down, paste0("up_DEGs_c",tag_b,"_c",tag_a,".csv"), row.names=FALSE) # 上下调互换  
}


#instersect
rm(list = ls())
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/DEG")

d12 <- read.csv("up_DEGs_c1_c2.csv", stringsAsFactors=FALSE,check.names=FALSE)
d13 <- read.csv("up_DEGs_c1_c3.csv", stringsAsFactors=FALSE,check.names=FALSE)
d14 <- read.csv("up_DEGs_c1_c4.csv", stringsAsFactors=FALSE,check.names=FALSE)

d12[1:3,1:3]
d13[1:3,1:3]
d14[1:3,1:3]

#加载dplyr包
library(dplyr)
#直接利用dplyr包里面的intersect函数对数据框取交集
result1 = Reduce(intersect,list(d12$gene_symbol,
                                d13$gene_symbol,
                                d14$gene_symbol,))

rownames(d12) <- d12[,1]
result1 <- d12[result1,]
result1[1:3,1:3]
#保存交集结果
write.csv(result1,"intersect1.csv")


###############2
d21 <- read.csv("up_DEGs_c2_c1.csv", stringsAsFactors=FALSE,check.names=FALSE)
d23 <- read.csv("up_DEGs_c2_c3.csv", stringsAsFactors=FALSE,check.names=FALSE)
d24 <- read.csv("up_DEGs_c2_c4.csv", stringsAsFactors=FALSE,check.names=FALSE)
#d25 <- read.csv("up_DEGs_c2_c5.csv", stringsAsFactors=FALSE,check.names=FALSE)

d21[1:3,1:3]
result2 = Reduce(intersect,list(d21$gene_symbol,
                                d23$gene_symbol,
                                d24$gene_symbol))
rownames(d21) <- d21[,1]
result2 <- d21[result2,]
result2[1:3,1:3]
#保存交集结果
write.csv(result2,"intersect2.csv")

######################3
d31 <- read.csv("up_DEGs_c3_c1.csv", stringsAsFactors=FALSE,check.names=FALSE)
d32 <- read.csv("up_DEGs_c3_c2.csv", stringsAsFactors=FALSE,check.names=FALSE)
d34 <- read.csv("up_DEGs_c3_c4.csv", stringsAsFactors=FALSE,check.names=FALSE)
#d35 <- read.csv("up_DEGs_c3_c5.csv", stringsAsFactors=FALSE,check.names=FALSE)

d31[1:3,1:3]
result3 = Reduce(intersect,list(d31$gene_symbol,
                                d32$gene_symbol,
                                d34$gene_symbol))
rownames(d31) <- d31[,1]
result3 <- d31[result3,]
result3[1:3,1:3]
#保存交集结果
write.csv(result3,"intersect3.csv")

#####################4
d41 <- read.csv("up_DEGs_c4_c1.csv", stringsAsFactors=FALSE,check.names=FALSE)
d42 <- read.csv("up_DEGs_c4_c2.csv", stringsAsFactors=FALSE,check.names=FALSE)
d43 <- read.csv("up_DEGs_c4_c3.csv", stringsAsFactors=FALSE,check.names=FALSE)
#d45 <- read.csv("up_DEGs_c4_c5.csv", stringsAsFactors=FALSE,check.names=FALSE)

result4 = Reduce(intersect,list(d41$gene_symbol,
                                d42$gene_symbol,
                                d43$gene_symbol))
rownames(d41) <- d41[,1]
result4 <- d41[result4,]
result4[1:3,1:3]
#保存交集结果
write.csv(result4,"intersect4.csv")



############################################7 Protein DEG#################################################

#######anova####
setwd("/Users/moonly/Desktop/pancan/protein/")
rppa = read.table("pancan_RPPA.csv", stringsAsFactors=FALSE, header=FALSE, check.names=FALSE, sep=",")
View(head(rppa))
dim(rppa)
# [1]  259 7755

rppa = na.omit(rppa)
dim(rppa)
# [1]  211 7755
#转置表达矩阵
rppa1 = t(rppa)
rppa1[1,1] = "sample" 
View(head(rppa1))
dim(rppa1)
# [1] 7755  211
rownames(rppa1)=NULL
colnames(rppa1)=rppa1[1,]
rppa2 = rppa1[-1,]
View(head(rppa2))
dim(rppa2)
# [1] 7754  211

group <- read.csv("/Users/moonly/Desktop/pancan/cluster.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

cordata = merge(group, rppa2, by.x="sample", by.y="sample")
View(head(cordata))
dim(cordata)
# [1] 7286  212

pv=c()  
meanresults = matrix()
for(f in 3:dim(cordata)[2]){
  aov<-aov(as.numeric(cordata[,f])~cordata[,2])
  r<-summary(aov)
  b<-r[[1]]["Pr(>F)"][[1]][1]
  b<-as.numeric(b)
  pv2<-cbind(colnames(cordata)[f],b)
  pv<-rbind(pv,pv2)
  #mean value
  meanresult <- aggregate(as.numeric(cordata[,f]),by=list(cordata$subtype),FUN=mean,na.rm=T) #cordata[,f]要保证是numeric,否则会返回NA
  colnames(meanresult) = c("",colnames(cordata)[f])
  meanresults=cbind(meanresults,meanresult[2])
}
result=t(meanresults)
res=cbind(pv, result[-1,])
colnames(res) = c("protein","pv_anova","St1","St2","St3","St4")
res = res[order(as.numeric(res[,2])),]
View(res)
write.csv(res,"protein_anova_test_new.csv",row.names=FALSE)



################################Student's t-Test 1v1##############################
setwd("/Users/moonly/Desktop/pan_cancer/pancan/protein/")
rppa = read.table("pancan_RPPA.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, row.names=1, sep=",")
View(head(rppa))
dim(rppa)
# [1]  258 7754

#rppa = rppa[-which(apply(rppa,1,function(x) all(is.na(x)))), ] #删除全为NA的行#没有全为NA的行
rppa = na.omit(rppa) #不去除NA无法进行student.t.test
dim(rppa)
# [1]  210 7754

group <- read.csv("/Users/moonly/Desktop/pancan/cluster.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

St = c("cluster 1","cluster 2","cluster 3","cluster 4")
for(j in 1:3){
  for(k in (j+1):4){
    #cluster$sample = substr(cluster$sample,1,15)
    c1 = subset(group[,1],group[,2]==cluster[j])
    c2 = subset(group[,1],group[,2]==cluster[k])
    
    pos_H=match(c1,colnames(rppa))
    pos_L=match(c2,colnames(rppa))
    pos_H=na.omit(pos_H)
    pos_L=na.omit(pos_L)
    
    rppa = rppa[rowSums(rppa)!=0,]
    p.value = apply(rppa, 1, function(x) {    #1表示行，2表示列
      t.test(x[pos_H], x[pos_L], var.equal=TRUE, na.rm=TRUE)$p.value
    })
    log2FC = apply(rppa, 1, function(x) {
      mean(x[pos_H], na.rm=TRUE) - mean(x[pos_L], na.rm=TRUE)
    })
    
    resC = cbind(rownames(rppa), p.value, log2FC)
    resC = resC[order(as.numeric(resC[,2])),]
    FDR = rep(1,dim(resC)[1]); for(i in 1:dim(resC)[1]) FDR[i] = as.numeric(resC[i,2]) * dim(resC)[1]/i
    resC = cbind(resC, FDR)
    colnames(resC) = c("protein","P_value", "log2FC","FDR")
    #以0.585为界的蛋白质太少，故以0.322(1.25倍)为界
    #fdr<0.05,log2FC</>0
    #fdr<0.05,log2FC</>0.138(1.1倍)
    sum_low=resC[which(as.numeric(resC[,3])<(-0.138 )&as.numeric(resC[,4])<0.05),]
    sum_low = sum_low[order(as.numeric(sum_low[,3])),] #按照log2FC升序排列
    write.csv(sum_low,paste0("result_1v1/fc1.1/",cluster[j],"vs",cluster[k],"_low_score.csv"),row.names=FALSE)
    
    sum_high=resC[which(as.numeric(resC[,3])> 0.138 &as.numeric(resC[,4])<0.05),]
    sum_high = sum_high[order(as.numeric(sum_high[,3]), decreasing = TRUE), ] #按照log2FC降序排列	
    write.csv(sum_high,paste0("result_1v1/fc1.1/",cluster[j],"vs",cluster[k],"_high_score.csv"),row.names=FALSE)
  }
}

##########对每一种亚型的差异蛋白进行overlap
setwd("/Users/moonly/Desktop/pancan/protein/result_1v1/fc1.1/")

rm(list=ls())
#StC1
gene1a = read.csv("cluster1vscluster2_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene2a = read.csv("cluster1vscluster3_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene3a = read.csv("cluster1vscluster4_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene1 = rbind(gene1a, gene2a, gene3a)
dim(gene1)
# [1] 265   4

View(head(gene1))
overlap1 = table(gene1[,1])
result1 = subset(overlap1, overlap1 == 3)
dim(result1)
# [1] 21

write.csv(result1, "St1_protein.csv", row.names = FALSE)

#StC2
gene1b = read.csv("cluster1vscluster2_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene2b = read.csv("cluster2vscluster3_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene3b = read.csv("cluster2vscluster4_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene2 = rbind(gene1b, gene2b, gene3b)
dim(gene2)
# [1] 244   4
View(head(gene2))
overlap2 = table(gene2[,1])
result2 = subset(overlap2, overlap2 == 3)
dim(result2)
# [1] 17

write.csv(result2, "cluster_protein.csv", row.names = FALSE)

#StC3
gene1c = read.csv("cluster1vscluster3_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene2c = read.csv("cluster2vscluster3_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene3c = read.csv("cluster3vscluster4_high_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene3 = rbind(gene1c, gene2c, gene3c) #gene1c, 
dim(gene3)
# [1] 294   4

overlap3 = table(gene3[,1])
result3 = subset(overlap3, overlap3 == 3)
dim(result3)
# [1] 70

write.csv(result3, "cluster3_protein.csv", row.names = FALSE)

#StC4
gene1d = read.csv("cluster1vscluster4_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene2d = read.csv("cluster2vscluster4_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene3d = read.csv("cluster3vscluster4_low_score.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
gene4 = rbind(gene1d, gene2d, gene3d)
dim(gene4)
# [1] 231   4

overlap4 = table(gene4[,1])
result4 = subset(overlap4, overlap4 == 3)
dim(result4)
# [1] 49

write.csv(result4, "cluster_protein.csv", row.names = FALSE)


################################Student's t-Test 1v4##############################
setwd("/Users/moonly/Desktop/pan_cancer/pancan/protein/")
rm(list=ls())
rppa = read.table("pancan_RPPA.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, row.names=1, sep=",")
View(head(rppa))
dim(rppa)

rppa = na.omit(rppa) #不去除NA无法进行student.t.test
dim(rppa)
# [1]  210 7754

group <- read.csv("/Users/moonly/Desktop/pan_cancer/pancan/cluster4.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

St = c("cluster1","cluster2","cluster3","cluster4")

for(j in 1:4){
  #cluster$sample = substr(cluster$sample,1,15)
  c1 = subset(group[,1],group[,2]==cluster[j])
  c2 = subset(group[,1],group[,2]!=cluster[j])
  
  pos_H=match(c1,colnames(rppa))
  pos_L=match(c2,colnames(rppa))
  pos_H=na.omit(pos_H)
  pos_L=na.omit(pos_L)
  
  rppa = rppa[rowSums(rppa)!=0,]
  p.value = apply(rppa, 1, function(x) {    #1表示行，2表示列
    t.test(x[pos_H], x[pos_L], var.equal=TRUE, na.rm=TRUE)$p.value
  })
  log2FC = apply(rppa, 1, function(x) {
    mean(x[pos_H], na.rm=TRUE) - mean(x[pos_L], na.rm=TRUE)
  })
  
  resC = cbind(rownames(rppa), p.value, log2FC)
  resC = resC[order(as.numeric(resC[,2])),]
  FDR = rep(1,dim(resC)[1]); for(i in 1:dim(resC)[1]) FDR[i] = as.numeric(resC[i,2]) * dim(resC)[1]/i
  resC = cbind(resC, FDR)
  colnames(resC) = c("protein","P_value", "log2FC","FDR")
  #以0.585为界的蛋白质太少，故以0.322(1.25倍)为界
  sum_low=resC[which(as.numeric(resC[,3])<(-0.322) &as.numeric(resC[,4])<0.05),]
  sum_low = sum_low[order(as.numeric(sum_low[,3])),] #按照log2FC升序排列
  write.csv(sum_low,paste0("result_1v4/",cluster[j],"_low_score.csv"),row.names=FALSE)
  
  sum_high=resC[which(as.numeric(resC[,3])>=0.322 &as.numeric(resC[,4])<0.05),]
  sum_high = sum_high[order(as.numeric(sum_high[,3]), decreasing = TRUE), ] #按照log2FC降序排列	
  write.csv(sum_high,paste0("result_1v4/",cluster[j],"_high_score.csv"),row.names=FALSE)
}

####################################################8 Pathway VS Pathway#########################################################
##########ratio: pathway vs pathway############
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS")
ssGSEA = read.csv("score_pathway.csv",header=T,check.names=FALSE, stringsAsFactor=FALSE, row.names = 1)
write.csv(ssGSEA, "score_col.csv")
ssGSEA = t(ssGSEA)
dim(ssGSEA)
# [1]     6 10323
View(head(ssGSEA))
write.csv(ssGSEA, "score_row.csv")



#############compare among subtypes#####################
###################合并##############

ki<-read.csv("score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(ki)[1]=c("sample")
group<-read.csv("cluster_2list.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

cordata<-merge(group,ki,by.x="sample",by.y="sample")
View(head(cordata))
dim(cordata)
# [1] 10323    8

write.csv(cordata,"cluster_combine_new.csv",row.names=FALSE)

##########KW test######################
rm(list=ls())
a = read.csv("cluster_combine_new.csv",check.names = FALSE, stringsAsFactors=FALSE,header=TRUE)
View(head(a))
colnames(a)[2] <- "subtype"
pv=c()
meanresults = matrix()
for (i in 3:dim(a)[2]){
  c1=subset(a[,i], a[,2]=="cluster 1") 
  c2=subset(a[,i], a[,2]=="cluster 2") 
  c3=subset(a[,i], a[,2]=="cluster 3") 
  c4=subset(a[,i], a[,2]=="cluster 4") 
  # c5=subset(a[,i], a[,2]=="St5")
  pv1<- kruskal.test(list(c1, c2, c3, c4))$p.value  #KW
  pv2<-cbind(colnames(a)[i],pv1)
  pv<-rbind(pv,pv2)
  #mean value
  meanresult <- aggregate(a[,i],by=list(a$subtype),FUN=mean,na.rm=T)
  colnames(meanresult) = c("",colnames(a)[i])
  meanresults=cbind(meanresults,meanresult[2])
}

result=t(meanresults)
res=cbind(pv, result[-1,])
colnames(res) = c("cluster_score","pv_kw","cluster 1","cluster 2","cluster 3","cluster 4")
View(res)

write.csv(res,"cluster_kw_new.csv",row.names = FALSE)


#########wilcox###################################
rm(list=ls())
cordata <- read.csv("cluster_combine_new.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)

St = c("cluster 1","cluster 2","cluster 3","cluster 4")
resj = c()
for(j in 1:3){
  resk = c()
  for(k in (j+1):4){
    pv=c() 
    result=c()
    for(f in 3:dim(cordata)[2]){ 
      high_pos = which(cordata[,2]==St[j])##高的组所对应样本位置
      low_pos = which(cordata[,2]==St[k])#低的组所对应样本位置
      
      pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
      pv2 = cbind(colnames(cordata)[f],pvG,pvL)
      pv<-rbind(pv,pv2)   
      rownames(pv)=NULL  
    }
    ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
    pvalue1<-pv[ord,]
    fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
    pv2<-cbind(pvalue1,fdr)
    colnames(pv2)<-c("St_score","pvG","pvL","fdr")
    result=rbind(result,pv2)
    resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
    rownames(resu1) = resu1[,1]
    resu1 = resu1[,-1]
    colnames(resu1)[1] = paste0("pvG_",St[j]) 
    colnames(resu1)[2] = paste0("pvL_",St[k])   
    resk = cbind(resk, resu1)   
  }
  resj = cbind(resj, resk) 
}
View(resj)
dim(resj)
# [1] 6 18
write.csv(resj,"cluster_score_wilcox_new.csv")

#################clinical###################################################
###################合并##############
rm(list = ls())
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS")
ki<-read.csv("score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
colnames(ki)[1]=c("sample")
group<-read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical/therapy/pan_therapy.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)
dim(group)
group = group[, c(1,5)]#tumor_stage
# [1] 9914    4 tide
# [1] 10295    26 stage
# [1] 12591     7 grade+treatment
# [1] 8472    3 metastatic 

cordata = group[, c(1, 5, 24:26)] #tumor_stage
cordata = group[, c(1, 7, 24:26)] #clinical_stage
group[,1] = substr(group[,1], 1, 15)
group = group[,c(-5:-2)] #grade+treatment

ki[,1] = substr(ki[,1], 1, 12) #metastatic 

cordata<-merge(group,ki,by.x="sample",by.y="sample")
View(head(cordata))
dim(cordata)
# [1] 9901    10
# [1] 10295     6 #grade+treatment
# [1] 10295     8 #tumor_stage
# [1] 8331    9 #metastatic
# dim(cordata[!duplicated(cordata[,1]),])
# [1] 8134    9
cordata = cordata[!duplicated(cordata[,1]),]
write.csv(cordata,"clinical/therapy/histological_grade_combine.csv",row.names=FALSE)

# 移除第二列中有空白值的行
#cordata <- cordata[!is.na(cordata[,2]) & cordata[,2] != "", ]

###########wilcox
pv=c() 
result=c()
for(f in 3:dim(cordata)[2]){ 
  high_pos = which(cordata[,2]== "Low Grade")##高的组所对应样本位置
  low_pos = which(cordata[,2]== "High Grade")#低的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pvG_","Low Grade") 
colnames(resu1)[2] = paste0("pvL_","High Grade")
View(resu1)

write.csv(resu1,"clinical/therapy/histological_grade_wilcox.csv")

###############gender(排除只有一种性别的癌症)
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
clin = read.csv("pan_clinical_10year.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)[,c(1,3)]
score <- read.csv("pancan_ratio_clinical.csv", stringsAsFactors = FALSE, check.names = FALSE, header=TRUE)[, c(1:3, 18:23)]

cordata = merge(score, clin, by.x = "Sample", by.y = "sample")
dim(cordata)
# [1] 10295     10
View(head(cordata))
result = cordata[which(cordata[,10]!="PRAD"& cordata[,10]!="OV"& cordata[,10]!="CESC"& cordata[,10]!="TGCT"& cordata[,10]!="UCEC"& cordata[,10]!="UCS"), ]
dim(result)
# [1] 8454    10
View(head(result))
write.csv(result, "gender_cancer_score.csv", row.names = F)

cordata = result[,-10]
table(cordata[,3])
#    0    1 
# 4072 4382 

pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,3]== 0)##female的组所对应样本位置
  low_pos = which(cordata[,3]== 1)#male的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","female") 
colnames(resu1)[2] = paste0("pv_","male")
View(resu1)

write.csv(resu1,"gender_wilcox_new.csv")
###############gender(不排除只有一种性别的癌症)
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
clin = read.csv("pan_clinical_10year.csv",stringsAsFactors = FALSE, check.names = FALSE,header=TRUE)[,c(1,3)]
score <- read.csv("pancan_ratio_clinical.csv", stringsAsFactors = FALSE, check.names = FALSE, header=TRUE)[, c(1:3, 18:23)]

cordata = merge(score, clin, by.x = "Sample", by.y = "sample")
dim(cordata)
# [1] 10295     10
View(head(cordata))
result = cordata
write.csv(result, "gender_cancer_score(不排除).csv", row.names = F)

cordata = result[,-10]
table(cordata[,3])
#    0    1 
# 5276 5019

pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,3]== 0)##female的组所对应样本位置
  low_pos = which(cordata[,3]== 1)#male的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","female") 
colnames(resu1)[2] = paste0("pv_","male")
View(resu1)

write.csv(resu1,"gender_wilcox_new(不排除).csv")

#########age
res1=c()
for(j in 4:9){
  p.value = cor.test(as.vector(as.matrix(cordata[,2])), as.vector(as.matrix(cordata[,j])), method="spearman")$p.value
  estimate = cor.test(as.vector(as.matrix(cordata[,2])), as.vector(as.matrix(cordata[,j])), method="spearman")$estimate
  res = cbind(p.value, estimate)
  rownames(res) = colnames(cordata)[j]
  res1 = rbind(res1, res)
}
View(res1)
write.csv(res1, "spearman_age.csv")

################smoke
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
cancer_name = read.csv("cancer_name.csv",stringsAsFactors = FALSE, check.names = FALSE,header=F)
smoke = c()
for(i in 1:10){
  cancer = read.csv(paste0("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical/clinical(TCGA)/",tolower(cancer_name[i,1]), "_tcga_clinical.csv"),stringsAsFactors = FALSE, check.names = FALSE,header=T)
  # pos1 = grep("PATIENT_ID",colnames(cancer))
  pos = grep("TOBACCO_SMOKING_HISTORY_INDICATOR",colnames(cancer))
  clinical = cancer[,c(2,pos)]
  smoke = rbind(smoke, clinical)
}
View(smoke)
dim(smoke)
# [1] 3584    2
write.csv(smoke, "smoke_10cancer.csv", row.names = F)

rm(list=ls())
smoke = read.csv("smoke_info.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score[,1] = substr(st_score[,1], 1, 12)
colnames(st_score)[1] = "sample"
cordata = merge(smoke, st_score, by.x = "PATIENT_ID", "sample")
View(head(cordata))
dim(cordata)
# [1] 3509    8
dim(cordata[duplicated(cordata[,1]),])
# [1] 10  8
cordata = cordata[!duplicated(cordata[,1]),]
dim(cordata)
# [1] 3499    8
#  NO  YES 
#  795 2063
result = na.omit(cordata)
# > dim(result)
# [1] 2858    8
write.csv(result, "smoke_cluster_combine.csv", row.names = F)
###wilcox
pv=c() 
result=c()
for(f in 3:dim(cordata)[2]){ 
  high_pos = which(cordata[,2]== "YES")##smoke的组所对应样本位置
  low_pos = which(cordata[,2]== "NO")#no smoke的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","YES") 
colnames(resu1)[2] = paste0("pv_","NO")
View(resu1)
write.csv(resu1,"smoke_wilcox.csv")

#############alcohol#####################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
alcohol = read.csv("TCGA_alcohol_history.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score[,1] = substr(st_score[,1], 1, 12)
colnames(st_score)[1] = "sample"
cordata = merge(alcohol, st_score, by.x = "patient_id", "sample")
View(head(cordata))
dim(cordata)
# [1] 1011    9
dim(cordata[duplicated(cordata[,1]),])
# [1] 4 9
cordata = cordata[!duplicated(cordata[,1]),]
dim(cordata)
# [1] 1007    9
table(cordata[,3])
#  No Yes 
# 349 658
write.csv(cordata, "alcohol_cluster_combine.csv", row.names = FALSE)
###wilcox
pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,3]== "Yes")##alcohol的组所对应样本位置
  low_pos = which(cordata[,3]== "No")#no alcohol的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","YES") 
colnames(resu1)[2] = paste0("pv_","NO")
View(resu1)
write.csv(resu1,"alcohol_wilcox.csv")

##########################cbioportal clinical info####################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
clinical_info = read.csv("clinical_riskfactor.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
colnames(st_score)[1] = "sample"
for(j in 1:5){
  cancer_name = clinical_info[,j]
  cancer_name = subset(cancer_name, cancer_name!= "")
  result = c()
  for(i in 1:length(cancer_name)[1]){
    cancer = read.csv(paste0("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical/clinical_cbioportal/",tolower(cancer_name[i]), "_tcga_pan_can_atlas_2018_clinical_data.csv"),stringsAsFactors = FALSE, check.names = FALSE,header=T)
    pos1 = grep("Sample ID",colnames(cancer))
    pos2 = grep(colnames(clinical_info)[j],colnames(cancer))
    clinical = cancer[,c(pos1,pos2)]
    result = rbind(result, clinical)
  }
  colnames(result)[1] = "Sample"
  cordata = merge(result, st_score, by.x = "Sample", "sample")
  write.csv(cordata, paste0("pan_", colnames(clinical_info)[j],"_combine.csv"), row.names = F)
}
###wilcox-radiation therapy
cordata = read.csv("pan_Radiation Therapy_combine.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
table(cordata[,2])
#   No  Yes 
# 5937 2528
pv=c() 
result=c()
for(f in 3:dim(cordata)[2]){ 
  high_pos = which(cordata[,2]== "Yes")##therapy的组所对应样本位置
  low_pos = which(cordata[,2]== "No")#no therapy的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","YES") 
colnames(resu1)[2] = paste0("pv_","NO")
View(resu1)
write.csv(resu1,"Radiation_Therapy_wilcox.csv")

###############spearman##############
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
rm(list = ls())
cordata = read.csv("pan_MSIsensor Score_combine.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
cordata = na.omit(cordata)
dim(cordata)
# [1] 7588    8 hypoxia
# [1] 9656    8 Aneuploidy
#[1] 9104    8 MSI MANTIS
#[1] 9403    8 MSIsensor
res1=c()
for(j in 3:8){
  p.value = cor.test(as.vector(as.matrix(cordata[,2])), as.vector(as.matrix(cordata[,j])), method="spearman")$p.value
  estimate = cor.test(as.vector(as.matrix(cordata[,2])), as.vector(as.matrix(cordata[,j])), method="spearman")$estimate
  res = cbind(p.value, estimate)
  rownames(res) = colnames(cordata)[j]
  res1 = rbind(res1, res)
}
View(res1)
write.csv(res1, "spearman_MSIsensor Score.csv")

###################with HPV####################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
cesc_hpv = read.csv("CESC_HPV.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
hnsc_hpv = read.csv("HNSC_HPV.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
hpv = rbind(cesc_hpv, hnsc_hpv)
dim(hpv)
# [1] 835   2
View(head(hpv))

st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score[,1] = substr(st_score[,1], 1, 12)
colnames(st_score)[1] = "sample"
cordata = merge(hpv, st_score, by.x = "patient", "sample")
View(head(cordata))
dim(cordata)
# [1] 828   8
dim(cordata[duplicated(cordata[,1]),])
# [1] 4 8
cordata = cordata[!duplicated(cordata[,1]),]
dim(cordata)
# [1] 824   8
table(cordata[,2])
# [Not Available]   indeterminate        Negative        Positive 
#             409               1              95             319
write.csv(cordata, "HPV_cluster_combine.csv", row.names = FALSE)
###wilcox
pv=c() 
result=c()
for(f in 3:dim(cordata)[2]){ 
  high_pos = which(cordata[,2]== "Positive") #Positive的组所对应样本位置
  low_pos = which(cordata[,2]== "Negative") #Negative的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[3:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","Positive") 
colnames(resu1)[2] = paste0("pv_","Negative")
View(resu1)
write.csv(resu1,"HPV_wilcox.csv")

###########################HBV#################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
hbv = read.csv("HBVHCV_LIHC.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)

dim(hbv)
# [1] 196   3
View(head(hbv))
hbv[,1] = substr(hbv[,1], 1, 15)
st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
# st_score[,1] = substr(st_score[,1], 1, 12)
colnames(st_score)[1] = "sample"
cordata = merge(hbv, st_score, by.x = "Barcode", "sample")
View(head(cordata))
dim(cordata)
#[1] 193   9
dim(cordata[duplicated(cordata[,1]),])
# [1] 0 9

table(cordata[,2]) #hbv
# neg pos 
# 149  44
table(cordata[,3]) #hcv
# neg pos 
# 158  35
write.csv(cordata, "HBVHCV_cluster_combine.csv", row.names = FALSE)
###wilcox
pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,3]== "pos") #Positive的组所对应样本位置
  low_pos = which(cordata[,3]== "neg") #Negative的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","Positive") 
colnames(resu1)[2] = paste0("pv_","Negative")
View(resu1)
write.csv(resu1,"HCV_wilcox.csv")
######################with HBV#######################
pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,2]== "pos") #Positive的组所对应样本位置
  low_pos = which(cordata[,2]== "neg") #Negative的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","Positive") 
colnames(resu1)[2] = paste0("pv_","Negative")
View(resu1)
write.csv(resu1,"HBV_wilcox.csv")
######################with EBV#######################
rm(list = ls())
gc()
setwd("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/clinical")
ebv = read.csv("EBV.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)

dim(ebv)
# [1] 390   3
View(head(ebv))
st_score = read.csv("/Users/moonly/Desktop/Pan-cancer/pan-results/VS/score_col.csv",stringsAsFactors = FALSE, check.names = FALSE,header=T)
st_score[,1] = substr(st_score[,1], 1, 12)
colnames(st_score)[1] = "sample"
cordata = merge(ebv, st_score, by.x = "sample", "sample")
View(head(cordata))
dim(cordata)
#[1] 390   9
dim(cordata[duplicated(cordata[,1]),])
# [1] 0 9

table(cordata[,3]) #ebv
# Negative Positive 
#      359       31
write.csv(cordata, "EBV_cluster_combine.csv", row.names = FALSE)
###wilcox
pv=c() 
result=c()
for(f in 4:dim(cordata)[2]){ 
  high_pos = which(cordata[,3]== "Positive") #Positive的组所对应样本位置
  low_pos = which(cordata[,3]== "Negative") #Negative的组所对应样本位置
  
  pvG = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="greater")$p.value
  pvL = wilcox.test(as.numeric(cordata[,f][high_pos]),as.numeric(cordata[,f][low_pos]),alternative="less")$p.value
  pv2 = cbind(colnames(cordata)[f],pvG,pvL)
  pv<-rbind(pv,pv2)   
  rownames(pv)=NULL  
}
ord<-order(as.numeric(pv[,2]),decreasing=FALSE)
pvalue1<-pv[ord,]
fdr<-rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,2])*dim(pvalue1)[1]/n
pv2<-cbind(pvalue1,fdr)
colnames(pv2)<-c("St_score","pvG","pvL","fdr")
result=rbind(result,pv2)
resu1 = result[sapply(colnames(cordata)[4:dim(cordata)[2]], function(x) {which(result[,1] == x)}), ] #reorder as the names
rownames(resu1) = resu1[,1]
resu1 = resu1[,-1]
colnames(resu1)[1] = paste0("pv_","Positive") 
colnames(resu1)[2] = paste0("pv_","Negative")
View(resu1)
write.csv(resu1,"EBV_wilcox.csv")










