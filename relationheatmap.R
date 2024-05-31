library(psych)
library(pheatmap)
library(reshape2)
library(magrittr)
library(scico)
library(scales)
phy <- read.csv(file = "top50.csv", header = T, row.names=1)  #读取微生物丰度信息表
met <- read.csv(file = "hepatic_49m34_forgenus.csv", header = T, row.names=1)  #读取代谢物信息表
cor <- corr.test(phy, met, method = "pearson")              #计算相关性矩阵、p值矩阵
cmt <- cor$r
pmt <- cor$p
cmt.out<-cbind(rownames(cmt),cmt)
pmt.out<-cbind(rownames(pmt),pmt)
write.table(cmt.out,file="cor.txt",sep="\t",row.names=F)       #输出相关系数表格
write.table(pmt.out,file="pvalue.txt",sep="\t",row.names=F)    #输出p值表格
df <- melt(cmt,value.name="cor")
df$pvalue <- as.vector(pmt)
head(df)
write.table(df,file="cor-p.txt",sep="\t")          #输出：相关系数、p值表格
if (!is.null(pmt)) {
  ssmt <- pmt < 0.01
  pmt[ssmt] <- '**'
  smt <- pmt > 0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt] <- ''
} else {
  pmt <- F
}
pmt
mycol <- scico(100, palette = "vik")


mycol <- scico(100, palette = "bam")

pheatmap(cmt, scale = "none", cluster_row = T, cluster_col = F,
         display_numbers = pmt, fontsize_number = 12, number_color = "white",
         cellwidth = 20, cellheight = 20,color = mycol,filename="heatmap2.pdf")               #输出：相关性热图pdf文件




pheatmap(cmt, scale = "none",cluster_row = F, cluster_col = F,
         display_numbers = pmt, fontsize_number = 12, number_color = "black",
         cellwidth = 20, cellheight = 20)

