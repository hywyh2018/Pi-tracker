#conda activate R44
library(data.table)
library(dplyr)
library(lessR)
library(ggplot2)
library(gridExtra)

args=commandArgs(TRUE)
gene = as.character(args[1])
#diam = as.character(args[2])
#thresh = as.character(args[2])

#candidates = data.table::fread(paste0("/mnt/sda/tsk/TWMR/Panc/",gene,"/gsmr/3.outcome.txt"),header=F,sep="\t")
gtf = data.table::fread("/mnt/sda/resource/ref/gtf/Homo_sapiens.GRCh38.100.gtf.h.geneANNO",header=T,sep="\t")
tissues = data.table::fread("/mnt/sdb/resource/xQTL/eQTL/catalogue/eQTL-Catalogue_info.txt",header=T,sep="\t")
Xin = data.table::fread("/mnt/sda/tsk/TWMR/Panc/pwas/UKB+deCODE.pwas.candidates9.txt",sep="\t")
geneID = gtf[gtf$GeneName==gene,"TargetID"]

index <- which(!is.na(Xin[, "TWAS.Z"]))  # 获取不包含NA的行索引
Xin <- Xin[index, ]  # 使用这些索引来过滤数据框
XinC <- Xin[Xin$ID %in% geneID & Xin$MODEL=="top1",]
XinC$TWAS.P = as.numeric(XinC$TWAS.P)
XinC$TWAS.Z = as.numeric(XinC$TWAS.Z)
XinC = XinC %>% mutate(Out = paste0(deCODE,"_",gene), bxy=TWAS.Z, log10P = -log10(TWAS.P))
XinC$ExpoID = gsub("QTD000", "", XinC$PANEL)

GTExID = tissues[tissues$study_label %in% c("GTEx") & tissues$quant_method=="ge","Phenotype"]
GTEx <- XinC[XinC$PANEL %in% GTExID$Phenotype,]
for (i in unique(GTEx$Out))
{
#Xin = Xin %>% mutate(log10P = -log10(p))
GTExi <- GTEx[GTEx$Out == i,]
ESsum <- ifelse(length(GTExi$bxy)<=1, GTExi$bxy, sum(abs(GTExi$bxy)))
GTExi$pct = round(abs(GTExi$bxy)/ESsum * 100,2)
FPATH=paste("/mnt/sda/tsk/TWMR/Panc/candidates9/pwas/pie/GTEx_",i,".donut.pdf",sep = "")
PieChart(x=ExpoID, y=pct, data = GTExi, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH)
}

TwinsUKID = tissues[tissues$study_label %in% c("TwinsUK") & tissues$quant_method=="ge","Phenotype"]
TwinsUK <- XinC[XinC$PANEL %in% TwinsUKID$Phenotype,]
for (i in unique(TwinsUK$Out))
{
#Xin = Xin %>% mutate(log10P = -log10(p))
TwinsUKi <- TwinsUK[TwinsUK$Out == i,]
ESsum <- ifelse(length(TwinsUKi$bxy)<=1, TwinsUKi$bxy, sum(abs(TwinsUKi$bxy)))
TwinsUKi$pct = round(abs(TwinsUKi$bxy)/ESsum * 100,2)
TwinsUKi$label <- paste(TwinsUKi$ExpoID, "\n", TwinsUKi$pct, "%", sep = "")  # 合并标签与比例
FPATH=paste("/mnt/sda/tsk/TWMR/Panc/candidates9/pwas/pie/TwinsUK_",i,".donut.pdf",sep = "")
#PieChart(x=label, y=pct, data = TwinsUKi, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH, values_position="label")
PieChart(x=ExpoID, y=pct, data = TwinsUKi, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH)
}
