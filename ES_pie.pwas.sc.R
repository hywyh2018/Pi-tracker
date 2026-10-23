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

Perez2022ID = tissues[tissues$study_label %in% c("Perez_2022") & tissues$quant_method=="ge","Phenotype"]
Perez2022 <- XinC[XinC$PANEL %in% Perez2022ID$Phenotype,]
for (i in unique(Perez2022$Out))
{
#Xin = Xin %>% mutate(log10P = -log10(p))
Perez2022i <- Perez2022[Perez2022$Out == i,]
ESsum <- ifelse(length(Perez2022i$bxy)<=1, Perez2022i$bxy, sum(abs(Perez2022i$bxy)))
Perez2022i$pct = round(abs(Perez2022i$bxy)/ESsum * 100,2)
FPATH=paste("/mnt/sda/tsk/TWMR/Panc/candidates9/pwas/pie/Perez2022_",i,".donut.pdf",sep = "")
PieChart(x=ExpoID, y=pct, data = Perez2022i, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH)
}

OneK1KID = tissues[tissues$study_label %in% c("OneK1K") & tissues$quant_method=="ge","Phenotype"]
OneK1K <- XinC[XinC$PANEL %in% OneK1KID$Phenotype,]
for (i in unique(OneK1K$Out))
{
#Xin = Xin %>% mutate(log10P = -log10(p))
OneK1Ki <- OneK1K[OneK1K$Out == i,]
ESsum <- ifelse(length(OneK1Ki$bxy)<=1, OneK1Ki$bxy, sum(abs(OneK1Ki$bxy)))
OneK1Ki$pct = round(abs(OneK1Ki$bxy)/ESsum * 100,2)
FPATH=paste("/mnt/sda/tsk/TWMR/Panc/candidates9/pwas/pie/OneK1K_",i,".donut.pdf",sep = "")
PieChart(x=ExpoID, y=pct, data = OneK1Ki, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH)
#PieChart(x=label, y=pct, data = OneK1Ki, main = i, hole =0.5,rows=(pct>0),trans=.6, values="%", values.color=rgb(.3,.3,.3),width=3.5,height=3.5,pdf.file=FPATH, values_position="label")
}

