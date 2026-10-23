#conda activate R412
library(data.table)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

args=commandArgs(TRUE)
gene = as.character(args[1])
diam = as.character(args[2])
outpath = as.character(args[3])
#thresh = as.character(args[4])

candidates = data.table::fread(paste0("/mnt/sda/tsk/TWMR/Panc/",gene,"/gsmr/3.outcome.txt"),header=F,sep="\t")
gtf = data.table::fread("/mnt/sda/resource/ref/gtf/Homo_sapiens.GRCh38.100.gtf.h.geneANNO",header=T,sep="\t")
tissues = data.table::fread("/mnt/sdb/resource/xQTL/eQTL/catalogue/eQTL-Catalogue_info.txt",header=T,sep="\t")
Xin = data.table::fread("/mnt/sda/tsk/TWMR/Panc/pwas/UKB+deCODE.pwas.candidates9.txt",sep="\t")
#chordIn = matrix(0,0,3)
#colnames(chordIn) = c("Taxa","Disease","log2FoldChange")
GTExID = tissues[tissues$study_label %in% c("GTEx","TwinsUK") & tissues$quant_method=="ge","Phenotype"]
geneID = gtf[gtf$GeneName==gene,"TargetID"]
#candidates = candidates[grep("INSRR", candidates$V1, invert = TRUE),]

#Xin$Expo = gsub(paste0(gene,"_"), "", Xin$Exposure)
index <- which(!is.na(Xin[, "TWAS.Z"]))  # 获取不包含NA的行索引
Xin <- Xin[index, ]  # 使用这些索引来过滤数据框
XinC <- Xin[Xin$ID %in% geneID & Xin$MODEL=="top1",]
XinC$TWAS.P = as.numeric(XinC$TWAS.P)
XinC$TWAS.Z = as.numeric(XinC$TWAS.Z)
XinC = XinC %>% mutate(Outcome = paste0(deCODE,"_",ID), bxy=TWAS.Z, log10P = -log10(TWAS.P))
XinC <- XinC[XinC$PANEL %in% GTExID$Phenotype,]
#XinC <- na.omit(XinC)
XinC$PANEL = gsub("QTD000", "", XinC$PANEL)
XinC$log10P[is.infinite(XinC$log10P)] <- 308  # 替换所有Inf为308
max_value <- max(XinC$bxy)
min_value <- min(XinC$bxy)
max_width <- max(XinC$log10P)
min_width <- min(XinC$log10P)
#max_width <- ifelse(is.infinite(max_width), 308, max_width)

chord <- XinC[,c("PANEL","Outcome","log10P")]
color_matrix <- XinC[,c("PANEL","Outcome","bxy")]
col_fun <- colorRamp2(c(min_value, 0, max_value), c("DeepSkyBlue","white","HotPink"))
color_links <- col_fun(color_matrix$bxy)

# Legend for colors
lgd_col <- Legend(
  title = "Effect size",
  col_fun = col_fun,
  at = c(min_value, 0, max_value),
  labels = c(min_value, 0, max_value)
)

# Legend for widths
lgd_width <- Legend(
  title = "-Log10(P)",
  labels = c(min_width, max_width),
  type = "lines",
  legend_gp = gpar(lwd = c(min_width, max_width), col = "black")
)

# Draw legends
pdf(paste0(outpath, "/", gene,".GTEx_TwinsUK.chord.pwas.pdf"), width=diam, height=diam)
chordDiagram(chord, directional = 1, transparency = 0.2, col=color_links)
draw(lgd_col, x = unit(0.95, "npc"), y = unit(0.9, "npc"))
draw(lgd_width, x = unit(0.95, "npc"), y = unit(0.2, "npc"))
dev.off()







