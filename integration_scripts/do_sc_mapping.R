#Integrate single-cell with spatial cluster data
load("Wk14NPC_SSB_resolvedv2.Robj")
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(plotly)
library(pdist)
#library(SAVER)

#Read single-cell data
genes <- c("CDH1", "CDH6",  "EMX2", "ERBB4", "FOXC2", "HNF1B", "JAG1", "LEF1", "LHX1", "MAFB", "MECOM", "PAPPA2", "PAX2", "POU3F3", "SALL1", "SIX2", "SOX9", "KRT8", "WT1")
saver <- readRDS("saver.rds") #gene expression imputed using SAVER for relevant genes
saver_scaled <- scale(t(saver$estimate))

#Map across only genes that have variable expression in the single-cell data
ssc_subset <- saver_scaled[,c("CDH6", "EMX2", "ERBB4",  "HNF1B", "JAG1","LHX1", "MAFB", "MECOM", "PAPPA2", "PAX2", "POU3F3", "SOX9", "KRT8", "WT1")]
tmp_names = sapply(strsplit(names(Wk14NPC_SSB_resolvedv2@ident), split="Wk14[A-Za-z]*_[SSB_]*"), "[[", 2)
ssc_subset = ssc_subset[tmp_names,]

#Read spatial cluster data
centers_14 <- read.table("centers/centers_14_nocldn5.txt", header=FALSE)
colnames(centers_14) <- c("CDH1", "CDH6", "EMX2", "ERBB4", "FOXC2", "HNF1B", "JAG1", "LEF1", "LHX1", "MAFB", "MECOM", "PAPPA2", "PAX2", "POU3F3", "SALL1", "SIX2", "SOX9", "KRT8", "WT1")

#Map across only genes that have variable expression in the single-cell data
centers_14_filt <- select(as.data.frame(centers_14), matches("EMX2|MECOM|POU3F3|SOX9|PAX2|PAPPA2|HNF1B|LHX1|ERBB4|KRT8|CDH6|JAG1|MAFB|WT1")) 

#Add a zero pattern to the spatial patterns
centers_14_filt_z <- rbind(centers_14_filt, rep(0,13))

#Find distace between each cell and each spatial pattern
m <- pdist(ssc_subset, centers_14_filt_z)
m2 <- as.matrix(m)

#Remove background patterns
m3 <- m2[,c(1,3,6,8,9,10,11,13,15)]
colnames(m3) <- c("p1", "p3", "p6", "p8", "p9", "p10", "p11", "p13", "z")
rownames(m3) <- rownames(ssc_subset)

#Map each cell to the closest spatial pattern
assignments <- colnames(m3)[apply(m3,1,which.min)]
names(assignments) <- names(Wk14NPC_SSB_resolvedv2@ident)
b <- cbind(assignments, as.numeric(as.character(Wk14NPC_SSB_resolvedv2@ident)))
colnames(b) <- c("assignments", "clust")
