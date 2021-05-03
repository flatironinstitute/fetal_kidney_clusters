#! /usr/bin/env Rscript

#Script for generating spatial pattern clusters

args = commandArgs(trailingOnly = TRUE) #k hyperparameter for k means clustering

library(MASS)
library(data.table)


genes = c('cdh1', 'cdh6', 'cldn5', 'emx2', 'erbb4', 'foxc2', 'hnf1b', 'jag1', 'lef1', 'lhx1', 'mafb', 'mecom', 'pax2', 'pou3f3', 'sall1', 'six2', 'sox9', 'troma1', 'wt1')
imgs_cdh1 = c('00000007', '00000009', '00000010', '00000011', '00000012', '00000013', '00000014', '00000060', '00000061', '00000062', '00000063', '00000064', '00000065', '00000066', '00000067', '00000068', '00000069', '00000115', '00000116', '00000117')
imgs_cdh6 = c('00000115', '00000116', '00000117')
imgs_cldn5 = c('00000129', '00000130', '00000131')
imgs_emx2 = c('00000107', '00000108', '00000109', '00000110')
imgs_erbb4  = c('00000107', '00000108', '00000109', '00000110')
imgs_foxc2 = c('00000067', '00000068', '00000069')
imgs_hnf1b = c('00000064', '00000065', '00000066')
imgs_jag1 = c('00000001', '00000002', '00000003', '00000004', '00000005', '00000006', '00000007', '00000009', '00000010', '00000011', '00000012', '00000013', '00000014', '00000064', '00000065', '00000066', '00000103', '00000104', '00000105', '00000106', '00000129', '00000130', '00000131')
imgs_lef1 = c('00000001', '00000002', '00000003')
imgs_lhx1 = c('00000060', '00000061', '00000062', '00000063')
imgs_mafb = c('00000103', '00000104', '00000105', '00000106')
imgs_mecom = c('00000129', '00000130', '00000131')
imgs_pax2 = c('00000060', '00000061', '00000062', '00000063')
imgs_pou3f3 = c('00000103', '00000104', '00000105', '00000106')
imgs_sall1 = c('00000001', '00000002', '00000003')
imgs_six2 = c('00000007', '00000009', '00000010')
imgs_sox9 = c('00000004', '00000005', '00000006')
imgs_troma1 = c('00000004', '00000005', '00000006')
imgs_wt1 = c('00000011', '00000012', '00000013', '00000014', '00000067', '00000068', '00000069', '00000107', '00000108', '00000109', '00000110')

n_clust = args[1]; #number of clusters

dir = "/Genomics/ogtr04/sealfon/fetal_kidney_usc/data/180730_RV_images/";

nslice = 224;

#Make number string for file name
make_number_str <- function(n) {
    number_str = ""             
    if (length(n) == 1) {
         number_str = paste("000", n, sep = "")
    }   
    else if (length(n) == 2) {
         number_str = paste("00", n, sep="")
    }
    else if (length(n) == 3) {
         number_str = paste("0", n, sep="")
    }
    number_str
}

for (i in c(1:length(genes))){
    files_path = paste(dir, genes[i], sep="")
    gene_imgs = get(paste("imgs_", genes[i], sep=""))
    for (j in c(1:length(gene_imgs))) {
        files = vector()
        for (k in 1:224) {
            files[k] = paste(dir, "RV_text_sequences/", "z:", k, "-224 - reg-", gene_imgs[j], "-scl-", genes[i], ".txt", sep="")
            #print(files[k])
        }       
        files.list <- lapply(files, read.table, header=FALSE)
        data_str = paste(genes[i], "_", gene_imgs[j], sep="")
        assign(data_str, rbindlist(files.list))
    }
    print(genes[i])
}

for (i in c(1:length(genes))){
    gene_imgs = get(paste("imgs_", genes[i], sep=""))
    assign(paste(genes[i], "_v", sep=""), rep(0, length(as.vector(as.matrix(get(paste(genes[1], "_",  imgs_cdh1[1], sep="")))))))       
    gene_v = get(paste(genes[i], "_v", sep=""))
    for (j in 1:length(gene_imgs)){
        gene_v = gene_v + as.vector(scale(as.vector(as.matrix(get(paste(genes[i], "_",  gene_imgs[j], sep=""))))))      
    }
    gene_v = gene_v / length(gene_imgs) 
    assign(paste(genes[i], "_v", sep=""), gene_v)
  
}


all <- cbind(get("cdh1_v"), get("cdh6_v"), get("cldn5_v"), get("emx2_v"), get("erbb4_v"), get("foxc2_v"), get("hnf1b_v"), get("jag1_v"), get("lef1_v"), get("lhx1_v"), get("mafb_v"), get("mecom_v"), get("pax2_v"), get("pou3f3_v"), get("sall1_v"), get("six2_v"), get("sox9_v"), get("troma1_v"), get("wt1_v"))
all_s <- scale(all)

clusters = kmeans(all_s, as.numeric(n_clust), iter.max = 100, nstart = 100)
clusters_image <- matrix(clusters$cluster, nrow = 36960, ncol = 151)

clust_name = paste(dir, "out/clusters_all_", n_clust, ".clust", sep="")
centers_name = paste(dir, "out/centers_all_", n_clust, ".clust", sep="")
tot_withinss_name = paste(dir, "out/tot_withinss_all_", n_clust,  ".clust", sep="")
write.table(clusters_image, clust_name, row.names=FALSE, col.names=FALSE)
write.table(clusters$centers, centers_name, row.names=FALSE, col.names=FALSE)
write.table(clusters$tot.withinss, tot_withinss_name, row.names=FALSE, col.names=FALSE)
