library(data.table)
library(spatstat)


genes = c('cdh1', 'cdh6', 'cldn5', 'emx2', 'erbb4', 'foxc2', 'hnf1b', 'jag1', 'lef1', 'lhx1', 'mafb', 'mecom', 'pappa2', 'pax2', 'pou3f3', 'sall1', 'six2', 'sox9', 'troma1', 'wt1')
imgs_cdh1 = c('00000035', '00000036', '00000037', '00000038', '00000039', '00000040', '00000041', '00000042', '00000080', '00000081', '00000082', '00000083', '00000084', '00000085', '00000086', '00000098', '00000118', '00000119', '00000120', '00000121')
imgs_cdh6 = c('00000118', '00000119', '00000120', '00000121')
imgs_cldn5 = c('00000126', '00000127', '00000128')
imgs_emx2 = c('00000111', '00000112', '00000113', '00000114')
imgs_erbb4  = c('00000111', '00000112', '00000113', '00000114', '00000122', '00000123', '00000124', '00000125')
imgs_foxc2 = c('00000085', '00000086')
imgs_hnf1b = c('00000082', '00000083', '00000084', '00000098')
imgs_jag1 = c('00000029', '00000031', '00000032', '00000033', '00000034', '00000035', '00000036', '00000037', '00000038', '00000039', '00000040', '00000041', '00000042', '00000082', '00000083', '00000084', '00000098', '00000099', '00000100', '00000101', '00000102', '00000126', '00000127', '00000128')
imgs_lef1 = c('00000029', '00000031')
imgs_lhx1 = c('00000080', '00000081')
imgs_mafb = c('00000099', '00000100', '00000101', '00000102')
imgs_mecom = c('00000126', '00000127', '00000128')
imgs_pappa2 = c('00000122', '00000123', '00000124', '00000125')
imgs_pax2 = c('00000080', '00000081')
imgs_pou3f3 = c('00000099', '00000100', '00000101', '00000102')
imgs_sall1 = c('00000029', '00000031')
imgs_six2 = c('00000035', '00000036', '00000037', '00000038', '00000039')
imgs_sox9 = c('00000032', '00000033', '00000034')
imgs_troma1 = c('00000032', '00000033', '00000034')
imgs_wt1 = c('00000040', '00000041', '00000042', '00000085', '00000086', '00000111', '00000112', '00000113', '00000114', '00000122', '00000123', '00000124', '00000125')

n_clust = args[1]; #number of clusters

dir = "/Genomics/ogtr04/sealfon/fetal_kidney_usc/data/180721_images/";

nslice = 286;

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

my_blur <- function(mat, s) {
        blur(im(as.matrix(mat)), sigma=s)
}

for (i in c(1:length(genes))){
    files_path = paste(dir, genes[i], sep="")
    gene_imgs = get(paste("imgs_", genes[i], sep=""))
    for (j in c(1:length(gene_imgs))) {
        files = vector()
        for (k in 1:286) {
            files[k] = paste(dir, "SSB_text_sequences/", "z:", k, "-286 - reg-", gene_imgs[j], "-scl-", genes[i], ".txt", sep="")
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


all <- cbind(get("cdh1_v"), get("cdh6_v"), get("cldn5_v"), get("emx2_v"), get("erbb4_v"), get("foxc2_v"), get("hnf1b_v"), get("jag1_v"), get("lef1_v"), get("lhx1_v"), get("mafb_v"), get("mecom_v"), get("pappa2_v"), get("pax2_v"), get("pou3f3_v"), get("sall1_v"), get("six2_v"), get("sox9_v"), get("troma1_v"), get("wt1_v"))
all_s <- scale(all)

clusters = kmeans(all_s, as.numeric(n_clust), iter.max = 100, nstart = 100)
clusters_image <- matrix(clusters$cluster, nrow = 62634, ncol = 190)

clust_name = paste(dir, "blur_out/clusters_all_", n_clust, "_", sigma, ".clust", sep="")
centers_name = paste(dir, "blur_out/centers_all_", n_clust, "_", sigma,  ".clust", sep="")
tot_withinss_name = paste(dir, "blur_out/tot_withinss_all_", n_clust, "_", sigma, ".clust", sep="")
write.table(clusters_image, clust_name, row.names=FALSE, col.names=FALSE)
write.table(clusters$centers, centers_name, row.names=FALSE, col.names=FALSE)
write.table(clusters$tot.withinss, tot_withinss_name, row.names=FALSE, col.names=FALSE)
