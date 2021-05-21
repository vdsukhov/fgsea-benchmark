library(ggplot2)
library(ggpubr)
source("methods/utils/type_I_plot_utils.R")
source("methods/utils/gseabenchmarker_colors.R")

path_to_kegg_res <- "outdata/downsampling/type_I_error/kegg/"
kegg_res_files <- list.files(path_to_kegg_res, pattern = ".rds")
kegg_data <- lapply(kegg_res_files, function(x){
    readRDS(paste0(path_to_kegg_res, x))
})
kegg_data <- do.call(rbind, kegg_data)
rownames(kegg_data) <- gsub(".rds", "", kegg_res_files)
kegg_data <- t(kegg_data)

plotTypeIError2(kegg_data) + labs(title="Type I error rate (KEGG)", subtitle = "Golub downsampled data")


path_to_gobp_res <- "outdata/downsampling/type_I_error/gobp/"
gobp_res_files <- list.files(path_to_gobp_res, pattern = ".rds")
gobp_data <- lapply(gobp_res_files, function(x){
    readRDS(paste0(path_to_gobp_res, x))
})
gobp_data <- do.call(rbind, gobp_data)
rownames(gobp_data) <- gsub(".rds", "", gobp_res_files)
gobp_data <- t(gobp_data)

plotTypeIError2(gobp_data) + labs(title="Type I error rate (GOBP)", subtitle = "Golub downsampled data")
