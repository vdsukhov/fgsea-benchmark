library(data.table)
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
library(fgsea)
library(BiocParallel)
library(parallel)
source("methods/utils/custom_ea_methods.R")
source("methods/utils/type_I_error_utils.R")

kegg.gs <- getGenesets(org="hsa", db="kegg")
gobp.gs <- getGenesets(org="hsa", db="go", onto="BP")

golub_se <- readRDS("inpdata/GSEABenchmarking/golub.rds")
nperm <- 1000
group_size <- 4

group_1_indxs <- which(golub_se$GROUP == 1)



set.seed(1)
col_ids <- lapply(seq_len(nperm), function(i) sample(group_1_indxs, size = 2 * group_size))
col_ids <- do.call(rbind, col_ids)


# ------------------------------------------------------------------------------------------------------------------------
# check with KEGG pathways
# ------------------------------------------------------------------------------------------------------------------------

camera_kegg_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    camera_res <- runEA(de_res, methods = "camera", gs = kegg.gs, perm = 0) # can't pass inter.gene.cor argument to camera
    camera_sigsets <- get_sig_fraction(camera_res)
    return(camera_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

camera_kegg_res <- summary(unlist(camera_kegg_res))


cameraPr_kegg_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    camera_res <- runEA(de_res, methods = cameraPrForSe, gs = kegg.gs, perm = 0)
    camera_sigsets <- get_sig_fraction(camera_res)
    return(camera_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

cameraPr_kegg_res <- summary(unlist(cameraPr_kegg_res))

path_prefix <- "outdata/downsampling/type_I_error/kegg/"
if (!dir.exists(path_prefix)){
    dir.create(path_prefix, showWarnings = FALSE)
}

saveRDS(camera_kegg_res, file = paste0(path_prefix, "camera.rds"))
saveRDS(cameraPr_kegg_res, file = paste0(path_prefix, "cameraPr.rds"))



# ------------------------------------------------------------------------------------------------------------------------
# check with GOBP pathways
# ------------------------------------------------------------------------------------------------------------------------

camera_kegg_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    camera_res <- runEA(de_res, methods = "camera", gs = gobp.gs, perm = 0) # can't pass inter.gene.cor argument to camera
    camera_sigsets <- get_sig_fraction(camera_res)
    return(camera_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

camera_kegg_res <- summary(unlist(camera_kegg_res))


cameraPr_kegg_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    camera_res <- runEA(de_res, methods = cameraPrForSe, gs = gobp.gs, perm = 0)
    camera_sigsets <- get_sig_fraction(camera_res)
    return(camera_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

cameraPr_kegg_res <- summary(unlist(cameraPr_kegg_res))

path_prefix <- "outdata/downsampling/type_I_error/gobp/"
if (!dir.exists(path_prefix)){
    dir.create(path_prefix, showWarnings = FALSE)
}

saveRDS(camera_kegg_res, file = paste0(path_prefix, "camera.rds"))
saveRDS(cameraPr_kegg_res, file = paste0(path_prefix, "cameraPr.rds"))


