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

kegg_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    fgsea_res <- runEA(de_res, fgseaForSe, gs = kegg.gs)
    fgsea_sigsets <- get_sig_fraction(fgsea_res)
    return(fgsea_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

kegg_res <- summary(unlist(kegg_res))

path_prefix <- "outdata/downsampling/type_I_error/kegg/"
if (!dir.exists(path_prefix)){
    dir.create(path_prefix, showWarnings = FALSE)
}

saveRDS(kegg_res, file = paste0(path_prefix, "fgsea.rds"))


# ------------------------------------------------------------------------------------------------------------------------
# check with GOBP pathways
# ------------------------------------------------------------------------------------------------------------------------

gobp_res <- bplapply(seq_len(nperm), function(i){
    se <- golub_se[, col_ids[i, ]]
    se$GROUP <- c(rep(0, group_size), rep(1, group_size))
    de_res <- runDE(list(se), de.method = "limma")
    fgsea_res <- runEA(de_res, fgseaForSe, gs = gobp.gs)
    fgsea_sigsets <- get_sig_fraction(fgsea_res)
    return(fgsea_sigsets)
}, BPPARAM = MulticoreParam(detectCores() - 1))

gobp_res <- summary(unlist(gobp_res))

path_prefix <- "outdata/downsampling/type_I_error/gobp/"
if (!dir.exists(path_prefix)){
    dir.create(path_prefix, showWarnings = FALSE)
}

saveRDS(gobp_res, file = paste0(path_prefix, "fgsea.rds"))

