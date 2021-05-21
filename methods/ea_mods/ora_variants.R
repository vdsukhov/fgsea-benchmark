library(data.table)
library(GSEABenchmarkeR)
library(EnrichmentBrowser)
source("methods/utils/custom_ea_methods.R")

# ----------------------------------------------------------------------------------------------------------------------------------
# Loading and processing gene sets collections
# ----------------------------------------------------------------------------------------------------------------------------------

kegg_gs <- getGenesets(org="hsa", db="kegg")
gobp_gs <- getGenesets(org="hsa", db="go", onto="BP")

gs_min_size <- EnrichmentBrowser::configEBrowser("GS.MIN.SIZE")
gs_max_size <- EnrichmentBrowser::configEBrowser("GS.MAX.SIZE")
kegg_to_keep <- (lengths(kegg_gs) >= gs_min_size) & (lengths(kegg_gs) <= gs_max_size)
gobp_to_keep <- (lengths(gobp_gs) >= gs_min_size) & (lengths(gobp_gs) <= gs_max_size)

kegg_gs <- kegg_gs[kegg_to_keep]
gobp_gs <- gobp_gs[gobp_to_keep]

# ----------------------------------------------------------------------------------------------------------------------------------
# Loading `SummarizedExperiment` files for microarray and rnaseq data
# ----------------------------------------------------------------------------------------------------------------------------------
path_to_tcga <- "inpdata/GSEABenchmarking/TCGA_preproc/GSE62944_matched_limmavoom/"

microarray_data <- loadEData("geo2kegg", preproc=TRUE)

rnaseq_data <- lapply(list.files(path_to_tcga, pattern = ".rds"), function(fnm){
    se <- readRDS(paste0(path_to_tcga, fnm))
    return(se)
})
names(rnaseq_data) <- gsub(".rds", "", list.files(path_to_tcga, pattern = ".rds"))


# ----------------------------------------------------------------------------------------------------------------------------------
# Run differential expression analysis for microarray data
# ----------------------------------------------------------------------------------------------------------------------------------
microarray_data <- runDE(microarray_data, de.method = "limma", padj.method = "flexible")

if (!dir.exists("outdata/ea_mods/")){
    dir.create("outdata/ea_mods/microarray/kegg/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/ea_mods/microarray/gobp/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/ea_mods/rnaseq/kegg/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/ea_mods/rnaseq/gobp/", recursive = TRUE, showWarnings = FALSE)
}

gene_sets_list <- list(kegg = kegg_gs, gobp = gobp_gs)
data_list <- list(microarray = microarray_data, rnaseq = rnaseq_data)

for (data_ind in seq_along(data_list)){
    for (gs_ind in seq_along(gene_sets_list)){
        out_path <- paste0("outdata/ea_mods/",
                           names(data_list)[data_ind], "/",
                           names(gene_sets_list[gs_ind]))
        ora_res <- runEA(data_list[[data_ind]], methods = "ora",
                         gs = gene_sets_list[[gs_ind]], perm = 0,
                         save2file = TRUE, out.dir = out_path)
    }
}


for (data_ind in seq_along(data_list)){
    for (gs_ind in seq_along(gene_sets_list)){
        out_path <- paste0("outdata/ea_mods/",
                           names(data_list)[data_ind], "/",
                           names(gene_sets_list[gs_ind]))

        ora_res <- runEA(data_list[[data_ind]], methods = list(oraUp = foraUpForSe),
                         gs = gene_sets_list[[gs_ind]], perm = 0,
                         save2file = TRUE, out.dir = out_path)
    }
}
