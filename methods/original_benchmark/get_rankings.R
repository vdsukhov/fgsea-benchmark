library(GSEABenchmarkeR)
library(data.table)
library(EnrichmentBrowser)
library(BiocParallel)
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

# ----------------------------------------------------------------------------------------------------------------------------------
# Check folders for saving output data
# ----------------------------------------------------------------------------------------------------------------------------------
if (!dir.exists("outdata/original_benchmark/")){
    dir.create("outdata/original_benchmark/microarray/kegg/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/original_benchmark/microarray/gobp/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/original_benchmark/rnaseq/kegg/", recursive = TRUE, showWarnings = FALSE)
    dir.create("outdata/original_benchmark/rnaseq/gobp/", recursive = TRUE, showWarnings = FALSE)
}

# ----------------------------------------------------------------------------------------------------------------------------------
# EA methods for microarray data
# ----------------------------------------------------------------------------------------------------------------------------------
m_kegg_out_path <- "outdata/original_benchmark/microarray/kegg/"
m_gobp_out_path <- "outdata/original_benchmark/microarray/gobp/"

register(SerialParam())
# KEGG
set.seed(1)
camera_res <- runEA(microarray_data, methods = "camera", gs = kegg_gs, perm = 0,
                    save2file = TRUE, out.dir = m_kegg_out_path)

set.seed(1)
cameraPr_res <- runEA(microarray_data, methods = list(cameraPr = cameraPrForSe), gs = kegg_gs, perm = 0,
                    save2file = TRUE, out.dir = m_kegg_out_path)

set.seed(1)
ora_res <- runEA(microarray_data, methods = "ora", gs = kegg_gs, perm = 0,
                 save2file = TRUE, out.dir = m_kegg_out_path)

set.seed(1)
fgsea_res <- runEA(microarray_data, methods = list(fgsea = fgseaForSe), gs = kegg_gs,
                   save2file = TRUE, out.dir = m_kegg_out_path)

# GOBP
set.seed(1)
camera_res <- runEA(microarray_data, methods = "camera", gs = gobp_gs, perm = 0,
                    save2file = TRUE, out.dir = m_gobp_out_path)

set.seed(1)
cameraPr_res <- runEA(microarray_data, methods = list(cameraPr = cameraPrForSe), gs = gobp_gs, perm = 0,
                      save2file = TRUE, out.dir = m_gobp_out_path)

set.seed(1)
ora_res <- runEA(microarray_data, methods = "ora", gs = gobp_gs, perm = 0,
                 save2file = TRUE, out.dir = m_gobp_out_path)

set.seed(1)
fgsea_res <- runEA(microarray_data, methods = list(fgsea = fgseaForSe), gs = gobp_gs,
                   save2file = TRUE, out.dir = m_gobp_out_path)


# ----------------------------------------------------------------------------------------------------------------------------------
# EA methods for rnaseq data
# ----------------------------------------------------------------------------------------------------------------------------------
r_kegg_out_path <- "outdata/original_benchmark/rnaseq/kegg/"
r_gobp_out_path <- "outdata/original_benchmark/rnaseq/gobp/"

# register(SerialParam())
# KEGG
set.seed(1)
camera_res <- runEA(rnaseq_data, methods = "camera", gs = kegg_gs, perm = 0,
                    save2file = TRUE, out.dir = r_kegg_out_path)


set.seed(1)
cameraPr_res <- runEA(rnaseq_data, methods = list(cameraPr = cameraPrForSe), gs = kegg_gs, perm = 0,
                      save2file = TRUE, out.dir = r_kegg_out_path)

set.seed(1)
ora_res <- runEA(rnaseq_data, methods = "ora", gs = kegg_gs, perm = 0,
                 save2file = TRUE, out.dir = r_kegg_out_path)

set.seed(1)
fgsea_res <- runEA(rnaseq_data, methods = list(fgsea = fgseaForSe), gs = kegg_gs,
                   save2file = TRUE, out.dir = r_kegg_out_path)

# GOBP
set.seed(1)
camera_res <- runEA(rnaseq_data, methods = "camera", gs = gobp_gs, perm = 0,
                    save2file = TRUE, out.dir = r_gobp_out_path)

register(SerialParam())
set.seed(1)
cameraPr_res <- runEA(rnaseq_data, methods = list(cameraPr = cameraPrForSe), gs = gobp_gs, perm = 0,
                      save2file = TRUE, out.dir = r_gobp_out_path)

set.seed(1)
ora_res <- runEA(rnaseq_data, methods = "ora", gs = gobp_gs, perm = 0,
                 save2file = TRUE, out.dir = r_gobp_out_path)

set.seed(1)
fgsea_res <- runEA(rnaseq_data, methods = list(fgsea = fgseaForSe), gs = gobp_gs,
                   save2file = TRUE, out.dir = r_gobp_out_path)
