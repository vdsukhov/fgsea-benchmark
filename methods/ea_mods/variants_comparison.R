library(GSEABenchmarkeR)
library(data.table)
source("methods/utils/facetplot.R")

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
# Load malacards based rankings
# ----------------------------------------------------------------------------------------------------------------------------------
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)

mala.kegg.file <- file.path(data.dir, "malacards", "KEGG.rds")
mala.kegg <- readRDS(mala.kegg.file)

mala.gobp.file <- file.path(data.dir, "malacards", "GO_BP.rds")
mala.gobp <- readRDS(mala.gobp.file)

tcga.ids <- names(rnaseq_data)
names(tcga.ids) <- names(rnaseq_data)

# ----------------------------------------------------------------------------------------------------------------------------------
# Load ea rankings
# ----------------------------------------------------------------------------------------------------------------------------------

m_kegg_out_path <- "outdata/ea_mods/microarray/kegg/"
m_gobp_out_path <- "outdata/ea_mods/microarray/gobp/"
r_kegg_out_path <- "outdata/ea_mods/rnaseq/kegg/"
r_gobp_out_path <- "outdata/ea_mods/rnaseq/gobp/"


methods <- c("ora", "oraUp")

m_kegg_ranks <- readResults(m_kegg_out_path, names(microarray_data),
                            methods = methods, type = "ranking")
m_gobp_ranks <- readResults(m_gobp_out_path, names(microarray_data),
                            methods = methods, type = "ranking")

r_kegg_ranks <- readResults(r_kegg_out_path, names(rnaseq_data),
                            methods=methods, type="ranking")
r_gobp_ranks <- readResults(r_gobp_out_path, names(rnaseq_data),
                            methods=methods, type="ranking")

# --------------------------------------------------------------------------------------------------------------------------------
# Original Relevance Score Plot
# --------------------------------------------------------------------------------------------------------------------------------
m_kegg_score <- evalRelevance(m_kegg_ranks, mala.kegg, d2d.map[names(microarray_data)])
m_gobp_score <- evalRelevance(m_gobp_ranks, mala.gobp, d2d.map[names(microarray_data)])
r_kegg_score <- evalRelevance(r_kegg_ranks, mala.kegg, tcga.ids)
r_gobp_score <- evalRelevance(r_gobp_ranks, mala.gobp, tcga.ids)


facetplot(m_kegg_score, m_gobp_score, r_kegg_score, r_gobp_score, ylab="% optimal relevance score") +
    ylim(0, 100) +
    labs(title = "Benchmark results",
         subtitle = "Optimal score")

# --------------------------------------------------------------------------------------------------------------------------------
# sig sets
# --------------------------------------------------------------------------------------------------------------------------------

m_kegg_sigsets <- evalNrSigSets(m_kegg_ranks, alpha = 0.05, padj = "BH")
m_gobp_sigsets <- evalNrSigSets(m_gobp_ranks, alpha = 0.05, padj = "BH")
r_kegg_sigsets <- evalNrSigSets(r_kegg_ranks, alpha = 0.05, padj = "BH")
r_gobp_sigsets <- evalNrSigSets(r_gobp_ranks, alpha = 0.05, padj = "BH")


facetplot(m_kegg_sigsets, m_gobp_sigsets, r_kegg_sigsets, r_gobp_sigsets, ylab="% sig sets") +
    labs(title = "Number of significant sets") +
    # ylim(0, 100) +
    NULL


