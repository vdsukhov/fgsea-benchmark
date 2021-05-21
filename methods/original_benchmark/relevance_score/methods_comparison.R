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

path_to_geo2kegg <- "inpdata/GSEABenchmarking/GEO2KEGG_preproc/"
path_to_tcga <- "inpdata/GSEABenchmarking/TCGA_preproc/GSE62944_matched_limmavoom/"

microarray_data <- lapply(list.files(path_to_geo2kegg, pattern = ".rds"), function(fnm){
    se <- readRDS(paste0(path_to_geo2kegg, fnm))
    return(se)
})
names(microarray_data) <- gsub(".rds", "", list.files(path_to_geo2kegg, pattern = ".rds"))

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

m_kegg_out_path <- "outdata/original_benchmark/microarray/kegg/"
m_gobp_out_path <- "outdata/original_benchmark/microarray/gobp/"
r_kegg_out_path <- "outdata/original_benchmark/rnaseq/kegg/"
r_gobp_out_path <- "outdata/original_benchmark/rnaseq/gobp/"


methods <- c("cameraPr", "fgsea", "ora", "camera")

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
# Original Relevance Score Plot for gene sets of different sizes
# --------------------------------------------------------------------------------------------------------------------------------
gs_cutoff <- 100

kegg_gs_small <- names(kegg_gs[lengths(kegg_gs) <= gs_cutoff])
kegg_gs_large <- names(kegg_gs[lengths(kegg_gs) > gs_cutoff])


gobp_gs_small <- names(gobp_gs[lengths(gobp_gs) <= gs_cutoff])
gobp_gs_large <- names(gobp_gs[lengths(gobp_gs) > gs_cutoff])


# results for small gene sets
m_kegg_ranks_small <- lapply(m_kegg_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% kegg_gs_small,])
    })
})
m_gobp_ranks_small <- lapply(m_gobp_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% gobp_gs_small,])
    })
})

r_kegg_ranks_small <- lapply(r_kegg_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% kegg_gs_small,])
    })
})
r_gobp_ranks_small <- lapply(r_gobp_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% gobp_gs_small,])
    })
})

# results for large gene sets
m_kegg_ranks_large <- lapply(m_kegg_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% kegg_gs_large,])
    })
})
m_gobp_ranks_large <- lapply(m_gobp_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% gobp_gs_large,])
    })
})

r_kegg_ranks_large <- lapply(r_kegg_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% kegg_gs_large,])
    })
})
r_gobp_ranks_large <- lapply(r_gobp_ranks, function(method_res){
    res <- lapply(method_res, function(df){
        return(df[df$GENE.SET %in% gobp_gs_large,])
    })
})


m_kegg_score_small <- evalRelevance(m_kegg_ranks_small, mala.kegg, d2d.map[names(microarray_data)])
m_gobp_score_small <- evalRelevance(m_gobp_ranks_small, mala.gobp, d2d.map[names(microarray_data)])
r_kegg_score_small <- evalRelevance(r_kegg_ranks_small, mala.kegg, tcga.ids)
r_gobp_score_small <- evalRelevance(r_gobp_ranks_small, mala.gobp, tcga.ids)


m_kegg_score_large <- evalRelevance(m_kegg_ranks_large, mala.kegg, d2d.map[names(microarray_data)])
m_gobp_score_large <- evalRelevance(m_gobp_ranks_large, mala.gobp, d2d.map[names(microarray_data)])
r_kegg_score_large <- evalRelevance(r_kegg_ranks_large, mala.kegg, tcga.ids)
r_gobp_score_large <- evalRelevance(r_gobp_ranks_large, mala.gobp, tcga.ids)


facetplot(m_kegg_score_small, m_gobp_score_small,
          r_kegg_score_small, r_gobp_score_small, ylab="% optimal relevance score") +
    ylim(0, 100) +
    labs(title = "Benchmark results",
         subtitle = "Optimal score (Small Genesets)")


facetplot(m_kegg_score_large, m_gobp_score_large,
          r_kegg_score_large, r_gobp_score_large, ylab="% optimal relevance score") +
    ylim(0, 100) +
    labs(title = "Benchmark results",
         subtitle = "Optimal score (Large Genesets)")
