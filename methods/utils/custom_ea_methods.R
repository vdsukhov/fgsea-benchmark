library(limma)
library(fgsea)

fgseaForSe <- function(se, gs){
    ranks <- setNames(rowData(se)[, "limma.STAT"], rownames(se))
    set.seed(1)
    fgseaRes <- fgsea(gs, ranks, eps=0.0, minSize=5, maxSize=500, BPPARAM = SerialParam())
    fgseaRes <- na.omit(fgseaRes)
    pvalRes <- setNames(fgseaRes$pval, fgseaRes$pathway)
    return(pvalRes)
}

fgseaAbsForSe <- function(se, gs){
    ranks <- setNames(rowData(se)[, "limma.STAT"], rownames(se))
    set.seed(1)
    fgseaRes <- fgsea(gs, abs(ranks), eps=0.0, minSize=5, maxSize=500,
                      scoreType = "pos", BPPARAM = SerialParam())
    fgseaRes <- na.omit(fgseaRes)
    pvalRes <- setNames(fgseaRes$pval, fgseaRes$pathway)
    return(pvalRes)
}


foraForSe <- function(se, gs, padj_trshld = 0.05){
    genes <- rownames(se)[which(rowData(se)[, "ADJ.PVAL"] < padj_trshld)]
    universe <- rownames(se)
    foraRes <- fora(gs, genes, universe, minSize=5, maxSize=500)
    pvalRes <- setNames(foraRes$pval, foraRes$pathway)
    return(pvalRes)
}

foraUpForSe <- function(se, gs, padj_trshld = 0.05){
    genes <- rownames(se)[which(rowData(se)[, "ADJ.PVAL"] < padj_trshld & rowData(se)[, "limma.STAT"] >= 0)]
    universe <- rownames(se)
    foraRes <- fora(gs, genes, universe, minSize=5, maxSize=500)
    pvalRes <- setNames(foraRes$pval, foraRes$pathway)
    return(pvalRes)
}


cameraPrForSe <- function(se, gs){
    ranks <- setNames(rowData(se)[, "limma.STAT"], rownames(se))
    universe <- rownames(se)
    gss <- unlist(lapply(gs, function(x) sum(x %in% universe)))
    gs_to_keep <- (gss >= 5) & (gss <= 500)
    set.seed(1)
    cameraRes <- limma::cameraPR(ranks, gs[gs_to_keep])
    pvalRes <- setNames(cameraRes$PValue, nm = rownames(cameraRes))
    return(pvalRes)
}
