library(ggplot2)
library(ggpubr)
source("methods/utils/gseabenchmarker_colors.R")

facetplot <- function(ma.kegg, ma.go, rseq.kegg, rseq.go,
                      ylab="% significant sets", vline=6.5, hline=NA, log=FALSE, orderf=median)
{
    l <- list(ma.kegg=ma.kegg, ma.go=ma.go, rseq.kegg=rseq.kegg, rseq.go=rseq.go)
    df <- reshape2::melt(l)
    gsc <- vapply(df$L1, function(x) unlist(strsplit(x,"\\."))[2],
                  character(1), USE.NAMES=FALSE)
    df <- cbind(df, gsc=gsc)
    df$gsc <- toupper(df$gsc)
    df$gsc <- vapply(df$gsc, function(n)
        ifelse(n == "GO", paste(n, "BP", sep="-"), n),
        character(1), USE.NAMES=FALSE)
    df$gsc <- factor(df$gsc, levels=c("KEGG", "GO-BP"))
    colnames(df)[1:2] <- c("dataset", "method")
    colnames(df)[4] <- "compendium"
    df$compendium <- sub("ma.kegg", "GEO2KEGG microarray", df$compendium)
    df$compendium <- sub("rseq.go", "TCGA RNA-seq", df$compendium)
    df$compendium <- sub("rseq.kegg", "TCGA RNA-seq", df$compendium)
    df$compendium <- sub("ma.go", "GEO2KEGG microarray", df$compendium)
    df$method <- substring(df$method, 1, 7)
    if(log) df$value <- log(df$value, base=10)
    o <- sort(vapply(split(df$value, df$method),
                     orderf, numeric(1), na.rm=TRUE))
    df$method <- factor(df$method, levels=names(o))
    p <- ggboxplot(df, x = "method", y = "value",
                   width = 0.8, ylab=ylab, xlab="", fill="method")
    p <- ggpar(p, x.text.angle=45, palette = "simpsons", legend="none")
    if(!is.na(vline))
        p <- p + geom_vline(xintercept=vline, linetype="dashed", color = cb.darkgrey)
    if(!is.na(hline))
        p <- p + geom_hline(yintercept=hline, linetype="dashed", color = cb.red)

    facet(p, facet.by=c("compendium", "gsc"))
}
