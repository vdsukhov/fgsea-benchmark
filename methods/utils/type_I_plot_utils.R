readTyI <- function(res.dir)
{
    res.files <- file.path(res.dir, paste0(c(sbeaMethods()[1:10],"camera_igcNA"), ".rds"))
    res.files <- res.files[file.exists(res.files)]
    res <- sapply(res.files, readRDS)
    colnames(res) <- basename(colnames(res))
    colnames(res) <- sub(".rds$", "", colnames(res))
    return(res)
}


plotTypeIError2 <- function(data, ylabel=0.4)
{
    data <- t(data)
    rownames(data) <- substring(rownames(data), 1, 7)
    data <- data[order(data[,"Max."] - data[,"Mean"]),]
    df <- data.frame(methods=rownames(data), data)
    colnames(df)[2:7] <- c("y0", "y25", "y50", "mean", "y75", "y100")
    df[,1] <- factor(df[,1], levels=df[,1])
    df.points <- data.frame(x=1:nrow(data), y=df$mean)
    p <- ggplot() +
        geom_boxplot(data=df, width = 0.8,
                     aes(x=df[,1], ymin=y0, lower=y25, middle=y50, upper=y75, ymax=y100),
                     stat="identity", fill="grey92") + theme_pubr()
    #get_palette(palette = "simpsons", 10)
    p <- ggpar(p, x.text.angle=45, legend="none") + xlab("") + ylab("type I error rate") +
        geom_point(data=df.points, aes(x=x, y=y), color=cb.blue) +
        geom_hline(yintercept=0.05, linetype="dashed", color = cb.red)
    return(p)
}
