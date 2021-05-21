get_sig_fraction <- function(ea_res, trshld = 0.05){
    ranking <- ea_res[[1]][[1]]$ranking
    sig_fraction <- sum(ranking$PVAL < trshld, na.rm = TRUE) / nrow(ranking)
    return(sig_fraction)
}
