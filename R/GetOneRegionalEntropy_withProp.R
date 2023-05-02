GetOneRegionalEntropy_withProp <- function(sce,
                                           regionOut,
                                           weight = NULL) {
    propmat <- sce@metadata$Proportions

    regionCT <- propmat[rownames(propmat) %in% regionOut$closeID, ]
    CTdist <- colSums(regionCT)
    proptab <- prop.table(CTdist)
    if (is.null(weight)) {
        thisent <- GetEntropy(proptab)
    } else {
        matchedweight <- weight$weight[stats::na.omit(match(names(proptab), weight$celltype))]
        thisent <- GetEntropy(proptab, weight = matchedweight)
    }
    return(thisent)
}
