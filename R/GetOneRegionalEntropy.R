GetOneRegionalEntropy <- function(sce,
                                  regionOut,
                                  enhanced = FALSE,
                                  weight = NULL,
                                  label = "celltype") {
  
    stopifnot(exprs = {
      is.character(label)
    })
  
    alllabel <- eval(parse(text = paste0("sce$", label)))

    if(enhanced) {
        regionCT <- alllabel[colnames(sce@assays@data$logcounts) %in% regionOut$closeID]
    } else {
        regionCT <- alllabel[sce$spot %in% regionOut$closeID]
    }
    proptab <- prop.table(table(regionCT))
    if (is.null(weight)) {
        thisent <- GetEntropy(proptab)
    } else {
        matchedweight <- weight$weight[stats::na.omit(match(names(proptab), weight$celltype))]
        thisent <- GetEntropy(proptab, weight = matchedweight)
    }
    return(thisent)
}
