FindUnitDist <- function(sce, avern = 5) {
    stopifnot(exprs = {
      is.numeric(avern)
    })
  
    xloc <- sce$imagecol
    yloc <- sce$imagerow

    ID <- sample(seq_len(length(xloc)), size = avern)
    rec_minDist <- rep(NA, avern)
    for(i in seq_len(avern)) {
        rec_minDist[i] <- min(sqrt((xloc[ID[i]] - xloc[-ID[i]])^2 + (yloc[ID[i]] - yloc[-ID[i]])^2))
    }
    return(mean(rec_minDist))
}
