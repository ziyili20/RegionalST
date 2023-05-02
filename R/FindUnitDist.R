FindUnitDist <- function(sce, avern = 5) {
    xloc <- sce$imagecol
    yloc <- sce$imagerow

    ID <- sample(1:length(xloc), size = avern)
    rec_minDist <- rep(NA, avern)
    for(i in seq_len(avern)) {
        rec_minDist[i] <- min(sqrt((xloc[ID[i]] - xloc[-ID[i]])^2 + (yloc[ID[i]] - yloc[-ID[i]])^2))
    }
    return(mean(rec_minDist))
}
