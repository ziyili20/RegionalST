FindRepelledID <- function(sce,
                           allID,
                           allEnt,
                           allGrp,
                           enhanced = FALSE,
                           min_radius = 5,
                           avern = 5,
                           MaxN = 10) {
    unitDist <- FindUnitDist(sce, avern = avern)

    idx <- order(allEnt, decreasing = TRUE)
    ranked_ID <- allID[idx]
    if (enhanced) {
        xvec <- sce$imagecol[match(ranked_ID, colnames(sce@assays@data$logcounts))]
        yvec <- sce$imagerow[match(ranked_ID, colnames(sce@assays@data$logcounts))]
    } else {
        xvec <- sce$imagecol[match(ranked_ID, sce$spot)]
        yvec <- sce$imagerow[match(ranked_ID, sce$spot)]
    }


    bestIdx <- rep(NA, MaxN)

    cumN <- 0
    while(cumN < MaxN & length(ranked_ID)>0) {
        cumN <- cumN + 1

        if(cumN == 1) {
            ranked_ID_new <- ranked_ID
            xvec_new <- xvec
            yvec_new <- yvec
        }

        bestIdx[cumN] <- ranked_ID_new[1]

        removeID <- which(sqrt((xvec_new[1] - xvec_new)^2 + (yvec_new[1] - yvec_new)^2)< min_radius*unitDist)
        ranked_ID_new <- ranked_ID_new[-removeID]
        xvec_new <- xvec_new[-removeID]
        yvec_new <- yvec_new[-removeID]
    }
    return(bestIdx)
}
