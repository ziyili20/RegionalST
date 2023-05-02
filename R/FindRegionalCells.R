#' Identify regional cells given centers and radiuses
#'
#' @param sce A single cell experiment object.
#' @param centerID One or a vector of spot IDs as centers of ROIs.
#' @param enhanced A logical variable for plotting enhanced plot or now. Default is FALSE.
#' @param radius A number of fixed ROI radius.
#' @param avern A number of the average sites used to compute unit distance, default is 5.
#' @param doPlot A logical variable to specify whether plot the figure or not.
#' @param returnPlot a logical variable to specify whether output the plot or not.
#'
#' @return A list including center spot ID and regional spot IDs.
#' @export
#' @examples
#' # FindRegionalCells(sce, centerID = "ACGCCTGACACGCGCT-1")
#'
FindRegionalCells <- function(sce,
                              centerID,
                              enhanced = FALSE,
                              radius = 10,
                              avern = 5,
                              doPlot = FALSE,
                              returnPlot = FALSE) {

    unitDist <- FindUnitDist(sce, avern = avern)

    xloc <- sce$imagecol
    yloc <- sce$imagerow

    if(enhanced) {
        sce.enhanced <- sce
        allspot <- colnames(sce.enhanced@assays@data$logcounts)
        center_x <- xloc[which(allspot == centerID)]
        center_y <- yloc[which(allspot == centerID)]

        alldist <- sqrt((xloc - center_x)^2 + (yloc - center_y)^2)
        closeID <- allspot[which(alldist <= (radius * unitDist+0.1))]

        featureDot <- rep(3, length(allspot))
        featureDot[allspot %in% closeID] <- 2
        featureDot[allspot == centerID] <- (-2)
    } else {
        center_x <- xloc[which(sce$spot == centerID)]
        center_y <- yloc[which(sce$spot == centerID)]

        alldist <- sqrt((xloc - center_x)^2 + (yloc - center_y)^2)
        closeID <- sce$spot[which(alldist <= (radius * unitDist+0.1))]

        featureDot <- rep(3, length(sce$spot))
        featureDot[sce$spot %in% closeID] <- 2
        featureDot[sce$spot == centerID] <- -2
    }

    sce$featureDot <- featureDot

    if (doPlot) {
        p1 <- BayesSpace::featurePlot(sce, feature = featureDot,
                                      low = "red", mid = "yellow",
                                      high = "gray") +
            ggplot2::theme(legend.position = "none")
        print(p1)
    }

    if (returnPlot) {
        p1 <- BayesSpace::featurePlot(sce, feature = featureDot,
                                     low = "red", mid = "yellow",
                                     high = "gray") +
            ggplot2::theme(legend.position = "none")
        return(p1)
    } else {
        return(list(centerID = centerID,
                    closeID = closeID))
    }

}
