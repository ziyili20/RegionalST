#' Automatically rank ROI centers based on entropy
#'
#' @param sce A single cell experiment object.
#' @param weight A data frame to specify the weights of all cell types.
#' @param enhanced A logical variable of whether using enhanced data.
#' @param selectN A total number for selected centers. Should be smaller than the total site number.
#' @param label A variable name that contains the cell type information.
#' @param topN A number to specify the total amount of top ranked ROIs.
#' @param min_radius The minimum repellent radius.
#' @param avern A number of the average sites used to compute unit distance, default is 5.
#' @param radius_vec A vector of numbers for candidate radiuses.
#' @param doPlot Logical variable about whether draw the plot.
#'
#' @return An sce object with selected ROI information.
#' @export
#' @examples
#' data("example_sce")
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#' weight <- data.frame(celltype = c("Cancer Epithelial", "CAFs", "T-cells", "Endothelial",
#'                                   "PVL", "Myeloid", "B-cells", "Normal Epithelial", "Plasmablasts"),
#'                      weight = c(0.25,0.05,
#'                                 0.25,0.05,
#'                                 0.025,0.05,
#'                                 0.25,0.05,0.025))
#' example_sce <- RankCenterByEntropy(example_sce, weight, label = "celltype",
#'                                   selectN = round(length(example_sce$spot)/10),
#'                                   topN = 3, min_radius = 10,
#'                                   radius_vec = c(10,15),
#'                                   doPlot = TRUE)
#'
RankCenterByEntropy <- function(sce,
                                weight,
                                enhanced = FALSE,
                                selectN = round(length(sce$spot)/10),
                                label = "celltype",
                                topN = 10,
                                min_radius = 10,
                                avern = 5,
                                radius_vec = c(10,15,20),
                                doPlot = TRUE) {
  
    stopifnot(exprs = {
      is.numeric(selectN)
      is.character(label)
      is.numeric(topN)
      is.numeric(min_radius)
      is.numeric(avern)
      is.numeric(radius_vec)
    })
  
    allID <- c()
    allEnt <- c()
    allGrp <- c()
    for(i in seq_len(length(radius_vec))) {
        radius <- radius_vec[i]
        oneradEnt <- GetOneRadiusEntropy(sce,
                                        selectN = selectN,
                                        enhanced = enhanced,
                                        weight = weight,
                                        label = label,
                                        radius = radius,
                                        doPlot = doPlot,
                                        mytitle = paste0("Radius:", radius))
        allID <- c(allID, oneradEnt$select_ID)
        allEnt <- c(allEnt, oneradEnt$select_ent)
        allGrp <- c(allGrp, rep(radius, length(oneradEnt$select_ent)))
    }
    allID2 <- allID[!duplicated(allID)]
    allEnt2 <- allEnt[!duplicated(allID)]
    allGrp2 <- allGrp[!duplicated(allID)]
    bestID <- FindRepelledID(sce,
                            allID = allID2,
                            allEnt = allEnt2,
                            allGrp = allGrp2,
                            enhanced = enhanced,
                            min_radius = min_radius,
                            avern = avern,
                            MaxN = topN)
    idx <- match(bestID, allID2)
    thisname <- paste0(seq_len(length(na.omit(bestID))),
                       paste0(":R(", allGrp2, ")",
                              "-Entp(", signif(allEnt2,4), ")")[idx])
    CellIDs <- allID2[idx]

    if(enhanced) {
        ThisCluster <- rep("Other", ncol(sce@assays@data$logcounts))
        tmpidx <- match(CellIDs, colnames(sce@assays@data$logcounts))
    } else {
        ThisCluster <- rep("Other", length(sce$spot))
        tmpidx <- match(CellIDs, sce$spot)
    }

    ThisCluster[na.omit(tmpidx)] <- thisname[!is.na(tmpidx)]
    thisname <- thisname[!duplicated(thisname)]
    sce$topRanked <- factor(ThisCluster, levels = c(thisname, "Other"))

    if (doPlot) {
        palette <- c("#0173b2", "#de8f05", "#029e73",
                     "#d55e00", "#cc78bc",
                     "#ca9161", "#fbafe4", "#949494",
                     "#ece133", "#56b4e9", "gray")
        BayesSpace::clusterPlot(sce, palette = palette,
                                label = "topRanked", size=0.1)
    }
    S4Vectors::metadata(sce)$selectCenters <- list(selectID = bestID,
                                       selectRadius = allGrp2[na.omit(idx)],
                                       selectEnt = allEnt2[na.omit(idx)],
                                       selectName = thisname)
    return(sce)
}


