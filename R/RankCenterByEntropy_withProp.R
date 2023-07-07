#' Automatically rank ROI centers based on entropy with proportions
#'
#' @param sce A single cell experiment object.
#' @param weight A data frame to specify the weights of all cell types.
#' @param selectN A total number for selected centers. Should be smaller than the total site number.
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
#' weight <- data.frame(celltype = c("Cancer Epithelial", "CAFs", "T-cells", "Endothelial",
#'                                 "PVL", "Myeloid", "B-cells", "Normal Epithelial", "Plasmablasts"),
#'                      weight = c(0.25,0.05,
#'                                 0.25,0.05,
#'                                 0.025,0.05,
#'                                 0.25,0.05,0.025))
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#' ## I set our min_raius as 10 and radius vector as 10 and 15 as the example dataset is very small
#' example_sce <- RankCenterByEntropy_withProp(example_sce, weight,
#'                                     selectN = round(length(example_sce$spot)/10),
#'                                     topN = 3, min_radius = 10,
#'                                     radius_vec = c(10,15),
#'                                     doPlot = TRUE)
#'
RankCenterByEntropy_withProp <- function(sce,
                                         weight,
                                         selectN = round(length(sce$spot)/10),
                                         topN = 10,
                                         min_radius = 10,
                                         avern = 5,
                                         radius_vec = c(10,15,20),
                                         doPlot = TRUE) {
    allID <- c()
    allEnt <- c()
    allGrp <- c()
    for(i in seq_len(length(radius_vec))) {
        radius <- radius_vec[i]
        oneradEnt <- GetOneRadiusEntropy_withProp(sce,
                                                 selectN = selectN,
                                                 weight = weight,
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
                            enhanced = FALSE,
                            min_radius = min_radius,
                            avern = avern,
                            MaxN = topN)
    idx <- match(bestID, allID2)
    thisname <- paste0(seq_len(length(stats::na.omit(bestID))), paste0(":R(", allGrp2, ")", "-Entp(", signif(allEnt2,4), ")")[idx])
    CellIDs <- allID2[idx]

    ThisCluster <- rep("Other", length(sce$spot))
    tmpidx <- match(CellIDs, sce$spot)

    ThisCluster[na.omit(tmpidx)] <- thisname[!is.na(tmpidx)]
    thisname <- thisname[!duplicated(thisname)]
    sce$topRanked <- factor(ThisCluster, levels = c(thisname, "Other"))

    if (doPlot) {
        palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
                     "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9", "gray")
        BayesSpace::clusterPlot(sce, palette = palette, label = "topRanked", size=0.1)
    }
    S4Vectors::metadata(sce)$selectCenters <- list(selectID = bestID,
                                       selectRadius = allGrp2[na.omit(idx)],
                                       selectEnt = allEnt2[na.omit(idx)],
                                       selectName = thisname)
    return(sce)
}
