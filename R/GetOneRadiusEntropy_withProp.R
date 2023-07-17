#' Computer the entropy for a fixed radius with cell type proportion
#'
#' @param sce A single cell experiment object.
#' @param selectN A total number for selected centers. Should be smaller than the total site number.
#' @param weight A data frame to specify the weights of all cell types.
#' @param label A variable name that contains the cell type information.
#' @param radius A number for fixed radius.
#' @param doPlot Logical variable about whether draw the plot.
#' @param mytitle A character string for the title of the plot.
#'
#' @return A list including the selected centers, computed entropies, radius.
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
#' GetOneRadiusEntropy_withProp(example_sce, selectN = round(length(example_sce$spot)/10),
#'                              weight = weight,
#'                              radius = 5,
#'                              doPlot = TRUE,
#'                              mytitle = "Radius 5 weighted entropy")
#'
GetOneRadiusEntropy_withProp <- function(sce,
                                         selectN,
                                         weight = NULL,
                                         label = "celltype",
                                         radius = 10,
                                         doPlot = FALSE,
                                         mytitle = NULL) {
    allspot <- sce$spot
    sel_center <- sample(allspot, selectN)

    thisRadiusEnt <- rep(NA, selectN)

    message(paste0("Processing to Radius = ", radius, "\n"))

    if (selectN > length(allspot)) {
        stop("selectN must be equal to or smaller than the number of the spots!")
    }

    # Initializes the progress bar
    pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = selectN, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 60,   # Progress bar width. Defaults to getOption("width")
                         char = "=")
    for(i in seq_len(selectN)) {
        utils::setTxtProgressBar(pb, i)

        regionOut <- FindRegionalCells(sce,
                                       centerID = sel_center[i],
                                       radius = radius,
                                       enhanced = FALSE,
                                       avern = 5,
                                       doPlot =  FALSE)
        thisRadiusEnt[i] <- GetOneRegionalEntropy_withProp(sce,
                                                           regionOut = regionOut,
                                                           weight = weight)
    }
    close(pb)

    eval(parse(text = paste0("Entro_Radius", radius, " <- rep(-", max(thisRadiusEnt)+0.3, ", length(allspot))")))
    eval(parse(text = paste0("Entro_Radius", radius, "[match(sel_center, allspot)] <- thisRadiusEnt")))
    eval(parse(text = paste0("sce$Entro_Radius", radius, " <- Entro_Radius", radius)))

    p1 <- c()
    if(doPlot) {
        eval(parse(text = paste0("p1 <- BayesSpace::featurePlot(sce, feature = Entro_Radius",
                                 radius,
                                 ", low = 'gray', mid = 'yellow', high = 'red') + ggplot2::ggtitle('", mytitle, "')")))
        print(p1)
    }
    return(list(select_ID = sel_center,
                select_ent = thisRadiusEnt,
                radius = radius))

}
