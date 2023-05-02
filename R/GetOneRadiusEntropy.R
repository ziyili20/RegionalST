#' Computer the entropy for a fixed radius
#'
#' @param sce A single cell experiment object.
#' @param selectN A total number for selected centers. Should be smaller than the total site number.
#' @param enhanced A logical variable of whether using enhanced data.
#' @param weight A data frame to specify the weights of all cell types.
#' @param label A variable name that contains the cell type information.
#' @param radius A number for fixed radius.
#' @param doPlot Logical variable about whether draw the plot.
#' @param mytitle A character string for the title of the plot.
#'
#' @importFrom ggplot2 ggtitle
#'
#' @return A list including the selected centers, computed entropies, radius.
#' @export
#' @examples
#' data("example_sce")
#' weight <- data.frame(celltype = c("Cancer Epithelial", "CAFs",
#'                                    "T-cells", "Endothelial",
#'                                    "PVL", "Myeloid", "B-cells",
#'                                    "Normal Epithelial", "Plasmablasts"),
#'                      weight = c(0.25,0.05,
#'                                 0.25,0.05,
#'                                 0.025,0.05,
#'                                 0.25,0.05,0.025))
#' example_sce <- BayesSpace::spatialPreprocess(example_sce, platform="Visium", log.normalize=TRUE)
#' GetOneRadiusEntropy(example_sce, selectN = round(length(example_sce$spot)/2),
#'                     weight = weight, radius = 5, doPlot = TRUE,
#'                     mytitle = "Radius 5 weighted entropy")
#'
GetOneRadiusEntropy <- function(sce,
                                selectN,
                                enhanced = FALSE,
                                weight = NULL,
                                label = "celltype",
                                radius = 10,
                                doPlot = FALSE,
                                mytitle = NULL) {

    if(enhanced) {
        allspot <- colnames(sce@assays@data$logcounts)
    } else {
        allspot <- sce$spot
    }
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
                                       enhanced = enhanced,
                                       avern = 5,
                                       doPlot =  FALSE)
        thisRadiusEnt[i] <- GetOneRegionalEntropy(sce,
                                                  regionOut = regionOut,
                                                  enhanced = enhanced,
                                                  weight = weight,
                                                  label = label)
    }
    close(pb)

    eval(parse(text = paste0("Entro_Radius", radius, " <- rep(-", max(thisRadiusEnt)+0.3, ", length(allspot))")))
    eval(parse(text = paste0("Entro_Radius", radius, "[match(sel_center, allspot)] <- thisRadiusEnt")))
    eval(parse(text = paste0("sce$Entro_Radius", radius, " <- Entro_Radius", radius)))

    if(doPlot) {
        eval(parse(text = paste0("p1 <- BayesSpace::featurePlot(sce, feature = Entro_Radius",
                                 radius,
                                 ", low = 'gray', mid = 'yellow', high = 'red') + ggtitle('", mytitle, "')")))
        print(p1)
    }
    return(list(select_ID = sel_center,
                select_ent = thisRadiusEnt,
                radius = radius))
}
