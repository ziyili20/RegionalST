#' Draw regional cell type distribution with cellular proportion information
#'
#' @param sce A single cell experiment object.
#' @param label A string character for the cell type variable.
#' @param selCenter A vector of the interested ROIs, e.g., 1:4.
#'
#' @return A plot object.
#' @export
#' @examples
#' data("example_sce")
#' DrawRegionProportion_withProp(example_sce,
#'                              label = "Proportions",
#'                              selCenter = 1:3)
#'
#'
DrawRegionProportion_withProp <- function(sce,
                                          label = "CARD_CellType",
                                          selCenter = seq_len(10)) {
    stopifnot(exprs = {
      is.character(label)
    })
  
    region <- c()
    prop <- c()
    eval(parse(text = paste0(label, " = c()")))

    for(thisCenter in selCenter) {
        thisID <- S4Vectors::metadata(sce)$selectCenters$selectID[thisCenter]
        thisRadius <- S4Vectors::metadata(sce)$selectCenters$selectRadius[thisCenter]

        OneRegOut <- FindRegionalCells(sce,
                                      centerID = thisID,
                                      radius = thisRadius,
                                      doPlot = FALSE,
                                      returnPlot = FALSE)
        tmp <- sce@metadata$Proportions[match(OneRegOut$closeID, sce$spot), ]
        ttmp <- colSums(tmp)
        region <- c(region, rep(thisCenter, length(ttmp)))
        eval(parse(text = paste0(label, " <- c(", label, ", names(ttmp))")))
        prop <- c(prop, c(ttmp)/sum(ttmp)*100)
    }

    eval(parse(text = paste0("dataframe <- data.frame(region = region,
                           celltype = ", label, ", prop = prop)")))
    p2 <- c()
    eval(parse(text = paste0("p2 <- ggplot2::ggplot(dataframe,  ggplot2::aes(fill = ", label, ",
                          y = prop, x = factor(region)))+
        ggplot2::geom_bar(position = 'stack', stat = 'identity')+
        ggplot2::ggtitle('Regional distribution of the cell types')+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
        ggplot2::theme_minimal() + ggplot2::ylab('Proportions') + ggplot2::xlab('Region') +
        ggplot2::theme(panel.background = ggplot2::element_blank()) +
        ggplot2::scale_fill_manual(values = colorspace::rainbow_hcl(length(unique(", label, "))))")))

    return(p2)
}
