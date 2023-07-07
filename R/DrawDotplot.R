#' Draw dot plot for GSEA results of cross-regional DE genes
#'
#' @param allCTres A list of GSEA results for all cell types.
#' @param CT A number of the interested cell type, e.g., 1, 2, 3.
#' @param angle A number of plotting parameter, angle of the x axis label.
#' @param vjust A number of vertical adjustment in plotting.
#' @param hjust A number of horizontal adjustment in plotting.
#' @param padj_cutoff A cutoff number of adjusted p value.
#' @param topN A number of the plotted top pathways.
#' @param chooseP A character string for the p value that used in plotting, e.g., "padj" or "pval".
#' @param eachN The maximum number of pathways in each cell type.
#'
#' @importFrom ggplot2 aes element_blank element_text geom_point labs guides
#' @importFrom grDevices rgb
#' @importFrom SummarizedExperiment colData
#'
#' @return A plot object
#' @export
#' @examples
#' data(exampleRes)
#' allCTres <- DoGSEA(exampleRes, whichDB = "hallmark", withProp = TRUE)
#' DrawDotplot(allCTres, CT = 1, angle = 15, vjust = 1, chooseP = "padj")
#'
DrawDotplot <- function(allCTres,
                        CT = 1,
                        angle = 20,
                        vjust = 0.9,
                        hjust = 1,
                        padj_cutoff = 1,
                        topN = 20,
                        chooseP = "padj",
                        eachN = NULL) {

    thisres <- allCTres[[CT]]
    thisres <- thisres[thisres$padj <= padj_cutoff, ]
    if (is.null(eachN)) {
        if (nrow(thisres) > topN) {
            rankedres <- thisres[order(thisres$pval, decreasing = FALSE),]
            thisres <- rankedres[seq_len(topN), ]
        }
    } else {
        chooseres <- c()
        uniqCT <- unique(thisres$group)
        for(q in seq_len(uniqCT)) {
            oneres <- thisres[thisres$group == uniqCT[q],]
            if (q == 1) {
                if(chooseP == "padj") {
                    chooseres <- oneres[order(oneres$padj, decreasing = FALSE)[seq_len(min(topN, nrow(oneres)))], ]
                } else if (chooseP == "pval") {
                    chooseres <- oneres[order(oneres$pval, decreasing = FALSE)[seq_len(min(topN, nrow(oneres)))], ]
                } else {
                    stop("chooseP need to be one of the following: padj, pval!")
                }
            } else {
                if(chooseP == "padj") {
                    chooseres <- rbind(chooseres,
                                       oneres[order(oneres$padj, decreasing = FALSE)[seq_len(min(topN, nrow(oneres)))], ])
                } else if (chooseP == "pval") {
                    chooseres <- rbind(chooseres,
                                       oneres[order(oneres$pval, decreasing = FALSE)[seq_len(min(topN, nrow(oneres)))], ])
                } else {
                    stop("chooseP need to be one of the following: padj, pval!")
                }
            }
        }
        thisres <- chooseres[!is.na(chooseres$pathway),]
    }

    if(chooseP == "padj") {
        S1 <- ggplot2::ggplot(thisres, aes(x= group, y=pathway, size=-log(padj, base = 10),
                                          color=NES, group=group)) + ggplot2::geom_point(alpha = 0.8) +
            ggplot2::theme_classic()+ ggplot2::labs(size='-log10(padj)') +
            ggplot2::theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                           axis.text = element_text(face="bold"),
                           title = element_text(face="bold"),
                           axis.text.x = element_text(angle = angle, vjust = vjust, hjust = hjust))+
            ggplot2::scale_color_gradient(high = "red2",  low = "mediumblue", space = "Lab")+
            ggplot2::ggtitle(names(allCTres[CT]))
    } else {
        S1 <- ggplot2::ggplot(thisres, aes(x= group, y=pathway, size=-log(pval, base = 10),
                                          color=NES, group=group)) + geom_point(alpha = 0.8) +
            ggplot2::theme_classic()+ labs(size='-log10(pval)') +
            ggplot2::theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                           axis.text = element_text(face="bold"),
                           title = element_text(face="bold"),
                           axis.text.x = element_text(angle = angle, vjust = vjust, hjust = hjust))+
            ggplot2::scale_color_gradient(high = "red2",  low = "mediumblue", space = "Lab")+
            ggplot2::ggtitle(names(allCTres[CT]))
    }

    return(S1)
}
