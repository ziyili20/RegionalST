#' Perform Preprocessing for spatial data (tailored from BayesSpace function)
#'
#' @param sce A SingleCellExperiment object.
#' @param platform Which platform the data are from, Visium or ST.
#' @param n.PCs Number of PCs used in the analysis.
#' @param n.HVGs Number of highly variable genes used in the analysis.
#' @param skip.PCA A boolean variable to choose whether skipping the PCA step or not.
#' @param assay.type Which assay to use, default is logcounts.
#'
#' @importFrom scater logNormCounts
#' @importFrom SingleCellExperiment logcounts
#' @importFrom scater runPCA
#' @importFrom stats sd
#'
#' @return A processed SingleCellExperiment object.
#'
#' @export
#' @examples
#' data(example_sce)
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#'
mySpatialPreprocess <- function (sce, 
                                 platform = c("Visium", "ST"), 
                                 n.PCs = 15, 
                                 n.HVGs = 2000,
                                 skip.PCA = FALSE, 
                                 assay.type = "logcounts")
{
    stopifnot(exprs = {
      platform %in% c("Visium", "ST")
      is.numeric(n.PCs)
      is.numeric(n.HVGs)
    })
  
    S4Vectors::metadata(sce)$BayesSpace.data <- list()
    S4Vectors::metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
    S4Vectors::metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    if (!skip.PCA) {
        sce <- scater::logNormCounts(sce)
        lmat <- SingleCellExperiment::logcounts(sce)
        lmat.sd <- apply(lmat, 1, sd)
        maxn <- min(n.HVGs, nrow(lmat))
        topGenes <- names(lmat.sd[order(lmat.sd, decreasing = TRUE)])[1:maxn]
        sce <- scater::runPCA(sce, subset_row = topGenes, ncomponents = n.PCs,
                      exprs_values = assay.type)
        SummarizedExperiment::rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% topGenes)
    }
    sce
}
