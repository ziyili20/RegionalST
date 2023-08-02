#' Perform GSEA analysis for cross-regional DE genes
#'
#' @param considerRes A list of cross-regional DE genes.
#' @param whichDB A character string to select the database names, e.g., "hallmark", "kegg", "reactome".
#' @param gmtdir Directory for external database gmt file location.
#' @param withProp Whether deconvolution proportion is used in previous steps.
#'
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @importFrom utils data
#'
#' @return A list including GSEA results for all cell types.
#'
#' @export
#' @examples
#' data(exampleRes)
#' allCTres <- DoGSEA(exampleRes, whichDB = "hallmark", withProp = TRUE)
#'
DoGSEA <- function(considerRes,
                   whichDB = "hallmark",
                   gmtdir = NULL,
                   withProp = FALSE) {

    stopifnot(exprs = {
        is.list(considerRes)
        is.character(whichDB)
        is.null(gmtdir) | is.character(gmtdir)
    })
  
    data_env <- new.env(parent = emptyenv())

    if (!is.null(gmtdir)) {
        thispathway <- fgsea::gmtPathways(gmtdir)
    } else {
        # utils::data(list=mydatabase)
        if(!whichDB %in% c("hallmark", "kegg", "reactome")) {
            message("If you want to customize database, provide through gmtdir.")
            stop("Otherwise,whichDB must be one of hallmark, kegg, reactome.")
        } else {
            if (whichDB == "hallmark") {
                data("pathways_hallmark", envir = data_env, package = "RegionalST")
                thispathway <- data_env[["pathways_hallmark"]]
            } else if (whichDB == "kegg") {
                data("pathways_kegg", envir = data_env, package = "RegionalST")
                thispathway <- data_env[["pathways_kegg"]]
            } else if (whichDB == "reactome") {
                data("pathways_reactome", envir = data_env, package = "RegionalST")
                thispathway <- data_env[["pathways_reactome"]]
            }
        }
    }

    allcelltype <- c()
    allpair <- c()
    for(i in seq_len(length(considerRes))) {
        tmp <- stats::na.omit(unique(considerRes[[i]]$allDE$Comparison))
        CellType <- unlist(lapply(strsplit(tmp, split = ": "), "[[", 1))
        allcelltype <- unique(c(allcelltype, CellType))
        thispair <- unlist(lapply(strsplit(tmp, split = ": "), "[[", 2))
        allpair <- unique(c(allpair, thispair))
    }
    allCTres <- list()
    for (i in seq_len(length(allcelltype))) {
        oneCTres <- c()
        for (j in seq_len(length(allpair))) {
            tmpcomb <- paste0(allcelltype[i], ": ", allpair[j])
            T1.mkers<- as.data.frame(considerRes[[j]]$allDE[considerRes[[j]]$allDE$Comparison == tmpcomb, ])
            if (withProp) {
                res.t1 <- T1.mkers %>%
                    dplyr::select(gene, effect_size) %>%
                    tibble::deframe()
            } else {
                res.t1 <- T1.mkers %>%
                    dplyr::select(gene, avg_log2FC) %>%
                    tibble::deframe()
            }
            if(length(res.t1) > 3) {
                fgseaRes.t1 <- fgsea::fgsea(pathways=thispathway, stats=na.omit(res.t1))
                fgseaResTidy.t1 <- fgseaRes.t1 %>%
                    dplyr::as_tibble() %>%
                    dplyr::arrange(pval)
                fgseaResTidy.t1$group <- allpair[j]
                if(j == 1) {
                    oneCTres <- fgseaResTidy.t1
                } else {
                    oneCTres <- rbind(oneCTres, fgseaResTidy.t1)
                }
            }
        }
        allCTres[[i]] <- oneCTres
    }
    names(allCTres) <- allcelltype

    return(allCTres)
}
