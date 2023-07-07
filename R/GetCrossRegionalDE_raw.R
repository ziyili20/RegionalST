#' Identify cross-regional differential analysis
#'
#' @param sce A single cell experiment object.
#' @param twoCenter A vector of two numbers for the interested ROI numbers.
#' @param enhanced A logical variable for using enhanced data or not.
#' @param label A variable name that contains the cell type information.
#' @param n_markers A number specifying the top DE gene number.
#' @param logfc.threshold A number for the cutoff threshold of log fold change.
#' @param angle A number for angle when plotting.
#' @param hjust A number for horizontal justification when plotting.
#' @param size A number for text font size.
#' @param min.pct A number of minimum percentage specified in the Seurat DE function.
#' @param padj_filter A number for filtering adjusted p values.
#' @param doHeatmap Logical variable for whether drawing the heatmap.
#'
#' @return A list including the top DE genes (topDE), and all DE genes (allDE).
#' @export
#' @examples
#' data("example_sce")
#' example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
#' # I used a very big padj filter here because this is just a toy data
#' GetCrossRegionalDE_raw(example_sce, twoCenter = c(1,2),
#'                        min.pct = 0.01, logfc.threshold = 0.01,
#'                        padj_filter = 0.5)
#'
GetCrossRegionalDE_raw <- function(sce,
                                   twoCenter = c(3,4),
                                   enhanced = FALSE,
                                   label = "celltype",
                                   n_markers = 10,
                                   logfc.threshold = 0.25,
                                   angle = 30,
                                   hjust = 0,
                                   size = 3,
                                   min.pct = 0.1,
                                   padj_filter = 0.05,
                                   doHeatmap = TRUE) {

    thisID1 <- S4Vectors::metadata(sce)$selectCenters$selectID[twoCenter[1]]
    thisRadius1 <- S4Vectors::metadata(sce)$selectCenters$selectRadius[twoCenter[1]]
    OneRegOut1 <- FindRegionalCells(sce,
                                   centerID = thisID1,
                                   radius = thisRadius1,
                                   enhanced = enhanced,
                                   doPlot = FALSE,
                                   returnPlot = FALSE)

    thisID2 <- S4Vectors::metadata(sce)$selectCenters$selectID[twoCenter[2]]
    thisRadius2 <- S4Vectors::metadata(sce)$selectCenters$selectRadius[twoCenter[2]]
    OneRegOut2 <- FindRegionalCells(sce,
                                   centerID = thisID2,
                                   radius = thisRadius2,
                                   enhanced = enhanced,
                                   doPlot = FALSE,
                                   returnPlot = FALSE)

    newcoldata <- SummarizedExperiment::colData(sce)
    newcoldata$TwoGroupInfo <- rep(NA, nrow(newcoldata))

    if (enhanced) {
        newcoldata$TwoGroupInfo[match(OneRegOut1$closeID, colnames(sce@assays@data$logcounts))] <- twoCenter[1]
        newcoldata$TwoGroupInfo[match(OneRegOut2$closeID, colnames(sce@assays@data$logcounts))] <- twoCenter[2]

        newcoldata$NewCellType <- paste0(newcoldata$TwoGroupInfo, "_", eval(parse(text = paste0("sce$", label))))
        CellType1 <- eval(parse(text = paste0("sce$", label)))[match(OneRegOut1$closeID, colnames(sce@assays@data$logcounts))]
        CellType2 <- eval(parse(text = paste0("sce$", label)))[match(OneRegOut2$closeID, colnames(sce@assays@data$logcounts))]

    } else {
        newcoldata$TwoGroupInfo[match(OneRegOut1$closeID, sce$spot)] <- twoCenter[1]
        newcoldata$TwoGroupInfo[match(OneRegOut2$closeID, sce$spot)] <- twoCenter[2]

        newcoldata$NewCellType <- paste0(newcoldata$TwoGroupInfo, "_", eval(parse(text = paste0("sce$", label))))
        CellType1 <- eval(parse(text = paste0("sce$", label)))[match(OneRegOut1$closeID, sce$spot)]
        CellType2 <- eval(parse(text = paste0("sce$", label)))[match(OneRegOut2$closeID, sce$spot)]
    }

    allcloseID <- c(OneRegOut1$closeID, OneRegOut2$closeID)
    seuratObj <- Seurat::CreateSeuratObject(counts=SingleCellExperiment::logcounts(sce)[, match(allcloseID, colnames(SingleCellExperiment::logcounts(sce)))],
                                            assay='Spatial',
                                            meta.data=as.data.frame(newcoldata)[match(allcloseID, colnames(SingleCellExperiment::logcounts(sce))), ])
    Seurat::Idents(seuratObj) <- seuratObj@meta.data$NewCellType

    allDE <- c()
    allCT <- unique(intersect(unique(CellType1), unique(CellType2)))
    for(i in seq_len(length(allCT))) {
        message(paste0("Processing to cell type:", allCT[i]))
        if(i == 1) {
            n1 <- sum(Seurat::Idents(seuratObj) == paste0(twoCenter[1], "_", allCT[i]))
            n2 <- sum(Seurat::Idents(seuratObj) == paste0(twoCenter[2], "_", allCT[i]))

            if (n1 > 3 & n2 > 3) {
                newDEres <- Seurat::FindMarkers(seuratObj, ident.1 = paste0(twoCenter[1], "_", allCT[i]),
                                                ident.2 = paste0(twoCenter[2], "_", allCT[i]),
                                                logfc.threshold = logfc.threshold,
                                                min.pct = min.pct)
                newDEres$`Comparison` <- paste0(allCT[i], ": Region", twoCenter[1], " vs Region", twoCenter[2])
                newDEres$`gene` <- rownames(newDEres)
                allDE <- newDEres
            }

        } else {
            n1 <- sum(Seurat::Idents(seuratObj) == paste0(twoCenter[1], "_", allCT[i]))
            n2 <- sum(Seurat::Idents(seuratObj) == paste0(twoCenter[2], "_", allCT[i]))

            if (n1 > 3 & n2 > 3) {
                newDEres <- Seurat::FindMarkers(seuratObj,
                                                ident.1 = paste0(twoCenter[1], "_", allCT[i]),
                                                ident.2 = paste0(twoCenter[2], "_", allCT[i]),
                                                logfc.threshold = logfc.threshold,
                                                min.pct = min.pct)
                newDEres$`Comparison` <- paste0(allCT[i], ": Region",
                                                twoCenter[1], " vs Region",
                                                twoCenter[2])
                newDEres$`gene` <- rownames(newDEres)
                allDE <- rbind(allDE, newDEres)
            }
        }
    }

    if(is.null(padj_filter)) {
        top_markers <- allDE %>%
            dplyr::group_by(Comparison) %>%
            dplyr::top_n(n_markers, avg_log2FC)
    } else {
        allDE <- allDE[allDE$p_val_adj <= padj_filter,]
        top_markers <- allDE %>%
            dplyr::group_by(Comparison) %>%
            dplyr::top_n(n_markers, avg_log2FC)
    }

    if (doHeatmap) {
        ## Scale data
        seuratObj@assays$Spatial@scale.data <-
            seuratObj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t

        palette <- RColorBrewer::brewer.pal(2*length(allCT), "Paired")
        Seurat::Idents(seuratObj) <- factor(Seurat::Idents(seuratObj), levels = paste0(twoCenter, "_", rep(allCT, each = 2)))
        ## Plot expression of markers
        p1 <- Seurat::DoHeatmap(seuratObj, features = top_markers$gene, slot='scale.data',
                                group.colors=palette,
                                angle=angle, size=size, hjust = hjust, label = TRUE, raster=FALSE) +
            guides(col = FALSE)
        print(p1)
    }
    return(list(allDE = allDE,
                topDE = top_markers))
}
