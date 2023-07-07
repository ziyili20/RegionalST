library(RegionalST)

test_that("GetCrossRegionalDE_raw works", {
    data("example_sce")
    example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
    res <- GetCrossRegionalDE_raw(example_sce, twoCenter = c(1,2),
                           min.pct = 0.001, logfc.threshold = 0,
                           padj_filter = 0.5,
                           doHeatmap = FALSE)

    expect_equal(ncol(res$allDE), 7)

})
#> Test passed
