library(RegionalST)

test_that("FindRegionalCells works", {
    data("example_sce")
    res <- FindRegionalCells(example_sce, centerID = "GTGCACGAAAGTGACT-1", radius = 3)

    expect_equal(length(res$closeID), 18)
})
#> Test passed

