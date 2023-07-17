library(RegionalST)

test_that("FindUnitDist works", {
    data("example_sce")
    res <- FindUnitDist(example_sce, avern = 5)

    expect_equal(round(res), 290)
})
#> Test passed

