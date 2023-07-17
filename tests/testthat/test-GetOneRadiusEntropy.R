library(RegionalST)

test_that("GetOneRadiusEntropy works", {
    data("example_sce")
    weight <- data.frame(celltype = c("Cancer Epithelial", "CAFs",
                                        "T-cells", "Endothelial",
                                        "PVL", "Myeloid", "B-cells",
                                        "Normal Epithelial", "Plasmablasts"),
                          weight = c(0.25,0.05,
                                     0.25,0.05,
                                     0.025,0.05,
                                     0.25,0.05,0.025))
    example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
    res <- GetOneRadiusEntropy(example_sce, selectN = round(length(example_sce$spot)/2),
                        weight = weight, radius = 5, doPlot = TRUE,
                        mytitle = "Radius 5 weighted entropy")
  expect_equal(res$radius, 5)
})
#> Test passed
