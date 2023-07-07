library(RegionalST)

data("example_sce")
example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
PlotOneSelectedCenter(example_sce, ploti = 1)
