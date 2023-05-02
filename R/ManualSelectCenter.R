#' Manually select top ROIs
#'
#' @param sce A single cell experiment object.
#'
#' @return An sce object with selected centers and radiuses.
#' @export
#' @examples
#' data("example_sce")
#' example_sce <- BayesSpace::spatialPreprocess(example_sce, platform="Visium", log.normalize=TRUE)
#' # I commented this out because the shiny app will get stuck without input.
#' # example_sce <- ManualSelectCenter(example_sce)
#'
ManualSelectCenter <- function(sce) {
    res <- ManualSelectCenter_app(sce)
    S4Vectors::metadata(sce)$selectCenters <- NULL
    S4Vectors::metadata(sce)$selectCenters <- list(selectCenters = res$centers,
                                        selectRadius = res$radius)
    return(sce)
}

ManualSelectCenter_app <- function(sce){

    outputVert <- GetVertices(sce, fill = "celltype", platform = "Visium")

    ui <- shiny::fluidPage(
        shiny::titlePanel("Manually choose your centers and radiuses"),

        # Sidebar layout with input and output definitions ----
        shiny::sidebarLayout(

            # Sidebar panel for inputs ----
            shiny::sidebarPanel(
                shiny::h3("Specify the number of centers:"),
                shiny::sliderInput("Ncenter",
                                   shiny::strong("Number of planned center points:"),
                                   value = 5,
                                   min = 2,
                                   max = 15),
                shiny::textInput("Radius", shiny::strong("Desired radiuses(separated by ',', below is an example):"),
                                 value = "10,10,5,10,10"),
                shiny::p("Make sure the provided radiuses have the equal number as your selected center points"),
                shiny::p("While you select centers, the order may be different from your selection order"),
                shiny::p("You can adjust the radius amount at the end to match your centers"),
                shiny::br(),
                shiny::h3("Instructions for selecting center points:"),
                shiny::tags$ol(
                    shiny::tags$li("Highlight nodes by clicking."),
                    shiny::tags$li("Reset the selection by double-clicking."),
                    shiny::tags$li(paste("Choose one dot each time. Repeat until you have",
                                         "chosen all your desired center nodes.")),
                    shiny::tags$li("Adjust the input radiuses to desirable amount."),
                    shiny::tags$li("Click 'Done'.")
                ),
                # clear button
                shiny::actionButton("reset", "Clear"),
                # done button
                shiny::actionButton("done", "Done")
            ),

            # Main panel for displaying outputs ----
            shiny::mainPanel(
                shiny::textOutput("selected_var"),
                shiny::plotOutput("plot", width = 500, height = 500,
                                  click = shiny::clickOpts(id = "plot_click"),
                                  dblclick = "plot_reset",
                                  brush = shiny::brushOpts(id = "plot1_brush"))
            )
        )
    )

    server <- function(input, output, session) {

        selected <- shiny::reactiveVal(rep(FALSE, nrow(outputVert)))

        # Toggle points that are clicked
        shiny::observeEvent(input$plot_click, {
            clicked <- shiny::nearPoints(outputVert,
                                         xvar = "x.pos",
                                         yvar = "y.pos",
                                         maxpoints = 1,
                                         input$plot_click, allRows = TRUE)$selected_
            selected(clicked | selected())
        })
        shiny::observeEvent(input$reset, {
            selected(rep(FALSE, nrow(outputVert)))
        })
        output$plot <- shiny::renderPlot({

            allR <- as.numeric(unlist(strsplit(input$Radius, split = ",")))
            selIdx <- which(selected())
            allneighbour <- c()

            output$selected_var <- shiny::renderText({
                if (sum(selected()) < input$Ncenter) {
                    paste0("Please select ", input$Ncenter - sum(selected()), " more centers!")
                } else {
                    paste0("Successfully selected ", input$Ncenter, " center!")
                }
            })

            for(i in seq_len(sum(selected()))) {
                reg <- FindRegionalCells(sce,
                                         centerID = outputVert$spot[selIdx[i]],
                                         radius = allR[i],
                                         avern = 5,
                                         doPlot = FALSE)
                allneighbour = c(allneighbour, reg$closeID)
            }

            close_df <- outputVert[match(allneighbour, outputVert$spot),]

            dotcolors <- c(RColorBrewer::brewer.pal(8, "Dark2"),
                           RColorBrewer::brewer.pal(9, "Set1"))
            mypalette <- makeTransparent(colorspace::qualitative_hcl(length(unique(outputVert$fill)),
                                                                     palette = "Set2"), 20)
            mypalette_ori <- mypalette

            alphavec <- rep(0.5, length(unique(outputVert$fill)))

            outputVert$newfill <- outputVert$fill

            if (sum(selected()) > 0) {
                outputVert$newfill[selected()] <- paste0("Center", seq_len(sum(selected())))
                outputVert$newfill <- factor(outputVert$newfill,
                                             levels = c(unique(outputVert$fill),
                                                        paste0("Center", seq_len(sum(selected())))))
                mypalette <- c(mypalette_ori,
                               dotcolors[seq_len(sum(selected()))])
                alphavec <- c(rep(0.5, length(unique(outputVert$fill))), rep(1, sum(selected())))
            } else {
                outputVert$newfill <- factor(outputVert$fill,
                                             levels = unique(outputVert$fill))}

            ggplot2::ggplot(outputVert, ggplot2::aes(x = x.pos, y = y.pos,
                                                     color = factor(newfill))) +
                ggplot2::geom_point(data = close_df,
                                    ggplot2::aes(x = x.pos,y = y.pos),
                                    color = "lightblue",
                                    size=2)+
                ggplot2::geom_point(size = 2) + ggplot2::labs(color = "Groups") +
                ggplot2::theme_void() + ggplot2::scale_color_manual(values = mypalette)+
                ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = alphavec)))

        }, res = 96, width = 800)

        shiny::observeEvent(input$done, {
            shiny::stopApp(list(centers = outputVert$spot[which(selected())],
                                radius = as.numeric(unlist(strsplit(input$Radius, split = ",")))))
        })

    }
    sel <- shiny::runApp(shiny::shinyApp(ui, server))

    return(sel)
}

