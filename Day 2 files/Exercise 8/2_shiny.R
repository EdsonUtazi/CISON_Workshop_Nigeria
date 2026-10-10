# # 
# # ===================================================
# # 10. Shiny
# # ===================================================
library(shiny)

ui <- fluidPage(

  titlePanel("Simple Iris Plot App"),

  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "plot_type",
        label = "Choose plot type:",
        choices = c("Scatterplot", "Histogram", "Boxplot")
      )
    ),

    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {

  output$plot <- renderPlot({

    if (input$plot_type == "Scatterplot") {
      plot(iris$Sepal.Length, iris$Petal.Length,
           main = "Scatterplot: Sepal vs Petal Length",
           xlab = "Sepal Length",
           ylab = "Petal Length",
           col = "blue", pch = 19)

    } else if (input$plot_type == "Histogram") {
      hist(iris$Sepal.Length,
           main = "Histogram of Sepal Length",
           xlab = "Sepal Length",
           col = "lightgreen")

    } else if (input$plot_type == "Boxplot") {
      boxplot(Sepal.Length ~ Species, data = iris,
              main = "Boxplot of Sepal Length by Species",
              xlab = "Species",
              ylab = "Sepal Length",
              col = "orange")
    }
  })
}

shinyApp(ui = ui, server = server)

# 
# 