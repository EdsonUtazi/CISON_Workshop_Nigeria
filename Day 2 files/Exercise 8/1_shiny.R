
# #===================================================
# 
# 
# 
library(shiny)

ui <- fluidPage(

  titlePanel("Random Scatterplot App"),

  sidebarLayout(
    sidebarPanel(

      sliderInput(
        inputId = "n_points",
        label = "Number of points:",
        min = 10,
        max = 500,
        value = 100,
        step = 10
      ),

      selectInput(
        inputId = "point_color",
        label = "Point color:",
        choices = c("blue", "red", "green", "purple", "black")
      ),

      checkboxInput(
        inputId = "show_lowess",
        label = "Add LOWESS curve",
        value = FALSE
      )
    ),

    mainPanel(
      plotOutput("scatter_plot", height = "600px")
    )
  )
)

server <- function(input, output) {

  output$scatter_plot <- renderPlot({

    # Create random data
    set.seed(123)
    df <- data.frame(
      x = rnorm(input$n_points),
      y = rnorm(input$n_points)
    )

    # Scatter plot
    plot(df$x, df$y,
         pch = 19,
         col = input$point_color,
         main = "Random Scatterplot",
         xlab = "X",
         ylab = "Y")

    # Add LOWESS if checked
    if (input$show_lowess) {
      lines(lowess(df$x, df$y), col = "red", lwd = 2)
    }
  })
}

shinyApp(ui = ui, server = server)