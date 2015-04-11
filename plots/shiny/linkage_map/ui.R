library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Kernel-smoothed Fst"),
  
  # Sidebar with a slider input for the number of bins
  sidebarPanel(
    numericInput("span", "alpha (smoothing parameter):", .2, min=.01, max =2.01, step=.01)
  ),
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("plot1")
  )
)
)
