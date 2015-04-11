library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Alelle-sharing distance"),
  
  # Sidebar with a slider input for the number of bins
    sidebarPanel(
      selectInput("xcol", "X axis:", c('dim_1', 'dim_2', 'dim_3', 'dim_4', 'dim_5', 'dim_6')),
      selectInput("ycol","Y axis:",  c('dim_1', 'dim_2', 'dim_3', 'dim_4', 'dim_5', 'dim_6')),
      selectInput("colorby","Color by:", c('Population', 'Region', 'Run Timing'))
      
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot1")
    )
  )
)
