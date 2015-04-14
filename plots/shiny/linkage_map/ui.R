library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Kernel-smoothed Fst along linkage groups"),
  
  # Sidebar with a slider input for the number of bins
  
  sidebarPanel(

    selectInput("plotfst", label = "Fst measure:", 
                list("PLINK" = 'Fst_plink', "BayeScan" = 'Fst_bayescan'), selected = 'Fst_plink'),
    selectInput("fdr", label = "False discovery rate:", 
                choices = c(0.05, 0.01, 0.005, 0.001, .0001), selected = 0.01),
    numericInput("span", label="alpha (kernel smoothing parameter):", .5, min=.01, max =2.01, step=.01), 
    h2('-----------------------------'), 
    h3('Loci are plotting along linkage groups, cM on the x-axis'),
    h3('Fst on the y-axis'),
    h3('grey dots are non-outlier loci'),
    h3('red dots are outlier loci at the selected FDR'), 
    h3('blue lines trace the kernel-smoothed Fst'),
    h3('shaded areas show the 95% CI')
    
  ),
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("plot1")
  )
)
)
