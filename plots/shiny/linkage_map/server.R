library(shiny)
library(ggplot2)
library(RColorBrewer)

#load base data
linkage_map = read.table('/home/ipseg/Desktop/waples/chum_populations/results/batch_4/pop_analysis/map_pop_fst.txt', 
                 header = TRUE,sep = "\t")

shinyServer(function(input, output) {
  # reactive data (not needed here)
  mydata.reac <- reactive({
    return(linkage_map)})
  
  # output
  output$plot1 <- renderPlot({
    #ggplot
    mydata <- mydata.reac()
    p = ggplot(mydata, aes(x=cM, y=FST)) + 
        facet_grid(LEP_LG ~ .) +
        stat_smooth(fill="blue", colour="darkblue", size=2, alpha = 0.1, method = 'loess', span=input$span)+
        geom_point(alpha = .5) +
        theme_bw()
    print(p)
    
  },height = 20000, width = 1000)
  
})
