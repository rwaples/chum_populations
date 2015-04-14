library(shiny)
library(ggplot2)
library(RColorBrewer)

#load base data
map_fst = read.table('map_pop_bayescan.txt', 
                 header = TRUE,sep = "\t")

shinyServer(function(input, output) {
  
  # reactive data (not needed here)
  mydata.reac <- reactive({
    map_fst_flagged = map_fst
    fst_to_plot = input$plotfst
    map_fst_flagged$plotfst = map_fst_flagged[[fst_to_plot]]
    map_fst_flagged$outlier = map_fst_flagged$qval < input$fdr
    return(map_fst_flagged)})
  
  # output
  output$plot1 <- renderPlot({
    #ggplot
    mydata <- mydata.reac()
    p = ggplot(mydata, aes(x=cM, y=plotfst)) + 
        facet_grid(LEP_LG ~ .) +
        geom_point(aes(color=outlier), alpha = .7, size = 3) +
        scale_color_manual(values=c("darkgrey","red")) +
        stat_smooth(fill="blue", colour="darkblue", size=2, alpha = 0.1, 
                    method = 'loess', span=input$span) +
        theme_bw() + theme(strip.text = element_text(size=25))
    print(p)
    
  },height = 20000, width = 1000)
  
})
