library(shiny)
library(ggplot2)
library(RColorBrewer)

#load base data
mds = read.table('non_paralogs.IND_md_poistions', 
                 header = TRUE,sep = "\t")

pop_colors <- c("#000000", brewer.pal(9, "Set1"))

shinyServer(function(input, output) {
  
  # reactive data
  mydata.reac <- reactive({
    mds['xx'] = mds[input$xcol]
    mds['yy'] = mds[input$ycol]
    #mds['colorby'] = mds[input$colorby]
    return(mds)})
  
  # output
  output$plot1 <- renderPlot({
    #ggplot
    mydata <- mydata.reac()
    
    if (input$colorby == "Population"){
      p = ggplot(mydata, aes(x=xx, y=yy)) + 
        geom_point(alpha = .5, size = 8, aes(color = POPNAME)) +
        scale_color_manual(values = pop_colors)

    } else if (input$colorby == "Region"){
      p = ggplot(mydata, aes(x=xx, y=yy)) + 
        geom_point(alpha = .5, size = 8, aes(color = REGION))+
        scale_color_brewer(type = 'qual', palette="Set1")
    } else if (input$colorby == "Run Timing"){
      p = ggplot(mydata, aes(x=xx, y=yy)) + 
        geom_point(alpha = .5, size = 8, aes(color = TIMING)) +
        scale_color_brewer(type = 'qual', palette="Set2")
    }
    p = p + geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      theme_bw() + coord_equal()  +
      xlab(input$xcol) + ylab(input$ycol) +
      ggtitle("Puget Sound Chum Salmon Populations")
    print(p)
    
  },height = 800, width = 800)
  
})
