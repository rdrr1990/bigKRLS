# this is an early version of a function ultimately to be included in the bigKRLS package.
# execute this code and then run 'shiny_bigKRLS(output)'
# to interact with the results.
# Set 'shiny_bigKRLS(output, export=T)' if you would like to be able 
# to use the output on a different machine or post it on a server.
# While we will not be able to fully anticipate everyone's needs
# we will be adding features and posting vignettes to facilitate customization. 
# In the meantime, feedback on features likely of interest to many is most welcome!
# 09.01.2016

# demo code
#library(bigKRLS)
#X <- matrix(runif(5000, -2*pi, 2*pi), ncol=5)
#y <- sin(X[,1])*X[,1] + X[,4:5] %*% 4:5 + rnorm(1000)
#out <- bigKRLS(y, X)
#summary(out)
#run code below here
#shiny_bigKRLS(out) 


library(shiny)
library(bigKRLS)


shiny_bigKRLS <- function(out, export=F){
  
  palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
            "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
  
  bigKRLS_server <- shinyServer(function(input, output, session) {
    
    selectedData <- reactive({
      return(list(cbind(out$derivatives[, input$dydxp], 
                        out$X[, c(input$xp)]), input$type))
    })
    
    output$graph <- renderPlot({
      
      if(selectedData()[[2]] == "Smooth"){
        
        L <- loess.smooth(x=selectedData()[[1]][,2], 
                          y=selectedData()[[1]][,1])
        
        plot(y=L$y, x=L$x, ylab=paste("Marginal Effect of",input$dydxp), pch = 19, bty = "n",
             main=main.label, 
             xlab=paste("Observed Level of", input$xp), cex=2, cex.axis=1.5,  cex.lab=1.4,
             col = colorRampPalette(c("blue", "red"))(length(L$y))[rank(L$y)])
        
      }else{
        plot(x=(selectedData()[[1]][,2]), xlab = paste("Observed Level of", input$xp),
             y=(selectedData()[[1]][,1]), ylab = paste("Marginal Effect of",input$dydxp), 
             pch = 4, bty = "n", cex=2, cex.axis=1.5,  cex.lab=1.4,
             main=main.label,
             col = colorRampPalette(c("green", "purple"))(nrow(out$X))[rank(out$coeffs^2)], 
             ylim = range(selectedData()[[1]][,1])*c(.8, 1.25), 
             xlim = range(selectedData()[[1]][,2])*c(.8, 1.25))
        
        fields::image.plot(legend.only = T, zlim=c(1/nrow(out$X), 1), 
                           legend.cex = 0.75,legend.shrink = .4,   
                           col = colorRampPalette(c("purple", "green"))(nrow(out$X)))
        text(x = 1.2*range(selectedData()[[1]][,2])[2], 
             y = .5*range(selectedData()[[1]][,1])[2], 
             "Relative Fit \nIn Full Model") 
      }
    })})
  
  bigKRLS_ui <- shinyUI(fluidPage(
    
    titlePanel(main.label),
    
    sidebarPanel(
      selectInput('dydxp', 'Local Derivatives (dy/dx)', colnames(out$derivatives)),
      selectInput('xp', 'x', colnames(out$X)), 
      selectInput('type', 'Plot Type', c("Smooth", "Scatter"))
    ),
    
    mainPanel(plotOutput('graph'))
    
  ))
  
  if(export){
    
    out <- out
    out$K <- tmp$vcov.c <- tmp$vcov.fitted <- NULL
    for(i in which(unlist(lapply(out, is.big.matrix)))){
      out[[i]] <- as.matrix(out[[i]])
    }
    
    save(out, file="shiny_out.rdata")
    
    cat("A re-formatted version of your output has been saved with file name \"shiny_out.rdata\" in your current working directory:\n", getwd(),
        "\nFor a few technical reasons, the big N * N matrices have been removed and the smaller ones converted back to base R;\nthis should make your output small enough for the free version of Shiny's server.\nTo access the Shiny app later or on a different machine, simply execute this script with the following commands:\n",
        "\nload(\"shiny_out.rdata\")\nNext, execute this script to make sure Shiny is initialized with current values. \nshiny_bigKRLS(out)")
  }else{
    shinyApp(ui = bigKRLS_ui, server = bigKRLS_server)
  }
}



  
