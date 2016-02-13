#### Load necessary packages and data ####
library(shiny)
library(networkD3)

data(MisLinks)
data(MisNodes)


shinyServer(function(input, output) {
  
  output$netGraph <- renderForceNetwork({
    forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
                 Target = "target", Value = "value", NodeID = "name",
                 Group = "group", opacity = input$opacity)
  })
  
})
