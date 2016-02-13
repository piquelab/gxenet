#### Load necessary packages and data ####
library(shiny)
library(networkD3)

load("./net.Rd")

names(comb) <- c("tr","ct","values")
names(cytlist)=row.names(comb)

trList <- levels(comb$tr)
names(trList) <- treatmentKEY[trList,"Common_Name"]


# cyt$nodeData$Group = "No_change"
# cyt$nodeData$Group[(cyt$nodeData$nodeAttribute.logFC>0) & (cyt$nodeData$nodeAttribute.padj<0.1)] = "Upregulated"
# cyt$nodeData$Group[(cyt$nodeData$nodeAttribute.logFC<0) & (cyt$nodeData$nodeAttribute.padj<0.1)] = "Downregulated"
prepCyt <- function(tr,ct,modnum){
  sel <- paste(tr,ct,sep=":")
  cyt <- cytlist[[sel]][[modnum]]
  rownames(cyt$nodeData) <- cyt$nodeData$altName
  cyt$nodeData$nodeNum=(1:nrow(cyt$nodeData))-1
  cyt$edgeData$Source = cyt$nodeData[cyt$edgeData$fromAltName,"nodeNum"]
  cyt$edgeData$Target = cyt$nodeData[cyt$edgeData$toAltName,"nodeNum"]
  cyt$nodeData$nodeAttribute.ModConnec =cyt$nodeData$nodeAttribute.ModConnec^2 
  cyt$nodeData$Group <- as.character(cut(cyt$nodeData$nodeAttribute.logFC,breaks=c(-Inf,-1,-0.5,-0.25,0,0.25,0.5,1,Inf)))
  cyt
}

##cyt <- prepCyt("T12C1:LCL",66)
##cyt <- prepCyt("T12C1","LCL",66)

shinyServer(function(input, output) {

  cyt <- reactive({
    prepCyt(input$tr,input$ct,input$modnum)
  })
  
  output$text1 <- renderText({ 
    paste("You have selected: ",input$tr,input$ct,input$modnum)
  })
  
  output$netGraph <- renderForceNetwork({
    forceNetwork(Links = cyt()$edgeData, Nodes = cyt()$nodeData, Source = "Source",
                 Target = "Target", Value = "weight", NodeID = "nodeName",
                 Group = "Group",
                 Nodesize = "nodeAttribute.ModConnec", 
                 legend = T,
                 colourScale = JS('d3.scale.category20c().domain(["(-Inf,-1]","(-1,-0.5]","(-0.5,-0.25]","(-0.25,0]","(1, Inf]","(0.5,1]","(0.25,0.5]","(0,0.25]"])'),
                 fontFamily = "helvetica",fontSize = 12,opacity=1,bounded=T)
    })
})

## simpleNetwork(cyt$edgeData, Source = "fromNode", Target = "toNode", fontFamily = "helvetica",fontSize = 12)


