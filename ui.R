
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(networkD3)

load("terms.Rd")

shinyUI(fluidPage(

  # Application title
  titlePanel("GxE Network Browser"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput('ct', 'Cell-type:', levels(comb$ct)),
      selectInput('treat', 'Treatment:', treatList,selected = "Dexamethasone"),
      numericInput('modnum', 'Module number:', 66,
                   min = 1, max = 87),
      p(textOutput("text1")),
      p("It may take a few seconds to load the data.")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      forceNetworkOutput("netGraph")
    )
  )
))
