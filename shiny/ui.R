library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Tumor Evolution Graph Maker"),

  sidebarPanel(
  	fileInput('file1', 'Choose TXT File',
  		accept=c('text/plain'))
  ),

  mainPanel()
))
