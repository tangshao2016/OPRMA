
library(shiny)



shinyUI(fluidPage(
  
  
  plotOutput('distPlot',height = "550px"),
  
  hr(),
  
  fluidRow(
    column(3,
           h2("Distribution Plot"),
           fileInput("file", label = h3("File input(experimental)"),accept=c('.txt')),
           fileInput("file2",label = h3("File input(control)"),accept = c('.txt'))
        
    ),
    column(4, offset = 1,
           sliderInput("range", h3("Range:"),
                       min = -1000, max = 1000, value = c(-1000,1000)), 
           
           selectInput("org", 
                       label = h3("Choose a org"),
                       choices = c("hg19", "mm10"),
                       selected = "mm10"),
           
           selectInput("analysis", 
                       label = h3("Choose an analysis method"),
                       choices = c("STOPCODON", "LASTEXON","53UTRCDS"),
                       selected = "LASTEXON")
           
    ),
    column(4,
           selectInput("type", 
                       label = h3("Choose a type"),
                       choices = c("histogram(only for experimental group)", "density"),
                       selected = "density"),
           radioButtons("color", label = h3("Color buttons(experimental)"),
                        choices = list("red" = "red", "blue" = "blue",
                                       "green" = "green"),selected = "red"),
           radioButtons("color2", label = h3("Color buttons(control)"),
                        choices = list("red" = "red", "blue" = "blue",
                                       "green" = "green"),selected = "blue")
           
    )
  )
))
