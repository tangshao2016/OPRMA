
library(shiny)



shinyUI(fluidPage(
  fluidRow(
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    h2("Online platform for peak information analysis in RNA editing using CLIP-Seq data", align = "center"),
    column(12,plotOutput('distPlot',height="600px"))
    
  ),
 
 
  hr(),
  
  fluidRow(
    column(1,
           h3("File input",align="center"),
           fileInput("file", label = h4("File input case"),accept=c('.txt')),
           fileInput("file2",label = h4("File input control(optional)"),accept = c('.txt'))
           
    ),
    column(2,offset=1,
           h3("Analysis type" ),
           selectInput("type", 
                       label = h4("Analysis method"),
                       choices = c("histogram(only for case group)", "density"),
                       selected = "density"),
           selectInput("org", 
                       label = h4("Choose model organism"),
                       choices = c("hg19", "mm10"),
                       selected = "mm10"),
           
           selectInput("analysis", 
                       label = h4("Choose an analysis method"),
                       choices = c("STOPCODON", "LASTEXON","53UTRCDS"),
                       selected = "LASTEXON")
           
    ),
    column(2,offset=1,
           h3("Display option"),
           sliderInput("range", h4("Range_X:"),
                       min = -1000, max = 1000, value = c(-1000,1000)), 
           sliderInput("rangey", h4("Range_Y:"),
                       min = 0, max = 0.0030, value = c(0,0.0020))
        
           
    ),
    column(2,offset=1,
           h3("Color option" ),
           radioButtons("color", label = h4("Color buttons(case)"),
                        choices = list("red" = "red", "blue" = "blue",
                                       "green" = "green"),selected = "red"),
           radioButtons("color2", label = h4("Color buttons(control)"),
                        choices = list("red" = "red", "blue" = "blue",
                                       "green" = "green"),selected = "blue")
           
    ),
    column(2,
           h3("Download",align="center"),
           selectInput("dataset", "Choose a dataset:", 
                       choices = c("merged_case", "peak_case","merged_con", "peak_con")),
           downloadButton('downloadData', 'Download'),
           h4("Explanation text"),
           helpText("merge(case/control):Summary table after merging the peak with gene annotation (case/control sample)", 
                    "peak(case/control):Distance metric between peak and stop-codon or last exon (case/control sample)")
           
    
    )
  ),
  fluidRow(
    column(12, DT::dataTableOutput('x1'))
  ),
  
  fluidRow(
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    h2("Plot peak distribution for selected gene", align = "center"),
    textInput("gene", label = h2("gene input"), 
              value = "URGCP"),
    column(12,plotOutput("rect1",height="300px"),plotOutput("rect2", height="300px"))
    
  )
  
  
  
))

