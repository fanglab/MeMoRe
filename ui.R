library(shiny)
library(plotly)

shinyUI(fluidPage(verticalLayout(
  titlePanel("SMRTMotif Plot"),
  wellPanel(fluidRow(
    column(
      3,
      fileInput('modfile', 'Choose CSV.GZ File', accept = c('.gz')),
      fileInput('genfile', 'Choose FASTA File', accept = c('.fasta'))
    ),
    column(
      5,
      offset = 1,
      div(
        style = "display: inline-block;vertical-align:top; width: 150px;",
        textInput(
          inputId = "motif",
          label = "Motif",
          value = "CAAAAA"
        )
      ),
      div(style = "display: inline-block;vertical-align:top; width: 150px;",
          textInput(
            inputId = "center",
            label = "Center",
            value = 5
          ))
    ),
    column(
      2,
      actionButton("submit", "Submit", style = "color: white; background-color: #337AB7"),
      actionButton("reset", "Clear")
    )
  )),
  tabsetPanel(
    type = "tabs",
    tabPanel("Score", plotlyOutput("score", width = "auto"), style = "overflow-y:scroll;"),
    tabPanel("ipdRatio", plotlyOutput("ipd", width = "auto"), style = "overflow-y:scroll;"),
    tabPanel("Coverage", plotlyOutput("coverage", width = "auto"), style = "overflow-y:scroll;")
  )
)))