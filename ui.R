library(shiny)
library(plotly)

shinyUI(fluidPage(verticalLayout(
  titlePanel("SMRTMotif Plot"),
  wellPanel(fluidRow(
    column(
      4,
      fileInput('modfile', 'modifications.csv.gz', accept = c('.gz')),
      fileInput('genfile', 'genome.fasta', accept = c('.fasta')),
      fileInput('motfile', 'motif_summary.csv', accept = c('.csv'))
    ),
    column(
      4,
      offset = 0,
      uiOutput('motifs'),
      h5("or", align = "center"),
      h5("specify a specific motif"),
      div(
        style = "display: inline-block;vertical-align:top; width: 180px;",
        textInput(
          inputId = "motif",
          label = NULL,
          value = "Motif (e.g. CAAAAA)"
        )
      ),
      div(style = "display: inline-block;vertical-align:top; width: 180px;",
          textInput(
            inputId = "center",
            label = NULL,
            value = "Center (e.g. 5)"
          )),
      hr(),
      radioButtons("radio", label = "Choose Motif Input Type:",
                   choices = list("Motif from File" = 1, 
                                  "Motif from Text Box" = 2), selected = 1)
    ),
    column(
      4, 
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