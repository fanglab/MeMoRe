library(shiny)

shinyUI(fluidPage(verticalLayout(
  titlePanel("SMRTMotif Plot"),
  
  tags$head(tags$style(
    type="text/css",
    "#combined img {max-width: auto; height: 900px;}",
    "#score img {max-width: auto; height: 300px;}",
    "#ipd img {max-width: auto; height: 300px;}",
    "#coverage img {max-width: auto; height: 300px;}"
  )),
  
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
    tabPanel("Combined", br(),
             downloadButton('dl_a'),
             imageOutput("combined"),
             style = "overflow-y:scroll;text-align: center;"),
    tabPanel("Score", br(),
             downloadButton('dl_s'),
             imageOutput("score"),
             style = "overflow-y:scroll;text-align: center;"),
    tabPanel("ipdRatio", br(),
             downloadButton('dl_i'),
             imageOutput("ipd"),
             style = "overflow-y:scroll;text-align: center;"),
    tabPanel("Coverage", br(),
             downloadButton('dl_c'),
             imageOutput("coverage"),
             style = "overflow-y:scroll;text-align: center;")
  )
)))