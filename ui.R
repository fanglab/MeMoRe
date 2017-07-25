library(shiny)
library(DT)

shinyUI(fluidPage(verticalLayout(
  titlePanel("SMRTMotif Plot"),
  
  tags$head(tags$style(
    type="text/css",
    "#combined img {max-width: auto; height: 900px;}",
    "#score img {max-width: auto; height: 300px;}",
    "#ipd img {max-width: auto; height: 300px;}",
    "#coverage img {max-width: auto; height: 300px;}",
    "#submit { margin-top: 100px;}",
    "#reset { margin-top: 100px;}"
  )),
  
  wellPanel(fluidRow(
    column(
      4,
      fileInput('modfile', 'modifications.csv.gz', accept = c('.gz'))
    ),
    column(
      4,
      fileInput('genfile', 'genome.fasta', accept = c('.fasta'))
    ),
    column(
      4,
      fileInput('motfile', 'motif_summary.csv', accept = c('.csv'))
    ),
    column(
      12, 
      uiOutput('motifs')
    ),
    column(
      4,
      offset = 4, align="center",
      h5("or"),
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
            value = "Modified Position (e.g. 5)"
          )),
      actionButton("addmotif", "Add motif"), actionButton("cleartable", "Clear table")
    ),
    column(
      4, align="right",
      actionButton("submit", "Generate plots for selected motif", style = "color: white; background-color: #337AB7"),
      actionButton("reset", "Clear plots")
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