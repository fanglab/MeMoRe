library(shiny)
library(DT)

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
      4,
      offset = 4, align="center",
      h5("Add a specific motif to table"),
      div(
        style = "display: inline-block;vertical-align:top; width: 180px;",
        textInput(
          inputId = "motif",
          label = NULL,
          value = "Motif (e.g. GATC)"
        )
      ),
      div(style = "display: inline-block;vertical-align:top; width: 180px;",
          textInput(
            inputId = "center",
            label = NULL,
            value = "Modified Position (e.g. 1)"
          )),
      actionButton("addmotif", "Add motif"), actionButton("cleartable", "Clear table")
    ),
    column(
      12, 
      uiOutput('motifs')
    ),
    column(
      4, offset = 4, align="center",
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