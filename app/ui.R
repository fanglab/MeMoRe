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
      fileInput('modfile', 'modifications.csv(.gz)', accept = c('.gz', '.csv'))
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
      3, align="center",
      h5("Add a specific motif to table"),
      div(
        style = "display: inline-block;vertical-align:top; width: 180px;",
        textInput(
          inputId = "motif",
          label = "Motif",
          value = "GATC"
        )
      ),
      div(
        style = "display: inline-block;vertical-align:top; width: 100px;",
        textInput(
          inputId = "center",
          label = "Modified Pos.",
          value = "1"
        )
      ),
      br(),
      actionButton("addmotif", "Add motif"), 
      actionButton("cleartable", "Clear table"),
      br(),br(),
      verbatimTextOutput("list_selected_motifs"),
      br(),
      actionButton("submit", "Generate plots for selected motifs", style = "color: white; background-color: #337AB7")
    ),
    column(
      9, 
      uiOutput('motifs')
    )
  )),
  h4("Generated Plots"),
  wellPanel(fluidRow(
    column(
      4, align="center",
      downloadButton('dl_everything', "Download all generated plots as .zip")
    ),
    column(
      8, 
      uiOutput('dl_dt')
    )
  )),
  uiOutput('images')
)))