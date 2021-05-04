library(shiny)
library(DT)
library(shinydashboard)
library(shinyjs)

shinyUI(fluidPage(verticalLayout(
  useShinyjs(),
  titlePanel("SMRTMotif Plot"),
  
  tags$head(tags$style(
    type="text/css",
    "[id^='combined_'] img {max-width: 100%; height: auto;}",
    "h4 {display: inline;}",
    "[id='download_button'] {display: inline;}"
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
          placeholder = "GATC"
        )
      ),
      div(
        style = "display: inline-block;vertical-align:top; width: 100px;",
        textInput(
          inputId = "center",
          label = "Modified Pos.",
          placeholder = "1"
        )
      ),
      br(),
      actionButton("addmotif", "Add motif"), 
      actionButton("cleartable", "Clear table"),
      br(),br(),
      verbatimTextOutput("list_selected_motifs"),
      br(),
      actionButton("submit", "Process selected motifs", style = "color: white; background-color: #337AB7"),
      br(),br(),
      actionButton("submit_all", "Process all motifs", style = "color: white; background-color: #337AB7")      
    ),
    column(
      9, 
      uiOutput('motifs')
    )
  )),

  # div(h4("Generated Plots"), HTML('&nbsp;'), downloadButton('dl_everything', "Download all")),
  div(h4("Generated Plots"), HTML('&nbsp;'), uiOutput("download_button")),
  tabsetPanel(id = "render_results", type = "tabs")
)))