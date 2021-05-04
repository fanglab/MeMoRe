library(shiny)
library(data.table)
library(R.utils)
library(Biostrings)
library(stringr)
library(parallel)
library(ggplot2)
library(gridExtra)
library(stringi)
library(DT)
library(pryr)
library(dplyr)
library(foreach)
library(GenomicRanges)
library(doMC)

# call functions necessary for analysis
source("global.R", keep.source=TRUE)

options(shiny.maxRequestSize=1000*1024^2) 

# TODO include ONT
# TODO plot type in table
# TODO add tabs colors
# TODO change motif name in tab + code to full name
function(input, output, session) {
  v <- reactiveValues(modFile=NULL, genFile=NULL, motFile=NULL, motiF=NULL, centeR=NULL, modType=NULL)

  debug_mode <<- TRUE
  processing_version <<- 2

  rendered_motif <<- NULL
  loaded_motif <<- NULL
  downloadable <- reactiveValues(status=FALSE)

  oldmodFile <<- NULL
  oldgenFile <<- NULL

  initialize.motif.summary(v, list_motif_summary_clean_cols)
  # v$dldf <- setNames(data.table(matrix(nrow = 0, ncol = length(list_motif_manual_clean_cols))), list_motif_manual_clean_cols)

  if(debug_mode){
    intiale <- reactiveValues(datapath=NULL)
    intiale$datapath <- "Clostridium_perfringens_ATCC13124.motif_summary.csv"
    observeEvent(intiale$datapath, {
      read.motif.summary(v, intiale, list_motif_summary_clean_cols, list_motif_summary_cols)
    })
  }

  # Upload motiffile to motiftable
  observeEvent(input$motfile, {
    inFile <- input$motfile
  
    # Read motif_summary.csv file and populate table
    read.motif.summary(v, inFile, list_motif_summary_clean_cols, list_motif_summary_cols)
  }) # End observe motfile

  # Clear table if no motif_summary.csv input, or reload from motif_summary.csv
  observeEvent(input$cleartable, {
    initialize.motif.summary(v, list_motif_summary_clean_cols, list_motif_summary_cols)

    inFile <- input$motfile

    # Read motif_summary.csv file and populate table
    read.motif.summary(v, inFile, list_motif_summary_clean_cols)
  }) # End observe cleartable

  # Render motif summary table
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable({v$df}, options=list(iDisplayLength=8, bLengthChange=0))
    dataTableOutput('motiftable')
  })

  observeEvent(input$addmotif, {
    # Update Motiftable
    num <- as.numeric(input$center)
    # modtypnum <- as.numeric(input$modtype)

    # check that motif contains only alpha
    if(check.valid.motif(input$motif)){
      motif_to_add <- toupper(input$motif)
    }else{
      showNotification("Enter valid motif.", type="error")

      return(NULL)
    }
    # check that center is numeric
    if(!is.na(num) & num > 0 & num <= nchar(input$motif)){
      center_to_add <- input$center
    }else{
      showNotification("Enter valid center position.", type="error")

      return(NULL)
    }
    
    # check if motif already exists
    if(motif_to_add %in% v$df$Motifs && center_to_add %in% v$df[v$df$"Motifs" == motif_to_add]$"Modified position"){
      showNotification("Motif already in table.", type="warning")

      return(NULL)
    }else{
      # if not add to v$df
      # assume modification occurs at specified center
      modtype <- paste0("x", toupper(substring(motif_to_add, as.numeric(center_to_add), as.numeric(center_to_add))))
      rowadd <- setNames(data.table(motif_to_add, center_to_add, modtype, "", "", "", "", "", "", "No"), list_motif_summary_clean_cols)
      isolate(v$df <- rbind(v$df, rowadd))    
    }
  }) # End observe addmotif
  
  # Submit all motif to processing function
  observeEvent(input$submit_all, {
    start_time <- proc.time()

    results <- wrapper.data.processing(input, v, input$motiftable_rows_all, list_motif_manual_clean_cols)

    if(is.null(results)){
      # At least an input is missing

      return(NULL)
    }
  
    print.runtime.message(start_time)
  })

  # Submit selected motif to processing function
  observeEvent(input$submit, {
    start_time <- proc.time()

    results <- wrapper.data.processing(input, v, input$motiftable_rows_selected, list_motif_manual_clean_cols)

    if(is.null(results)){
      # At least an input is missing

      return(NULL)
    }
  
    print.runtime.message(start_time)
  })

  # Populate new tab and render new plots
  observeEvent(v$df, {
    processed_motifs <- v$df[v$df$"Plots Generated"=="Yes",]
    nb_processed_motifs <- nrow(processed_motifs)

    if(nb_processed_motifs==0){
      # No motif processed

      return(NULL)
    }else{
      for(idx_processed_motifs in 1:nb_processed_motifs){
        print_db(nb_processed_motifs)
        print_db(idx_processed_motifs)
        print_db(rendered_motif)
        current_motif <- processed_motifs$"Motifs"[idx_processed_motifs]
        current_modpos <- processed_motifs$"Modified position"[idx_processed_motifs]
        current_motif_detail <- paste0(current_motif,"_",current_modpos)
        print_db(current_motif_detail)

        first_render <- FALSE
        to_render <- FALSE
        if(is.null(rendered_motif)){ # Create tabsetPanel
          print_db(paste0("First render"))
          rendered_motif <<- current_motif_detail
          downloadable$status <- TRUE
          first_render <- TRUE
          to_render <- TRUE
        }else if(current_motif_detail %in% rendered_motif){
          # Motif already rendered
          print_db(paste0("Already render"))
        }else{
          print_db(paste0("New render"))
          rendered_motif <<- c(rendered_motif, current_motif_detail)
          to_render <- TRUE
        }
        if(to_render){
          appendTab(
            inputId = "render_results",
            tab = tabPanel(title=current_motif_detail, value=current_motif_detail, downloadButton(paste0('dl_',current_motif_detail)), br(), imageOutput(paste0('combined_',current_motif_detail)), style="overflow-y:scroll;"),
            select = first_render
          )
        }
        print_db(paste0(current_motif_detail," is done."))
      }
      print_db("All done!")
    }
  })

  observeEvent(input$render_results, {
    selected_motif <- input$render_results

    if(!selected_motif %in% loaded_motif){
      output[[paste0('combined_',selected_motif)]] <- renderImage({  
        print_db(paste0("Rendering combine: ", selected_motif))      

        path_graph_data <- paste0(selected_motif, ".gpa")
        gpa <- readRDS(file=path_graph_data)
      
        outfile <- tempfile(fileext='.png')
        ggsave(outfile, gpa, width=nchar(selected_motif)*2.8, height=9) # TODO remove added char
        
        list(src = outfile, contentType = 'image/png')
      }, deleteFile = TRUE)

      loaded_motif <<- c(loaded_motif, selected_motif)
    }
  })
  
  observeEvent(input$render_results ,{
    selected_motif <- input$render_results

    output[[paste0('dl_',selected_motif)]] <- downloadHandler(
      filename = function(){
        selected_motif <- input$render_results

        return(paste0(selected_motif, '_combined.pdf'))
      },
      content = function(file){
        selected_motif <- input$render_results

        path_graph_data <- paste0(selected_motif, ".gpa")
        gpa <- readRDS(file=path_graph_data)
      
        ggsave(file, gpa, width=nchar(selected_motif)*2.8, height=9, device="pdf")
      } # TODO remove added char
    )
  })

  output$list_selected_motifs <- renderPrint({
    list_selected_motifs <- input$motiftable_rows_selected
    if(!is.null(list_selected_motifs)){
      cat("Selected motif(s):\n")
      cat(paste0(v$df$Motifs[list_selected_motifs], collapse="\n"))
    }
  })

  output$download_button <- renderUI({
    if(downloadable$status){
      downloadButton("dl_everything", "Download all")
    }
  })

  # Generate .zip file with graphs from all processed motifs 
  output$dl_everything <- downloadHandler(
    filename = function(){
      return(paste0(Sys.Date(),"_motifplots.zip"))
    },
    content = function(file) {
      graph_output_dir <- getwd()
      owd <- setwd(dir=tempdir()) # Maybe not needed
      on.exit(setwd(dir=owd)) # Maybe not needed
      list_graph_files <- NULL
      
      processed_motifs <- v$df[v$df$"Plots Generated"=="Yes",]      
      if(nrow(processed_motifs)==0){
        showNotification(paste0("No motifs were processed."), type="warning")

        return(NULL)
      }

      withProgress(message='Zipping files', detail="Initialization", value=0, {
        for(idx_processed_motifs in 1:nrow(processed_motifs)){
          current_motif <- processed_motifs$"Motifs"[idx_processed_motifs]
          current_modpos <- processed_motifs$"Modified position"[idx_processed_motifs]
          current_motif_detail <- paste0(current_motif,"_",current_modpos)
          
          path_graph_data <- paste0(graph_output_dir, "/", current_motif_detail)

          graph_width <- nchar(current_motif_detail)*2.8
          graph_height_single <- 3
          graph_height_combined <- 8.5

          gpa_f <- paste0(current_motif_detail, '_combined.pdf')
          gp <- readRDS(file=paste0(path_graph_data, ".gpa"))
          ggsave(filename=gpa_f, plot=gp, width=graph_width, height=graph_height_combined, limitsize=FALSE, device="pdf")

          gps_f <- paste0(current_motif_detail, '_score.pdf')
          gp <- readRDS(file=paste0(path_graph_data, ".gps"))
          ggsave(filename=gps_f, plot=gp, width=graph_width, height=graph_height_single, limitsize=FALSE, device="pdf")

          gpi_f <- paste0(current_motif_detail, '_ipdRatio.pdf')
          gp <- readRDS(file=paste0(path_graph_data, ".gpi"))
          ggsave(filename=gpi_f, plot=gp, width=graph_width, height=graph_height_single, limitsize=FALSE, device="pdf")

          gpc_f <- paste0(current_motif_detail, '_coverage.pdf')
          gp <- readRDS(file=paste0(path_graph_data, ".gpc"))
          ggsave(filename=gpc_f, plot=gp, width=graph_width, height=graph_height_single, limitsize=FALSE, device="pdf")

          list_graph_files <- c(list_graph_files, gpa_f, gps_f, gpi_f, gpc_f)
          
          incProgress((1)/length(processed_motifs), detail=paste("Downloading Plots for Motif", idx_processed_motifs))
        }
      })
      
      # Create the zip from graph files
      zip(file, list_graph_files)
    }
  )

  session$onSessionEnded(function(){
    cat("Removing temporary files\n")
    list_graph_files <- list.files(pattern="*.gp[a|i|c|s]")
    print_db(list_graph_files)
    unlink(list_graph_files)
    cat("Session ended\n")
  })
}