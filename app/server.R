# Run locally: shiny::runApp('./app', host='0.0.0.0', port=3838)
# Deploy: library(rsconnect); rsconnect::deployApp(appName="MeMoRe", appDir='app/')
# Track ressources: ~/Library/Python/2.7/bin/psrecord $(pgrep -x R) --include-children --interval 0.1  --plot plot.png

app_version <<- "0.1.0_dev"

library(shiny)
# options(shiny.trace=TRUE)
options(shiny.maxRequestSize=1000*1024^2) 

# call functions necessary for analysis
source("global.R", keep.source=TRUE)

load.libraries()

# TODO add like tutorial
function(input, output, session) {
  v <- reactiveValues(modFile=NULL, genFile=NULL, motFile=NULL, motiF=NULL, centeR=NULL, modType=NULL)

  debug_mode <<- FALSE
  testing_mode <<- FALSE
  processing_version <<- 2

  initial <- reactiveValues(datapath=NULL)
  rendered_motif <<- NULL
  loaded_motif <<- NULL
  downloadable_motif <<- NULL
  downloadable <- reactiveValues(status=FALSE)
  process_all <<- 0
  load_testing_data <<- 0

  oldmodFile <<- NULL
  oldgenFile <<- NULL

  initialize.motif.summary(v, list_motif_summary_clean_SMRT_cols)

  if(debug_mode){
    # Use testing data directly
    initial$datapath <- "data/motif_summary.tsv"
  }
  
  # If testing with SMRT data
  observeEvent(input$smrt_test, {
    if(!file.exists("../iam_a_container")){
      testing_mode <<- TRUE # Override input checking
      initial$datapath <- NULL

      # Reset motif summary table
      initialize.motif.summary(v, list_motif_summary_clean_SMRT_cols)

      v$modFile <- "data/modification.smrt.csv.gz"
      v$genFile <- "data/reference.fasta"
      v$motFile <- "data/motif_summary.tsv"

      initial$datapath <- v$motFile

      load_testing_data <<- load_testing_data + 1 # Keep track of actionButton status
      # load_testing_data <<- 1 # Keep track of actionButton status

      # Hide inputs
      updateActionButton(session, "input_toggle", label="Show")
      shinyjs::hide(id="input_subpanel")
    }else{
      showNotification("Testing datasets are not available for local execution, please use the shiny web application.", type="warning")
    }
  })
  
  # If testing with ONT data
  observeEvent(input$ont_test, {
    if(!file.exists("../iam_a_container")){
      testing_mode <<- TRUE # Override input checking
      initial$datapath <- NULL

      # Reset motif summary table
      initialize.motif.summary(v, list_motif_summary_clean_SMRT_cols)

      v$modFile <- "data/modification.ont.rds"
      v$genFile <- "data/reference.fasta"
      v$motFile <- "data/motif_summary.tsv"
      
      initial$datapath <- v$motFile

      load_testing_data <<- load_testing_data + 1 # Keep track of actionButton status
      # load_testing_data <<- 1 # Keep track of actionButton status

      # Hide inputs
      updateActionButton(session, "input_toggle", label="Show")
      shinyjs::hide(id="input_subpanel")
    }else{
      showNotification("Testing datasets are not available for local execution, please use the shiny web application.", type="warning")
    }
  })

  # Initialize automated input processes
  observeEvent(initial$datapath, {
    # Read motif summary
    read.motif.summary(v, initial)
  })

  # Upload motiffile to motiftable
  observeEvent(input$motfile, {
    inFile <- input$motfile
  
    # Read motif_summary.csv file and populate table
    read.motif.summary(v, inFile)
  }) # End observe motfile

  # Clear table if no motif_summary.csv input, or reload from motif_summary.csv
  observeEvent(input$cleartable, {
    initialize.motif.summary(v, list_motif_summary_clean_SMRT_cols)

    inFile <- input$motfile

    # Read motif_summary.csv file and populate table
    read.motif.summary(v, inFile)
  })

  # Render motif summary table
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable({v$df}, options=list(iDisplayLength=8, bLengthChange=0, pageLength=8, dom='tip'))
    dataTableOutput('motiftable')
  })

  observeEvent(input$addmotif, {
    # Update Motiftable
    mod_pos <- as.integer(input$center)

    # check that motif contains only alpha
    if(check.valid.motif(input$motif)){
      motif_to_add <- toupper(input$motif)
    }else{
      showNotification("Enter valid motif.", type="error")

      return(NULL)
    }
    # check that center is numeric
    if(!is.na(mod_pos) & mod_pos > 0 & mod_pos <= nchar(input$motif)){
      center_to_add <- mod_pos
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
      # rowadd <- setNames(data.table(motif_to_add, center_to_add, modtype, "", "", "", "", "", "", "No"), list_motif_summary_clean_SMRT_cols)
      rowadd <- setNames(data.table(motif_to_add, center_to_add, modtype, "No"), list_motif_summary_clean_SMRT_cols)
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
  
    process_all <<- process_all + 1 # Keep track of actionButton status
    load_testing_data <<- 0 # Reset value to keep track of actionButton status

    # Hide inputs
    updateActionButton(session, "input_toggle", label="Show")
    shinyjs::hide(id="input_subpanel")

    # Hide motif summary
    updateActionButton(session, "motif_toggle", label="Show")
    shinyjs::hide(id="motif_subpanel")

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
    processed_motifs <- v$df[v$df$"Plots Generated"!="No",]
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
        current_mod_pos <- processed_motifs$"Modified position"[idx_processed_motifs]
        current_mod_type <- processed_motifs$"Type"[idx_processed_motifs]
        current_data_type <- processed_motifs$"Plots Generated"[idx_processed_motifs]
        current_motif_detail <- paste0(current_motif,"_",current_mod_pos)
        print_db(current_motif_detail)

        current_motif_clean <- paste0(substr(current_motif, 1, current_mod_pos-1),current_mod_type,substr(current_motif, current_mod_pos+1, nchar(current_motif)))

        if(current_data_type=="Both"){
          current_data_type <- c("ONT","SMRT")
        }

        first_render <- FALSE
        to_render <- FALSE
        if(is.null(rendered_motif)){ # Create tabsetPanel
          print_db(paste0("First render"))
          downloadable$status <- TRUE
          first_render <- TRUE
          to_render <- TRUE
        }else if(all(paste0(current_motif_detail,"_",current_data_type) %in% rendered_motif)){
          # All motifs already rendered
          print_db(paste0("Already render"))
        }else{
          print_db(paste0("New render"))
          to_render <- TRUE
        }
        if(to_render){
          for(selected_data_type in current_data_type){
            # Only render if not already done 
            if(!paste0(current_motif_detail,"_",selected_data_type) %in% rendered_motif){
              appendTab(
                inputId="render_results",
                tab=tabPanel(title=current_motif_clean, value=paste0(current_motif_detail,'_',selected_data_type),
                  downloadButton(paste0('dl_',current_motif_detail,'_',selected_data_type), label="", style='padding:4px; font-size:80%; margin-bottom:0.4em'),
                  br(),
                  imageOutput(paste0('combined_',current_motif_detail,'_',selected_data_type)),
                  style="overflow-y: visible;"
                ),
                select=first_render
              )

              if(is.null(rendered_motif)){
                rendered_motif <<- paste0(current_motif_detail,"_",selected_data_type)
              }else{
                rendered_motif <<- c(rendered_motif, paste0(current_motif_detail,"_",selected_data_type))
              }
            }
          }
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

        selected_motif_info <- str_split(selected_motif, pattern="_", simplify=TRUE)
        clean_selected_motif <- paste0(selected_motif_info[,c(1,2)], collapse="_")
        if(selected_motif_info[,3]=="SMRT"){
          path_graph_data <- paste0(clean_selected_motif, ".gpa")
        }else{
          path_graph_data <- paste0(clean_selected_motif, ".gpo")
        }

        withProgress(message=paste0("Rendering ",selected_motif_info[,1],"."), value=0, {
          incProgress(.05, detail="Read ggplot data.")
          # Read ggplot data
          gp <- readRDS(file=path_graph_data)

          incProgress(.5, detail="Create png file.")
          outfile <- tempfile(fileext='.png')
          ggsave(outfile, gp, width=nchar(str_split(clean_selected_motif, pattern="_", simplify=TRUE)[,1])*2.8, height=9) # TODO remove added char

          incProgress(.45, detail="Done.")
        })

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

        clean_selected_motif <- paste0(str_split(selected_motif, pattern="_", simplify=TRUE)[,c(1,2)], collapse="_")
        path_graph_data <- paste0(clean_selected_motif, ".gpa")
        if(file.exists(path_graph_data)){
          gp <- readRDS(file=path_graph_data)
        }else{
          path_graph_data <- paste0(clean_selected_motif, ".gpo")
          gp <- readRDS(file=path_graph_data)
        }

        ggsave(file, gp, width=nchar(str_split(clean_selected_motif, pattern="_", simplify=TRUE)[,1])*2.8, height=9, device="pdf")
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
      downloadButton("dl_everything", "All", style='padding:4px; font-size:80%; margin-bottom:0.4em')
    }
  })

  observeEvent(input$input_toggle, {
    if( ((input$input_toggle + load_testing_data + process_all) %% 2 == 0) ){
      updateActionButton(session, "input_toggle", label="Hide")
      shinyjs::show(id="input_subpanel")
    }else{
      updateActionButton(session, "input_toggle", label="Show")
      shinyjs::hide(id="input_subpanel")
    }
  })

  observeEvent(input$motif_toggle, {
    if(((input$motif_toggle + process_all) %% 2) == 0){
      updateActionButton(session, "motif_toggle", label="Hide")
      shinyjs::show(id="motif_subpanel")
      shinyjs::runjs("$('#motiftable table.dataTable[id]').DataTable().draw();")
    }else{
      updateActionButton(session, "motif_toggle", label="Show")
      shinyjs::hide(id="motif_subpanel")
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
      
      processed_motifs <- v$df[v$df$"Plots Generated"!="No",]      
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

          # ONT plot
          if(file.exists(paste0(path_graph_data, ".gpo"))){
            gpo_f <- paste0(current_motif_detail, '_ont.pdf')
            gp <- readRDS(file=paste0(path_graph_data, ".gpo"))
            ggsave(filename=gpo_f, plot=gp, width=graph_width, height=graph_height_combined, limitsize=FALSE, device="pdf")

            list_graph_files <- c(list_graph_files, gpo_f)
          }

          # SMRT plots
          if(file.exists(paste0(path_graph_data, ".gpa"))){
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
          }
          
          incProgress((1)/length(processed_motifs), detail=paste("Downloading Plots for Motif", idx_processed_motifs))
        }
      })
      
      # Create the zip from graph files
      zip(file, list_graph_files)
    }
  )

  output$app_version <- renderUI({
    tags$a(href="https://github.com/fanglab/MeMoRe", app_version, style="color:#333")
  })

  session$onSessionEnded(function(){
    cat("Removing temporary files\n")
    list_graph_files <- list.files(pattern="*.gp[a|i|c|s|o]")
    print_db(list_graph_files)
    unlink(list_graph_files)
    cat("Session ended\n")
  })
}
