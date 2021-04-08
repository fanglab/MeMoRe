library("shiny")
library("data.table")
library("R.utils")
library("Biostrings")
library("stringr")
library("parallel")
library("ggplot2")
library("gridExtra")
library("stringi")
library("DT")
library("pryr")
library("dplyr")

# call functions necessary for analysis
source("global.R", keep.source=TRUE)

options(shiny.maxRequestSize=1000*1024^2) 

oldmodFile <<- NULL
oldgenFile <<- NULL

# TODO clean up empty raw
# TODO properly pass variable
# TODO try moving function out again
function(input, output, session) {
  v <- reactiveValues(modFile=NULL, genFile=NULL, motFile=NULL, motiF=NULL, centeR=NULL, modType=NULL)

  newcols <- c("Motifs","Modified position","Type","% motifs detected","# motifs in genome","Partner motif","Mean Score","Mean IPD ratio","Mean Coverage","Plots Generated")
  newcols_dl <- c("Motifs","Modified position","Type")
  
  initialize.motif.summary <- function(){
    # Initialize motif summary
    v$df <- setNames(data.table(matrix(nrow = 0, ncol = length(newcols))), newcols)
  }

  reformat.motif.summary <- function(motif_summary){
    motif_summary <- motif_summary %>%
      mutate(`Modified position`=`Modified position` + 1) %>%
      mutate(Type=ifelse(is.na(Type),paste0("x",substr(Motifs,`Modified position`,`Modified position`)),Type)) %>%
      mutate(`% motifs detected`=format(round(as.numeric(`% motifs detected`),2), nsmall=2)) %>%
      mutate(`Mean Score`=round(as.numeric(`Mean Score`))) %>%
      mutate(`Mean IPD ratio`=format(round(as.numeric(`Mean IPD ratio`),2), nsmall=2)) %>%
      mutate(`Mean Coverage`=round(as.numeric(`Mean Coverage`)))

    return(motif_summary)
  }
  
  read.motif.summary <- function(inFile){
    if(!is.null(inFile)){
      list_motif_summary_cols <- c("motifString","centerPos","modificationType","fraction","nGenome","partnerMotifString","meanScore","meanIpdRatio","meanCoverage")
      motif_summary <- read.table(inFile$datapath, sep=",", header=TRUE) # Read motif_summary.csv file
      motif_summary <- motif_summary[,colnames(motif_summary) %in% list_motif_summary_cols] # Select useful columns

      # Add missing columns      
      missing_columns <- setdiff(list_motif_summary_cols, names(motif_summary))
      # Check that necessary columns are here
      if(any(missing_columns %in% list_motif_summary_cols[c(1,2)])){
        showNotification("At least one column in motifString & centerPos is missing.", type="error")

        return(NULL)
      }
      motif_summary[missing_columns] <- NA
      motif_summary <- motif_summary[list_motif_summary_cols] # Reorder columns

      # Add plotting status column
      motif_summary$plotted <- "No"
      motif_summary$plotted[which(motif_summary$motifString %in% v$dldf$Motifs)] <- "Yes"

      colnames(motif_summary) <- newcols # Rename columns with clean names

      motif_summary <- reformat.motif.summary(motif_summary)

      # Remove duplicate entries
      new_motifs <- paste0(motif_summary$`Motifs`,"_",motif_summary$`Modified position`)
      isolate(current_motifs <- paste0(v$df$`Motifs`,"_",v$df$`Modified position`))
      motif_summary <- motif_summary[!new_motifs %in% current_motifs,]

      # Update main reactive object without re-execution
      isolate(v$df <- rbind(v$df, motif_summary))
    }
  }

  initialize.motif.summary()
  v$dldf <- setNames(data.table(matrix(nrow = 0, ncol = length(newcols_dl))), newcols_dl)

  # upload motiffile to motiftable
  observeEvent(input$motfile, {
    inFile <- input$motfile
  
    # Read motif_summary.csv file and populate table
    read.motif.summary(inFile)
  })
  
  # clear table if motif_summary.csv not selected, else replace table with motif_summary.csv
  observeEvent(input$cleartable, {
    initialize.motif.summary()

    inFile <- input$motfile

    # Read motif_summary.csv file and populate table
    read.motif.summary(inFile)
  })
  # div(DT::dataTableOutput("motiftable"), style = "font-size: 75%; width: 75%")

  # render DT
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable({v$df}, options=list(iDisplayLength=8, bLengthChange=0))
    dataTableOutput('motiftable')
  })
  
  # render dl DT
  output$dl_dt <- renderUI({
    output$dl_dt2 <- renderDataTable({v$dldf}, selection='single', options=list(iDisplayLength=5, bLengthChange=0, bFilter=0, bInfo=0, bAutoWidth=0))
    dataTableOutput('dl_dt2')
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
      rowadd <- setNames(data.table(motif_to_add, center_to_add, modtype, "", "", "", "", "", "", "No"), newcols)
      isolate(v$df <- rbind(v$df, rowadd))    
    }
  })
  
  observeEvent(input$submit, {
    # modifications.csv.gz check
    if(is.null(input$modfile)){
      showNotification("Upload modifications.csv(.gz)", type="error")

      return(NULL)
    }
    # genome.fasta check
    if(is.null(input$genfile)){
      showNotification("Upload genome.fasta", type="error")

      return(NULL)
    }
    # modifications.csv.gz check
    if(is.null(input$motfile)){
      showNotification("Upload motif_summary.csv", type="error")

      return(NULL)
    }
    # row motif selected check
    if(is.null(input$motiftable_rows_selected)){
      showNotification("Select a motif.", type="error")

      return(NULL)
    }
    
    # DATA INPUT
    v$modFile <- input$modfile$datapath
    v$genFile <- input$genfile$datapath
    v$motFile <- input$motfile$datapath
    
    v$motiF <- NULL
    v$centeR <- NULL
    v$modType <- NULL
    
    for(k in 1:length(input$motiftable_rows_selected)){
      v$motiF[k] <- toString(v$df$Motifs[input$motiftable_rows_selected[k]])
      v$centeR[k] <- as.numeric(v$df$"Modified position"[input$motiftable_rows_selected[k]])
      v$modType[k] <- toString(v$df$Type[input$motiftable_rows_selected[k]])
    }
    
    if(any(v$motiF %in% v$dldf$Motifs)){
      showNotification("Plots for at least one of the selected motifs were already generated.", type="warning")
    }else{ # if not add to v$dldf
      
      for(j in 1:length(input$motiftable_rows_selected)){
        # Progress Bar (most time intensive)
        print("Beginning Motif Loop")
        print(mem_used())
        withProgress(message = 'Making plots', value = 0, {
          # test if uploaded genome and modifications files are changed
          isemptyorchanged = (is.null(oldmodFile) & is.null(oldgenFile)) || !(oldmodFile == v$modFile & oldgenFile == v$genFile)
          
          # If changed, reupload
          if(isemptyorchanged) {
            incProgress(.2, detail = "Uploading Files")
            uploaddat(v$modFile, v$genFile, v$motiF[j], v$centeR[j] - 1)
            memuse("After upload dat")          
          }
          
          # Continue to process files
          incProgress(.2, detail = "Processing Files")
          graphs <- processdat(v$motiF[j], v$centeR[j] - 1, v$modType[j])
          memuse("After process dat") 
        })
        
        # GRAPHS
        assign(paste0(v$motiF[j], ".gpa"), graphs$ga, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".gps"), graphs$gs, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".gpi"), graphs$gi, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".gpc"), graphs$gc, envir = .GlobalEnv)
        
        # Below irrelevant rn
        assign(paste0(v$motiF[j], ".mcount"), graphs$mc, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".mscore"), graphs$ms, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".mipd"), graphs$mi, envir = .GlobalEnv)
        assign(paste0(v$motiF[j], ".mcov"), graphs$mco, envir = .GlobalEnv)
        
        rowadd <- setNames(data.table(v$motiF[j], v$centeR[j], v$modType[j]), newcols_dl)
        isolate(v$dldf <- rbind(v$dldf, rowadd))
        
        isolate(v$df$"Plots Generated" <- as.character(v$df$"Plots Generated"))
        isolate(v$df$"Plots Generated"[input$motiftable_rows_selected[j]] <- "Yes")
        isolate(v$df$"Plots Generated" <- as.factor(v$df$"Plots Generated"))
        
        oldmodFile <<- v$modFile
        oldgenFile <<- v$genFile

        memuse("After assigning variables")

        print(lsos())
      }
    }

  }) #end observe submit
  
  # observe selection in dl_dt
  observeEvent(input$dl_dt2_rows_selected, {
    motif_selected <- v$dldf$Motifs[input$dl_dt2_rows_selected]
    gpa <- get(paste0(motif_selected, ".gpa"))
    gps <- get(paste0(motif_selected, ".gps"))
    gpi <- get(paste0(motif_selected, ".gpi"))
    gpc <- get(paste0(motif_selected, ".gpc"))
    
    output$combined <- renderImage({
      if ((is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) | length(input$dl_dt2_rows_selected) == 0) return()
      
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gpa, width=nchar(motif_selected)*2.8, height=9, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$score <- renderImage({
      if ((is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) | length(input$dl_dt2_rows_selected) == 0) return()
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gps, width=nchar(motif_selected)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$ipd <- renderImage({
      if ((is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) | length(input$dl_dt2_rows_selected) == 0) return()
      
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gpi, width=nchar(motif_selected)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$coverage <- renderImage({
      if ((is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) | length(input$dl_dt2_rows_selected) == 0) return()
      
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gpc, width=nchar(motif_selected)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    # Download buttons
    output$dl_a = downloadHandler(
      filename = function() { paste0(motif_selected, '_combined.pdf') },
      content = function(file) {
        ggsave(file, gpa, width=nchar(motif_selected)*2.8, height=8.5, limitsize = F, device="pdf")
      })
    output$dl_s = downloadHandler(
      filename = function() { paste0(motif_selected, '_score.pdf') },
      content = function(file) {
        ggsave(file, gps, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
      })
    output$dl_i = downloadHandler(
      filename = function() { paste0(motif_selected, '_ipdRatio.pdf') },
      content = function(file) {
        ggsave(file, gpi, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
      })
    output$dl_c = downloadHandler(
      filename = function() { paste0(motif_selected, '_coverage.pdf') },
      content = function(file) {
        ggsave(file, gpc, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
      })
    
    output$images <- renderUI({
      tabsetPanel(
        type = "tabs",
        tabPanel("Combined", br(),
                 downloadButton('dl_a'),br(),
                 imageOutput("combined"),
                 style = "overflow-y:scroll;"),
        tabPanel("Score", br(),
                 downloadButton('dl_s'),br(),
                 imageOutput("score"),
                 style = "overflow-y:scroll;"),
        tabPanel("ipdRatio", br(),
                 downloadButton('dl_i'),br(),
                 imageOutput("ipd"),
                 style = "overflow-y:scroll;"),
        tabPanel("Coverage", br(),
                 downloadButton('dl_c'),br(),
                 imageOutput("coverage"),
                 style = "overflow-y:scroll;")
      )
    })
  })
  
  output$updatenom = renderPrint({
    s = input$motiftable_rows_selected
    if (length(s)) {
      cat("These motifs were selected: \n\n")
      cat(paste(v$df$Motifs[s], collapse="\n "))
    }
  })
  
  output$dl_everything <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_motifplots.zip")
    },
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;
      
      withProgress(message = 'Zipping files', value = 0, {
        #loop through the plots
        for (i in 1:nrow(v$dldf)){
          motif_selected <- v$dldf$Motifs[i]
          
          gpa <- get(paste0(motif_selected, ".gpa"))
          gps <- get(paste0(motif_selected, ".gps"))
          gpi <- get(paste0(motif_selected, ".gpi"))
          gpc <- get(paste0(motif_selected, ".gpc"))
          gpa_f <- paste0(motif_selected, '_combined.pdf')
          gps_f <- paste0(motif_selected, '_score.pdf')
          gpi_f <- paste0(motif_selected, '_ipdRatio.pdf')
          gpc_f <- paste0(motif_selected, '_coverage.pdf')
          
          ggsave(gpa_f, gpa, width=nchar(motif_selected)*2.8, height=8.5, limitsize = F, device="pdf")
          ggsave(gps_f, gps, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
          ggsave(gpi_f, gpi, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
          ggsave(gpc_f, gpc, width=nchar(motif_selected)*2.8, height=3, limitsize = F, device="pdf")
          files <- c(gpa_f,gps_f,gpi_f,gpc_f,files)
          
          incProgress((1)/nrow(v$dldf), detail = paste("Downloading Plots for Motif", i))
        }
      })
      
      #create the zip file
      zip(file,files)
    }
  )
}