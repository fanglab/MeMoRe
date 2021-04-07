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

function(input, output, session) {
  v <- reactiveValues(modFile=NULL, genFile=NULL, motFile=NULL, motiF=NULL, centeR=NULL, modType=NULL)

  newcols <- c("Motifs","Modified position","Type","% motifs detected","# motifs in genome","Partner motif","Mean Score","Mean IPD ratio","Mean Coverage","Plots Generated")
  newcols_dl <- c("Motifs","Modified position","Type")
  
  v$df <- setNames(data.table(matrix(nrow = 0, ncol = length(newcols))), newcols)
  v$dldf <- setNames(data.table(matrix(nrow = 0, ncol = length(newcols_dl))), newcols_dl)
  
  refinedat <- function(x) {
    methnum <- gsub("[^0-9]","",x[3])
    modletter <- ifelse(methnum=="", "x", paste0("m",methnum))
    o_modnum <- as.numeric(x[2]) + 1
    o_type <- paste0(modletter, substring(x[1], o_modnum, o_modnum))
    o_fraction <- format(round(as.numeric(x[4]),2), nsmall=2)
    o_score <- round(as.numeric(x[7]))
    o_ipd <- format(round(as.numeric(x[8]),2), nsmall=2) 
    o_cov <- round(as.numeric(x[9]))
    list(o_modnum, o_type, o_fraction, o_score, o_ipd, o_cov)
  }
  
  updatetable <- function() {
    source <- t(as.data.table(apply(v$df, 1, refinedat)))
    isolate(v$df$"Modified position" <- unlist(source[,1]))
    isolate(v$df$Type <- unlist(source[,2]))
    isolate(v$df$"% motifs detected" <- unlist(source[,3]))
    isolate(v$df$"Mean Score" <- unlist(source[,4]))
    isolate(v$df$"Mean IPD ratio" <- unlist(source[,5]))
    isolate(v$df$"Mean Coverage" <- unlist(source[,6]))
  }

  read.motif.summary <- function(inFile){
    if (!is.null(inFile)){
      list_motif_summary_cols <- c("motifString","centerPos","modificationType","fraction","nGenome","partnerMotifString","meanScore","meanIpdRatio","meanCoverage")
      motif_summary <- read.table(inFile$datapath, sep=",", header=TRUE) # Read motif_summary.csv file
      motif_summary <- motif_summary[,colnames(motif_summary) %in% list_motif_summary_cols] # Select useful columns

      # Add missing columns      
      missing_columns <- setdiff(list_motif_summary_cols, names(motif_summary))
      # TODO Error if missing "motifString","centerPos"
      motif_summary[missing_columns] <- NA
      motif_summary <- motif_summary[list_motif_summary_cols] # Reorder columns

      motif_summary$plotted <- "No" # Add plotting status column
      colnames(motif_summary) <- newcols # Rename columns with clean names

      # Update main reactive object without re-execution
      isolate(v$df <- rbind(v$df, motif_summary) %>% distinct(across("Motifs","Modified position"), .keep_all=TRUE))
    }
  }

  # upload motiffile to motiftable
  outMotifs <- observeEvent(input$motfile, {
    inFile <- input$motfile
  
    # Read motif_summary.csv file
    read.motif.summary(inFile)

    # Update the table in the UI
    updatetable()
  })
  
  # render DT
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable(v$df, options=list(iDisplayLength=8, bLengthChange=0))
    dataTableOutput('motiftable')
  })
  
  # render dl DT
  output$dl_dt <- renderUI({
    output$dl_dt2 <- renderDataTable({v$dldf}, selection = 'single', options=list(iDisplayLength=5, bLengthChange=0, bFilter=0, bInfo=0, bAutoWidth=0))
    dataTableOutput('dl_dt2')
  })
  
  # clear table if motif_summary.csv not selected, else replace table with motif_summary.csv
  observeEvent(input$cleartable, {
    inFile = input$motfile
    if (!is.null(inFile)){
      d <- cbind(read.table(inFile$datapath, sep = ",", header = TRUE)[c(1:4,6,8:11)], "Plots Generated" = "No", stringsAsFactors = FALSE)
      d$"Plots Generated"[which(d$motifString %in% v$dldf$Motifs)] <- "Yes"
      
      
      isolate(v$df <- setnames(d, old = colnames(d), new = newcols))
      updatetable()
    }else{
      v$df <- setNames(data.table(matrix(nrow = 0, ncol = length(newcols))), newcols)
    }
  })
  
  observeEvent(input$addmotif, {
    # Update Motiftablet
    num <- as.numeric(input$center)
    modtypnum <- as.numeric(input$modtype)
    
    # check that motif contains only alpha
    if(grepl('^[A-Za-z]+$', input$motif)){
      motif_to_add <- input$motif
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
      showNotification("Motif already in table.", type="error")
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
      showNotification("Plots for one of the selected motifs were already generated.", type="error")
    }
    else{ # if not add to v$dldf
      
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