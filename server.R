# PICK UP FROM LINES 33-48 TRYING TO AUTOMATE MODIFICATIOn TYPE WITHOUT HAVING USER TYPE IT IN

library(shiny)
library(data.table)
library(R.utils)
library("Biostrings")
library(stringr)
library(parallel)
library(ggplot2)
library(gridExtra)
library(stringi)
library(DT)

# call functions necessary for analysis
source("data.R")

options(shiny.maxRequestSize=200*1024^2) 

oldmodFile <<- NULL
oldgenFile <<- NULL

function(input, output, session) {
  v <- reactiveValues(modFile=NULL, 
                      genFile=NULL, 
                      motFile=NULL,
                      motiF=NULL, 
                      centeR=NULL,
                      modType=NULL)
  newcols = c("Motifs","Modified position","Type","% motifs detected","# motifs in genome","Partner motif","Mean Score","Mean IPD ratio","Mean Coverage")
  
  v$df <- setNames(data.table(matrix(nrow = 0, ncol = 9)), newcols)
  
  outMotifs <- observeEvent(input$motfile, {
    inFile = input$motfile
    if (!is.null(inFile)){
      d <- read.table(inFile$datapath, sep = ",", header = TRUE)[c(1:4,6,8:11)]
      for (i in 1:nrow(d)) {
        print(i)
        modiftype <- d$modificationType[i]
        checkC <- ifelse(grepl("\\d", modiftype), gsub('.*-([0-9]+).*','\\1',modiftype), "X")
        checkpos <- as.numeric(d$centerPos[i]) + 1
        print(paste0(checkC,"m", substring(d$motifString[i], checkpos, checkpos)))
        d$modificationType[i] <- paste0(checkC,"m", substring(d$motifString[i], checkpos, checkpos))
      }
      isolate(v$df <- rbind(v$df, setnames(d, old = colnames(d), new = newcols)))
      print(v$df)
    } 
  })

  # render DT
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable(v$df, selection = 'single')
    dataTableOutput('motiftable')
  })
  
  # clear table if motif_summary.csv not selected, else replace table with motif_summary.csv
  observeEvent(input$cleartable, {
    inFile = input$motfile
    if (!is.null(inFile) & isuploaded == F){
      d <- read.table(inFile$datapath, sep = ",", header = TRUE)[c(1:4,6,8:11)]
      isolate(v$df <- rbind(v$df, setnames(d, old = colnames(d), new = newcols)))
      print(v$df)
    }else{
      v$df <- setNames(data.table(matrix(nrow = 0, ncol = 9)), newcols)
    } 
  })
  
  observeEvent(input$addmotif, {
    # Update Motiftablet
    num = as.numeric(input$center)
    
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
    if(motif_to_add %in% v$df$Motifs){
      showNotification("Motif already in table.", type="error")
    }
    else{ # if not add to v$df
      # assume modification occurs at specified center
      

      rowadd <- setNames(data.table(motif_to_add, center_to_add, "", "", "", "", "", "", ""), newcols)
      print(rowadd)
      isolate(v$df <- rbind(v$df, rowadd))
    }
  })

  observeEvent(input$submit, {
    # modifications.csv.gz check
    if(is.null(input$modfile)){
      showNotification("Upload modifications.csv.gz", type="error")
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
    v$motiF <- toString(v$df$Motifs[input$motiftable_rows_selected])
    v$centeR <- as.numeric(v$df$"Modified position"[input$motiftable_rows_selected])
    v$modType <- toString(v$df$Type[input$motiftable_rows_selected])
    
    # Progress Bar (most time intensive)
    withProgress(message = 'Making plots', value = 0, {
      # test if uploaded genome and modifications files are changed
      testifemptyorchanged = (is.null(oldmodFile) & is.null(oldgenFile)) || !(oldmodFile == v$modFile & oldgenFile == v$genFile)
      
      # If changed, reupload
      if(testifemptyorchanged) {
        incProgress(.2, detail = "Uploading Files")
        uploaddat(v$modFile, v$genFile, v$motiF, v$centeR)
      }
      
      # Continue to process files
      incProgress(.2, detail = "Processing Files")
      graphs <- processdat(v$motiF, v$centeR, v$modType)
    })
    
    # GRAPHS
    gpa <<- graphs$ga
    gps <<- graphs$gs
    gpi <<- graphs$gi
    gpc <<- graphs$gc
    
    # Below irrelevant rn
    mcount <<- graphs$mc
    mscore <<- graphs$ms
    mipd <<- graphs$mi
    mcov <<- graphs$mco
    
    oldmodFile <<- v$modFile
    oldgenFile <<- v$genFile
    
    # Update Motiftablet
    # ## newcols = c("Motifs","Modified position","Type","% motifs detected","# motifs in genome","Partner motif","Mean Score","Mean IPD ratio","Mean Coverage")
    # rowadd <- setNames(data.table(v$motiF, v$centeR, v$modType, "", mcount, "", mscore, mipd, mcov), newcols)
    # print(rowadd)
    # outMotifs = rbind(outMotifs(), rowadd)
    
    output$combined <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) return()
      
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gpa, width=nchar(v$motiF)*2.8, height=9, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$score <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) return()
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gps, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    print("WHAT")
    
    output$ipd <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) return()
      
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gpi, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$coverage <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile) && is.null(v$motFile)) return()
      
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gpc, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
  }) #end observe submit
  
  observeEvent(input$reset, {
    v$modFile <- NULL
    v$genFile <- NULL
    v$motFile <- NULL
  })  
  
  # Download buttons
  output$dl_a = downloadHandler(
    filename = function() { paste0(v$motiF, '_combined.pdf') },
    content = function(file) {
      ggsave(file, gpa, width=nchar(v$motiF)*2.8, height=8.5, limitsize = F, device="pdf")
    })
  output$dl_s = downloadHandler(
    filename = function() { paste0(v$motiF, '_score.pdf') },
    content = function(file) {
      ggsave(file, gps, width=nchar(v$motiF)*2.8, height=3, limitsize = F, device="pdf")
    })
  output$dl_i = downloadHandler(
    filename = function() { paste0(v$motiF, '_ipdRatio.pdf') },
    content = function(file) {
      ggsave(file, gpi, width=nchar(v$motiF)*2.8, height=3, limitsize = F, device="pdf")
    })
  output$dl_c = downloadHandler(
    filename = function() { paste0(v$motiF, '_coverage.pdf') },
    content = function(file) {
      ggsave(file, gpc, width=nchar(v$motiF)*2.8, height=3, limitsize = F, device="pdf")
    })
}