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

options(shiny.maxRequestSize=200*1024^2) 

uploaddat <- function(rmodfile, rgenfile, motif, center) {
  csv_temp <- gunzip(rmodfile, remove = T, temporary = T, overwrite = T)
  csv <<- fread(csv_temp, sep = ",", header = T, verbose = F, drop = c(6, 7, 8, 11, 12, 13))  # Set csv parameters
  gene <<- readDNAStringSet(rgenfile)
  print(motif)
  print("BLAH")
  print(center)
  print(strsplit(motif,""))
  print(strsplit(motif,"")[[1]])
  print(as.integer(center)+1)
  print(strsplit(motif,"")[[1]][as.integer(center)+1])
  motiftable <<- data.table(motifString = motif, centerPos = center, modificationType =strsplit(motif,"")[[1]][as.integer(center)+1])
}

processdat <- function(motif, center, modificationtype) {
  print(motif)
  print(center)
  motiftable <<- data.table(motifString = motif, centerPos = center, modificationType =strsplit(motif,"")[[1]][as.integer(center)+1])
  motif_ref <- data.frame(sym = c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"), bases = c("(A|T)", "(C|G)", "(A|C)", "(G|T)", "(A|G)", "(C|T)", "(C|G|T)", "(A|G|T)", "(A|C|T)", "(A|C|G)", "(A|C|G|T)"))
  cl_max <- 8
  modtype <- as.vector(motiftable$modificationType)
  methvec <- gsub("\\d|[[:lower:]]","",modtype)
  motif_f <- as.vector(motiftable$motifString)
  cov_cut_spec <- 10
  size_spec <- 100
  
  gag <- which(width(gene) == max(width(gene)))
  csv2 <- csv[csv$refName == names(gene)[gag]] #csv2 = csv
  
  # R and F prep
  csv2 <- csv2[complete.cases(csv2), ]
  csv2 <- csv2[order(csv2$tpl),]
  param <- c(csv2$tpl[1], csv2$tpl[nrow(csv2)])  # create reverse genome
  genome_f <- toString(gene[[gag]])
  genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
  csv_r <- csv2[(csv2$strand == 1), ]  # separate csv r by strandness
  csv_f <- csv2[(csv2$strand == 0), ]  # separate csv f by strandness
  # cluster prep
  csv_sp <- csv[csv$coverage >= cov_cut_spec, ]
  csv_sp_f <- csv_f[csv_f$coverage >= cov_cut_spec, ]
  csv_sp_r <- csv_r[csv_r$coverage >= cov_cut_spec, ]
  
  #exportme2 <- list(data.table(csv_sp_f), data.table(csv_sp_r), genome_f)
  
  # create cluster
  #cat(paste0("  + Making cluster of ", cl_max, " cores"), "\n")
  #cat("   + Exporting global variables/libraries", "\r")
  # cl <- makeCluster(cl_max)
  # clusterExport(cl, c("exportme2"))
  # clusterEvalQ(cl, {
  #   library(stringi)
  #   library(IRanges)
  #   library(data.table)
  # })
  
  for (j in 1:length(motif_f)){ 
    new_mots <- nchar(motif_f[j])
    nom_matrix2 <- matrix(seq(1, 4*new_mots), nrow = 4, ncol = new_mots)
    meth <- methvec[j]
    transl_r <- transl_f <- as.integer(motiftable$centerPos[j])
    
    #Get new_list_o
    new_list_o <- list()
    new_list_o2 <- list()
    new_mot_ref <- rbind(motif_ref, data.frame(sym = c("A", "C", "G", "T"), bases = c("A", "C", "G", "T")))
    for (q in 1:nchar(motif_f[j])) {
      cat <- 1:nchar(motif_f[j])
      x <- motif_f[j]
      pat2 <- paste0("(", paste(as.vector(new_mot_ref$sym), collapse = "|"), ")")
      rep <- as.vector(new_mot_ref$bases)
      g <- (1:nchar(x))[cat[!cat %in% q]]  #gregexpr(pat2, x)[[1]][cat[!cat %in% q]]
      cat2 <- as.vector(new_mot_ref$sym)
      k <- unlist(strsplit(x, ""))
      for (a in 1:length(g)) {
        k[g[a]] <- rep[match(substr(x, g[a], g[a]), cat2)]
      }
      new_mots2 <- paste0(k, collapse = "")
      num_deg2 <- c("A", "C", "G", "T")  #unlist(strsplit(num_deg[q],'[|]'))
      for (s in 1:length(num_deg2)) {
        new_mots2.a <- strsplit(new_mots2, "(?<=\\))|(?<=[[:alpha:]])(?=[[:alpha:]\\(])", perl = TRUE)[[1]]
        new_mots2.a[q] <- num_deg2[s]
        new_mots3 <- paste0(new_mots2.a, collapse = "")
        mottleme <- unlist(strsplit(motif_f[j],""))
        mottleme[q] <- num_deg2[s]
        mottleme2 <- paste(mottleme,collapse="")
        nom_brak <- str_count(new_mots3, "[\\(]")
        nom_brak_ref <- str_count(motif_f[j], paste0("(", paste(as.vector(motif_ref$sym), collapse = "|"), ")"))
        # NEW DATA DATA FOR LISTED DATAFRAME new_listo <- subset(csv_mot_listo, grepl(new_mots3, motif))
        o_motif <- new_mots3
        o_motif_n <- motif_f[j]
        o_pos <- gregexpr(pat2, x)[[1]][cat[cat %in% q]]
        o_deg_y.or.n <- ifelse(nom_brak < nom_brak_ref, "Y", "N")
        new_list_o[[nom_matrix2[s, q]]] <- data.frame(motif_ = mottleme2, motif = o_motif, motif_n = o_motif_n,
                                                      pos = o_pos, deg_y.or.n = o_deg_y.or.n)
      }
    }
    new_list_o <- lapply(new_list_o, function(z) replace(z, is.na(z), 0))  # replace NaN's with 0
    
    #PREPARING DATA
    dog1.me <- do.call(rbind, new_list_o)  # all checked at this point # count>=5 & mean_r<1.5 & -log10 p <3
    dog1.me2 <- dog1.me[1]
    dog1 <- dog1.me[2:ncol(dog1.me)]          ## NEW DATA prep CSV
    csv_sp_2 <- csv_sp[order(csv_sp$tpl)][complete.cases(csv_sp[order(csv_sp$tpl)]), ]
    param_sp <- c(csv_sp_2$tpl[1], csv_sp_2$tpl[nrow(csv_sp_2)])
    dat2.A.1 <- csv_sp[csv_sp$base == meth, ]
    dat2.A <- dat2.A.1[sample(nrow(dat2.A.1), size_spec), ]
    
    exportme <- list(transl_f, size_spec, param_sp)
    ## PARALLELIZE
    #clstart1 <- Sys.time()
    #clusterExport(cl, c("exportme"))
    #clend.a <- Sys.time()
    #cat("   + Merging csv and dat", "\r")
    
    ######
    dat2 <- list()
    for(z.count in 1:length(dog1$motif)) {
      z <- dog1$motif[z.count]
      ar <- paste(rev(stri_extract_all(z, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
      arya <- unlist(stri_extract_all_regex(genome_f, z))[1]
      f_count <- stri_count_regex(genome_f, z)
      r_count <- stri_count_regex(genome_r, ar)
      locations_f2 <- data.table(tpl = as.data.table(stri_locate_all_regex(genome_f, z))$start + transl_f,
                                 motif = if (f_count == 0) 0 else rep(arya, f_count))
      locations_r2 <- data.table(tpl = as.data.table(stri_locate_all_regex(genome_r, ar))$end - transl_r,
                                 motif = if (r_count == 0) 0 else rep(arya, r_count))
      loc_f_s2 <- locations_f2[!(locations_f2$tpl < param_sp[1] | locations_f2$tpl > param_sp[2]), ]
      loc_r_s2 <- locations_r2[!(locations_r2$tpl < param_sp[1] | locations_r2$tpl > param_sp[2]), ]
      merged_f2 <- merge(loc_f_s2, csv_sp_f, all = F)
      merged_r2 <- merge(loc_r_s2, csv_sp_r, all = F)
      dat2[[z.count]] <- data.table(rbind(merged_f2, merged_r2))
      incProgress(.6/length(dog1$motif), detail = paste("Analyzing Motif", z.count))
      print(z.count)
    }
    
    # mcount, mscore, mipd, mcov
    mcount = NA
    mscore = NA
    mipd = NA
    mcov = NA
    
    ##################################### BOTTOM IN R/RSHINY, TOP IN C ################################33
    
    # WRITING TO DEBUG_GRAPH 2 NEW
    bases <- c("A","C","G","T")
    scorelist <- list()
    ipdlist <- list()
    coveragelist <- list()
    countlist <- list()
    
    for(newy in 1:nchar(motif_f[j])){
      for(newy2 in 1:4){
        namepos <- dat2[[((newy-1)*4+newy2)]] #paste0(dir,"/data/pos", sprintf("%02d",newy), locpos, ".csv")
        dataworky10 <- namepos #read.csv(namepos, sep = ",", header = T)
        scorelist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                      pos = newy,
                                                      mod = 0, #Added
                                                      values = c(dataworky10$score),
                                                      type = c(rep("score", nrow(dataworky10))))
        ipdlist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                    pos = newy,
                                                    mod = 0, #Added
                                                    values = c(dataworky10$ipdRatio),
                                                    type = c(rep("ipdRatio", nrow(dataworky10))))
        coveragelist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                         pos = newy,
                                                         mod = 0, #Added
                                                         values = c(dataworky10$coverage),
                                                         type = c(rep("coverage", nrow(dataworky10))))
        countlist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                      pos = newy,
                                                      labels = nrow(dataworky10))
      }
    }
    
    datascore <- do.call(rbind, scorelist)
    dataipd <- do.call(rbind, ipdlist)
    datacoverage <- do.call(rbind, coveragelist)
    datacount <- do.call(rbind, countlist)
    
    prettify_base <- theme(
      panel.background=element_rect(fill = NA,color="gray"), 
      panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(size=.1, color="black",linetype="dotted"), 
      panel.grid.minor.y=element_blank(),
      panel.grid.minor.x=element_line(size=.1, color="black"),
      legend.position="none"
    )
    
    prettify_top <- prettify_base +
      theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
      )
    
    prettify_mid <- prettify_top +
      theme(
        strip.background=element_blank(),
        strip.text.x=element_blank()
      )
    
    prettify_btm <- prettify_base +
      theme(
        strip.background=element_blank(),
        strip.text.x=element_blank()
      )
    
    datascore <<- datascore
    dataipd <- dataipd
    datacoverage <- datacoverage
    graphtitle <- paste0(substr(motif_f[j], 1, center), modificationtype, substr(motif_f[j], center+2, nchar(motif_f[j]))) 
    print(graphtitle)
    print(datascore)
    print(dataipd)
    print(datacoverage)
    
    # Mark expected nucleotide 
    motif_split <- unlist(strsplit(motif_f[j],""))
    indx <- match(motif_split, motif_ref$sym, nomatch=0) # Added no match --> as is
    list_pattern <- as.character(motif_ref$bases[match(motif_split,motif_ref$sym)])
    list_pattern[indx == 0] <- motif_split[indx == 0]
    for(idx_motif in 1:length(list_pattern)){
      opt_nuc <- unlist(strsplit(gsub(" ","",chartr("(|)", "   ", list_pattern[idx_motif])),""))
      #separated for faster data processing
      datascore$mod[datascore$pos==idx_motif & datascore$base %in% opt_nuc] <- 1
      dataipd$mod[dataipd$pos==idx_motif & dataipd$base %in% opt_nuc] <- 1
      datacoverage$mod[datacoverage$pos==idx_motif & datacoverage$base %in% opt_nuc] <- 1
    }
    
    # Combined Graphs Below
    
    gp_cov <- ggplot(datacoverage) + 
      geom_violin(aes(x=base, y=values, color="red", group=base)) +
      geom_text(data=datacount, aes(x=base, y=10, label=labels), size=5, colour="grey25") + 
      facet_grid(. ~ pos) +
      scale_colour_manual(values=c("red"="#F8766D")) +
      labs(y="Coverage") +
      prettify_top
    
    gp_ipd <- ggplot(dataipd) + 
      geom_violin(aes(x=base, y=values, color="green", group=base, fill=as.factor(mod))) +
      coord_cartesian(ylim=c(0,10)) +
      facet_grid(. ~ pos) +
      scale_colour_manual(values=c("green"="#00BA38")) +
      scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
      labs(y="IPD ratio") +
      prettify_mid
    
    gp_scr <- ggplot(datascore) + 
      geom_violin(aes(x=base, y=values, color="blue", group=base, fill=as.factor(mod))) +
      facet_grid(. ~ pos) +
      scale_colour_manual(values=c("blue"="#619CCF")) +
      scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="Score") +
      prettify_btm
    
    gp_cov <- ggplotGrob(gp_cov)
    gp_ipd <- ggplotGrob(gp_ipd) # Margin too large, see rbind below
    gp_scr <- ggplotGrob(gp_scr)
    
    graphcombined <- arrangeGrob(rbind(gp_cov, gp_ipd, gp_scr), ncol=1, top=graphtitle)
    
    # Separate Graphs Below
    print("Separate Graphs Below")
    
    graphscore <- ggplot(datascore) + 
      geom_violin(aes(x=base, y=values,color = "blue", group=base, fill=as.factor(mod))) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="score") +
      scale_colour_manual(values=c("blue"="#619CCF")) +
      scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
      prettify_base
    
    graphipd <- ggplot(dataipd) + 
      geom_violin(aes(x=base, y=values,color = "green", group=base, fill=as.factor(mod))) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="ipdRatio") +
      scale_colour_manual(values=c("green"="#00BA38")) +
      scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
      prettify_base
    
    graphcoverage <- ggplot(datacoverage) + 
      geom_violin(aes(x=base, y=values,color = "red", group=base)) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="coverage") +
      scale_colour_manual(values=c("red"="#F8766D")) +
      prettify_base
  
    return(list(ga = graphcombined, gs = graphscore, gi = graphipd, gc = graphcoverage, 
                mc = mcount, ms = mscore, mi = mipd, mco = mcov))
    
  }
}

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
    isuploaded = F
    if (!is.null(inFile) & isuploaded == F){
      isuploaded = T
      print("YES")
      d <- read.table(inFile$datapath, sep = ",", header = TRUE)[c(1:4,6,8:11)]
      isolate(v$df <- rbind(v$df, setnames(d, old = colnames(d), new = newcols)))
      print(v$df)
    } 
  })

  #for x values
  output$motifs <- renderUI({
    output$motiftable <- renderDataTable(v$df, selection = 'single')
    dataTableOutput('motiftable')
    #selectInput("motifdropdown", "Choose a motif from motif_summary.csv", outMotifs()$motifString)
  })
  
  observeEvent(input$cleartable, {
    v$df <- setNames(data.table(matrix(nrow = 0, ncol = 9)), newcols)
  })
  
  observeEvent(input$addmotif, {
    # Update Motiftablet
    # ## newcols = c("Motifs","Modified position","Type","% motifs detected","# motifs in genome","Partner motif","Mean Score","Mean IPD ratio","Mean Coverage")
    num = as.numeric(input$center)
    if(grepl('^[A-Za-z]+$', input$motif)){
      motif_to_add <- input$motif
    }else{
      showNotification("Enter valid motif.", type="error")
      return(NULL)
    }
    if(!is.na(num) & num > 0 & num <= nchar(input$motif)){
      center_to_add <- input$center
    }else{
      showNotification("Enter valid center position.", type="error")
      return(NULL)
    }
    
    # check if it already exists
    if(motif_to_add %in% v$df$Motifs){
      showNotification("Motif already in table.", type="error")
    }
    else{
      rowadd <- setNames(data.table(motif_to_add, center_to_add, "", "", "", "", "", "", ""), newcols)
      print(rowadd)
      isolate(v$df <- rbind(v$df, rowadd))
    }
  })

  observeEvent(input$submit, {
    # DATA INPUT
    v$modFile <- input$modfile$datapath
    v$genFile <- input$genfile$datapath
    v$motFile <- input$motfile$datapath
    v$motiF <- toString(v$df$Motifs[input$motiftable_rows_selected])
    v$centeR <- as.numeric(v$df$"Modified position"[input$motiftable_rows_selected])
    v$modType <- toString(v$df$Type[input$motiftable_rows_selected])
    
    print(v$modFile)
    print(v$genFile)
    print(v$motiF)
    print(v$centeR)
    print(v$modType)
    # print(oldmodFile)
    # print(oldgenFile)
    
    withProgress(message = 'Making plots', value = 0, {
      # UPLOAD
      testifemptyorchanged = (
        is.null(oldmodFile) & is.null(oldgenFile)) || 
        !(oldmodFile == v$modFile & oldgenFile == v$genFile)
      print(testifemptyorchanged)
      if(testifemptyorchanged) {
        print("preupload")
        incProgress(.2, detail = "Uploading Files")
        uploaddat(v$modFile, v$genFile, v$motiF, v$centeR)
        print(csv)
        print(gene)
        print(motiftable)
      }
      
      # PROCESSING FILES
      incProgress(.2, detail = "Processing Files")
      graphs <- processdat(v$motiF, v$centeR, v$modType)
    })
    
    # GRAPHS
    gpa <<- graphs$ga
    gps <<- graphs$gs
    gpi <<- graphs$gi
    gpc <<- graphs$gc
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
      if (is.null(v$modFile) && is.null(v$genFile)) return()
      
      outfile <- tempfile(fileext='.png')
      ggsave(outfile, gpa, width=nchar(v$motiF)*2.8, height=9, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$score <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile)) return()
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gps, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    print("WHAT")
    
    output$ipd <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile)) return()
      
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gpi, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
    
    output$coverage <- renderImage({
      if (is.null(v$modFile) && is.null(v$genFile)) return()
      
      outfile <- tempfile(fileext='.png')
        ggsave(outfile, gpc, width=nchar(v$motiF)*2.8, height=3, limitsize = F)
      
      list(src = outfile,
           contentType = 'image/png')
    }, deleteFile = TRUE)
  }) #end observe submit
  
  observeEvent(input$reset, {
    v$modFile <- NULL
    v$genFile <- NULL
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