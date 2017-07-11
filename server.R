library(shiny)
library(data.table)
library(R.utils)
library("Biostrings")
library(stringr)
library(parallel)
library(ggplot2)
library(gridExtra)
library(plotly)
library(stringi)

options(shiny.maxRequestSize=200*1024^2) 

uploaddat <- function(rmodfile, rgenfile, motif, center) {
  csv_temp <- gunzip(rmodfile, remove = T, temporary = T, overwrite = T)
  csv <<- fread(csv_temp, sep = ",", header = T, verbose = F, drop = c(6, 7, 8, 11, 12, 13))  # Set csv parameters
  gene <<- readDNAStringSet(rgenfile)
  motiftable <<- data.table(motifString = motif, centerPos = center, modificationType =strsplit(motif,"")[[1]][as.integer(center)+1])
}

processdat <- function(motif, center) {
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
    myreverse <- function(x,...){unlist(lapply(strsplit(as.vector(x),""),function(z)paste(rev(z),collapse="")))}
    dat2 <- list()
    for(z.count in 1:length(dog1$motif)) {
      z <- dog1$motif[z.count]
      ar <- paste(rev(stri_extract_all(z, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
      arya <- unlist(stri_extract_all_regex(genome_f, z))[1]
      locations_f2 <- data.table(tpl = as.data.table(stri_locate_all_regex(genome_f, z))$start + transl_f,
                                 motif = if (stri_count_regex(genome_f, z) == 0) 0 else rep(arya, stri_count_regex(genome_f, z)))
      locations_r2 <- data.table(tpl = as.data.table(stri_locate_all_regex(genome_r, ar))$end - transl_r,
                                 motif = if (stri_count_regex(genome_r, ar) == 0) 0 else rep(arya, stri_count_regex(genome_r, ar)))
      loc_f_s2 <- locations_f2[!(locations_f2$tpl < param_sp[1] | locations_f2$tpl > param_sp[2]), ]
      loc_r_s2 <- locations_r2[!(locations_r2$tpl < param_sp[1] | locations_r2$tpl > param_sp[2]), ]
      merged_f2 <- merge(loc_f_s2, csv_sp_f, all = F)
      merged_r2 <- merge(loc_r_s2, csv_sp_r, all = F)
      dat2[[z.count]] <- data.table(rbind(merged_f2, merged_r2))
      print(z.count)
    }

    # dat2 <- clusterApply(cl, dog1$motif, function(z) {
    #   # VARS
    #   csv_sp_f <- exportme2[[1]]
    #   csv_sp_r <- exportme2[[2]]
    #   genome_f <- exportme2[[3]]
    #   transl_f <- exportme[[1]]
    #   size_spec <- exportme[[2]]
    #   param_sp <- exportme[[3]]
    #   transl_r <- transl_f
    #   genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
    #   myreverse <- function(x,...){unlist(lapply(strsplit(as.vector(x),""),function(z)paste(rev(z),collapse="")))}
    #   ar <- paste(rev(stri_extract_all(z, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")
    #   arya <- unlist(stri_extract_all_regex(genome_f, z))[1]
    #   locations_f2 <- data.frame(tpl = as.data.table(stri_locate_all_regex(genome_f, z))$start + transl_f,
    #                             motif = rep(arya, stri_count_fixed(genome_f, z)))
    #   locations_r2 <- data.table(tpl = as.data.table(stri_locate_all_regex(genome_r, ar))$end - transl_r,
    #                             motif = rep(arya, stri_count_fixed(genome_r, ar)))
    #   loc_f_s2 <- locations_f2[!(locations_f2$tpl < param_sp[1] | locations_f2$tpl > param_sp[2]), ]
    #   loc_r_s2 <- locations_r2[!(locations_r2$tpl < param_sp[1] | locations_r2$tpl > param_sp[2]), ]
    #   merged_f2 <- merge(loc_f_s2, csv_sp_f, all = F)
    #   merged_r2 <- merge(loc_r_s2, csv_sp_r, all = F)
    #   dat2[[z.count]] <- data.table(rbind(merged_f2, merged_r2))
    # })
    
    ##################################### BOTTOM IN R/RSHINY, TOP IN C ################################33
    
    # WRITING TO DEBUG_GRAPH 2 NEW
    bases <- c("A","C","G","T")
    scorelist <- list()
    ipdlist <- list()
    coveragelist <- list()
    counts <- vector()
    
    for(newy in 1:nchar(motif_f[j])){
      for(newy2 in 1:4){
        namepos <- dat2[[((newy-1)*4+newy2)]] #paste0(dir,"/data/pos", sprintf("%02d",newy), locpos, ".csv")
        dataworky10 <- namepos #read.csv(namepos, sep = ",", header = T)
        scorelist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                      pos = newy,
                                                      values = c(dataworky10$score),
                                                      type = c(rep("score", nrow(dataworky10))))
        ipdlist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                    pos = newy,
                                                    values = c(dataworky10$ipdRatio),
                                                    type = c(rep("ipdRatio", nrow(dataworky10))))
        coveragelist[[((newy-1)*4+newy2)]] <- data.table(base = bases[newy2], 
                                                         pos = newy,
                                                         values = c(dataworky10$coverage),
                                                         type = c(rep("coverage", nrow(dataworky10))))
        counts[((newy-1)*4+newy2)] <- nrow(dataworky10)
      }
    }
    
    datascore <- do.call(rbind, scorelist)
    dataipd <- do.call(rbind, ipdlist)
    datacoverage <- do.call(rbind, coveragelist)
    
    #write.csv(datascore, paste0(dir,"/datascore.csv"), row.names = F)
    #write.csv(dataipd, paste0(dir,"/dataipd.csv"), row.names = F)
    #write.csv(datacoverage, paste0(dir,"/datacoverage.csv"), row.names = F)
    #write(paste0(dir,"/dataw.csv"), file = "graph_dirs.txt",append = TRUE)
    #write(motif_f[j], file = "graph_dirs.txt",append = TRUE)
    
    prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_line(size=.1, color="black"),
                      legend.position="none")
    
    datascore <- datascore
    dataipd <- dataipd
    datacoverage <- datacoverage
    
    graphscore <- ggplot(datascore, aes(x=base, y=values,color = base)) + 
      ggtitle(paste0(motif_f[j])) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="score") +
      geom_violin() + 
      facet_grid(.~pos) + 
      theme_gray() %+replace% prettify
    graphipd <- ggplot(dataipd, aes(x=base, y=values,color = base)) + 
      ggtitle(paste0(motif_f[j])) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="ipdRatio") +
      geom_violin() + 
      facet_grid(.~pos) + 
      theme_gray() %+replace% prettify
    graphcoverage <- ggplot(datacoverage, aes(x=base, y=values,color = base)) + 
      ggtitle(paste0(motif_f[j])) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="coverage") +
      geom_violin() + 
      facet_grid(.~pos) + 
      theme_gray() %+replace% prettify 
    
    gp_cov <- ggplotGrob(graphscore)
    gp_ipd <- ggplotGrob(graphipd) # Margin too large, see rbind below
    gp_scr <- ggplotGrob(graphcoverage)
    
    gp <- arrangeGrob(rbind(gp_cov, gp_ipd, gp_scr), ncol=1, top=motif_f[j])

    return(list(gs = graphscore, gi = graphipd, gc = graphcoverage))
    
  }
}

oldmodFile <<- NULL
oldgenFile <<- NULL

function(input, output, session) {
  v <- reactiveValues(modFile=NULL, genFile=NULL, motiF=NULL, centeR=NULL)

  
  observeEvent(input$submit, {
    v$modFile <- input$modfile$datapath
    v$genFile <- input$genfile$datapath
    v$motiF <- input$motif
    v$centeR <- input$center
    print(v$modFile)
    print(v$genFile)
    # print(v$motiF)
    # print(v$centeR)
    # print(oldmodFile)
    # print(oldgenFile)
    print((is.null(oldmodFile) & is.null(oldgenFile)) || !(oldmodFile == v$modFile & oldgenFile == v$genFile))
    if((is.null(oldmodFile) & is.null(oldgenFile)) || !(oldmodFile == v$modFile & oldgenFile == v$genFile)) {
      print("preupload")
      uploaddat(v$modFile, v$genFile, v$motiF, v$centeR)
      print(csv)
      print(gene)
      print(motiftable)
    }
    graphs <- processdat(v$motiF, v$centeR)
    gps <<- graphs$gs
    gpi <<- graphs$gi
    gpc <<- graphs$gc
    
    oldmodFile <<- v$modFile
    oldgenFile <<- v$genFile
  })
  
  observeEvent(input$reset, {
    v$modFile <- NULL
    v$genFile <- NULL
  })  
  
  output$score <- renderPlotly({
    if (is.null(v$modFile) && is.null(v$genFile)) return()
    ggplotly(gps, width = nchar(v$motiF)*300, height = 350)
  })
  
  output$ipd <- renderPlotly({
    if (is.null(v$modFile) && is.null(v$genFile)) return()
    ggplotly(gpi, width = nchar(v$motiF)*300, height = 350)
  })
  
  output$coverage <- renderPlotly({
    if (is.null(v$modFile) && is.null(v$genFile)) return()
    ggplotly(gpc, width = nchar(v$motiF)*300, height = 350)
  })

}