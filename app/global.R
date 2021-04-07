# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=200) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

# print mem usage
memuse <- function(string) {
  print(string)
  print(mem_used())
}

uploaddat <- function(rmodfile, rgenfile, motif, center) {  
  if (summary(file(rmodfile))$class == "gzfile") {
    csv_temp <- gunzip(rmodfile, remove = T, temporary = T, overwrite = T)
  } else {
    csv_temp <- rmodfile
  }
  csv_t <- fread(csv_temp, sep = ",", header = T, verbose = F, drop = c(6, 7, 8, 11, 12, 13))  # Set csv parameters
  memuse("- after csv_t")
  gene <<- readDNAStringSet(rgenfile)
  memuse("- after gene")
  csv2 <- csv_t[csv_t$refName == names(gene)[which(width(gene) == max(width(gene)))]]
  csv2 <- csv2[complete.cases(csv2), ]
  csv2 <<- csv2[order(csv2$tpl),]
  memuse("- after csv2")
  rm(csv_t)
  memuse("- after rm(csv_t)")
}

processdat <- function(motif, center, modificationtype) {
  motiftable <- data.table(motifString = motif, centerPos = center, modificationType =strsplit(motif,"")[[1]][as.integer(center)+1])
  motif_ref <- data.frame(sym = c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"), bases = c("(A|T)", "(C|G)", "(A|C)", "(G|T)", "(A|G)", "(C|T)", "(C|G|T)", "(A|G|T)", "(A|C|T)", "(A|C|G)", "(A|C|G|T)"))
  cl_max <- 8
  modtype <- as.vector(motiftable$modificationType)
  methvec <- gsub("\\d|[[:lower:]]","",modtype)
  motif_f <- as.vector(motiftable$motifString)
  cov_cut_spec <- 10
  size_spec <- 100
  
  gag <- which(width(gene) == max(width(gene)))
  # csv2 <<- csv[csv$refName == names(gene)[gag]] #csv2 = csv
  
  # # R and F prep
  # csv2 <- csv2[complete.cases(csv2), ]
  # csv2 <- csv2[order(csv2$tpl),]
  param <- c(csv2$tpl[1], csv2$tpl[nrow(csv2)])  # create reverse genome
  memuse("- after param")
  genome_f <- toString(gene[[gag]])
  memuse("- after genome_f")
  genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
  memuse("- after genome_r")
  # cluster prep
  csv_sp <- csv2[csv2$coverage >= cov_cut_spec, ]
  memuse("- after csv_sp")
  csv_sp_f <- csv2[(csv2$strand == 0), ][csv2[(csv2$strand == 0), ]$coverage >= cov_cut_spec, ]
  memuse("- after csv_sp_f")
  csv_sp_r <- csv2[(csv2$strand == 1), ][csv2[(csv2$strand == 1), ]$coverage >= cov_cut_spec, ]
  memuse("- after csv_sp_r")
  
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
    dog1 <- dog1.me[2:ncol(dog1.me)]          ## NEW DATA prep CSV
    csv_sp_2 <- csv_sp[order(csv_sp$tpl)][complete.cases(csv_sp[order(csv_sp$tpl)]), ]
    param_sp <- c(csv_sp_2$tpl[1], csv_sp_2$tpl[nrow(csv_sp_2)])
    dat2.A.1 <- csv_sp[csv_sp$base == meth, ]
    dat2.A <- dat2.A.1[sample(nrow(dat2.A.1), size_spec), ]
    
    exportme <- list(transl_f, size_spec, param_sp)
    memuse("- preparing data for moteity loop")
    
    ## Future Improvement: Parallelize

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
      #print(z.count)
    }
    memuse("- after loop through motif moteities")
    
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
    memuse("- after create graphing data")
    
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
    
    dataipd <- dataipd
    datacoverage <- datacoverage
    graphtitle <- paste0(substr(motif_f[j], 1, center), "(",modificationtype,")", substr(motif_f[j], center+2, nchar(motif_f[j]))) 
    
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

    memuse("- combined graphs")
    
    # Separate Graphs Below
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
    
    memuse("- separate graphs")

    # garbage collection to reduce memory
    gc()
    memuse("- garbage collection")
    rm(list=setdiff(ls(), c("graphcombined","graphscore","graphipd","graphcoverage",
                            "mcount","mscore","mipd","mcov", lsf.str(), "gene","csv2")))
    memuse("- manual garbage collection")
  
    return(list(ga = graphcombined, gs = graphscore, gi = graphipd, gc = graphcoverage, 
                mc = mcount, ms = mscore, mi = mipd, mco = mcov))
    
  }
}