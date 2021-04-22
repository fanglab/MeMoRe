iupac_nc <<- data.frame(
  code=c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N"),
  pattern=c("A","C","G","T","[AG]","[CT]","[CG]","[AT]","[GT]","[AC]","[CGT]","[AGT]","[ACT]","[ACG]","[ACGT]"),
  choice=c("A","C","G","T","AG","CT","CG","AT","GT","AC","CGT","ATG","ACT","ACG","ACGT")
) # [^] do not rev.comp easily

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
  cat(paste0(string,"\n"))
  cat(paste0(print(mem_used()),"\n"))
}

check.valid.motif <- function(motif){
  results <- grepl('^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$', motif)

  return(results)
}

read.modification.file <- function(modFile, genFile){
  # modFile <- "~/Desktop/Clostridium_perfringens_ATCC13124.modifications.csv.gz"
  # genFile <- "~/Desktop/Clostridium_perfringens_ATCC13124.fasta"

  # Retrieve modFile information
  modFile_con <- file(modFile)
  modFile_info <- summary(modFile_con)
  close(modFile_con)

  if(modFile_info$class == "gzfile"){
    path_modification_file <- gunzip(modFile, remove = FALSE, temporary = TRUE, overwrite = TRUE)
  }else{
    path_modification_file <- modFile
  }

  modification_file <<- fread(path_modification_file, sep=",", header=TRUE, verbose=FALSE, drop=c(6, 7, 8, 11, 12, 13), stringsAsFactors=TRUE) # Read modification file
  memuse("- after modification_file")
  g_seq <<- readDNAStringSet(genFile)
  memuse("- after genome")
  if(!all(levels(modification_file$refName) %in% names(g_seq))){
    showNotification("At least one contig from modifications.csv(.gz) don't match the ones from genome.fasta.", type="error")

    return(c("Not matching"))
  }else{
   
    return(c("Matching"))
  }
}

generate.mutated.motif <- function(motifs_summary){
  # Generate sets of mutated motifs
  # print(paste0("  Generating set of mutated motifs."))
  mutated_motifs <- foreach(idx_motif=seq_along(motifs_summary$motifString), .combine=rbind) %do% {
    len_motif <- nchar(motifs_summary$motifString[idx_motif])
    motif <- motifs_summary$motifString[idx_motif]
    mod_pos <- motifs_summary$centerPos[idx_motif]

    splited_motif <- strsplit(motif,"")[[1]]
    list_mutated_motif <- foreach(idx=seq(1,len_motif), .combine=rbind) %do% {
      tmp_motif <- splited_motif
      sublist_mutated_motif <- foreach(nucleotide=c("A","C","G","T"), .combine=rbind) %do% {
        tmp_motif[idx] <- nucleotide

        return(data.frame(mutated_motif=paste0(tmp_motif,collapse=""), mutation_type=nucleotide))
      }
      sublist_mutated_motif$mutated_motif <- as.character(sublist_mutated_motif$mutated_motif)
      sublist_mutated_motif$pos_mutation <- idx

      return(sublist_mutated_motif)
    }
    list_mutated_motif$original_motif <- motif
    list_mutated_motif$mod_pos <- mod_pos

    return(list_mutated_motif)
  }

  return(mutated_motifs)
}

define.ambiguous.genomic.ranges <- function(g_seq){
  gr_ambiguous_position <- foreach(direction=c("fwd","rev"), .combine=c) %do% {
    tmp_ambiguous_position <- vmatchPattern("N", g_seq, fixed=TRUE)
    tmp_ambiguous_position <- as(tmp_ambiguous_position, "GRanges")
    strand(tmp_ambiguous_position) <- ifelse(direction=="fwd","+","-")

    tmp_ambiguous_position <- reduce(tmp_ambiguous_position)

    return(tmp_ambiguous_position)
  }

  return(gr_ambiguous_position)
}

hide.ambiguous.matches <- function(gr_ambiguous_position, tmp_motif_matches, direction){
  gr_motif_matches <- as(tmp_motif_matches, "GRanges")
  strand(gr_motif_matches) <- ifelse(direction=="fwd","+","-")

  overlaps <- findOverlaps(gr_ambiguous_position, gr_motif_matches, type="any", select="all")
  gr_motif_matches <- gr_motif_matches[setdiff(seq_along(gr_motif_matches), overlaps@to)] # Select non overalpping motif matches

  return(gr_motif_matches)
}

find.motifs.sub <- function(mutated_motifs, idx_motif, direction, g_seq, gr_ambiguous_position){
  original_motif <- mutated_motifs$mutated_motif[idx_motif]
  mod_pos <- mutated_motifs$mod_pos[idx_motif]
  len_motif <- nchar(original_motif)

  if(direction=="rev"){
    motif <- reverseComplement(DNAString(original_motif)) # Double checked
    mod_pos <- (len_motif - mod_pos) + 1
  }else{
    motif <- original_motif
  }

  tmp_motif_matches <- vmatchPattern(motif, g_seq, fixed=FALSE)
  tmp_motif_matches <- hide.ambiguous.matches(gr_ambiguous_position, tmp_motif_matches, direction)

  motif_matches <- data.frame(contig_name=seqnames(tmp_motif_matches), contig_pos_motif=start(tmp_motif_matches) + (mod_pos - 1)) # Mark mod_pos
  if(nrow(motif_matches)>0){
    motif_matches$motif <- as.factor(as.character(original_motif)) # Add motif
  }

  return(motif_matches)
}

find.motifs <- function(g_seq, mutated_motifs, iupac_nc, left_signal, right_signal, error_margin, nb_threads, verbose=TRUE){
  gr_ambiguous_position <- define.ambiguous.genomic.ranges(g_seq) # TODO move out

  pos_isolated_motifs <- foreach(direction=c("fwd","rev"), .combine=rbind) %do% { # Process both strand
    if(verbose){
      print(paste0("Processing ",direction," strand."))
      print(paste0("    Researching motifs from mutated_motifs"))
    }

    # List all modified bases and keep duplicate
    if(nb_threads>1){
      registerDoMC(nb_threads)
      motifs_matches <- foreach(idx_motif=seq(1,nrow(mutated_motifs)), .combine=rbind) %dopar% {
        motif_matches <- find.motifs.sub(mutated_motifs, idx_motif, direction, g_seq, gr_ambiguous_position)

        return(motif_matches)
      }
      registerDoSEQ()
    }else{
      motifs_matches <- foreach(idx_motif=seq(1,nrow(mutated_motifs)), .combine=rbind) %do% {
        motif_matches <- find.motifs.sub(mutated_motifs, idx_motif, direction, g_seq, gr_ambiguous_position)

        return(motif_matches)
      }
    }

    left_free <- left_signal + right_signal + error_margin  # same as right_free
    right_free <- left_signal + right_signal + error_margin # same as left_free
    if(left_signal!=-1 | right_signal!=-1 | error_margin!=-1){
      if(verbose){
        print(paste0("    Removing overlapping motifs"))
      }
      dir_pos_isolated_motifs <- motifs_matches %>%
        arrange(contig_name, contig_pos_motif) %>%
        mutate(dist_next_contig_motif=lead(contig_pos_motif,1) - contig_pos_motif, dist_prev_contig_motif=lag(dist_next_contig_motif,1)) %>%
        mutate(is_contig_end=ifelse(lead(contig_name,1)!=contig_name | is.na(dist_next_contig_motif), TRUE, FALSE)) %>% # Find contigs ends
        mutate(dist_next_contig_motif=ifelse(is_contig_end, nchar(g_seq)[match(contig_name, names(g_seq))] - contig_pos_motif, dist_next_contig_motif)) %>%
        mutate(is_contig_start=ifelse(lag(contig_name,1)!=contig_name | is.na(dist_prev_contig_motif), TRUE, FALSE)) %>% # Find contigs starts
        mutate(dist_prev_contig_motif=ifelse(is_contig_start, contig_pos_motif - 1, dist_prev_contig_motif)) %>%
        dplyr::select(-c(is_contig_start,is_contig_end)) %>%
        filter(dist_next_contig_motif>right_free & dist_prev_contig_motif>left_free) %>%
        mutate(dir=as.factor(direction)) %>%
        dplyr::select(-c(dist_next_contig_motif,dist_prev_contig_motif))
    }else{
      dir_pos_isolated_motifs <- motifs_matches %>%
        mutate(dir=as.factor(direction))
    }

    return(dir_pos_isolated_motifs)
  }

  return(pos_isolated_motifs)
}

extract.motifs.signal <- function(modification_file, g_seq, mutated_motifs, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nb_threads){

  # Creat GRanges of modification information
  gr_modification <- GRanges(
    seqnames=modification_file$refName,
    ranges=IRanges(modification_file$tpl, modification_file$tpl),
    strand=ifelse(modification_file$strand==0,"+","-")
  )

  memuse("- After gr_modification")

  motifs <- find.motifs(g_seq, mutated_motifs, iupac_nc, left_signal, right_signal, error_margin, nb_threads, FALSE)
  expected_signal <- motifs %>%
    mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
    mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)
  gr_signal_motifs <- GRanges(
    seqnames=expected_signal$contig_name,
    ranges=IRanges(expected_signal$left_side, expected_signal$right_side),
    strand=as.factor(ifelse(expected_signal$dir=="fwd","+","-"))
  )

  memuse("- After find.motifs")

  if(filter_iso){
    overlapping_motifs <- overlapsAny(gr_signal_motifs, type="any", drop.self=TRUE)
    gr_signal_motifs <- gr_signal_motifs[!overlapping_motifs]
    expected_signal <- expected_signal[!overlapping_motifs,]
  }

  overlaps_motifs_modification <- findOverlaps(gr_signal_motifs, gr_modification, type="any", select="all")
  modification_at_motifs <- data.frame(
    contig=expected_signal$contig_name[overlaps_motifs_modification@from],
    pos_motif=expected_signal$contig_pos_motif[overlaps_motifs_modification@from],
    motif=expected_signal$motif[overlaps_motifs_modification@from],
    pos_signal=modification_file$tpl[overlaps_motifs_modification@to],
    dir=modification_file$strand[overlaps_motifs_modification@to],
    ipdRatio=modification_file$ipdRatio[overlaps_motifs_modification@to],
    coverage=modification_file$coverage[overlaps_motifs_modification@to],
    score=modification_file$score[overlaps_motifs_modification@to]
  )

  memuse("- After findOverlaps")

  modification_at_motifs <- modification_at_motifs %>% 
    mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
    filter(coverage>=min_cov)
  
  modification_at_motifs <- merge(modification_at_motifs, mutated_motifs, by.x=c("motif"), by.y=c("mutated_motif")) %>%
    mutate(expected_base=ifelse(mutation_type==substr(original_motif,pos_mutation,pos_mutation),1,0))

  return(modification_at_motifs)
}

generate.plots <- function(modification_at_motifs, modificationtype){
  motif <- unique(modification_at_motifs$original_motif)
  mod_pos <- unique(modification_at_motifs$mod_pos)
  graphtitle <- paste0(substr(motif, 1, mod_pos),modificationtype , substr(motif, mod_pos+2, nchar(motif))) 

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
  
  # Separate Graphs Below
  gp_score <- ggplot(modification_at_motifs) + 
    geom_violin(aes(x=mutation_type, y=score, color="blue", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=graphtitle) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="Score") +
    scale_colour_manual(values=c("blue"="#619CCF")) +
    scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
    prettify_base
  
  gp_ipd <- ggplot(modification_at_motifs) + 
    geom_violin(aes(x=mutation_type, y=ipdRatio, color="green", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=graphtitle) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="IPD ratio") +
    scale_colour_manual(values=c("green"="#00BA38")) +
    scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
    prettify_base
  
  gp_cov <- ggplot(modification_at_motifs) + 
    geom_violin(aes(x=mutation_type, y=coverage, color="red", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=graphtitle) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="Coverage") +
    scale_colour_manual(values=c("red"="#F8766D")) +
    scale_fill_manual(values=c("0"="white", "1"="#F8766D")) +
    prettify_base
  
  list_plots <- list(ga=NA, gs=gp_score, gi=gp_ipd, gc=gp_cov)
  
  memuse("- Separate graphs")

  pdf(file=NULL)
  gp_cov <- ggplotGrob(gp_cov + prettify_top) # + labs(title=NULL)
  gp_ipd <- ggplotGrob(gp_ipd + prettify_mid + labs(title=NULL)) # Margin too large, see rbind below
  gp_score <- ggplotGrob(gp_score + prettify_btm + labs(title=NULL))
  
  graphcombined <- arrangeGrob(rbind(gp_cov, gp_ipd, gp_score), ncol=1) #, top=graphtitle
  list_plots$ga <- graphcombined
  dev.off()

  memuse("- Sombined graphs")

  return(list_plots) # , mc=mcount, ms=mscore, mi=mipd, mco=mcov
}

generate.process.data <- function(motif, center, modificationtype){
 # motif <- "GATC"
 # center <- 2
 # modificationtype <- "6mA"
 left_signal <- -1 # No overlap filtering
 right_signal <- -1 # No overlap filtering
 error_margin <- -1 # No overlap filtering
 expected_signal_left <- 0 # Conserve only one position
 expected_signal_right <- 0 # Conserve only one position
 signal_margin <- 0 # Conserve only one position
 filter_iso <- FALSE # Remove overlapping motifs
 min_cov <- 0 # No coverage threshold
 nb_threads <- 1

  memuse("- Start generate.process.data")

  motifs_summary <- data.table(motifString=motif, centerPos=center, modificationType=modificationtype)

  # Generate mutated motifs
  mutated_motifs <- generate.mutated.motif(motifs_summary)
  
  memuse("- After generate.mutated.motif")

  modification_at_motifs <- extract.motifs.signal(modification_file, g_seq, mutated_motifs, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nb_threads)

  memuse("- After extract.motifs.signal")

  list_plots <- generate.plots(modification_at_motifs, modificationtype)

  memuse("- After generate.plots")

  return(list_plots)
}

processdat <- function(motif, center, modificationtype) {
 # motif <- "VGACAT"
 # center <- "2"
 # modificationtype <- "6mA"

  motiftable <- data.table(motifString = motif, centerPos = center, modificationType =strsplit(motif,"")[[1]][as.integer(center)+1])
  motif_ref <- data.frame(sym = c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"), bases = c("(A|T)", "(C|G)", "(A|C)", "(G|T)", "(A|G)", "(C|T)", "(C|G|T)", "(A|G|T)", "(A|C|T)", "(A|C|G)", "(A|C|G|T)"))
  cl_max <- 8
  modtype <- as.vector(motiftable$modificationType)
  methvec <- gsub("\\d|[[:lower:]]","",modtype)
  motif_f <- as.vector(motiftable$motifString)
  cov_cut_spec <- 10
  size_spec <- 100
  
  gag <- which(width(g_seq) == max(width(g_seq)))
  csv2 <- modification_file
  # csv2 <<- csv[csv$refName == names(genome)[gag]] #csv2 = csv

  # # R and F prep
  # csv2 <- csv2[complete.cases(csv2), ]
  # csv2 <- csv2[order(csv2$tpl),]
  param <- c(csv2$tpl[1], csv2$tpl[nrow(csv2)])  # create reverse genome
  memuse("- after param")
  genome_f <- toString(g_seq[[gag]])
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
    
    # gp_cov <- ggplot(datacoverage) + 
    #   geom_violin(aes(x=base, y=values, color="red", group=base)) +
    #   geom_text(data=datacount, aes(x=base, y=10, label=labels), size=5, colour="grey25") + 
    #   facet_grid(. ~ pos) +
    #   scale_colour_manual(values=c("red"="#F8766D")) +
    #   labs(y="Coverage") +
    #   prettify_top
    
    # gp_ipd <- ggplot(dataipd) + 
    #   geom_violin(aes(x=base, y=values, color="green", group=base, fill=as.factor(mod))) +
    #   coord_cartesian(ylim=c(0,10)) +
    #   facet_grid(. ~ pos) +
    #   scale_colour_manual(values=c("green"="#00BA38")) +
    #   scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
    #   labs(y="IPD ratio") +
    #   prettify_mid
    
    # gp_scr <- ggplot(datascore) + 
    #   geom_violin(aes(x=base, y=values, color="blue", group=base, fill=as.factor(mod))) +
    #   facet_grid(. ~ pos) +
    #   scale_colour_manual(values=c("blue"="#619CCF")) +
    #   scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
    #   labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
    #   labs(y="Score") +
    #   prettify_btm
    
    # Separate Graphs Below
    gp_scr <- ggplot(datascore) + 
      geom_violin(aes(x=base, y=values,color = "blue", group=base, fill=as.factor(mod))) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="score") +
      scale_colour_manual(values=c("blue"="#619CCF")) +
      scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
      prettify_base
    
    gp_ipd <- ggplot(dataipd) + 
      geom_violin(aes(x=base, y=values,color = "green", group=base, fill=as.factor(mod))) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="ipdRatio") +
      scale_colour_manual(values=c("green"="#00BA38")) +
      scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
      prettify_base
    
    gp_cov <- ggplot(datacoverage) + 
      geom_violin(aes(x=base, y=values,color = "red", group=base)) + 
      facet_grid(.~pos) + 
      ggtitle(graphtitle) +
      labs(x=paste0("Alternate base in ",motif_f[j]," motif")) +
      labs(y="coverage") +
      scale_colour_manual(values=c("red"="#F8766D")) +
      prettify_base
    
    memuse("- separate graphs")

    gp_cov <- ggplotGrob(gp_cov + prettify_top)
    gp_ipd <- ggplotGrob(gp_ipd + prettify_mid) # Margin too large, see rbind below
    gp_scr <- ggplotGrob(gp_scr + prettify_btm)
    
    graphcombined <- arrangeGrob(rbind(gp_cov, gp_ipd, gp_scr), ncol=1, top=graphtitle)

    memuse("- combined graphs")
    
    # garbage collection to reduce memory
    gc()
    memuse("- garbage collection")
    rm(list=setdiff(ls(), c("graphcombined","gp_scr","gp_ipd","gp_cov",
                            "mcount","mscore","mipd","mcov", lsf.str(), "g_seq","csv2")))
    memuse("- manual garbage collection")
  
    return(list(ga = graphcombined, gs = gp_scr, gi = gp_ipd, gc = gp_cov)) # , mc = mcount, ms = mscore, mi = mipd, mco = mcov
  }
}
