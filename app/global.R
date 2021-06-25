options(repos=BiocManager::repositories()) # Needed for shinyapp.io

iupac_nc <<- data.frame(
  code=c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N"),
  pattern=c("A","C","G","T","[AG]","[CT]","[CG]","[AT]","[GT]","[AC]","[CGT]","[AGT]","[ACT]","[ACG]","[ACGT]"),
  choice=c("A","C","G","T","AG","CT","CG","AT","GT","AC","CGT","ATG","ACT","ACG","ACGT")
) # [^] do not rev.comp easily
# Not used

# list_motif_summary_SMRT_cols <<- c("motifString", "centerPos", "modificationType", "fraction", "nGenome", "partnerMotifString", "meanScore", "meanIpdRatio", "meanCoverage")
# list_motif_summary_clean_SMRT_cols <<- c("Motifs", "Modified position", "Type", "% motifs detected", "# motifs in genome", "Partner motif", "Mean Score", "Mean IPD ratio", "Mean Coverage", "Plots Generated")
list_motif_summary_SMRT_cols <<- c("motifString", "centerPos", "modificationType")
list_motif_summary_clean_SMRT_cols <<- c("Motifs", "Modified position", "Type", "Plots Generated")
list_motif_summary_ONT_cols <<- c("Motif", "Predicted_position", "Predicted_type")
list_motif_summary_clean_ONT_cols <<- c("Motifs", "Modified position", "Type", "Plots Generated")
list_motif_manual_clean_cols <<- list_motif_summary_clean_SMRT_cols[c(1,2,3)]

load.libraries <- function(){
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
  library(RColorBrewer)
  library(egg)
  library(cowplot)
  library(rhdf5)
}

# Improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5){
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))

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

# Shorthand .ls.objects
lsos <- function(..., n=200){
  results <- .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)

  return(results)
}

# Print current memory usage
memuse <- function(string){
  cat(paste0(string,"\n"))
  cat(paste0(print(mem_used()),"\n"))
}

print_db <- function(message){
  if(debug_mode){
    print(message)
  }
}

initialize.motif.summary <- function(v, list_motif_summary_clean_cols){
  # Initialize motif summary
  v$df <- setNames(data.table(matrix(nrow = 0, ncol = length(list_motif_summary_clean_cols))), list_motif_summary_clean_cols)
}

reformat.motif.summary <- function(motif_summary){
  # Correct modified base position depending on summary type
  if(any(motif_summary$`Modified position`==0)){
    # It's definitely 0-based
    mod_pos_offset <- 1
  }else{
    motif_summary_version <- motif_summary %>%
      mutate(mod_base=substr(Type, nchar(Type), nchar(Type))) %>%
      mutate(mod_base=ifelse(mod_base %in% c("A","C"), mod_base, "x")) %>%
      mutate(zero_based=substr(Motifs,`Modified position`+1,`Modified position`+1)) %>%
      mutate(one_based=substr(Motifs,`Modified position`,`Modified position`)) %>%
      mutate(zero_based_match=ifelse(mod_base=="x",ifelse(zero_based %in% c("A","C"),1,0),ifelse(zero_based==mod_base,1,0))) %>%
      mutate(one_based_match=ifelse(mod_base=="x",ifelse(one_based %in% c("A","C"),1,0),ifelse(one_based==mod_base,1,0))) %>%
      summarize(zero_based_res=sum(zero_based_match)/n(), one_based_res=sum(one_based_match)/n())
    if(motif_summary_version$zero_based_res > motif_summary_version$one_based_res){
      # It's likely 0-based
      mod_pos_offset <- 1
    }else{
      # It's likely 1-based
      mod_pos_offset <- 0
    }
  }

  motif_summary <- motif_summary %>%
    mutate(`Modified position`=`Modified position` + mod_pos_offset) %>%
    mutate(Type=ifelse(is.na(Type),paste0("x",substr(Motifs,`Modified position`,`Modified position`)),Type)) %>%
    mutate(Type=ifelse(Type=="modified_base",paste0("x",substr(Motifs,`Modified position`,`Modified position`)),Type)) # From SMRTLink
  if("% motifs detected" %in% colnames(motif_summary)){
    motif_summary <- motif_summary %>%
    mutate(`% motifs detected`=round(as.numeric(`% motifs detected`),2)) %>%
    mutate(`Mean Score`=round(as.numeric(`Mean Score`))) %>%
    mutate(`Mean IPD ratio`=round(as.numeric(`Mean IPD ratio`),2)) %>%
    mutate(`Mean Coverage`=round(as.numeric(`Mean Coverage`)))
  }

  return(motif_summary)
}

read.motif.summary <- function(v, inFile){
  if(!is.null(inFile)){
    if(grepl(inFile$datapath, pattern=".tsv$")){
      data_type <- "ONT"
      columns_sep <- '\t'
      list_motif_summary_cols <- list_motif_summary_ONT_cols
      list_motif_summary_clean_cols <- list_motif_summary_clean_ONT_cols
    }else if(grepl(inFile$datapath, pattern=".csv$")){
      data_type <- "SMRT"
      columns_sep <- ','
      list_motif_summary_cols <- list_motif_summary_SMRT_cols
      list_motif_summary_clean_cols <- list_motif_summary_clean_SMRT_cols
    }else{
      showNotification("Motif summary file format not recognized.", type="warning")
      Sys.sleep(1)
      showNotification("SMRT Link (.csv) or Nanodisco (.tsv) output expected.", type="error")

      return(NULL)
    }

    motif_summary <- read.table(inFile$datapath, sep=columns_sep, header=TRUE) # Read motif_summary.csv file
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
    # motif_summary$plotted[which(motif_summary$motifString %in% v$df$Motifs)] <- "Yes"

    colnames(motif_summary) <- list_motif_summary_clean_cols # Rename columns with clean names

    if(data_type=="SMRT"){
      motif_summary <- reformat.motif.summary(motif_summary)
    }

    # Remove duplicate entries
    new_motifs <- paste0(motif_summary$`Motifs`,"_",motif_summary$`Modified position`)
    isolate(current_motifs <- paste0(v$df$`Motifs`,"_",v$df$`Modified position`))
    motif_summary <- motif_summary[!new_motifs %in% current_motifs,]

    # Update main reactive object without re-execution
    isolate(v$df <- base::rbind(v$df, motif_summary, fill=TRUE))
  }
}

check.valid.motif <- function(motif){
  results <- grepl('^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$', motif)

  return(results)
}

check.inputs  <- function(input){
  # Check modifications.csv.gz 
  if(!testing_mode & !debug_mode & is.null(input$modfile)){
    showNotification("Upload modifications (.csv.gz, .h5)", type="error")

    return(NULL)
  }
  # Check genome.fasta
  if(!testing_mode & !debug_mode & is.null(input$genfile)){
    showNotification("Upload genome (.fa, .fasta)", type="error")

    return(NULL)
  }
  # Check modifications.csv.gz
  if(!testing_mode & !debug_mode & is.null(input$motfile)){
    showNotification("Upload motif_summary (.csv, .tsv)", type="error")

    return(NULL)
  }

  return("OK.")
}

wrapper.data.processing <- function(input, v, list_selected_motifs, list_motif_manual_clean_cols){
  if(is.null(check.inputs(input))){
    # At least an input is missing

    return(NULL)
  }

  # Checking selected motifs or all motifs
  if(is.null(list_selected_motifs)){
    showNotification("Select a motif.", type="error")

    return(NULL)
  }

  # DATA INPUT
  if(!testing_mode){ # Use users inputs if not in testing_mode
    v$modFile <- input$modfile$datapath
    v$genFile <- input$genfile$datapath
    v$motFile <- input$motfile$datapath
  }

  if(debug_mode){ # TODO remove
    if(is.null(rendered_motif) | is.null(v$modFile)){
      v$modFile <- "data/modification.smrt.csv.gz"
      # v$modFile <- "data/modification.ont.rds"
    }
    v$genFile <- "data/reference.fasta"
    v$motFile <- "data/motif_summary.tsv"
  }
  
  v$motiF <- NULL
  v$centeR <- NULL
  v$modType <- NULL

  data_type <- find.data.type(v$modFile)
  if(is.null(data_type)){
    showNotification("Modification file format not recognized.", type="warning")
    Sys.sleep(1)
    showNotification("SMRT Link (.csv.gz or .h5) or Nanodisco (.RDS) output expected.", type="error")

    return(NULL)
  }
  # Stop if all selected motifs were already processed
  if(all(v$df$"Plots Generated"[list_selected_motifs] %in% c("Both", data_type))){
    showNotification(paste0("Plots for all motifs were already generated."), type="warning")
    
    return(NULL)
  }

  nb_selected_motifs <- length(list_selected_motifs)
  for(k in 1:nb_selected_motifs){
    v$motiF[k] <- toString(v$df$Motifs[list_selected_motifs[k]])
    v$centeR[k] <- as.numeric(v$df$"Modified position"[list_selected_motifs[k]])
    v$modType[k] <- toString(v$df$Type[list_selected_motifs[k]])
  }

  for(j in 1:nb_selected_motifs){
    # Skip motifs already processed with a warning.
    if(v$df$"Plots Generated"[list_selected_motifs[j]] %in% c("Both", data_type)){
      showNotification(paste0("Plots for ",v$motiF[j]," was already generated."), type="warning")
      next
    }

    # Progress Bar (most time intensive)
    print_db("Start processing motif(s).")
    memuse("After start processing motif(s).")
    if(nb_selected_motifs>1){
      progress_message <- paste0('Processing motifs (',j,'/',nb_selected_motifs,').')
    }else{
      progress_message <- paste0('Processing motif.')
    }
    withProgress(message = progress_message, value = 0, {
      incProgress(.05, detail = "Checking motif(s).")
      # test if uploaded genome and modifications files are changed
      isemptyorchanged <- (is.null(oldmodFile) & is.null(oldgenFile)) || !(oldmodFile == v$modFile & oldgenFile == v$genFile)
      
      # If changed, reupload
      if(isemptyorchanged) {
        incProgress(.05, detail = "Reading input files.")
        results <- read.modification.file(v$modFile, v$genFile)
        memuse("After read.modification.file")          
        if(results=="Not matching"){
          break
        }
      }else{
        incProgress(.05, detail = "Input files already loaded.")
      }
      
      # Continue to process files
      # incProgress(.2, detail = paste0("Extracting ",v$motiF[j]," information."))
      if(processing_version==1){
        graphs <- processdat(v$motiF[j], v$centeR[j] - 1, v$modType[j]) # 137s and peak mem usage at 1.2 GB

        # Save graphs
        motif_detail <- paste0(v$motiF[j], "_", v$centeR[j])
        saveRDS(graphs$ga, file=paste0(motif_detail, ".gpa"))
        saveRDS(graphs$gs, file=paste0(motif_detail, ".gps"))
        saveRDS(graphs$gi, file=paste0(motif_detail, ".gpi"))
        saveRDS(graphs$gc, file=paste0(motif_detail, ".gpc"))
        rm(graphs)
      }else if(processing_version==2){
        gc()
        generate.process.data(v$motiF[j], v$centeR[j], v$modType[j], v$df) # 82s and peak mem usage at 750 MB
        gc()
      }
      memuse("After process dat")

      incProgress(.05, detail = paste0("Update tables with ",v$motiF[j]," information."))

      isolate({
        v$df$"Plots Generated" <- as.character(v$df$"Plots Generated")
        if(v$df$"Plots Generated"[list_selected_motifs[j]] == "No"){
          v$df$"Plots Generated"[list_selected_motifs[j]] <- data_type
        }else{
          v$df$"Plots Generated"[list_selected_motifs[j]] <- "Both"
        }
        v$df$"Plots Generated" <- as.factor(v$df$"Plots Generated")
      })

      oldmodFile <<- v$modFile
      oldgenFile <<- v$genFile

      memuse("After assigning variables")

      incProgress(.35, detail = paste0("Done."))
    })
  }

  return("Done.")
}

print.runtime.message <- function(start_time){
  run_time <- proc.time() - start_time
  run_time_message <- paste0("Run time was ",round(run_time[[3]],0),"s.")
  print_db(run_time_message)
  showNotification(run_time_message, type="message")
}

find.data.type <- function(path_modification_file){
  if(grepl(path_modification_file, pattern="*.RDS$|*.rds$")){
    data_type <- "ONT"
  }else if(grepl(path_modification_file, pattern="*.csv.gz$|*.csv$|*.h5$|*/0.gz$|*/0$")){ # If already in tmp */0.gz$ */0$; Not clean
    data_type <- "SMRT"
  }else{

    return(NULL)
  }

  return(data_type)
}

read.hdf.file <- function(path_modification_file){
  modification_info_tmp <- foreach(refName=subset(h5ls(path_modification_file), group=="/")$name, .combine=rbind) %do% {
    subset_modification_info <- data.frame(
      refName=as.factor(refName),
      tpl=h5read(file=path_modification_file, name=paste0(refName, "/tpl")),
      strand=as.integer(h5read(file=path_modification_file, name=paste0(refName, "/strand"))),
      score=h5read(file=path_modification_file, name=paste0(refName, "/score")),
      ipdRatio=h5read(file=path_modification_file, name=paste0(refName, "/ipdRatio")),
      coverage=h5read(file=path_modification_file, name=paste0(refName, "/coverage"))
    )

    return(subset_modification_info)
  }
  h5closeAll()

  return(modification_info_tmp)
}

read.modification.file <- function(modFile, genFile){
  is_h5 <- FALSE

  # Retrieve modFile information
  modFile_con <- file(modFile)
  modFile_info <- summary(modFile_con)
  close(modFile_con)

  if(modFile_info$class == "gzfile"){
    path_modification_file <- gunzip(modFile, temporary=TRUE, remove=FALSE, overwrite=TRUE)
  }else if(grepl("*.zip$", modFile)){
    path_modification_file <- unzip(modFile, overwrite=TRUE, junkpaths=TRUE)
  }else if(grepl("*.h5$", modFile)){
    is_h5 <- TRUE
    path_modification_file <- modFile
  }else{
    path_modification_file <- modFile
  }

  data_type <- find.data.type(path_modification_file)
  if(is.null(data_type)){
    showNotification("Modification file format not recognized.", type="warning")
    Sys.sleep(1)
    showNotification("SMRT Link (.csv.gz or .h5) or Nanodisco (.RDS) output expected.", type="error")

    return(NULL)
  }

  # Reading dataset
  if(data_type=="SMRT"){
    if(is_h5){
      modification_info <<- read.hdf.file(path_modification_file)
    }else{
      modification_info <<- fread(path_modification_file, sep=",", header=TRUE, verbose=FALSE, drop=c(4, 6, 7, 8, 11, 12, 13), stringsAsFactors=TRUE) # Read modification file
    }

    list_contigs <- str_split(levels(modification_info$refName), pattern=" ", simplify=TRUE)[,1]
  }else if(data_type=="ONT"){
    modification_info <<- readRDS(path_modification_file)
    list_contigs <- levels(modification_info$contig)
  }

  memuse("- after modification_info")
  
  g_seq <<- readDNAStringSet(genFile)
  
  memuse("- after genome")

  if(!all(list_contigs %in% str_split(names(g_seq), pattern=" ", simplify=TRUE)[,1])){
    showNotification("At least one contig from the modification file don't match the ones from the genome.fasta.", type="error")

    return(c("Not matching"))
  }else{
   
    return(c("Matching"))
  }
}

generate.mutated.motif <- function(motifs_summary){
  # Generate sets of mutated motifs
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

  # Annotated expected base from queried motif
  mutated_motifs <- mutated_motifs %>%
    rowwise() %>%
    mutate(expected_base=ifelse(mutation_type %in% as.vector(str_split(iupac_nc$choice[iupac_nc$code==substr(original_motif, pos_mutation, pos_mutation)], "", simplify=TRUE)), 1, 0))

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
  original_motif <- mutated_motifs$motif[idx_motif]
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

extract.motifs.signal <- function(modification_info, data_type, g_seq, mutated_motifs, motifs_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nb_threads){

  incProgress(.04, detail=paste0("Localize ",unique(mutated_motifs$original_motif)," motif sites."))
  mutated_motifs$motif <- mutated_motifs$mutated_motif # Adapt for detection of mutated motif
  motifs <- find.motifs(g_seq, mutated_motifs, iupac_nc, left_signal, right_signal, error_margin, nb_threads, FALSE)
  expected_signal <- motifs %>%
    mutate(contig_name=str_split(contig_name, " ", simplify=TRUE)[,1]) %>%
    mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
    mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)
  gr_signal_motifs <- GRanges(
    seqnames=expected_signal$contig_name,
    ranges=IRanges(expected_signal$left_side, expected_signal$right_side),
    strand=as.factor(ifelse(expected_signal$dir=="fwd","+","-"))
  )

  memuse("- After find.motifs")

  if(filter_iso){
    memuse("- Filtering isolated motifs")

    incProgress(.04, detail=paste0("Localize background motif sites."))
    background_motifs_summary <- subset(motifs_summary, ! Motifs %in% unique(mutated_motifs$original_motif))
    # Adapt for detection of background motifs
    background_motifs_summary$motif <- background_motifs_summary$Motifs
    background_motifs_summary$mod_pos <- background_motifs_summary$'Modified position'
    background_motifs <- find.motifs(g_seq, background_motifs_summary, iupac_nc, left_signal, right_signal, error_margin, nb_threads, FALSE)
    expected_background_signal <- background_motifs %>%
      mutate(contig_name=str_split(contig_name, " ", simplify=TRUE)[,1]) %>%
      mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
      mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)
    gr_background_signal_motifs <- GRanges(
      seqnames=expected_background_signal$contig_name,
      ranges=IRanges(expected_background_signal$left_side, expected_background_signal$right_side),
      strand=as.factor(ifelse(expected_background_signal$dir=="fwd","+","-"))
    )

    incProgress(.03, detail=paste0("Remove overlapping motif sites."))
    withCallingHandlers(expr={
      overlapping_background_motifs <- findOverlaps(gr_signal_motifs, gr_background_signal_motifs, type="any", select="all")
    }, warning=function(w){
      if(grepl("Each of the 2 combined objects has sequence levels not in the other", w$message)){
        invokeRestart("muffleWarning")
      }
    })
    if(length(overlapping_background_motifs)>0){
      gr_signal_motifs <- gr_signal_motifs[-overlapping_background_motifs@from]
      expected_signal <- expected_signal[-overlapping_background_motifs@from,]
    }
    memuse("- After isolated motifs")
  }else{
    incProgress(.07, detail="")
  }

  if(data_type=="SMRT"){
    # Creat GRanges of modification information
    incProgress(.04, detail=paste0("Integrate modification information."))
    gr_modification <- GRanges(
      seqnames=str_split(modification_info$refName, " ", simplify=TRUE)[,1],
      ranges=IRanges(modification_info$tpl, modification_info$tpl),
      strand=ifelse(modification_info$strand==0,"+","-")
    )

    incProgress(.03, detail=paste0("Filter modification information."))
    overlaps_motifs_modification <- findOverlaps(gr_signal_motifs, gr_modification, type="any", select="all")
    rm(gr_signal_motifs, gr_modification)
    gc()
    modification_at_motifs <- data.frame(
      contig=expected_signal$contig_name[overlaps_motifs_modification@from],
      pos_motif=expected_signal$contig_pos_motif[overlaps_motifs_modification@from],
      motif=expected_signal$motif[overlaps_motifs_modification@from],
      pos_signal=modification_info$tpl[overlaps_motifs_modification@to],
      dir=modification_info$strand[overlaps_motifs_modification@to],
      ipdRatio=modification_info$ipdRatio[overlaps_motifs_modification@to],
      coverage=modification_info$coverage[overlaps_motifs_modification@to],
      score=modification_info$score[overlaps_motifs_modification@to]
    )
    rm(expected_signal, overlaps_motifs_modification)
    gc()

    modification_at_motifs <- modification_at_motifs %>% 
      mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
      filter(coverage>=min_cov)
  }else if(data_type=="ONT"){
    # Creat GRanges of modification information
    incProgress(.04, detail=paste0("Integrate modification information."))
    gr_modification <- GRanges(
      seqnames=modification_info$contig,
      ranges=IRanges(modification_info$position, modification_info$position),
      strand=ifelse(modification_info$dir=="fwd","+","-")
    )

    incProgress(.03, detail=paste0("Filter modification information."))
    overlaps_motifs_modification <- findOverlaps(gr_signal_motifs, gr_modification, type="any", select="all")
    rm(gr_signal_motifs, gr_modification)
    gc()
    modification_at_motifs <- data.frame(
      contig=expected_signal$contig_name[overlaps_motifs_modification@from],
      pos_motif=expected_signal$contig_pos_motif[overlaps_motifs_modification@from],
      motif=expected_signal$motif[overlaps_motifs_modification@from],
      pos_signal=modification_info$position[overlaps_motifs_modification@to],
      dir=modification_info$dir[overlaps_motifs_modification@to],
      strand=modification_info$strand[overlaps_motifs_modification@to],
      N_wga=modification_info$N_wga[overlaps_motifs_modification@to],
      N_nat=modification_info$N_nat[overlaps_motifs_modification@to],
      mean_diff=modification_info$mean_diff[overlaps_motifs_modification@to]
    )
    rm(expected_signal, overlaps_motifs_modification)
    gc()

    modification_at_motifs <- modification_at_motifs %>%
      mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
      filter(N_wga>=min_cov & N_nat>=min_cov)
  }

  memuse("- After findOverlaps")
  modification_at_motifs <- merge(modification_at_motifs, mutated_motifs, by.x=c("motif"), by.y=c("mutated_motif"))

  return(modification_at_motifs)
}

generate.process.data <- function(motif, center, modification_type, motifs_summary){
  # motif <- "GATC"
  # center <- 3
  # modification_type <- "5mC"
  # read.modification.file("app/Clostridium_perfringens_ATCC13124.RDS","app/Clostridium_perfringens_ATCC13124.fasta")

  left_signal <- -1 # No overlap filtering
  right_signal <- -1 # No overlap filtering
  error_margin <- -1 # No overlap filtering

  filter_iso <- FALSE # Remove overlapping motifs
  min_cov <- 0 # No coverage threshold
  nb_threads <- 1

  if(c("tpl") %in% colnames(modification_info)){
    data_type <- "SMRT"
    expected_signal_left <- 0 # Conserve only one position
    expected_signal_right <- 0 # Conserve only one position
    signal_margin <- 0 # Conserve only one position
  }else if(c("mean_diff") %in% colnames(modification_info)){
    data_type <- "ONT"
    expected_signal_left <- -6
    expected_signal_right <- -1
    signal_margin <- 4
  }

  memuse("- Start generate.process.data")

  motif_to_process <- data.table(motifString=motif, centerPos=center, modificationType=modification_type)

  # Generate mutated motifs
  incProgress(.01, detail = paste0("Generating mutated motifs for ",motif,"."))
  mutated_motifs <- generate.mutated.motif(motif_to_process)

  memuse("- After generate.mutated.motif")

  incProgress(.01, detail = paste0("Extracting ",motif," information."))
  modification_at_motifs <- extract.motifs.signal(modification_info, data_type, g_seq, mutated_motifs, motifs_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nb_threads)

  memuse("- After extract.motifs.signal")

  incProgress(.3, detail = paste0("Saving plot for ",motif,"."))
  generate.plots(modification_at_motifs, data_type, modification_type, filter_iso)

  memuse("- After generate.plots")
}

generate.smrt.plots <- function(modification_at_motifs, modification_type){
  motif <- unique(modification_at_motifs$original_motif)
  mod_pos <- unique(modification_at_motifs$mod_pos)
  motif_detail <- paste0(motif, "_", mod_pos)
  clean_motif <- paste0(substr(motif, 1, mod_pos-1), modification_type, substr(motif, mod_pos+1, nchar(motif))) 

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
  gp_score <- ggplot(subset(modification_at_motifs, select=c("mutation_type","score","expected_base","pos_mutation"))) + 
    geom_violin(aes(x=mutation_type, y=score, color="blue", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=clean_motif) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="Score") +
    scale_colour_manual(values=c("blue"="#619CCF")) +
    scale_fill_manual(values=c("0"="white", "1"="#619CCF")) +
    prettify_base
  saveRDS(gp_score, file=paste0(motif_detail, ".gps"))

  gp_ipd <- ggplot(subset(modification_at_motifs, select=c("mutation_type","ipdRatio","expected_base","pos_mutation"))) + 
    geom_violin(aes(x=mutation_type, y=ipdRatio, color="green", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=clean_motif) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="IPD ratio") +
    scale_colour_manual(values=c("green"="#00BA38")) +
    scale_fill_manual(values=c("0"="white", "1"="#00BA38")) +
    prettify_base
  saveRDS(gp_ipd, file=paste0(motif_detail, ".gpi"))

  gp_cov <- ggplot(subset(modification_at_motifs, select=c("mutation_type","coverage","expected_base","pos_mutation"))) + 
    geom_violin(aes(x=mutation_type, y=coverage, color="red", group=mutation_type, fill=as.factor(expected_base))) + 
    facet_grid(.~pos_mutation) + 
    labs(title=clean_motif) +
    labs(x=paste0("Alternate base in ",motif," motif")) +
    labs(y="Coverage") +
    scale_colour_manual(values=c("red"="#F8766D")) +
    scale_fill_manual(values=c("0"="white", "1"="#F8766D")) +
    prettify_base
  saveRDS(gp_cov, file=paste0(motif_detail, ".gpc"))

  memuse("- Separate graphs")

  pdf(file=NULL)
  gp_ipd <- ggplotGrob(gp_ipd + prettify_top)
  gp_score <- ggplotGrob(gp_score + prettify_mid + labs(title=NULL)) # Margin too large, see rbind below
  gp_cov <- ggplotGrob(gp_cov + prettify_btm + labs(title=NULL))
  
  gp_all <- arrangeGrob(rbind(gp_ipd, gp_score, gp_cov), ncol=1) #, top=clean_motif
  saveRDS(gp_all, file=paste0(motif_detail, ".gpa"))
  dev.off()

  memuse("- Combined graphs")
}

score.mutated.motifs <- function(modification_at_motifs){
  mutated_motif_score <- modification_at_motifs %>% # TODO try to improve summary function
    group_by(mutation_type, pos_mutation, distance) %>%
    summarize(score2=abs(mean(mean_diff, na.rm=TRUE)), .groups="drop_last") %>%
    group_by(mutation_type, pos_mutation) %>%
    summarize(score=sum(score2), .groups="drop_last")
  mutated_motif_score$mutation_type <- ordered(mutated_motif_score$mutation_type, levels=c("T","G","C","A")) # Same order as facet_grid

  return(mutated_motif_score)
}

generate.ont.plots <- function(modification_at_motifs, modification_type, filter_iso){
  xmin_value <- -7.5
  xmax_value <- 6.5

  motif <- unique(modification_at_motifs$original_motif)
  mod_pos <- unique(modification_at_motifs$mod_pos)
  motif_detail <- paste0(motif, "_", mod_pos)
  clean_motif <- paste0(substr(motif, 1, mod_pos-1), modification_type, substr(motif, mod_pos+1, nchar(motif))) 

  modification_at_motifs <- subset(modification_at_motifs, select=c("mutation_type", "pos_mutation", "distance", "mean_diff"), !is.na(mean_diff)) %>%
    mutate(distance=distance + 3)

  if(nrow(modification_at_motifs)>15){ # If motif found and more than 15 data point left TODO arbitrary threshold
    # Compute mutated motif scores
    mutated_motif_score <- score.mutated.motifs(modification_at_motifs)

    # Plot mutated motif signal
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    ymin_value <- max(-10,min(modification_at_motifs$mean_diff))
    ymax_value <- min(10,max(modification_at_motifs$mean_diff))
    gp_detail <- ggplot(modification_at_motifs) + # modification_at_motifs[sample(nrow(modification_at_motifs), 100000), ]
      geom_rect(data=mutated_motif_score, aes(xmin=xmin_value, xmax=xmax_value, ymin=ymin_value, ymax=ymin_value + 2, fill=score)) +
      geom_jitter(aes(x=distance, y=mean_diff), pch=46, height=0) +
      geom_violin(aes(x=distance, y=mean_diff, group=distance), alpha=0.6) + # Throw warning if < one motif per strand|dir
      geom_hline(yintercept=0, col="red") +
      facet_grid(mutation_type~pos_mutation) +
      scale_fill_gradientn(colours=myPalette(100), guide=FALSE) +
      labs(title=paste0("Refinement plot for ",clean_motif," motifs")) +
      labs(y="Mean current differences (pA)") +
      coord_cartesian(xlim=c(xmin_value, xmax_value), ylim=c(ymin_value, ymax_value), expand=FALSE) +
      theme_bw()

    if(filter_iso){
      gp_detail <- gp_detail + labs(x=paste0("Distance from modified base at isolated ",clean_motif," motifs"))
    }else{
      gp_detail <- gp_detail + labs(x=paste0("Distance from modified base at ",clean_motif," motifs"))
    }
    # saveRDS(gp_detail, file=paste0(motif_detail, ".gpa"))

    # Plot mutated motif scores
    gp_score <- ggplot(mutated_motif_score) +
      geom_tile(aes(x=pos_mutation, y=mutation_type, fill=score)) +
      geom_text(aes(x=pos_mutation, y=mutation_type, label=mutation_type)) +
      scale_fill_gradientn(colours=myPalette(100)) +
      # labs(title=paste0(clean_motif," motif scores")) +
      labs(x="Mutated position", y="Mutated base", fill="Score") +
      coord_cartesian(expand=FALSE) +
      theme_bw() +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(legend.position="right", legend.direction="horizontal") +
      guides(fill=guide_colourbar(title.vjust=0.5, title.position="top"))
    # saveRDS(gp_score, file=paste0(motif_detail, ".gps"))
    gp_empty <- ggplot() +
      geom_blank() +
      theme_bw() +
      theme(panel.border=element_blank())

    pdf(file=NULL) # Avoid opening empty window
    gp_combine <- ggarrange(gp_detail, ggdraw(ggarrange(gp_score, gp_empty, ncol=2, widths=c(1,2))), ncol=1, heights=c(5,1.4))
    dev.off()
    saveRDS(gp_combine, file=paste0(motif_detail, ".gpo"))
  }
}

generate.plots <- function(modification_at_motifs, data_type, modification_type, filter_iso){
  if(data_type=="ONT"){
    generate.ont.plots(modification_at_motifs, modification_type, filter_iso)
  }else if(data_type=="SMRT"){
    generate.smrt.plots(modification_at_motifs, modification_type)
  }
}

processdat <- function(motif, center, modificationtype){
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
  csv2 <- modification_info
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

# Not used
render.refine.plot <- function(selected_motif){
  renderImage({  
    print_db(paste0("Rendering combine: ", selected_motif))      

    path_graph_data <- paste0(selected_motif, ".gpa")
    gpa <- readRDS(file=path_graph_data)
  
    outfile <- tempfile(fileext='.png')
    ggsave(outfile, gpa, width=nchar(selected_motif)*2.8, height=9) # TODO remove added char
    
    list(src = outfile, contentType = 'image/png')
  }, deleteFile = TRUE)
}

prepare.download <- function(input, selected_motif){
  downloadHandler(
    filename = function(){
      selected_motif <- input$render_results

      return(paste0(selected_motif, '_combined.pdf'))
    },
    content = function(file){
      selected_motif <- input$render_results

      path_graph_data <- paste0(selected_motif, ".gpa")
      gpa <- readRDS(file=path_graph_data)
    
      ggsave(file, gpa, width=nchar(selected_motif)*2.8, height=8.5, device="pdf")
    } # TODO remove added char
  )
}

