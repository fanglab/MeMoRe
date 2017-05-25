cat("Prerequisite libraries: data.table, R.utils, Biostrings, stringr, parallel, qpcR, ggplot2")

if (length(args)==0) {
  stop("parameters.txt argument must be supplied", call.=FALSE)
}

library(data.table)
library(R.utils)
library("Biostrings")
library(stringr)
library(parallel)
library(qpcR)
library(ggplot2)

motif_ref <- data.frame(sym = c("W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"),
                              bases = c("(A|T)", "(C|G)", "(A|C)", "(G|T)", "(A|G)", "(C|T)", "(C|G|T)", "(A|G|T)", "(A|C|T)", "(A|C|G)", "(A|C|G|T)"))
parameters = args[1]
lines <- readLines(parameters)
cl_max <- as.numeric(gsub("^.*\\|","", lines[1]))
fasta_file <-gsub("^.*\\|","", lines[2])
csv_file <- gsub("^.*\\|","", lines[3])
destination <- gsub("^.*\\|","", lines[5])
tab1 <- read.table(gsub("^.*\\|","", lines[4]), sep = ",", header = TRUE)[1:4]
tab2 <-  read.table(textConnection(lines[7:length(lines)]), header = F,sep=",",
                    col.names = colnames(tab1))
motiftable <- rbind(tab1,tab2)
modtype <- as.vector(motiftable$modificationType)
methvec <- gsub("\\d|[[:lower:]]","",modtype)
motif_f <- as.vector(motiftable$motifString)
cov_cut_spec <- 10
size_spec <- 100
unlink("graph_dirs.txt")

#Get csv_sp
  csv_temp <- gunzip(csv_file, remove = F, temporary = F, overwrite = T)
  csv <- fread(csv_temp, sep = ",", header = T, verbose = F, drop = c(6, 7, 8, 11, 12, 13))  # Set csv parameters
  gene <- readDNAStringSet(fasta_file)
  gag <- which(width(gene) == max(width(gene)))
  csv2 <- csv[csv$refName == names(gene)[gag]]

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

  exportme2 <- list(data.table(csv_sp_f), data.table(csv_sp_r), genome_f)

  # create cluster
  #cat(paste0("  + Making cluster of ", cl_max, " cores"), "\n")
  #cat("   + Exporting global variables/libraries", "\r")
  cl <- makeCluster(cl_max)
  clusterExport(cl, c("exportme2"))
  clusterEvalQ(cl, {
	library(stringi)
	library(IRanges)
	library(data.table)
  })


  for (j in 1:length(motif_f)){
    dir <- paste0(gsub("\\\\", "/", destination), "/", motif_f[j]) #, "/", dirname(fasta_file),
    unlink(dir, recursive = T)
    ifelse(!dir.exists(dir), dir.create(dir), FALSE)

    dir2 <- paste0(dir, "/data/")
    unlink(dir2, recursive = T)
    ifelse(!dir.exists(dir2), dir.create(dir2), FALSE)

    new_mots <- unlist(strsplit(motif_f[j], ""))
    nom3 <- rep(4, length(new_mots))
    nom_matrix2 <- do.call(qpcR:::cbind.na, split(1:sum(nom3), rep(1:length(nom3), nom3)))
    colnames(nom_matrix2) <- NULL
    meth <- methvec[j]
    transl_r <- transl_f <- motiftable$centerPos[j]

    #Get new_list_o
    new_list_o <- list()
    new_list_o2 <- list()
    new_mot_ref <- rbind(motif_ref, data.frame(sym = c("A", "C", "G", "T"), bases = c("A", "C", "G", "T")))
    #cat(paste0("   + Methy Checking child motifs from pos. 1-", nchar(motif_f[j])), "\n")
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
    clstart1 <- Sys.time()
    clusterExport(cl, c("exportme"))
    clend.a <- Sys.time()
    #cat("   + Merging csv and dat", "\r")
    dat2 <- clusterApply(cl, dog1$motif, function(z) {
      # VARS
      csv_sp_f <- exportme2[[1]]
      csv_sp_r <- exportme2[[2]]
      genome_f <- exportme2[[3]]
      transl_f <- exportme[[1]]
      size_spec <- exportme[[2]]
      param_sp <- exportme[[3]]
      transl_r <- transl_f
      genome_r <- chartr("GATC", "CTAG", genome_f)  # chart
      myreverse <- function(x,...){unlist(lapply(strsplit(as.vector(x),""),function(z)paste(rev(z),collapse="")))}
      ar <- paste(rev(stri_extract_all(z, regex = "\\([^)]+\\)|.")[[1]]), collapse = "")

      locations_f2 <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_f, z)), ncol = 2))[1] + transl_f,
                                 motif = unlist(stri_extract_all_regex(genome_f, z)))
      locations_r2 <- data.frame(data.frame(matrix(unlist(stri_locate_all_regex(genome_r, ar)), ncol = 2))[2] - transl_r,
                                 motif = myreverse(unlist(stri_extract_all_regex(genome_r, ar))))
      loc_f_s2 <- locations_f2[!(locations_f2$X1 < param_sp[1] | locations_f2$X1 > param_sp[2]), ]
      loc_r_s2 <- locations_r2[!(locations_r2$X2 < param_sp[1] | locations_r2$X2 > param_sp[2]), ]
      merged_f2 <- merge(data.frame(tpl = loc_f_s2$X1, motif = loc_f_s2$motif), csv_sp_f, all = F)
      merged_r2 <- merge(data.frame(tpl = loc_r_s2$X2, motif = loc_r_s2$motif), csv_sp_r, all = F)
      # merge forward and reverse motif locations
      data.table(rbind(merged_f2, merged_r2))
    })

    #WRITING POS CSVS NEEDED FOR DEBUG_GRAPH2
    for(dutu2 in 1:NROW(dat2)){
      filenamer2 <- paste0("pos", formatC(ceiling(dutu2/4), width = 2, format = "d", flag = "0"),
                           unlist(strsplit(as.character(dat2[[dutu2]]$motif[1]),""))[ceiling(dutu2/4)],".csv")
      write.csv(dat2[[dutu2]], paste0(dir2,filenamer2), row.names = F)
    }


    ##################################### BOTTOM IN R/RSHINY, TOP IN C ################################33

     # WRITING TO DEBUG_GRAPH 2 NEW
     bases <- c("A","C","G","T")
     deglist<- list()

     for(newy in 1:nchar(motif_f[j])){
       for(newy2 in 1:4){
         locpos <- bases[newy2]
         namepos <- paste0(dir,"/data/pos", sprintf("%02d",newy), locpos, ".csv")
         dataworky10 <- read.csv(namepos, sep = ",", header = T)
         deglist[[((newy-1)*4+newy2)]] <- data.frame(base = locpos, pos = newy,
                                                     values = c(dataworky10$score,
                                                                dataworky10$ipdRatio,
                                                                dataworky10$coverage),
                                                     type = c(rep("score", nrow(dataworky10)),
                                                              rep("ipdRatio", nrow(dataworky10)),
                                                              rep("coverage", nrow(dataworky10))),
                                                     count = c(rep("", nrow(dataworky10)),
                                                               rep("", nrow(dataworky10)),
                                                               rep(nrow(dataworky10), nrow(dataworky10))),
                                                     county = c(rep(max(dataworky10$score), nrow(dataworky10)),
                                                                rep(max(dataworky10$ipdRatio), nrow(dataworky10)),
                                                                rep(max(dataworky10$ipdRatio), nrow(dataworky10))))
       }
     }
     dataw <- do.call(rbind, deglist)
     valueme2 <- data.frame(valueme2 = 0.1 + max(as.numeric(as.vector(sapply(dataw[dataw$type=="ipdRatio",]$values, function (x) rep(x,3))))))
     dataworky <- cbind(dataw, valueme2)
     dataworky2 <- data.frame()
     dataworky2 <- dataworky[dataworky$base=="A",]

     write.csv(dataw, paste0(dir,"/dataw.csv"), row.names = F)
     write(paste0(dir,"/dataw.csv"), file = "graph_dirs.txt",append = TRUE)
     write(motif_f[j], file = "graph_dirs.txt",append = TRUE)
     
     prettify <- theme(panel.background = element_rect(fill = NA,color="gray"), 
                       panel.grid.major.y = element_blank(),
                       panel.grid.major.x = element_line(size=.1, color="black",linetype="dotted"), 
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_line(size=.1, color="black"),
                       legend.position="bottom")
     
    graph <- ggplot(dataw, aes (x = base, y = values, color = type, group = base)) +
      geom_violin() +
      facet_grid(type ~ pos, scales = 'free') +
      theme_gray() %+replace% prettify + 
      geom_text(aes(y=county, label=count),size = 5, colour = "grey25") + 
      ggtitle(motif_f[j])
    ggsave(filename=paste0(dir, "/", motif_f[j], ".png"), graph, width=length(unique(dataw$pos))*2.8, height=8, limitsize = F)
     
  }
