#' @include utilities.R
#' @title Detect splicing events and generate event id
#' @description Detect splicing events based on Pseudobulk bam
#' in work_path/data and generate event id saved to work_path/events

#'
#' @param paras list readSCSESconfig(paras_file)
#' @param star_ref_path path to STAR reference.
#' Default: NULL, the STAR reference will be generated.
#' Providing STAR reference can speed up the function.

#' @return Splicing event directory
#' @export

detectEvents <- function(paras,star_ref_path=NULL) {
  options("scipen" = 100)
  # input
  print("Checking cells...")
  bam_path <- paste0(paras$Basic$work_path, "/data/")
  cells <- list.files(bam_path, "*bam$", recursive = F)
  if (length(cells) == 0) {
    msg <- paste("There is no bam file in", bam_path, "!")
    stop(msg)
  }
  print("Checking events...")
  event_types <- paras$Task$event$event_type
  event_types <- unlist(strsplit(event_types, ";"))
  event_types <- check.valid(x = event_types, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
  print(paste0("event_type=", paste(event_types, collapse = ";"), "  checked"))
  SS <- ifelse(length(intersect(event_types, c("A3SS", "A5SS", "AL"))) > 0, TRUE, FALSE)
  SEMX <- ifelse(length(intersect(event_types, c("SE", "MXE"))) > 0, TRUE, FALSE)
  RI <- ifelse(length(intersect(event_types, "RI")) > 0, TRUE, FALSE)
  # parameter
  gff <- paras$Basic$refgenome$gff_path
  gtf <- paras$Basic$refgenome$gtf_path
  ref <- paras$Basic$refgenome$ref_path
  genome_name <- paras$Basic$refgenome$genome_name
  readlength <- paras$Basic$readlength
  paired <- paras$Basic$paired
  remove.chr <- paras$Task$event$remove_chr
  core <- paras$Basic$core
  ExonToIntronReads <- paras$Basic$filter_merged_bam$ExonToIntronReads
  junctionReads <- paras$Basic$filter_merged_bam$junctionReads

  # output
  work_path <- paste0(paras$Basic$work_path, "/events/")
  print(paste0("Output: ", work_path))
  msg <- paste0("[", Sys.time(), "] ", "Creating events directory...")
  print(msg)
  dir.create(path = work_path)
  log_file <- paste0(work_path, "/detectEvents.log")
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  # script
  dir_shell <- system.file("shell", package = "SCSES")
  res.file <- list.files(paste0(work_path,'/majiq/'),
                         pattern = "alt3prime.tsv",
                         recursive = TRUE, full.names = TRUE)
  if (SS & length(res.file)==0) {
    majiq.work.path <- paste0(work_path, "/majiq/")
    msg <- paste0("[", Sys.time(), "] ", "MAJIQ...")
    print(msg)
    dir.create(path = majiq.work.path, recursive = T)
    # script
    filename <- "run_majiq.sh"
    script_majiq <- file.path(dir_shell, filename)
    license_file <-paras$Task$event$majiq_license_file
    condabin_path <- paras$Basic$conda_binpath
    majiq_env <- paras$Basic$MAJIQ_env
    res <- getA35SSALevent(
      majiq.work.path, bam_path, gff, script_majiq,
      core, genome_name, junctionReads, log_file,
      majiq_env,condabin_path,license_file
    )
  }
  res.file <- list.files(paste0(work_path,'/IRFinder/'),
                         pattern = "IRFinder-IR-nondir.txt",
                         recursive = TRUE, full.names = TRUE)
  if (RI & length(res.file)==0) {
    irfinder.work.path <- paste0(work_path, "/IRFinder/")
    msg <- paste0("[", Sys.time(), "] ", "IRFinder...")
    print(msg)
    dir.create(path = irfinder.work.path, recursive = T)
    # script
    filename <- "run_irfinder.sh"
    script_irfinder <- file.path(dir_shell, filename)
    irfinder_path <- paras$Basic$IRFinder_path
    samtools_path <- paras$Basic$samtools_path
    star_path <- paras$Basic$STAR_path
    res <- getRIevent(
      irfinder.work.path, bam_path, gtf, ref, core,readlength,
      script_irfinder, irfinder_path, samtools_path,
      star_path, log_file
    )
  }
  res.file <- list.files(paste0(work_path,'/rMats/'),
                         pattern = "SE.MATS.JCEC.txt",
                         recursive = TRUE, full.names = TRUE)
  if (SEMX & length(res.file)==0) {
    rmats.work.path <- paste0(work_path, "/rMats/")
    msg <- paste0("[", Sys.time(), "] ", "rMats...")
    print(msg)
    dir.create(path = rmats.work.path, recursive = T)
    # script
    filename <- "run_rmats.sh"
    script_rmats <- file.path(dir_shell, filename)
    rmats_path <- paras$Basic$rMATS_path
    py_path <- paras$Basic$python_path
    res <- getSEevent(
      rmats.work.path, bam_path, gtf, paired, readlength,
      rmats_path, py_path, script_rmats, core, log_file
    )
  }
  event_num_list <- lapply(event_types, function(type) {
    fun <- get(paste0("get", type, "id"))
    msg <- paste0("[", Sys.time(), "] ", "Generating ", type, " event id")
    print(msg)
    if (type == "RI") {
      event_num <- fun(
        file = paste0(work_path, "/IRFinder/"),
        outfile = paste0(work_path, "/", type, ".txt"),
        remove.chr = as.logical(remove.chr),
        base = 5,
        junctionReads = junctionReads,
        ExonToIntronReads = ExonToIntronReads,
        gtf = gtf
      )
    } else if (type %in% c("SE", "MXE")) {
      event_num <- fun(
        file = paste0(work_path, "/rMats/"),
        outfile = paste0(work_path, "/", type, ".txt"),
        remove.chr = as.logical(remove.chr),
        junctionReads = junctionReads
      )
    } else if (type %in% c("A3SS", "A5SS", "AL")) {
      event_num <- fun(
        file = paste0(work_path, "/majiq/"),
        outfile = paste0(work_path, "/", type, ".txt"),
        remove.chr = as.logical(remove.chr),
        junctionReads = junctionReads
      )
    }
    msg <- paste0("[", Sys.time(), "] ", "Generating ", type, " event id Finish.")
    print(msg)
    return(event_num)
  })
  event_num_list <- unlist(event_num_list)
  names(event_num_list) <- event_types
  msg <- paste(c("Total splicing event:", paste(names(event_num_list), event_num_list, sep = "=")),
    collapse = " "
  )
  print(msg)
  return(work_path)
}


#' @title Detect a specific type of splicing events
#' @description Detect a specific type of splicing events based on Pseudobulk
#' bam in work_path/data and generate event id saved to work_path/events

#' @param paras list readSCSESconfig(paras_file)
#' @param event_type splicing event type (A3SS,A5SS,AL,SE,MXE,RI)

#' @return Splicing event directory
#' @export
#'
getEvent <- function(paras, event_type) {
  options("scipen" = 100)
  bam_path <- paste0(paras$Basic$work_path, "/data/")
  work_path <- paste0(paras$Basic$work_path, "/events/")
  msg <- paste0("[", Sys.time(), "] ", "Creating events directory...")
  print(msg)
  dir.create(path = work_path)
  log_file <- paste0(work_path, "/detectEvents_", event_type, ".log")
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  gff <- paras$Basic$refgenome$gff_path
  gtf <- paras$Basic$refgenome$gtf_path
  ref <- paras$Basic$refgenome$ref_path
  genome_name <- paras$Basic$refgenome$genome_name
  readlength <- paras$Basic$readlength
  paired <- paras$Basic$paired
  core <- paras$Basic$core
  junctionReads <- paras$Basic$filter_merged_bam$junctionReads
  license_file <-paras$Task$event$majiq_license_file
  if (!event_type %in% c("A3SS", "A5SS", "AL", "SE", "MXE", "RI")) {
    print(paste(c("Invalid type:", event_type, "Supported splicing event types: A3SS,A5SS,AL,SE,MXE,RI."), collapse = " "))
  } else {
    msg <- paste(c("Detecting splicing event types:", event_type), collapse = " ")
    print(msg)
    SS <- ifelse(length(intersect(event_type, c("A3SS", "A5SS", "AL"))) > 0, TRUE, FALSE)
    SEMX <- ifelse(length(intersect(event_type, c("SE", "MXE"))) > 0, TRUE, FALSE)
    RI <- ifelse(length(intersect(event_type, "RI")) > 0, TRUE, FALSE)
    # script
    dir_shell <- system.file("shell", package = "SCSES")
    if (SS) {
      majiq.work.path <- paste0(work_path, "/majiq/")
      msg <- paste0("[", Sys.time(), "] ", "MAJIQ...")
      print(msg)
      dir.create(path = majiq.work.path, recursive = T)
      # script
      filename <- "run_majiq.sh"
      script_majiq <- file.path(dir_shell, filename)
      condabin_path <- paras$Basic$conda_binpath
      majiq_env <- paras$Basic$MAJIQ_env
      res <- getA35SSALevent(
        majiq.work.path, bam_path, gff, script_majiq,
        core, genome_name, junctionReads, log_file,
        majiq_env,condabin_path,license_file
      )
    }
    if (RI) {
      irfinder.work.path <- paste0(work_path, "/IRFinder/")
      msg <- paste0("[", Sys.time(), "] ", "IRFinder...")
      print(msg)
      dir.create(path = irfinder.work.path, recursive = T)
      # script
      filename <- "run_irfinder.sh"
      script_irfinder <- file.path(dir_shell, filename)
      irfinder_path <- paras$Basic$IRFinder_path
      samtools_path <- paras$Basic$samtools_path
      star_path <- paras$Basic$STAR_path
      res <- getRIevent(
        irfinder.work.path, bam_path, gtf, ref, core,readlength,
        script_irfinder, irfinder_path, samtools_path,
        star_path, log_file
      )
    }
    if (SEMX) {
      rmats.work.path <- paste0(work_path, "/rMats/")
      msg <- paste0("[", Sys.time(), "] ", "rMats...")
      print(msg)
      dir.create(path = rmats.work.path, recursive = T)
      # script
      filename <- "run_rmats.sh"
      script_rmats <- file.path(dir_shell, filename)
      rmats_path <- paras$Basic$rMATS_path
      py_path <- paras$Basic$python_path
      res <- getSEevent(
        rmats.work.path, bam_path, gtf, paired, readlength,
        rmats_path, py_path, script_rmats, core, log_file
      )
    }
  }
  return(work_path)
}

#' @title Detect SE events
#'
#' @param work_path directory to deposite all SCSES results in
#' @param bam_path directory to Pseudobulk bam file
#' @param gtf the gene annotation in gtf format
#' @param readlength the length of each read
#' @param paired type of read: either "paired" for paired-end  or "single" for single-end
#' @param core the number of threads
#' @param script_rmats directory to script for running rMats
#' @param rmats_path directory to rmats.py
#' @param py_path directory to python for running rMats
#' @param log_file file saving stdout and stderr information
#'
#' @return Splicing event directory
#'
#' @export

getSEevent <- function(work_path, bam_path, gtf, paired, readlength, rmats_path, py_path, script_rmats, core, log_file) {
  cmd <- paste(
    "bash", script_rmats,
    work_path,
    bam_path,
    gtf,
    paired,
    readlength,
    core,
    rmats_path,
    py_path, ">>", log_file, "2>&1"
  )
  msg <- paste0("[", Sys.time(), "] ", "Run rMats: ", cmd)
  print(msg)
  res <- system(command = cmd, intern = T, wait = T)
  if (!is.null(attributes(res)) && attributes(res)$status == 1) {
    msg <- paste0("[", Sys.time(), "] ", "Run rMats Error.")
    stop(msg)
  } else {
    msg <- paste0("[", Sys.time(), "] ", "Run rMats Finish.")
    print(msg)
    return(rmats_path)
  }
}

#' @title Detect RI events
#'
#' @param work_path directory to deposite all SCSES results in
#' @param bam_path directory to Pseudobulk bam file
#' @param gtf the gene annotation in gtf format
#' @param ref fasta file
#' @param core the number of threads
#' @param readlength the length of each read
#' @param script_irfinder directory to script for running IRFinder
#' @param irfinder_path directory to executable file of IRFinder
#' @param samtools_path directory to samtools
#' @param star_path directory to STAR
#' @param log_file file saving stdout and stderr information
#'
#' @return Splicing event directory
#'
#' @export

getRIevent <- function(
    work_path, bam_path, gtf, ref, core,readlength,
    script_irfinder, irfinder_path, samtools_path,
    star_path, log_file, star_ref_path=NULL) {
  old.path=Sys.getenv("PATH")
  Sys.setenv(PATH=paste0(dirname(samtools_path),":",old.path))
  if(dir.exist(paste0(work_path,'/REF'))){
    stop(paste0(work_path,'/REF should not yet exist!'))
  }
  cmd <- paste(
    "bash", script_irfinder,
    work_path,
    ref,
    gtf,
    bam_path,
    core,
    readlength,
    irfinder_path,
    samtools_path,
    star_path, star_ref_path, ">>", log_file, "2>&1"
  )
  msg <- paste0("[", Sys.time(), "] ", "Run IRFinder: ", cmd)
  print(msg)
  res <- system(command = cmd, intern = T, wait = T)
  if (!is.null(attributes(res)) && attributes(res)$status == 1) {
    msg <- paste0("[", Sys.time(), "] ", "Run IRFinder Error.")
    stop(msg)
  } else {
    res.file <- list.files(work_path,
      pattern = "IRFinder-IR-nondir.txt",
      recursive = TRUE, full.names = TRUE
    )
    print(res.file)
    if(file.exists(res.file)){
      msg <- paste0("[", Sys.time(), "] ", "Run IRFinder Finish.")
      print(msg)
    } else{
      irfinder.stderr <- list.files(work_path,
        pattern = "irfinder.stderr",
        recursive = TRUE, full.names = TRUE
      )
      file_content <- readLines(irfinder.stderr)
      cat(file_content, sep = "\n")
      msg <- paste0("[", Sys.time(), "] ", paste("Run IRFinder Error. See", irfinder.stderr, "."))
      stop(msg)
    }
    return(work_path)
  }
}

#' @title Detect A3SS,A5SS,AL events
#'
#' @param work_path directory to deposite all SCSES results in
#' @param bam_path directory to Pseudobulk bam file
#' @param gff the gene annotation in gff format
#' @param genome_name Genome name required by MAJIQ
#' @param core the number of threads
#' @param junctionReads the minimum total number of reads for any junction
#' @param script_majiq directory to script for running MAJIQ
#' @param log_file file saving stdout and stderr information
#' @param license_file the license file of MAJIQ
#' @param majiq_env the conda environment name of MAJIQ
#' @param condabin_path directory to The bin directory in conda contains the executable binary files
#'
#' @return Splicing event directory
#' @export

getA35SSALevent <- function(work_path, bam_path, gff, script_majiq,
                            core, genome_name, junctionReads, log_file,
                            majiq_env,condabin_path,license_file) {
  cmd <- paste(
    "bash", script_majiq,
    work_path,
    bam_path,
    gff,
    core,
    genome_name,
    junctionReads,
    majiq_env,
    condabin_path,
    license_file, ">>", log_file, "2>&1"
  )
  msg <- paste0("[", Sys.time(), "] ", "Run MAJIQ: ", cmd)
  print(msg)
  res <- system(command = cmd, intern = T, wait = T)
  if (!is.null(attributes(res)) && attributes(res)$status == 1) {
    msg <- paste0("[", Sys.time(), "] ", "Run MAJIQ Error.")
    stop(msg)
  } else {
    msg <- paste0("[", Sys.time(), "] ", "Run MAJIQ Finish.")
    print(msg)
    return(work_path)
  }
}


#' @title Generate SE event id
#'
#' @param file directory to rMats output
#' @param outfile directory to output file SE.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction

#' @return The number of SE events
#' @export
#' @importFrom utils read.table write.table

getSEid <- function(file, outfile, remove.chr = F, junctionReads) {
  file2 <- system(paste0("find ", file, " -name SE.MATS.JCEC.txt"), show.output.on.console = F, intern = T)
  data <- read.table(file2, header = T, sep = "\t")
  if (remove.chr) {
    data$chr <- sub(pattern = "chr", replacement = "", x = data$chr)
  }
  data <- data[which(data$IJC_SAMPLE_1 > junctionReads | data$SJC_SAMPLE_1 > junctionReads), ]
  if (nrow(data) > 0) {
    chr <- data$chr
    strand <- data$strand
    exon1 <- paste0(data$upstreamES + 1, "-", data$upstreamEE)
    exon2 <- paste0(data$exonStart_0base + 1, "-", data$exonEnd)
    exon3 <- paste0(data$downstreamES + 1, "-", data$downstreamEE)
    iso1 <- paste0(data$upstreamEE + 1, "-", data$downstreamES)
    iso2_1 <- paste0(data$upstreamEE + 1, "-", data$exonStart_0base)
    iso2_2 <- paste0(data$exonEnd + 1, "-", data$downstreamES)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon3, ":", strand)
    iso2_id <- paste0("isoform2=junction:", chr, ":", iso2_1, ":", strand, "@exon:", chr, ":", exon2, ":", strand, "@junction:", chr, ":", iso2_2, ":", strand)
    event_id <- paste(iso1_id, iso2_id, sep = "|")
    event_id <- paste(event_id, data$geneSymbol, "SE", sep = "|")
    write.table(x = unique(event_id), file = outfile, append = F, quote = F, col.names = F, row.names = F)
  } else {
    print(paste0("NO SE events under the cutoff: junctionReads = ", junctionReads))
  }
  return(nrow(data))
}

#' @title Generate RI event id
#'
#' @param file directory to IRFinder output
#' @param base default=5
#' @param outfile directory to output file RI.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction
#' @param ExonToIntronReads the minimum total number of reads that overlap the 3'/5' flanking exon and the intron
#' @param gtf the gene annotation in gtf format,used to get exon position

#' @return The number of RI events
#' @export
#' @import rtracklayer

getRIid <- function(file, base = 5, outfile, remove.chr = F, junctionReads, ExonToIntronReads, gtf) {
  file2 <- system(paste0("find ", file, " -name IRFinder-IR-nondir.txt"), show.output.on.console = F, intern = T)
  data <- read.table(file2, header = T, sep = "\t")
  data <- data[which((data$ExonToIntronReadsLeft > ifelse(ExonToIntronReads / 10 > 10, ExonToIntronReads / 10, 10) | data$ExonToIntronReadsRight > ifelse(ExonToIntronReads / 10 > 10, ExonToIntronReads / 10, 10)) & data$SpliceExact > junctionReads & data$IntronDepth > 0 & data$Warnings == "-"), ]
  if (nrow(data) > 0) {
    if (remove.chr) {
      data$Chr <- sub(pattern = "chr", replacement = "", x = data$Chr)
    }
    gtf <- rtracklayer::import(gtf)
    gtf <- as.data.frame(gtf)
    exon <- gtf[which(gtf[, "type"] == "exon"), 1:3]
    exon$exon_id1 <- paste(exon[, 1], exon[, 2], sep = ":")
    exon$exon_id2 <- paste(exon[, 1], exon[, 3], sep = ":")
    data$exon_id1 <- paste(data[, 1], data[, 3] + 1, sep = ":")
    data$exon_id2 <- paste(data[, 1], data[, 2], sep = ":")
    data$exon1 <- paste(exon[match(data$exon_id1, exon$exon_id1), 2], exon[match(data$exon_id1, exon$exon_id1), 3], sep = "-")
    data$exon2 <- paste(exon[match(data$exon_id2, exon$exon_id2), 2], exon[match(data$exon_id2, exon$exon_id2), 3], sep = "-")
    junction <- paste0("junction:", data$Chr, ":", data$Start + 1, "-", data$End, ":", data$Strand)
    retention1 <- paste0("retention:", data$Chr, ":", data$Start - base + 1, "-", data$Start + base, ":", data$Strand)
    retention2 <- paste0("retention:", data$Chr, ":", data$End - base + 1, "-", data$End + base, ":", data$Strand)
    exon1 <- paste0("exon:", data$Chr, ":", data$exon2, ":", data$Strand)
    exon2 <- paste0("exon:", data$Chr, ":", data$exon1, ":", data$Strand)
    geneSymbol <- sapply(data$Name, function(x) {
      unlist(strsplit(x, "[/]"))[1]
    })
    events <- paste0("isoform1=", exon1, "@", junction, "@", exon2, "|", "isoform2=", retention1, "@", retention2, "|", geneSymbol, "|", "RI")
    write.table(x = unique(events), file = outfile, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  } else {
    print(paste0("NO RI events under the cutoff: junctionReads = ", junctionReads, "; ExonToIntronReads", ExonToIntronReads))
  }
  return(nrow(data))
}

#' @title Generate A3SS event id
#'
#' @param file directory to MAJIQ output
#' @param outfile directory to output file A3SS.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction

#' @return The number of A3SS events
#' @export
#'
getA3SSid <- function(file, outfile, remove.chr = F, junctionReads) {
  file2 <- system(paste0("find ", file, " -name alt3prime.tsv"), intern = T)
  data <- read.table(file2, header = T, sep = "\t", comment.char = "#")
  data <- data[, 1:15]
  if (nrow(data) > 0) {
    data <- data[which(data$junction_name == "Distal"), ]
    if (remove.chr) {
      data$seqid <- sub(pattern = "chr", replacement = "", x = data$seqid)
    }
    id_exon_uniq <- apply(data, 1, function(x) {
      id <- paste(c(x["reference_exon_coord"], x["spliced_with_coord"], x["junction_coord"]), collapse = "-")
      id_ordered <- sort(as.numeric(unique(unlist(strsplit(id, "-")))))
      if (length(id_ordered) == 5 & length(which(id_ordered == 1)) == 0) {
        return(data.frame(
          event_id = x["event_id"], gene_name = x["gene_name"], seqid = x["seqid"], strand = x["strand"],
          A = id_ordered[1], B = id_ordered[2], C = id_ordered[3], D = id_ordered[4], E = id_ordered[5]
        ))
      } else {
        return(NULL)
      }
    })
    id_exon_uniq <- id_exon_uniq[which(sapply(id_exon_uniq, function(x) !is.null(x)))]
    id_exon_uniq <- do.call(what = rbind, args = id_exon_uniq)
    id_exon_uniq1 <- subset(id_exon_uniq, strand == "+")
    chr <- id_exon_uniq1$seqid
    strand <- id_exon_uniq1$strand
    exon1 <- paste0(id_exon_uniq1$A, "-", id_exon_uniq1$B)
    exon2 <- paste0(id_exon_uniq1$D, "-", id_exon_uniq1$E)
    exon3 <- paste0(id_exon_uniq1$C, "-", id_exon_uniq1$E)
    iso1 <- paste0(id_exon_uniq1$B + 1, "-", id_exon_uniq1$D - 1)
    iso2 <- paste0(id_exon_uniq1$B + 1, "-", id_exon_uniq1$C - 1)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    iso2_id <- paste0("isoform2=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon3, ":", strand)
    event_id1 <- paste(iso1_id, iso2_id, sep = "|")
    event_id1 <- paste(event_id1, id_exon_uniq1$gene_name, "A3SS", sep = "|")
    id_exon_uniq2 <- subset(id_exon_uniq, strand == "-")
    chr <- id_exon_uniq2$seqid
    strand <- id_exon_uniq2$strand
    exon1 <- paste0(id_exon_uniq2$A, "-", id_exon_uniq2$B)
    exon2 <- paste0(id_exon_uniq2$D, "-", id_exon_uniq2$E)
    exon3 <- paste0(id_exon_uniq2$A, "-", id_exon_uniq2$C)
    iso1 <- paste0(id_exon_uniq2$B + 1, "-", id_exon_uniq2$D - 1)
    iso2 <- paste0(id_exon_uniq2$C + 1, "-", id_exon_uniq2$D - 1)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    iso2_id <- paste0("isoform2=exon:", chr, ":", exon3, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    event_id2 <- paste(iso1_id, iso2_id, sep = "|")
    event_id2 <- paste(event_id2, id_exon_uniq2$gene_name, "A3SS", sep = "|")
    event_id <- c(event_id1, event_id2)
    write.table(x = unique(event_id), file = outfile, append = F, quote = F, col.names = F, row.names = F)
  } else {
    print(paste0("NO A3SS events under the cutoff: junctionReads = ", junctionReads))
  }
  return(nrow(data))
}

#' @title Generate AL event id
#'
#' @param file directory to MAJIQ output
#' @param outfile directory to output file AL.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction

#' @return The number of AL events
#' @export
#'

getALid <- function(file, outfile, remove.chr = F, junctionReads) {
  file2 <- system(paste0("find ", file, " -name alternate_last_exon.tsv"), intern = T)
  data <- read.table(file2, header = T, sep = "\t", comment.char = "#")
  data <- data[, 1:15]
  if (nrow(data) > 0) {
    if (remove.chr) {
      data$seqid <- sub(pattern = "chr", replacement = "", x = data$seqid)
    }
    event_id_all <- lapply(unique(data$event_id), function(x) {
      tmp <- data[which(data$event_id == x), ]
      chr <- tmp$seqid[1]
      strand <- tmp$strand[1]
      gene_name <- tmp$gene_name[1]
      exon1 <- tmp[1, "reference_exon_coord"]
      exon2 <- tmp[2, "spliced_with_coord"]
      exon3 <- tmp[1, "spliced_with_coord"]
      iso1 <- tmp[2, "junction_coord"]
      iso2 <- tmp[1, "junction_coord"]
      id <- paste(c(exon1, exon2, exon3, iso1, iso2), collapse = "-")
      id_ordered <- sort(as.numeric(unique(unlist(strsplit(id, "-")))))
      if (sum(id_ordered == 1 | id_ordered == (-1)) == 0) {
        iso1 <- unlist(strsplit(iso1, "-"))
        iso1 <- paste(as.numeric(iso1[1]) + 1, as.numeric(iso1[2]) - 1, sep = "-")
        iso2 <- unlist(strsplit(iso2, "-"))
        iso2 <- paste(as.numeric(iso2[1]) + 1, as.numeric(iso2[2]) - 1, sep = "-")
        if (tmp[1, 5] == "+") {
          iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
          iso2_id <- paste0("isoform2=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon3, ":", strand)
          event_id <- paste(iso1_id, iso2_id, sep = "|")
          event_id <- paste(event_id, gene_name, "AL", sep = "|")
        } else {
          iso1_id <- paste0("isoform1=exon:", chr, ":", exon2, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon1, ":", strand)
          iso2_id <- paste0("isoform2=exon:", chr, ":", exon3, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon1, ":", strand)
          event_id <- paste(iso1_id, iso2_id, sep = "|")
          event_id <- paste(event_id, gene_name, "AL", sep = "|")
        }
        return(event_id)
      }
    })
    event_id_all <- unlist(event_id_all)
    write.table(x = unique(event_id_all), file = outfile, append = F, quote = F, col.names = F, row.names = F)
  } else {
    print(paste0("NO AL events under the cutoff: junctionReads = ", junctionReads))
  }
  return(nrow(data))
}

#' @title Generate A5SS event id
#'
#' @param file directory to MAJIQ output
#' @param outfile directory to output file A5SS.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction

#' @return The number of A5SS events
#' @export
#'
getA5SSid <- function(file, outfile, remove.chr = F, junctionReads) {
  file2 <- system(paste0("find ", file, " -name alt5prime.tsv"), intern = T)
  data <- read.table(file2, header = T, sep = "\t", comment.char = "#")
  if (nrow(data) > 0) {
    data <- data[which(data$junction_name == "Distal"), ]
    if (remove.chr) {
      data$seqid <- sub(pattern = "chr", replacement = "", x = data$seqid)
    }
    id_exon_uniq <- apply(data, 1, function(x) {
      id <- paste(c(x["reference_exon_coord"], x["spliced_with_coord"], x["junction_coord"]), collapse = "-")
      id_ordered <- sort(as.numeric(unique(unlist(strsplit(id, "-")))))
      if (length(id_ordered) == 5 & length(which(id_ordered == 1)) == 0) {
        return(data.frame(
          event_id = x["event_id"], gene_name = x["gene_name"], seqid = x["seqid"], strand = x["strand"],
          A = id_ordered[1], B = id_ordered[2], C = id_ordered[3], D = id_ordered[4], E = id_ordered[5]
        ))
      } else {
        return(NULL)
      }
    })
    id_exon_uniq <- id_exon_uniq[which(sapply(id_exon_uniq, function(x) !is.null(x)))]
    id_exon_uniq <- do.call(what = rbind, args = id_exon_uniq)
    id_exon_uniq1 <- subset(id_exon_uniq, strand == "-")
    chr <- id_exon_uniq1$seqid
    strand <- id_exon_uniq1$strand
    exon1 <- paste0(id_exon_uniq1$A, "-", id_exon_uniq1$B)
    exon2 <- paste0(id_exon_uniq1$D, "-", id_exon_uniq1$E)
    exon3 <- paste0(id_exon_uniq1$C, "-", id_exon_uniq1$E)
    iso1 <- paste0(id_exon_uniq1$B + 1, "-", id_exon_uniq1$D - 1)
    iso2 <- paste0(id_exon_uniq1$B + 1, "-", id_exon_uniq1$C - 1)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    iso2_id <- paste0("isoform2=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon3, ":", strand)
    event_id1 <- paste(iso1_id, iso2_id, sep = "|")
    event_id1 <- paste(event_id1, id_exon_uniq1$gene_name, "A5SS", sep = "|")
    id_exon_uniq2 <- subset(id_exon_uniq, strand == "+")
    chr <- id_exon_uniq2$seqid
    strand <- id_exon_uniq2$strand
    exon1 <- paste0(id_exon_uniq2$A, "-", id_exon_uniq2$B)
    exon2 <- paste0(id_exon_uniq2$D, "-", id_exon_uniq2$E)
    exon3 <- paste0(id_exon_uniq2$A, "-", id_exon_uniq2$C)
    iso1 <- paste0(id_exon_uniq2$B + 1, "-", id_exon_uniq2$D - 1)
    iso2 <- paste0(id_exon_uniq2$C + 1, "-", id_exon_uniq2$D - 1)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    iso2_id <- paste0("isoform2=exon:", chr, ":", exon3, ":", strand, "@junction:", chr, ":", iso2, ":", strand, "@exon:", chr, ":", exon2, ":", strand)
    event_id2 <- paste(iso1_id, iso2_id, sep = "|")
    event_id2 <- paste(event_id2, id_exon_uniq2$gene_name, "A5SS", sep = "|")
    event_id <- c(event_id1, event_id2)
    write.table(x = unique(event_id), file = outfile, append = F, quote = F, col.names = F, row.names = F)
  } else {
    print(paste0("NO A5SS events under the cutoff: junctionReads = ", junctionReads))
  }
  return(nrow(data))
}

#' @title Generate MXE event id
#'
#' @param file directory to rMats output
#' @param outfile directory to output file MXE.txt
#' @param remove.chr if TRUE,remove the "chr" character of event id
#' @param junctionReads the minimum total number of reads for any junction

#' @return The number of MXE events
#' @export
#'

getMXEid <- function(file, outfile, remove.chr = F, junctionReads) {
  file2 <- system(paste0("find ", file, " -name MXE.MATS.JCEC.txt"), intern = T)
  data <- read.table(file2, header = T, sep = "\t")
  if (remove.chr) {
    data$chr <- sub(pattern = "chr", replacement = "", x = data$chr)
  }
  data <- data[which(data$IJC_SAMPLE_1 > junctionReads | data$SJC_SAMPLE_1 > junctionReads), ]
  if (nrow(data) > 0) {
    chr <- data$chr
    strand <- data$strand
    exon1 <- paste0(data$upstreamES + 1, "-", data$upstreamEE)
    exon2 <- paste0(data$X1stExonStart_0base + 1, "-", data$X1stExonEnd)
    exon3 <- paste0(data$X2ndExonStart_0base + 1, "-", data$X2ndExonEnd)
    exon4 <- paste0(data$downstreamES + 1, "-", data$downstreamEE)
    iso1_1 <- paste0(data$upstreamEE + 1, "-", data$X2ndExonStart_0base)
    iso1_2 <- paste0(data$X2ndExonEnd + 1, "-", data$downstreamES)
    iso2_1 <- paste0(data$upstreamEE + 1, "-", data$X1stExonStart_0base)
    iso2_2 <- paste0(data$X1stExonEnd + 1, "-", data$downstreamES)
    iso1_id <- paste0("isoform1=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso1_1, ":", strand, "@exon:", chr, ":", exon3, ":", strand, "@junction:", chr, ":", iso1_2, ":", strand, "@exon:", chr, ":", exon4, ":", strand)
    iso2_id <- paste0("isoform2=exon:", chr, ":", exon1, ":", strand, "@junction:", chr, ":", iso2_1, ":", strand, "@exon:", chr, ":", exon2, ":", strand, "@junction:", chr, ":", iso2_2, ":", strand, "@exon:", chr, ":", exon4, ":", strand)
    event_id <- paste(iso1_id, iso2_id, sep = "|")
    event_id <- paste(event_id, data$geneSymbol, "MXE", sep = "|")
    write.table(x = unique(event_id), file = outfile, append = F, quote = F, col.names = F, row.names = F)
  } else {
    print(paste0("NO MXE events under the cutoff: junctionReads = ", junctionReads))
  }
  return(nrow(data))
}
