#' @title Split SE event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#'
SE.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  return(data.frame(
    event = event, exon1 = iso1.info[2], exon2 = iso2.info[3],
    exon3 = iso1.info[4], junction1 = iso2.info[2], junction2 = iso2.info[4],
    gene = gene, strand = strand
  ))
}
#' @title Split RI event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#'
RI.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  return(data.frame(
    event = event, exon1 = iso1.info[2], exon2 = iso1.info[4],
    junction1 = iso1.info[3], gene = gene, strand = strand
  ))
}
#' @title Split A3SS event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#' @importFrom GenomicRanges GRanges
#'
A3SS.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  if (strand == "+") {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[4]
    exon2.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon2))
    exon3.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon3))
    exon2 <- paste0(
      "exon:", seqnames(exon2.GR), ":", start(exon2.GR),
      "-", start(exon3.GR) - 1, ":", strand(exon2.GR)
    )
  } else {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[2]
    exon2.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon2))
    exon1.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon1))
    exon2 <- paste0("exon:", seqnames(exon2.GR), ":", end(exon1.GR) + 1, "-", end(exon2.GR), ":", strand(exon2.GR))
  }

  return(data.frame(
    event = event, exon1 = exon1, exon2 = exon2,
    exon3 = exon3, junction1 = junction,
    gene = gene, strand = strand
  ))
}
#' @title Split AL event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#'
AL.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  if (strand == "+") {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[4]
  } else {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[2]
  }

  return(data.frame(
    event = event, exon1 = exon1, exon2 = exon2,
    exon3 = exon3, junction1 = junction,
    gene = gene, strand = strand
  ))
}
#' @title Split A5SS event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#' @importFrom GenomicRanges GRanges
#'
A5SS.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  if (strand == "+") {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[2]
    exon2.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon2))
    exon1.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon1))
    exon2 <- paste0("exon:", seqnames(exon2.GR), ":", end(exon1.GR) + 1, "-", end(exon2.GR), ":", strand(exon2.GR))
  } else {
    exon1 <- iso1.info[2]
    exon3 <- iso1.info[4]
    junction <- iso2.info[3]

    exon2 <- iso2.info[4]
    exon2.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon2))
    exon3.GR <- GRanges(sub(pattern = "exon:", replacement = "", x = exon3))
    exon2 <- paste0("exon:", seqnames(exon2.GR), ":", start(exon2.GR), "-", start(exon3.GR) - 1, ":", strand(exon2.GR))
  }

  return(data.frame(
    event = event, exon1 = exon1, exon2 = exon2,
    exon3 = exon3, junction1 = junction,
    gene = gene, strand = strand
  ))
}
#' @title Split MXE event information
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
#'
MXE.info.seperation <- function(event) {
  infos <- unlist(strsplit(x = event, split = "\\|"))
  gene <- infos[[3]]
  iso1.info <- unlist(strsplit(x = infos[1], split = "isoform1=|@"))
  iso2.info <- unlist(strsplit(x = infos[2], split = "isoform2=|@"))
  strand <- ifelse(grepl(pattern = "\\+", x = iso1.info[2]), "+", "-")
  return(data.frame(
    event = event,
    exon1 = iso1.info[2],
    exon2 = iso2.info[4],
    exon3 = iso1.info[4],
    exon4 = iso2.info[6],
    junction1 = iso1.info[3],
    junction2 = iso1.info[5],
    junction3 = iso2.info[3],
    junction4 = iso2.info[5],
    gene = gene, strand = strand
  ))
}

# SE event associated functions----

#' @title Extract splice sequence for SE event
#' @param events.info data.frame SE event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings

SE.splice.seq.extraction <- function(events.info, bs.genome,core) {
  cluster <- makeCluster(spec = core)
  # 1. I1-5----
  I1.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I1.5.region <- do.call(what = rbind, I1.5.region)
  I1.5.region <- I1.5.region[I1.5.region$chr %in% names(bs.genome), ]
  I1.5.seq <- getSeq(
    x = bs.genome, names = I1.5.region$chr,
    start = I1.5.region$start,
    end = I1.5.region$end,
    strand = I1.5.region$strand
  )
  I1.5.region$seq <- as.character(I1.5.seq)
  # 2. I1-3----
  I1.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I1.3.region <- do.call(what = rbind, I1.3.region)
  I1.3.region <- I1.3.region[I1.3.region$chr %in% names(bs.genome), ]
  I1.3.seq <- getSeq(
    x = bs.genome, names = I1.3.region$chr,
    start = I1.3.region$start,
    end = I1.3.region$end,
    strand = I1.3.region$strand
  )
  I1.3.region$seq <- as.character(I1.3.seq)

  # 3. I2-5----
  I2.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I2.5.region <- do.call(what = rbind, I2.5.region)
  I2.5.region <- I2.5.region[I2.5.region$chr %in% names(bs.genome), ]
  I2.5.seq <- getSeq(
    x = bs.genome, names = I2.5.region$chr,
    start = I2.5.region$start,
    end = I2.5.region$end,
    strand = I2.5.region$strand
  )
  I2.5.region$seq <- as.character(I2.5.seq)
  # 4. I2-3----
  I2.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I2.3.region <- do.call(what = rbind, I2.3.region)
  I2.3.region <- I2.3.region[I2.3.region$chr %in% names(bs.genome), ]
  I2.3.seq <- getSeq(
    x = bs.genome, names = I2.3.region$chr,
    start = I2.3.region$start,
    end = I2.3.region$end,
    strand = I2.3.region$strand
  )
  I2.3.region$seq <- as.character(I2.3.seq)
  stopCluster(cluster)
  return(list(I1.5 = I1.5.region, I1.3 = I1.3.region, I2.5 = I2.5.region, I2.3 = I2.3.region))
}
#' @title Extract exon sequence for SE event
#' @param events.info data.frame SE event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
SE.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. exon3----
  exon3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon3.region <- do.call(what = rbind, exon3.region)
  exon3.region <- exon3.region[exon3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon3.region$chr,
    start = exon3.region$start,
    end = exon3.region$end,
    strand = exon3.region$strand
  )
  exon3.region$seq <- as.character(seq)

  stopCluster(cluster)
  return(list(exon1 = exon1.region, exon2 = exon2.region, exon3 = exon3.region))
}
#' @title get length of exons and introns associated with SE event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export

SE.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  exon3 <- unlist(strsplit(x = event$exon3, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction1, split = ":|-"))
  jc2 <- unlist(strsplit(x = event$junction2, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  exon2.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  exon3.length <- as.numeric(exon3[4]) - as.numeric(exon3[3]) + 1
  jc1.length <- as.numeric(jc1[4]) - as.numeric(jc1[3]) + 1
  jc2.length <- as.numeric(jc2[4]) - as.numeric(jc2[3]) + 1
  e2.jc1.ratio <- exon2.length / jc1.length
  e2.jc2.ratio <- exon2.length / jc2.length
  jc1.jc2.ratio <- jc1.length / jc2.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      exon2 = exon2.length,
      exon3 = exon3.length,
      jc1 = jc1.length,
      jc2 = jc2.length,
      e2.jc1 = e2.jc1.ratio,
      e2.jc2 = e2.jc2.ratio,
      jc1.jc2 = jc1.jc2.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon3.length,
      exon2 = exon2.length,
      exon3 = exon1.length,
      jc1 = jc2.length,
      jc2 = jc1.length,
      e2.jc1 = e2.jc2.ratio,
      e2.jc2 = e2.jc1.ratio,
      jc1.jc2 = 1 / jc1.jc2.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}

#' @title get kmer of exons and splice region associated with SE event
#' @param splice.region data.frame with sequence of splice region, \code{SE.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{SE.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings
SE.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # I1.5----
  features <- getKmer(splice.region[["I1.5"]], "I1.5", 2)
  kmer <- c(kmer, features)
  # I2.3----
  features <- getKmer(splice.region[["I2.3"]], "I2.3", 2)
  kmer <- c(kmer, features)
  # I1.3----
  features <- getKmer(splice.region[["I1.3"]], "I1.3", 3)
  kmer <- c(kmer, features)
  # I2.5----
  features <- getKmer(splice.region[["I2.5"]], "I2.5", 3)
  kmer <- c(kmer, features)
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 3)
  kmer <- c(kmer, features)
  # exon3----
  features <- getKmer(exon.region[["exon3"]], "exon3", 3)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["exon2"]], "exon2", 4)
  kmer <- c(kmer, features)

  return(kmer)
}
# RI event associated functions----

#' @title Extract splice sequence for RI event
#' @param events.info data.frame RI event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
RI.splice.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. I1-5----
  A.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  A.5.region <- do.call(what = rbind, A.5.region)
  A.5.region <- A.5.region[A.5.region$chr %in% names(bs.genome), ]
  A.5.seq <- getSeq(
    x = bs.genome, names = A.5.region$chr,
    start = A.5.region$start,
    end = A.5.region$end,
    strand = A.5.region$strand
  )
  A.5.region$seq <- as.character(A.5.seq)
  # 2. I1-3----
  A.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  A.3.region <- do.call(what = rbind, A.3.region)
  A.3.region <- A.3.region[A.3.region$chr %in% names(bs.genome), ]
  A.3.seq <- getSeq(
    x = bs.genome, names = A.3.region$chr,
    start = A.3.region$start,
    end = A.3.region$end,
    strand = A.3.region$strand
  )
  A.3.region$seq <- as.character(A.3.seq)

  stopCluster(cluster)
  return(list(A.5.region = A.5.region, A.3.region = A.3.region))
}
#' @title Extract exon sequence for RI event
#' @param events.info data.frame RI event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
RI.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. junction----
  junction.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  junction.region <- do.call(what = rbind, junction.region)
  junction.region <- junction.region[junction.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = junction.region$chr,
    start = junction.region$start,
    end = junction.region$end,
    strand = junction.region$strand
  )
  junction.region$seq <- as.character(seq)

  stopCluster(cluster)
  return(list(exon1 = exon1.region, exon2 = exon2.region, junction = junction.region))
}

#' @title get length of exons and introns associated with RI event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export

RI.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  exon2.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  jc1.length <- as.numeric(jc1[4]) - as.numeric(jc1[3]) + 1
  e1.jc1.ratio <- exon1.length / jc1.length
  e2.jc1.ratio <- exon1.length / jc1.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      exon2 = exon2.length,
      jc1 = jc1.length,
      e1.jc1 = e1.jc1.ratio,
      e2.jc1 = e2.jc1.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon2.length,
      exon2 = exon1.length,
      jc1 = jc1.length,
      e1.jc1 = e2.jc1.ratio,
      e2.jc1 = e1.jc1.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}

#' @title get kmer of exons and splice region associated with RI event
#' @param splice.region data.frame with sequence of splice region, \code{RI.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{RI.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings
RI.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 4)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["exon2"]], "exon2", 4)
  kmer <- c(kmer, features)
  # junction----
  features <- getKmer(exon.region[["junction"]], "junction", 4)
  kmer <- c(kmer, features)
  return(kmer)
}
# A3SS event associated functions----
#' @title Extract splice sequence for A3SS event
#' @param events.info data.frame A3SS event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
A3SS.splice.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(spec = core)
  # 1. I-5----
  I.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.5.region <- do.call(what = rbind, I.5.region)
  I.5.region <- I.5.region[I.5.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.5.region$chr,
    start = I.5.region$start,
    end = I.5.region$end,
    strand = I.5.region$strand
  )
  I.5.region$seq <- as.character(seq)
  # 2. I-3----
  I.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.3.region <- do.call(what = rbind, I.3.region)
  I.3.region <- I.3.region[I.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.3.region$chr,
    start = I.3.region$start,
    end = I.3.region$end,
    strand = I.3.region$strand
  )
  I.3.region$seq <- as.character(seq)


  # 3. A-3----
  A.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  A.3.region <- do.call(what = rbind, A.3.region)
  A.3.region <- A.3.region[A.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = A.3.region$chr,
    start = A.3.region$start,
    end = A.3.region$end,
    strand = A.3.region$strand
  )
  A.3.region$seq <- as.character(seq)
  stopCluster(cluster)
  return(list(I.5 = I.5.region, I.3 = I.3.region, A.3 = A.3.region))
}
#' @title Extract exon sequence for A3SS event
#' @param events.info data.frame A3SS event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
A3SS.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. exon3----
  exon3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon3.region <- do.call(what = rbind, exon3.region)
  exon3.region <- exon3.region[exon3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon3.region$chr,
    start = exon3.region$start,
    end = exon3.region$end,
    strand = exon3.region$strand
  )
  exon3.region$seq <- as.character(seq)

  stopCluster(cluster)
  return(list(exon1 = exon1.region, A = exon2.region, exon3 = exon3.region))
}
#' @title get length of exons and introns associated with A3SS event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
A3SS.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  exon3 <- unlist(strsplit(x = event$exon3, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction1, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  A.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  exon3.length <- as.numeric(exon3[4]) - as.numeric(exon3[3]) + 1
  jc1.length <- as.numeric(jc1[4]) - as.numeric(jc1[3]) + 1
  A.exon1.ratio <- A.length / exon1.length
  A.exon3.ratio <- A.length / exon3.length
  A.I.ratio <- A.length / jc1.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      A = A.length,
      exon3 = exon3.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon3.length,
      A = A.length,
      exon3 = exon1.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}

#' @title get kmer of exons and splice region associated with A3SS event
#' @param splice.region data.frame with sequence of splice region, \code{A3SS.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{A3SS.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings
A3SS.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # I-5----
  features <- getKmer(splice.region[["I.5"]], "I.5", 3)
  kmer <- c(kmer, features)
  # I-3----
  features <- getKmer(splice.region[["I.3"]], "I.3", 3)
  kmer <- c(kmer, features)
  # A-3----
  features <- getKmer(splice.region[["A.3"]], "A.3", 3)
  kmer <- c(kmer, features)
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 4)
  kmer <- c(kmer, features)
  # exon3----
  features <- getKmer(exon.region[["exon3"]], "exon3", 4)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["A"]], "A", 4)
  kmer <- c(kmer, features)

  return(kmer)
}
# AL event associated functions----
#' @title Extract splice sequence for AL event
#' @param events.info data.frame AL event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
AL.splice.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. I-5----
  I.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.5.region <- do.call(what = rbind, I.5.region)
  I.5.region <- I.5.region[I.5.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.5.region$chr,
    start = I.5.region$start,
    end = I.5.region$end,
    strand = I.5.region$strand
  )
  I.5.region$seq <- as.character(seq)
  # 2. I-3----
  I.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.3.region <- do.call(what = rbind, I.3.region)
  I.3.region <- I.3.region[I.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.3.region$chr,
    start = I.3.region$start,
    end = I.3.region$end,
    strand = I.3.region$strand
  )
  I.3.region$seq <- as.character(seq)


  # 3. A-3----
  A.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  A.3.region <- do.call(what = rbind, A.3.region)
  A.3.region <- A.3.region[A.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = A.3.region$chr,
    start = A.3.region$start,
    end = A.3.region$end,
    strand = A.3.region$strand
  )
  A.3.region$seq <- as.character(seq)

  stopCluster(cluster)
  return(list(I.5 = I.5.region, I.3 = I.3.region, A.3 = A.3.region))
}
#' @title Extract exon sequence for AL event
#' @param events.info data.frame AL event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
AL.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. exon3----
  exon3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon3.region <- do.call(what = rbind, exon3.region)
  exon3.region <- exon3.region[exon3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon3.region$chr,
    start = exon3.region$start,
    end = exon3.region$end,
    strand = exon3.region$strand
  )
  exon3.region$seq <- as.character(seq)
  stopCluster(cluster)
  return(list(exon1 = exon1.region, A = exon2.region, exon3 = exon3.region))
}
#' @title get length of exons and introns associated with AL event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
AL.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  exon3 <- unlist(strsplit(x = event$exon3, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction1, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  A.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  exon3.length <- as.numeric(exon3[4]) - as.numeric(exon3[3]) + 1
  jc1.length <- as.numeric(jc1[4]) - as.numeric(jc1[3]) + 1
  A.exon1.ratio <- A.length / exon1.length
  A.exon3.ratio <- A.length / exon3.length
  A.I.ratio <- A.length / jc1.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      A = A.length,
      exon3 = exon3.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon3.length,
      A = A.length,
      exon3 = exon1.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}
#' @title get kmer of exons and splice region associated with AL event
#' @param splice.region data.frame with sequence of splice region, \code{AL.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{AL.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings
AL.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # I-5----
  features <- getKmer(splice.region[["I.5"]], "I.5", 3)
  kmer <- c(kmer, features)
  # I-3----
  features <- getKmer(splice.region[["I.3"]], "I.3", 3)
  kmer <- c(kmer, features)
  # A-3----
  features <- getKmer(splice.region[["A.3"]], "A.3", 3)
  kmer <- c(kmer, features)
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 4)
  kmer <- c(kmer, features)
  # exon3----
  features <- getKmer(exon.region[["exon3"]], "exon3", 4)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["A"]], "A", 4)
  kmer <- c(kmer, features)

  return(kmer)
}
# A5SS event associated functions----
#' @title Extract splice sequence for A5SS event
#' @param events.info data.frame A5SS event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
A5SS.splice.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. I-5----
  I.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.5.region <- do.call(what = rbind, I.5.region)
  I.5.region <- I.5.region[I.5.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.5.region$chr,
    start = I.5.region$start,
    end = I.5.region$end,
    strand = I.5.region$strand
  )
  I.5.region$seq <- as.character(seq)
  # 2. I-3----
  I.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I.3.region <- do.call(what = rbind, I.3.region)
  I.3.region <- I.3.region[I.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I.3.region$chr,
    start = I.3.region$start,
    end = I.3.region$end,
    strand = I.3.region$strand
  )
  I.3.region$seq <- as.character(seq)


  # 3. A-5----
  A.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  A.5.region <- do.call(what = rbind, A.5.region)
  A.5.region <- A.5.region[A.5.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = A.5.region$chr,
    start = A.5.region$start,
    end = A.5.region$end,
    strand = A.5.region$strand
  )
  A.5.region$seq <- as.character(seq)

  stopCluster(cluster)
  return(list(I.5 = I.5.region, I.3 = I.3.region, A.5 = A.5.region))
}
#' @title Extract exon sequence for A5SS event
#' @param events.info data.frame A5SS event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
A5SS.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. exon3----
  exon3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon3.region <- do.call(what = rbind, exon3.region)
  exon3.region <- exon3.region[exon3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon3.region$chr,
    start = exon3.region$start,
    end = exon3.region$end,
    strand = exon3.region$strand
  )
  exon3.region$seq <- as.character(seq)
  stopCluster(cluster)
  return(list(exon1 = exon1.region, A = exon2.region, exon3 = exon3.region))
}
#' @title get length of exons and introns associated with A5SS event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
A5SS.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  exon3 <- unlist(strsplit(x = event$exon3, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction1, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  A.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  exon3.length <- as.numeric(exon3[4]) - as.numeric(exon3[3]) + 1
  jc1.length <- as.numeric(jc1[4]) - as.numeric(jc1[3]) + 1
  A.exon1.ratio <- A.length / exon1.length
  A.exon3.ratio <- A.length / exon3.length
  A.I.ratio <- A.length / jc1.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      A = A.length,
      exon3 = exon3.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon3.length,
      A = A.length,
      exon3 = exon1.length,
      I = jc1.length,
      A.exon1 = A.exon1.ratio,
      A.exon3 = A.exon3.ratio,
      A.I = A.I.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}
#' @title get kmer of exons and splice region associated with A5SS event
#' @param splice.region data.frame with sequence of splice region, \code{A5SS.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{A5SS.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings
A5SS.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # I-5----
  features <- getKmer(splice.region[["I.5"]], "I.5", 3)
  kmer <- c(kmer, features)
  # I-3----
  features <- getKmer(splice.region[["I.3"]], "I.3", 3)
  kmer <- c(kmer, features)
  # A-5----
  features <- getKmer(splice.region[["A.5"]], "A.5", 3)
  kmer <- c(kmer, features)
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 4)
  kmer <- c(kmer, features)
  # exon3----
  features <- getKmer(exon.region[["exon3"]], "exon3", 4)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["A"]], "A", 4)
  kmer <- c(kmer, features)

  return(kmer)
}


# MXE event associated functions----
#' @title Extract splice sequence for MXE event
#' @param events.info data.frame MXE event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
MXE.splice.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. I1-5----
  I1.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I1.5.region <- do.call(what = rbind, I1.5.region)
  I1.5.region <- I1.5.region[I1.5.region$chr %in% names(bs.genome), ]
  I1.5.seq <- getSeq(
    x = bs.genome, names = I1.5.region$chr,
    start = I1.5.region$start,
    end = I1.5.region$end,
    strand = I1.5.region$strand
  )
  I1.5.region$seq <- as.character(I1.5.seq)
  # 2. I1-3----
  I1.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I1.3.region <- do.call(what = rbind, I1.3.region)
  I1.3.region <- I1.3.region[I1.3.region$chr %in% names(bs.genome), ]
  I1.3.seq <- getSeq(
    x = bs.genome, names = I1.3.region$chr,
    start = I1.3.region$start,
    end = I1.3.region$end,
    strand = I1.3.region$strand
  )
  I1.3.region$seq <- as.character(I1.3.seq)

  # 3. I2-5----
  I2.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction4[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I2.5.region <- do.call(what = rbind, I2.5.region)
  I2.5.region <- I2.5.region[I2.5.region$chr %in% names(bs.genome), ]
  I2.5.seq <- getSeq(
    x = bs.genome, names = I2.5.region$chr,
    start = I2.5.region$start,
    end = I2.5.region$end,
    strand = I2.5.region$strand
  )
  I2.5.region$seq <- as.character(I2.5.seq)
  # 4. I2-3----
  I2.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction4[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I2.3.region <- do.call(what = rbind, I2.3.region)
  I2.3.region <- I2.3.region[I2.3.region$chr %in% names(bs.genome), ]
  I2.3.seq <- getSeq(
    x = bs.genome, names = I2.3.region$chr,
    start = I2.3.region$start,
    end = I2.3.region$end,
    strand = I2.3.region$strand
  )
  I2.3.region$seq <- as.character(I2.3.seq)

  # 5. I3-5----
  I3.5.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[3])
      start <- splice.site - 4
      end <- splice.site + 5
      strand <- "+"
    } else {
      jc <- events.info$junction3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[4])
      chr <- jc[2]
      start <- splice.site - 5
      end <- splice.site + 4
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I3.5.region <- do.call(what = rbind, I3.5.region)
  I3.5.region <- I3.5.region[I3.5.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I3.5.region$chr,
    start = I3.5.region$start,
    end = I3.5.region$end,
    strand = I3.5.region$strand
  )
  I3.5.region$seq <- as.character(seq)
  # 6. I3-3----
  I3.3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$junction2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      splice.site <- as.numeric(jc[4])
      start <- splice.site - 15
      end <- splice.site + 4
      strand <- "+"
    } else {
      jc <- events.info$junction3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      splice.site <- as.numeric(jc[3])
      chr <- jc[2]
      start <- splice.site - 4
      end <- splice.site + 15
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  I3.3.region <- do.call(what = rbind, I3.3.region)
  I3.3.region <- I3.3.region[I3.3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = I3.3.region$chr,
    start = I3.3.region$start,
    end = I3.3.region$end,
    strand = I3.3.region$strand
  )
  I3.3.region$seq <- as.character(seq)
  stopCluster(cluster)
  return(list(I1.5 = I1.5.region, I1.3 = I1.3.region, I2.5 = I2.5.region, I2.3 = I2.3.region, I3.5 = I3.5.region, I3.3 = I3.3.region))
}
#' @title Extract exon sequence for MXE event
#' @param events.info data.frame MXE event information
#' @param bs.genome \code{get(pkg)}, see \code{createBSgenome} function
#' @param core the number of threads
#' @return list
#'
#' @keywords internal
#' @export
#' @import parallel
#' @import BSgenome
#' @import Biostrings
MXE.exon.seq.extraction <- function(events.info, bs.genome, core) {
  cluster <- makeCluster(core)
  # 1. exon1----
  exon1.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon1[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon4[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon1.region <- do.call(what = rbind, exon1.region)
  exon1.region <- exon1.region[exon1.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon1.region$chr,
    start = exon1.region$start,
    end = exon1.region$end,
    strand = exon1.region$strand
  )
  exon1.region$seq <- as.character(seq)
  # 2. exon2----
  exon2.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon2[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon3[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon2.region <- do.call(what = rbind, exon2.region)
  exon2.region <- exon2.region[exon2.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon2.region$chr,
    start = exon2.region$start,
    end = exon2.region$end,
    strand = exon2.region$strand
  )
  exon2.region$seq <- as.character(seq)
  # 3. exon3----
  exon3.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon3[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon2[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon3.region <- do.call(what = rbind, exon3.region)
  exon3.region <- exon3.region[exon3.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon3.region$chr,
    start = exon3.region$start,
    end = exon3.region$end,
    strand = exon3.region$strand
  )
  exon3.region$seq <- as.character(seq)
  # 4. exon4----
  exon4.region <- parLapply(cl = cluster, X = as.list(seq(1, nrow(events.info))), function(index) {
    if (events.info$strand[index] == "+") {
      jc <- events.info$exon4[index]

      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "+"
    } else {
      jc <- events.info$exon1[index]
      jc <- unlist(strsplit(x = jc, split = ":|-"))
      chr <- jc[2]
      start <- as.numeric(jc[3])
      end <- as.numeric(jc[4])
      strand <- "-"
    }
    return(data.frame(event = events.info$event[index], chr = chr, start = start, end = end, strand = strand))
  })
  exon4.region <- do.call(what = rbind, exon4.region)
  exon4.region <- exon4.region[exon4.region$chr %in% names(bs.genome), ]
  seq <- getSeq(
    x = bs.genome, names = exon4.region$chr,
    start = exon4.region$start,
    end = exon4.region$end,
    strand = exon4.region$strand
  )
  exon4.region$seq <- as.character(seq)
  stopCluster(cluster)
  return(list(exon1 = exon1.region, exon2 = exon2.region, exon3 = exon3.region, exon4 = exon4.region))
}
#' @title get length of exons and introns associated with MXE event
#' @param event a event id
#' @return data.frame
#'
#' @keywords internal
#' @export
MXE.length.feature <- function(event) {
  exon1 <- unlist(strsplit(x = event$exon1, split = ":|-"))
  exon2 <- unlist(strsplit(x = event$exon2, split = ":|-"))
  exon3 <- unlist(strsplit(x = event$exon3, split = ":|-"))
  exon4 <- unlist(strsplit(x = event$exon4, split = ":|-"))
  jc1 <- unlist(strsplit(x = event$junction1, split = ":|-"))
  jc2 <- unlist(strsplit(x = event$junction2, split = ":|-"))
  jc3 <- unlist(strsplit(x = event$junction3, split = ":|-"))
  jc4 <- unlist(strsplit(x = event$junction4, split = ":|-"))

  exon1.length <- as.numeric(exon1[4]) - as.numeric(exon1[3]) + 1
  exon2.length <- as.numeric(exon2[4]) - as.numeric(exon2[3]) + 1
  exon3.length <- as.numeric(exon3[4]) - as.numeric(exon3[3]) + 1
  i1.length <- as.numeric(jc3[4]) - as.numeric(jc3[3]) + 1
  i2.length <- as.numeric(jc1[4]) - as.numeric(jc4[3]) + 1
  i3.length <- as.numeric(jc2[4]) - as.numeric(jc2[3]) + 1
  e2.i1.ratio <- exon2.length / i1.length
  e2.i2.ratio <- exon2.length / i2.length
  e3.i2.ratio <- exon3.length / i1.length
  e3.i3.ratio <- exon3.length / i2.length
  i1.i2.ratio <- i1.length / i2.length
  i3.i2.ratio <- i3.length / i2.length
  i1.i3.ratio <- i1.length / i3.length
  if (event$strand == "+") {
    length.feature <- data.frame(
      exon1 = exon1.length,
      exon2 = exon2.length,
      exon3 = exon3.length,
      i1 = i1.length,
      i2 = i2.length,
      i3 = i3.length,
      e2.i1 = e2.i1.ratio,
      e2.i2 = e2.i2.ratio,
      e3.i2 = e3.i2.ratio,
      e3.i3 = e3.i3.ratio,
      i1.i2 = i1.i2.ratio,
      i3.i2 = i3.i2.ratio,
      i1.i3 = i1.i3.ratio
    )
  } else {
    length.feature <- data.frame(
      exon1 = exon3.length,
      exon2 = exon2.length,
      exon3 = exon1.length,
      i1 = i3.length,
      i2 = i2.length,
      i3 = i1.length,
      e2.i1 = e3.i3.ratio,
      e2.i2 = e3.i2.ratio,
      e3.i2 = e2.i2.ratio,
      e3.i3 = e2.i1.ratio,
      i1.i2 = i3.i2.ratio,
      i3.i2 = i1.i2.ratio,
      i1.i3 = 1 / i1.i3.ratio
    )
  }
  length.feature$event <- event$event
  return(length.feature)
}

#' @title get kmer of exons and splice region associated with MXE event
#' @param splice.region data.frame with sequence of splice region, \code{MXE.splice.seq.extraction}
#' @param exon.region data.frame with sequence of exon region, \code{MXE.exon.seq.extraction}
#' @return vector of sequence kmer
#'
#' @keywords internal
#' @export
#'
#' @import Biostrings

MXE.kmer.extraction <- function(splice.region, exon.region) {
  kmer <- list()
  getKmer <- function(region, name, k) {
    seq <- DNAStringSet(region$seq)
    features <- matrix(NA, nrow = nrow(region), ncol = 0)
    for (i in seq(1, k))
    {
      tmp <- oligonucleotideFrequency(x = seq, width = i, as.prob = T)
      features <- cbind(features, tmp)
    }
    rownames(features) <- region$event
    features <- list(features)
    names(features) <- name
    return(features)
  }
  # I1.5----
  features <- getKmer(splice.region[["I1.5"]], "I1.5", 3)
  kmer <- c(kmer, features)
  # I2.3----
  features <- getKmer(splice.region[["I2.3"]], "I2.3", 3)
  kmer <- c(kmer, features)
  # I1.3----
  features <- getKmer(splice.region[["I1.3"]], "I1.3", 3)
  kmer <- c(kmer, features)
  # I2.5----
  features <- getKmer(splice.region[["I2.5"]], "I2.5", 3)
  kmer <- c(kmer, features)
  # I3.3----
  features <- getKmer(splice.region[["I3.3"]], "I3.3", 3)
  kmer <- c(kmer, features)
  # I3.5----
  features <- getKmer(splice.region[["I3.5"]], "I3.5", 3)
  kmer <- c(kmer, features)
  # exon1----
  features <- getKmer(exon.region[["exon1"]], "exon1", 3)
  kmer <- c(kmer, features)
  # exon3----
  features <- getKmer(exon.region[["exon3"]], "exon3", 3)
  kmer <- c(kmer, features)
  # exon2----
  features <- getKmer(exon.region[["exon2"]], "exon2", 4)
  kmer <- c(kmer, features)
  # exon4----
  features <- getKmer(exon.region[["exon4"]], "exon4", 4)
  kmer <- c(kmer, features)
  return(kmer)
}

#' @title extract conservation scores of genome region
#' @param regionlist list of genome region
#' @param bwf \code{BigWigFile(phast.path)},
#' phast.path can be downloaded from ucsc database,
#' e.g. https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/
#' @param chr.prefix prefix added to chromosome id, "chr" or ""
#' @param core the number of threads
#'
#' @return list
#'
#' @keywords internal
#' @export
#'
#' @import rtracklayer
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
phastScore.extraction <- function(regionlist, bwf, chr.prefix, core) {
  regionlist <- lapply(X = regionlist, function(region) {
    range <- GRanges(
      seqnames = paste0(chr.prefix, region$chr),
      ranges = IRanges(start = region$start, end = region$end),
      strand = region$strand
    )
    range <- summary(bwf, range, type = "mean", )
    region$phastscore <- score(unlist(range))
    region$phastscore[is.na(region$phastscore)] <- 0
    return(region)
  })
  return(regionlist)
}

#' @title Calculate alternative branch point selection
#' Specifically, for an alternative region in an event, we divide the region
#' into \bold{seperate} bins with equal length, and calculate the cumulative adenine ratio
#' for each bin as the adenine ratio features.
#' @param region data.frame containing column event id and sequence
#' @param seperate the number of bins, default: 100
#' @param core the number of threads
#'
#' @return list
#'
#' @keywords internal
#' @export
#'
#' @import parallel
#'
aPercentage <- function(region, seperate = 100, core) {
  cluster <- makeCluster(core)
  ratio <- parLapply(cl = cluster, X = as.list(region$seq), function(seq) {
    seq <- unlist(strsplit(seq, split = ""))
    index <- ceiling(seq(0, length(seq), length.out = seperate + 1))[-1]
    ratio <- c()
    for (i in index) {
      ratio <- c(ratio, sum(seq[1:i] == "A") / i)
    }
    return(ratio)
  })
  ratio <- do.call(what = rbind, args = ratio)
  rownames(ratio) <- region$event
  stopCluster(cluster)
  return(ratio)
}
