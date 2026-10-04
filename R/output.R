#' @title Generate plots and tables from MMAPPR2 data
#' 
#' @name outputMmapprData
#' 
#' @param mmapprData The \linkS4class{MmapprData} object to be output
#'
#' @return A \linkS4class{MmapprData} object after writing output files
#'   to the folder specified in the \code{outputFolder} slot of the
#'   \code{\link{MmapprParam}} used.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
#' 
#' \dontrun{
#' md <- mmapprData(mmappr_param)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#'
#' outputMmapprData(postCandidatesMD)
#' }
#' 
NULL


outputMmapprData <- function(mmapprData) {
  stopifnot(is(mmapprData, "MmapprData"))
  
  if (!dir.exists(outputFolder(param(mmapprData)))) {
    dir.create(outputFolder(param(mmapprData)), recursive = TRUE, showWarnings = FALSE)
  }
  current_devs = dev.list() # get open graphics devices to help cleanup
  if (length(snpDistance(mmapprData)) > 0) {
    tryCatch({
      .plotGenomeDistance(mmapprData)
      .plotPeaks(mmapprData)
      if (isTRUE(exportAiccPlots(param(mmapprData)))) .plotAicc(mmapprData)
    }, finally = {
      now <- dev.list()
      opened_devs <- if (is.null(now)) integer() else now[!(now %in% current_devs)]
      if (length(opened_devs) > 0L) {
        for (d in rev(opened_devs)) try(grDevices::dev.off(which = d), silent = TRUE)
      }
    })
  }
  
  if (length(candidates(mmapprData)) > 0) {
    .writeCandidateTables(mmapprData@candidates,
                          outputFolder(param(mmapprData)))
  }
  
  # Return the object documented by the public API, invisibly so scripts
  # can chain output without noisy printing.
  invisible(mmapprData)
}


# Export the AICc evaluations retained by loessFit(). One page is written per
# successfully fitted chromosome, with all finite evaluated spans shown and the
# selected span highlighted. This is output-only: it never reruns a LOESS fit or
# changes the optimized result stored in the MmapprData object.
.plotAicc <- function(mmapprData) {
  usable <- names(mmapprData@snpDistance)[vapply(
    mmapprData@snpDistance,
    function(x) is.list(x) && is.data.frame(x$aicc) &&
      all(c("spans", "aiccValues") %in% names(x$aicc)) &&
      any(is.finite(x$aicc$aiccValues)),
    logical(1)
  )]
  if (!length(usable)) {
    warning("AICc plot export was requested, but no finite AICc search results are available")
    return(invisible(NULL))
  }

  grDevices::pdf(file.path(outputFolder(param(mmapprData)), "aicc_plots.pdf"),
                 width = 8.5, height = 6.5)
  deviceOpen <- TRUE
  on.exit(if (deviceOpen) try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  for (seqname in usable) {
    result <- mmapprData@snpDistance[[seqname]]
    tab <- result$aicc
    ok <- is.finite(tab$spans) & is.finite(tab$aiccValues)
    tab <- tab[ok, , drop = FALSE]
    tab <- tab[order(tab$spans), , drop = FALSE]
    if (!nrow(tab)) next

    graphics::plot(tab$spans, tab$aiccValues, type = "b", pch = 16,
                   xlab = "LOESS span", ylab = "AICc",
                   main = paste(seqname, "LOESS span optimization"),
                   xlim = .safeXLim(tab$spans), ylim = .safeYLim(tab$aiccValues))

    best <- result$bestSpan
    if (length(best) == 1L && is.finite(best)) {
      nearest <- which.min(abs(tab$spans - best))
      graphics::abline(v = best, lty = 2)
      graphics::points(tab$spans[nearest], tab$aiccValues[nearest],
                       pch = 19, cex = 1.4)
      graphics::legend("topright",
                       legend = sprintf("selected span = %s", format(best)),
                       lty = 2, pch = 19, bty = "n")
    }
  }

  grDevices::dev.off()
  deviceOpen <- FALSE
  invisible(NULL)
}


.atomicSaveRDS <- function(object, file, .renameFile = file.rename) {
  dir <- dirname(file)
  if (!dir.exists(dir)) stop("RDS destination directory does not exist: ", dir)
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dir)
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)

  saveRDS(object, tmp)
  info <- file.info(tmp)
  if (!file.exists(tmp) || is.na(info$size) || info$size <= 0)
    stop("Temporary RDS serialization failed: ", tmp)

  if (!.renameFile(tmp, file)) {
    if (!file.exists(file))
      stop("Could not move completed RDS into place: ", file)

    backup <- tempfile(pattern = paste0(".", basename(file), ".previous."), tmpdir = dir)
    if (file.exists(backup)) unlink(backup, force = TRUE)
    if (!.renameFile(file, backup))
      stop("Could not move existing RDS aside for replacement: ", file)

    installed <- .renameFile(tmp, file)
    if (!installed) {
      restored <- .renameFile(backup, file)
      if (!restored) {
        stop("Could not install completed RDS and could not restore the previous file. ",
             "Previous state remains at: ", backup)
      }
      stop("Could not install completed RDS; previous file was restored: ", file)
    }

    if (file.exists(backup)) {
      rc <- unlink(backup, force = TRUE)
      if (!identical(as.integer(rc), 0L) || file.exists(backup))
        warning("Replacement succeeded but previous RDS backup could not be removed: ", backup)
    }
  }
  invisible(file)
}


.defaultOutputFolder <- function() {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # tempfile() adds a collision-resistant suffix without consuming R's
  # random-number stream, so two analyses started in the same second do not clash.
  basename(tempfile(pattern = paste0("mmappr2_", stamp, "_"), tmpdir = getwd()))
}


#' @title Generate temporary output folder
#' @name tempOutputFolder
#' @description
#' Returns a unique temporary output-directory path. The directory is created
#' by \code{mmapprParam()} when the parameter object is constructed.
#'
#' @return The path to the temporary directory
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'     mmappr_param <- mmapprParam(refFasta = MMAPPR2data::goldenFasta(),
#'                                 wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#' }
NULL


tempOutputFolder <- function() {
  # Return a unique path; mmapprParam() creates it. This avoids same-second name
  # collisions in tests and parallel workflows.
  tempfile(pattern = paste0("mmappr2_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_"),
           tmpdir = tempdir())
}


.prepareOutputFolder <- function(outputFolder, overwrite = FALSE) {
  if (!is.character(outputFolder) || length(outputFolder) != 1L || !nzchar(outputFolder))
    stop("outputFolder must be one non-empty path")

  if (dir.exists(outputFolder)) {
    existingPath <- normalizePath(outputFolder, mustWork = TRUE)

    if (isTRUE(overwrite)) {
      protected <- unique(normalizePath(c(path.expand("~"), getwd(), tempdir()),
                                        mustWork = TRUE))
      isRoot <- identical(dirname(existingPath), existingPath)
      parent <- dirname(existingPath)
      isTopLevel <- !isRoot && identical(dirname(parent), parent)
      if (isRoot || isTopLevel || existingPath %in% protected)
        stop("Refusing to overwrite protected directory: ", existingPath)
    }

    contents <- list.files(existingPath, all.files = TRUE, no.. = TRUE)
    if (length(contents) > 0L) {
      if (!isTRUE(overwrite)) {
        stop("Output folder already exists and is not empty: ", outputFolder,
             ". Choose a new folder or set overwrite=TRUE explicitly.")
      }
      unlink(file.path(existingPath, contents), recursive = TRUE, force = TRUE)
      stillThere <- list.files(existingPath, all.files = TRUE, no.. = TRUE)
      if (length(stillThere) > 0L)
        stop("Could not fully clear output folder: ", existingPath)
    }
  } else {
    ok <- dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(outputFolder)) stop("Could not create output folder: ", outputFolder)
  }

  outputFolder <- normalizePath(outputFolder, mustWork = TRUE)
  logPath <- file.path(outputFolder, "mmappr2.log")
  if (!file.exists(logPath) && !file.create(logPath))
    stop("Could not create log file in output folder: ", outputFolder)
  outputFolder
}



.safeYLim <- function(x, upperPad = 0.10) {
  # Base plot() rejects zero-width or non-finite limits. Constant/degenerate
  # fits can occur on small datasets, so construct a finite visible range.
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  lo <- min(x); hi <- max(x)
  if (lo == hi) {
    pad <- if (lo == 0) 1 else max(abs(lo) * upperPad, .Machine$double.eps^0.5)
    return(c(lo - pad, hi + pad))
  }
  c(lo, hi + (hi - lo) * upperPad)
}

.safeXLim <- function(x, padFraction = 0.01) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  lo <- min(x); hi <- max(x)
  if (lo == hi) {
    pad <- if (lo == 0) 1 else max(abs(lo) * padFraction, 1)
    return(c(lo - pad, hi + pad))
  }
  c(lo, hi)
}

# Genome-level plotting uses safe limits, normalized numeric LOESS coordinates, and
# namespaced sequence ordering. Peak-cutoff/refined-region diagnostics are drawn by
# .plotPeaks(); device handling here avoids closing unrelated user devices.
.plotGenomeDistance <- function(mmapprData, savePdf = TRUE) {
  #generate one big dataframe for plots, along with break and label points
  tailPos <- 0
  plotDf <- NULL
  breaks <- tailPos
  labelpos <- NULL
  for (i in GenomeInfoDb::orderSeqlevels(names(mmapprData@snpDistance))) {
    if (!is(mmapprData@snpDistance[[i]], 'list'))
      next
    else if (!("loess" %in% names(mmapprData@snpDistance[[i]])))
      stop(sprintf("Distance list for sequence %s missing loess fit data",
                   names(mmapprData@snpDistance)[i]))
    
    chrLoess <- mmapprData@snpDistance[[i]]$loess
    chrX <- as.numeric(chrLoess$x)
    chrDf <- data.frame(pos = chrX + tailPos,
                        seqname = names(mmapprData@snpDistance)[i],
                        fitted = as.numeric(chrLoess$fitted),
                        unfitted = as.numeric(chrLoess$y))
    plotDf <- rbind(plotDf, chrDf)
    tailPos <- chrDf$pos[nrow(chrDf)]
    breaks <- c(breaks, tailPos)
    labelpos <-
      c(labelpos, (breaks[length(breaks)]+breaks[length(breaks)-1])/2)
  }

  if (is.null(plotDf) || nrow(plotDf) == 0L)
    stop("No successfully fitted chromosome data are available for genome plotting")

  openedPdf <- FALSE
  if (savePdf) {
    pdf(file.path(mmapprData@param@outputFolder, "genome_plots.pdf"),
        width = 11, height = 8.5)
    openedPdf <- TRUE
    on.exit(if (openedPdf) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  }
  par(mfrow = c(2,1))
  
  plot(x = plotDf$pos, y = plotDf$fitted, type = 'l',
       ylim = .safeYLim(plotDf$fitted),
       ylab = NA,
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome")
  abline(v = breaks[seq_len(length(breaks) - 1L)], col = "grey")
  mtext(unique(sub("^chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p~ ~"(Loess fit)",
                   list(p = mmapprData@param@distancePower)), 
        side = 2, line = 2)
  
  plot(x = plotDf$pos, y = plotDf$unfitted, 
       pch = 16, cex = .8, col = "#999999AA",
       ylim = .safeYLim(plotDf$unfitted),
       ylab = NA,
       xaxt = 'n', xaxs = 'i', xlab = "Chromosome" )
  abline(v = (breaks), col = "grey")
  mtext(unique(sub("^chr", "", plotDf$seqname[!is.na(plotDf$seqname)])),
        at = labelpos, side = 1, cex = .6)
  mtext(substitute("ED"^p,
                   list(p = mmapprData@param@distancePower)),
        side = 2, line = 2)
  
  if (openedPdf) {
    grDevices::dev.off()
    openedPdf <- FALSE
  }
}


.plotPeaks <- function(mmapprData) {
  if (length(mmapprData@peaks) == 0L) return(invisible(NULL))
  pdf(file.path(mmapprData@param@outputFolder, "peak_plots.pdf"),
      width = 11, height = 8.5)
  deviceOpen <- TRUE
  on.exit(if (deviceOpen) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  par(mfrow = c(2,1), mar = c(5, 4, 4, 4))

  for (seqname in names(mmapprData@peaks)) {
    chrLoess <- mmapprData@snpDistance[[seqname]]$loess
    chrX <- as.numeric(chrLoess$x)
    fitted <- as.numeric(chrLoess$fitted)
    rawY <- as.numeric(chrLoess$y)
    start <- mmapprData@peaks[[seqname]]$start
    end <- mmapprData@peaks[[seqname]]$end
    densityData <- mmapprData@peaks[[seqname]]$densityData
    yLim <- .safeYLim(fitted)
    xLim <- .safeXLim(chrX) * 1E-6

    plot(chrX/1000000, fitted, type = 'l',
         ylim = yLim, xlim = xLim,
         xlab = paste(seqname, "Base Position (MB)"),
         ylab = NA, xaxs = 'i')
    mtext(substitute("ED"^p~ ~"(Loess fit)",
                     list(p = mmapprData@param@distancePower)),
          side = 2, line = 2)

    inside <- chrX >= start & chrX <= end & is.finite(fitted)
    if (any(inside)) {
      # Shade to the actual plot baseline rather than an arbitrary -5,
      # which could distort clipping when the plotted range is small/positive.
      shadeX <- c(start, chrX[inside], end)
      shadeY <- c(yLim[1], fitted[inside], yLim[1])
      polygon(shadeX/1000000, shadeY, col = '#2ecc71', border = NA)
    }
    cutoff <- mmapprData@peaks[[seqname]]$cutoff
    if (length(cutoff) == 1L && is.finite(cutoff))
      abline(h = cutoff, lty = 2, col = "grey")

    # Overlay the resampling density only when finite support remains inside the
    # plotted chromosome range. Degenerate one-point densities use a point rather
    # than an invisible zero-length line.
    if (length(densityData) > 0L && length(densityData$x) > 0L) {
      dx <- as.numeric(densityData$x)
      dy <- as.numeric(densityData$y)
      posMatch <- is.finite(dx) & is.finite(dy) & dx >= min(chrX, na.rm = TRUE) &
                  dx <= max(chrX, na.rm = TRUE)
      if (any(posMatch)) {
        dx <- dx[posMatch]; dy <- dy[posMatch]
        par(new = TRUE)
        plot(dx, dy, type = if (length(dx) > 1L) 'l' else 'p',
             ylim = .safeYLim(dy), xlim = .safeXLim(chrX),
             ann = FALSE, xaxs = 'i', xaxt = 'n', yaxt = 'n', col = '#502ecc')
        axis(side = 4, col = '#502ecc')
        mtext(side = 4, line = 2, 'Density', col = '#502ecc')
        legend('topright',
               legend = c("Fitted Distance Curve", "Peak Resampling Distribution"),
               col = c("black", "#502ecc"), lty = c(1, 1), cex = 0.8, lwd = 3)
      }
    }

    if (length(chrX) > 0L) {
      plot(chrX/1000000, rawY, pch = 16, cex = .6,
           ylim = .safeYLim(rawY), xlim = xLim,
           ylab = NA, xlab = paste(seqname, "Base Position (MB)"), xaxs = 'i')
      mtext(substitute("ED"^p, list(p = mmapprData@param@distancePower)),
            side = 2, line = 2)
    }
  }
  grDevices::dev.off()
  deviceOpen <- FALSE
  invisible(NULL)
}



.tsvSafeDataFrame <- function(x) {
  df <- as.data.frame(x)
  # Even a zero-row table can contain a list/List column whose class makes
  # write.table() fail before it writes the header. Only the truly zero-column
  # case can bypass column sanitization safely.
  if (ncol(df) == 0L) return(df)

  sanitizeText <- function(x) {
    # A literal tab/newline inside an annotation value would
    # corrupt an unquoted TSV row. GTF-derived labels should not contain these,
    # but sanitize defensively so output remains one record per line.
    out <- x
    ok <- !is.na(out)
    out[ok] <- gsub("[\t\r\n]+", " ", out[ok])
    out
  }

  collapseCell <- function(x) {
    if (length(x) == 0L) return("")
    vals <- tryCatch(as.character(x), error = function(e) character())
    if (length(vals) == 0L) return("")
    if (all(is.na(vals))) return(NA_character_)
    paste(sanitizeText(vals[!is.na(vals)]), collapse = ",")
  }

  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.character(col)) {
      df[[nm]] <- sanitizeText(col)
      next
    }
    # Leave ordinary scalar vectors alone; only flatten columns that base
    # write.table() cannot safely serialize. The [[i]] path is important for
    # CharacterList/IntegerList/list columns because one table row can contain
    # multiple annotation values.
    listLike <- is.list(col) || methods::is(col, "List") || !is.atomic(col)
    if (!listLike) next
    df[[nm]] <- vapply(seq_len(nrow(df)), function(i) {
      cell <- tryCatch(col[[i]], error = function(e) col[i])
      collapseCell(cell)
    }, character(1))
  }
  df
}

.writeTsv <- function(x, file) {
  utils::write.table(.tsvSafeDataFrame(x), file = file, sep = "\t",
                     quote = FALSE, row.names = FALSE, na = "NA")
  invisible(file)
}


.originalDetectedMutationTable <- function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  originalCols <- c("seqnames", "start", "end", "width", "strand",
                    "ref", "alt", "refDepth", "altDepth", "peakDensity",
                    "GENEID", "TXID", "PROTEINLOC", "CONSEQUENCE",
                    "REFAA", "VARAA")

  if (!"refDepth" %in% names(df) && "mutRefDepth" %in% names(df))
    df$refDepth <- df$mutRefDepth
  if (!"altDepth" %in% names(df) && "mutAltDepth" %in% names(df))
    df$altDepth <- df$mutAltDepth

  if ("PROTEINLOC" %in% names(df)) {
    proteinLoc <- df$PROTEINLOC
    listLike <- is.list(proteinLoc) || methods::is(proteinLoc, "List") || !is.atomic(proteinLoc)
    if (listLike) {
      df$PROTEINLOC <- vapply(seq_len(nrow(df)), function(i) {
        value <- tryCatch(proteinLoc[[i]], error = function(e) proteinLoc[i])
        value <- tryCatch(as.character(value), error = function(e) character())
        if (length(value)) value[[1L]] else NA_character_
      }, character(1))
    }
  }

  # Assign a zero-length vector when x has zero rows. Assigning scalar NA to a
  # zero-row data.frame can otherwise create a replacement-length error instead
  # of producing the desired header-only old-format table.
  for (nm in setdiff(originalCols, names(df))) df[[nm]] <- rep(NA, nrow(df))
  df[, originalCols, drop = FALSE]
}


.writeCandidateTables <- function(candList, outputFolder){
  seqnames <- unique(c(names(candList$snps), names(candList$effects), names(candList$diff)))
  for (seqname in seqnames) {
    # Always write the full candidate SNV set, including noncoding
    # variants and the new WT-vs-mutant allele-frequency evidence.
    if (!is.null(candList$snps[[seqname]])) {
      .writeTsv(candList$snps[[seqname]],
                file.path(outputFolder, paste0("AllCandidateVariantsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$effects[[seqname]])) {
      .writeTsv(.originalDetectedMutationTable(candList$effects[[seqname]]),
                file.path(outputFolder, paste0("DetectedMutationsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$diff[[seqname]])) {
      .writeTsv(candList$diff[[seqname]],
                file.path(outputFolder, paste0("DifferentiallyExpressedGenesFor", seqname, ".tsv")))
    }
  }
}

