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


# [CHANGE — OUTPUT/RECOVERY CONTRACT FIXES]
# The new implementation repairs the missing-output-directory recovery branch (the old implementation passed a
# MmapprData object where a path was expected), closes only graphics devices opened
# by this function, writes the richer candidate tables through safe serializers,
# and returns the documented MmapprData object invisibly rather than integer 1.
outputMmapprData <- function(mmapprData) {
  stopifnot(is(mmapprData, "MmapprData"))
  
  if (!dir.exists(outputFolder(param(mmapprData)))) {
    # [FIX] The old recovery branch passed the entire MmapprData object to a
    # function expecting a path, then overwrote mmapprData with the returned path.
    dir.create(outputFolder(param(mmapprData)), recursive = TRUE, showWarnings = FALSE)
  }
  current_devs = dev.list() # get open graphics devices to help cleanup
  if (length(snpDistance(mmapprData)) > 0) {
    tryCatch({
      .plotGenomeDistance(mmapprData)
      .plotPeaks(mmapprData)
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
  
  # [FIX] Return the object documented by the public API, invisibly so scripts
  # can chain output without noisy printing.
  invisible(mmapprData)
}


# [CHANGE — ATOMIC RDS CHECKPOINT HELPER]
# The old checkpoint path wrote directly to the target RDS, so interruption could leave a partial file.
# The new helper writes and verifies a same-directory temporary RDS before replacement; when direct
# overwrite-by-rename is unavailable, it preserves the old target and restores it on failure.
.atomicSaveRDS <- function(object, file, .renameFile = file.rename) {
  dir <- dirname(file)
  if (!dir.exists(dir)) stop("RDS destination directory does not exist: ", dir)
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dir)
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)

  saveRDS(object, tmp)
  info <- file.info(tmp)
  if (!file.exists(tmp) || is.na(info$size) || info$size <= 0)
    stop("Temporary RDS serialization failed: ", tmp)

  # Same-directory rename replaces atomically on normal POSIX filesystems. Some
  # platforms (notably Windows) refuse to rename over an existing destination.
  # In that case move the old readable file aside first, then restore it if the
  # completed temporary file still cannot be installed. Never delete the only
  # readable prior state before the replacement succeeds.
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


# [CHANGE — COLLISION-RESISTANT OUTPUT NAMES]
# The old implementation used a timestamp alone, so runs started in the same second could
# collide. The new implementation uses tempfile-generated unique suffixes for default and temporary
# output paths without consuming the analysis RNG stream.
.defaultOutputFolder <- function() {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # [IMPROVE] tempfile() adds a collision-resistant suffix without consuming R's
  # random-number stream, so two analyses started in the same second do not clash.
  basename(tempfile(pattern = paste0("mmappr2_", stamp, "_"), tmpdir = getwd()))
}


#' @title Generate temporary output folder
#' @name tempOutputFolder
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


# [CHANGE — NONINTERACTIVE, PROTECTED OVERWRITE]
# The old implementation prompted interactively and could recursively clear an existing
# directory after a single response. The new implementation requires explicit overwrite=TRUE, refuses
# root/top-level/home/current/temp-session directories, verifies clearing succeeds,
# creates paths recursively, normalizes the accepted path, and creates the log file
# without interactive input.
.prepareOutputFolder <- function(outputFolder, overwrite = FALSE) {
  if (!is.character(outputFolder) || length(outputFolder) != 1L || !nzchar(outputFolder))
    stop("outputFolder must be one non-empty path")

  if (dir.exists(outputFolder)) {
    existingPath <- normalizePath(outputFolder, mustWork = TRUE)

    # [FIX] Evaluate overwrite safety before checking directory contents. The
    # old implementation had no protected-path rule; the new rule applies even
    # when the target is empty, because safety is a property of the path itself.
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



# [CHANGE — DEGENERATE PLOT LIMIT GUARDS]
# The old plotting code constructed limits directly from min/max values, which can fail for constant,
# empty, or non-finite data. The new helpers return finite non-zero-width limits for those edge cases.
.safeYLim <- function(x, upperPad = 0.10) {
  # [FIX] Base plot() rejects zero-width or non-finite limits. Constant/degenerate
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

# [CHANGE — GENOME-PLOT HARDENING/METADATA]
# The old genome plot used raw min/max limits and did not display peak-cutoff/refined-region metadata.
# The new plot uses safe limits, normalized numeric LOESS coordinates, namespaced sequence ordering, and
# explicit cutoff/region overlays; device handling also avoids closing unrelated user devices.
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


# [CHANGE — REFINED-PEAK DIAGNOSTICS]
# The old peak plot shaded to a hard-coded -5 baseline, omitted the initial cutoff, and labeled the KDE
# overlay as probability while relying on raw ranges. The new plot shades from the actual plotting
# baseline, shows the stored cutoff, handles degenerate density support safely, and labels the overlay
# as density.
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
      # [FIX] Shade to the actual plot baseline rather than an arbitrary -5,
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



# Convert Bioconductor-rich table columns into values that base write.table() can
# serialize reliably. VariantAnnotation/GenomicRanges metadata can contain List,
# DNAStringSet, IRanges, or other vector-like columns; base write.table() errors on
# genuine list columns ("unimplemented type 'list'"). Atomic columns are preserved
# as-is so numeric depths/frequencies stay numeric in the TSV.
# [CHANGE — ROBUST TSV SERIALIZATION]
# Old write.table() calls can fail on VariantAnnotation/GenomicRanges list-like
# metadata. The new implementation flattens only non-atomic cells deterministically, sanitizes embedded
# tabs/newlines, preserves ordinary numeric columns, and handles zero-row list
# columns so header-only result files still serialize correctly.
.tsvSafeDataFrame <- function(x) {
  df <- as.data.frame(x)
  # Even a zero-row table can contain a list/List column whose class makes
  # write.table() fail before it writes the header. Only the truly zero-column
  # case can bypass column sanitization safely.
  if (ncol(df) == 0L) return(df)

  sanitizeText <- function(x) {
    # [FIX/ROBUSTNESS] A literal tab/newline inside an annotation value would
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
    # [FIX] Leave ordinary scalar vectors alone; only flatten columns that base
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


# [CHANGE — PRESERVE OLD MUTATION TABLE WHILE ADDING METADATA]
# The old implementation wrote a fixed DetectedMutationsFor*.tsv schema. The new implementation
# reconstructs that same schema/column order from the richer coding-effect object, using explicit mutant
# depth metadata when predictCoding() does not carry the standard depth columns. WT/delta-AF evidence
# is written separately in AllCandidateVariantsFor*.tsv rather than changing the old table.
.originalDetectedMutationTable <- function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  originalCols <- c("seqnames", "start", "end", "width", "strand",
                    "ref", "alt", "refDepth", "altDepth", "peakDensity",
                    "GENEID", "TXID", "PROTEINLOC", "CONSEQUENCE",
                    "REFAA", "VARAA")

  # Current predictCoding versions do not consistently propagate the query's
  # standard depth columns, but the current object carries equivalent mutant
  # evidence explicitly. Use it only to reconstruct the old table fields.
  if (!"refDepth" %in% names(df) && "mutRefDepth" %in% names(df))
    df$refDepth <- df$mutRefDepth
  if (!"altDepth" %in% names(df) && "mutAltDepth" %in% names(df))
    df$altDepth <- df$mutAltDepth

  # The old table emitted the first protein-location value for this
  # table. VariantAnnotation commonly stores PROTEINLOC as an IRanges/S4Vectors
  # List derivative rather than a base list, so handle both without flattening the
  # richer annotation retained in the in-memory object.
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


# [CHANGE — FULL SNV OUTPUT IN ADDITION TO OLD FILES]
# The old implementation wrote coding-effect rows and expression output only, so noncoding
# candidate SNVs disappeared from tabular output. The new implementation adds AllCandidateVariantsFor*
# with the complete SNV/WT evidence while retaining the old mutation and
# expression filenames/schema through the safe TSV writer.
.writeCandidateTables <- function(candList, outputFolder){
  seqnames <- unique(c(names(candList$snps), names(candList$effects), names(candList$diff)))
  for (seqname in seqnames) {
    # [IMPROVE] Always write the full candidate SNV set, including noncoding
    # variants and the new WT-vs-mutant allele-frequency evidence.
    if (!is.null(candList$snps[[seqname]])) {
      .writeTsv(candList$snps[[seqname]],
                file.path(outputFolder, paste0("AllCandidateVariantsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$effects[[seqname]])) {
      # Preserve the old DetectedMutations table schema and column order.
      # The richer SNV/WT metadata stays available
      # in AllCandidateVariantsFor*.tsv and in the returned MmapprData object.
      .writeTsv(.originalDetectedMutationTable(candList$effects[[seqname]]),
                file.path(outputFolder, paste0("DetectedMutationsFor", seqname, ".tsv")))
    }

    if (!is.null(candList$diff[[seqname]])) {
      .writeTsv(candList$diff[[seqname]],
                file.path(outputFolder, paste0("DifferentiallyExpressedGenesFor", seqname, ".tsv")))
    }
  }
}

