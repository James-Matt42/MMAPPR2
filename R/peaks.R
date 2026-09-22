#' 
#' @title Characterize Euclidean distance peaks using resampling simulation
#'
#' @name peakRefinement
#'
#' Follows the \code{\link{prePeak}} step and precedes
#' \code{\link{generateCandidates}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{peaks}
#'   slot filled and populated.
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
#'
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' }
#' 
#' @import BiocParallel
#' 
NULL


# [CHANGE — PARAMETERIZED, REPRODUCIBLE PEAK REFINEMENT]
# The old implementation hard-coded the peak cutoff multiplier and 1,000 resamples, used one interval
# rule, and depended on the ambient RNG state. The new implementation records those choices in
# MmapprParam, supports an explicit interval method, and uses deterministic chromosome-specific seeds.
peakRefinement <- function(mmapprData){
    mmapprData@peaks <- BiocParallel::bplapply(
        mmapprData@peaks,
        .peakRefinementChr,
        mmapprData = mmapprData
    )
    mmapprData
}


# [CHANGE — FLAT-APEX POSITION FIX]
# When multiple markers share the maximum LOESS value, the old implementation effectively
# selected the first/left-most maximum. The new implementation reports the median genomic position of
# the plateau, removing that arbitrary positional bias.
.peakApexPosition <- function(x, fitted) {
    valid <- is.finite(x) & is.finite(fitted)
    if (!any(valid)) return(NA_real_)
    x <- as.numeric(x[valid])
    fitted <- as.numeric(fitted[valid])
    maxFit <- max(fitted)
    # [FIX] Flat-topped fits should not be biased toward the first/left-most point.
    # The median coordinate is stable and symmetric across a plateau.
    stats::median(x[fitted == maxFit])
}


# [CHANGE — RESAMPLE FAILURE TOLERANCE]
# In the old implementation, an error in any half-marker LOESS resample could abort the
# entire refinement. The new implementation converts a failed resample to NA, excludes it from
# the density estimate, and records the resulting resample success rate.
.getSubsampleLoessMax <- function(rawData, loessSpan) {
    n <- nrow(rawData)
    if (n < 5L) return(NA_real_)
    sampleSize <- max(5L, floor(n * 0.5))
    idx <- sample.int(n, size = sampleSize, replace = FALSE)
    tempData <- rawData[idx, , drop = FALSE]
    tempData <- tempData[order(tempData$pos), , drop = FALSE]

    fit <- suppressWarnings(try(
        stats::loess(euclideanDistance ~ pos,
                     data = tempData,
                     span = loessSpan,
                     degree = 1,
                     family = "symmetric"),
        silent = TRUE
    ))
    # [FIX] One bad half-sample should not abort all 1,000 resamples.
    if (inherits(fit, "try-error") || !inherits(fit, "loess")) return(NA_real_)
    .peakApexPosition(as.numeric(fit$x), fit$fitted)
}


# [CHANGE — KDE INTEGRATED AS PROBABILITY MASS]
# The old implementation's interval selection treated sampled KDE heights as if they were
# directly comparable probability mass. The new implementation multiplies density by local grid-cell
# widths before selecting the requested mass, which is correct on nonuniform grids.
.densityPointMass <- function(x, y) {
    stopifnot(length(x) == length(y), length(x) >= 1L)
    if (length(x) == 1L) return(1)
    dx <- diff(x)
    if (any(!is.finite(dx)) || any(dx <= 0)) stop("Density x coordinates must be strictly increasing")
    # Approximate the integration cell represented by each KDE grid point. This
    # makes the mass calculation correct even if a future density grid is not
    # perfectly equally spaced.
    widths <- c(dx[1] / 2, (head(dx, -1L) + tail(dx, -1L)) / 2, tail(dx, 1L) / 2)
    rawMass <- pmax(y, 0) * widths
    if (!is.finite(sum(rawMass)) || sum(rawMass) <= 0) stop("Peak-density mass is empty")
    rawMass / sum(rawMass)
}


# [CHANGE — EXPLICIT REFINED-INTERVAL METHODS]
# The old implementation had one density-threshold rule and returned the span from the leftmost
# to rightmost selected position. The new implementation names that envelope-style choice `hpd_span`
# (now using integrated KDE mass) and adds `shortest_contiguous` as an opt-in shortest single
# interval; both target peakIntervalWidth probability mass.
.getPeakFromTopP <- function(data, topP, method = c("hpd_span", "shortest_contiguous")) {
    stopifnot(ncol(data) == 2L, topP > 0, topP <= 1)
    method <- match.arg(method)
    names(data) <- c("x", "y")
    data <- data[is.finite(data$x) & is.finite(data$y) & data$y >= 0, , drop = FALSE]
    data <- data[order(data$x), , drop = FALSE]
    if (nrow(data) == 0L || sum(data$y) <= 0) stop("Peak-density data are empty")

    mass <- .densityPointMass(data$x, data$y)
    peakPos <- .peakApexPosition(data$x, data$y)

    if (method == "hpd_span") {
        # [FIX/IMPROVE] Preserve the old interval rule without constructing one row
        # per genomic base. Rank the
        # compact KDE grid by density, accumulate probability mass, then span all
        # grid points at or above the resulting density threshold. If the density
        # is multimodal this intentionally keeps both plausible modes inside the
        # single interval supported by the current data model.
        ranked <- order(data$y, decreasing = TRUE)
        cumulative <- cumsum(mass[ranked])
        firstEnough <- which(cumulative >= topP)[1]
        cutoffDensity <- data$y[ranked[firstEnough]]
        selected <- data$y >= cutoffDensity
        return(list(minPos = min(data$x[selected]),
                    maxPos = max(data$x[selected]),
                    peakPos = peakPos))
    }

    # Optional alternative: the shortest ONE-PIECE interval containing the
    # requested mass. This is narrower, but can discard a secondary plausible
    # mode, so it is deliberately not the default that preserves the old interval behavior.
    cs <- c(0, cumsum(mass))
    bestLeft <- 1L
    bestRight <- nrow(data)
    bestWidth <- Inf
    for (left in seq_len(nrow(data))) {
        target <- cs[left] + topP
        right <- which(cs[-1L] >= target)[1]
        if (is.na(right) || right < left) next
        width <- data$x[right] - data$x[left]
        if (width < bestWidth) {
            bestWidth <- width
            bestLeft <- left
            bestRight <- right
        }
    }
    list(minPos = data$x[bestLeft],
         maxPos = data$x[bestRight],
         peakPos = peakPos)
}


# [CHANGE — REPRODUCIBLE CHROMOSOME RNG]
# The old implementation sampled from the ambient RNG stream, so results could depend on prior RNG use
# and chromosome processing order. The new implementation derives a stable chromosome-specific seed
# from the recorded base seed.
.derivedSeed <- function(baseSeed, seqname) {
    # Deterministic per chromosome and stable across serial/parallel execution.
    chars <- utf8ToInt(as.character(seqname))
    offset <- if (length(chars)) sum(chars * seq_along(chars)) else 0
    as.integer((as.double(baseSeed) + offset) %% (.Machine$integer.max - 1L) + 1L)
}


# [CHANGE — PEAK-REFINEMENT HARDENING/METADATA]
# The old implementation returned a minimal interval/density record, used the caller RNG stream, and
# assumed resampling/density coordinates were directly usable. The new implementation normalizes LOESS
# predictors, restores caller RNG state, tolerates failed resamples, clips density support to chromosome
# bounds, returns valid integer coordinates, and records LOESS/density apices, success rate, seed, and
# interval-construction metadata.
.peakRefinementChr <- function(inputList, mmapprData) {
    stopifnot("seqname" %in% names(inputList))
    seqname <- inputList$seqname

    loessObj <- mmapprData@snpDistance[[seqname]]$loess
    loessSpan <- loessObj$pars$span
    # [IMPROVE] loess stores predictors in a matrix-like `x` component on some R
    # versions. Normalize the single genomic predictor to a plain numeric vector so
    # downstream data.frame/subsetting behavior is unambiguous.
    rawData <- data.frame(pos = as.numeric(loessObj$x),
                          euclideanDistance = as.numeric(loessObj$y))
    rawData <- rawData[is.finite(rawData$pos) & is.finite(rawData$euclideanDistance), , drop = FALSE]
    if (nrow(rawData) < 10L) stop("Too few finite markers to refine peak on ", seqname)

    # [FIX] Peak refinement used to depend on the ambient RNG state, so identical
    # inputs could yield slightly different intervals. Use a deterministic seed
    # derived from the user-recorded base seed and chromosome name.
    seed <- .derivedSeed(randomSeed(mmapprData@param), seqname)
    oldSeedExists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (oldSeedExists) oldSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
        if (oldSeedExists) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)

    nIter <- as.integer(peakResampleIterations(mmapprData@param))
    maxValues <- replicate(nIter,
                           .getSubsampleLoessMax(rawData = rawData,
                                                loessSpan = loessSpan))
    success <- is.finite(maxValues)
    successRate <- mean(success)
    maxValues <- maxValues[success]
    if (length(maxValues) < max(5L, ceiling(0.5 * nIter)))
        stop("Too many failed peak-resampling LOESS fits on ", seqname,
             " (successful: ", length(maxValues), "/", nIter, ")")

    # A degenerate distribution can occur when every resample chooses exactly the
    # same marker. In that case density() cannot estimate bandwidth sensibly.
    if (length(unique(maxValues)) == 1L) {
        densityData <- list(x = maxValues[[1]], y = 1)
        peak <- list(minPos = maxValues[[1]], maxPos = maxValues[[1]],
                     peakPos = maxValues[[1]])
        densityFunction <- function(x) as.numeric(x == maxValues[[1]])
    } else {
        densityData <- stats::density(maxValues)
        densityDf <- data.frame(x = densityData$x, y = densityData$y)

        # [FIX] KDE support can extend beyond the chromosome. Clip it before using
        # the density to define a physical genomic interval.
        fa_si <- .faSeqinfo(mmapprData@param@refGenome)
        chrLen <- GenomeInfoDb::seqlengths(fa_si)[seqname]
        if (length(chrLen) == 1L && is.finite(chrLen)) {
            densityDf <- densityDf[densityDf$x >= 1 & densityDf$x <= chrLen, , drop = FALSE]
        }
        if (nrow(densityDf) == 0L) stop("No valid KDE support remained on ", seqname)

        peak <- .getPeakFromTopP(densityDf, peakIntervalWidth(mmapprData@param),
                                 method = peakIntervalMethod(mmapprData@param))
        densityData$x <- densityDf$x
        densityData$y <- densityDf$y
        densityFunction <- stats::approxfun(x = densityDf$x, y = densityDf$y,
                                            yleft = 0, yright = 0, rule = 1)
    }

    # Keep integer, valid coordinates for downstream GRanges queries.
    fa_si <- .faSeqinfo(mmapprData@param@refGenome)
    chrLen <- GenomeInfoDb::seqlengths(fa_si)[seqname]
    lower <- max(1L, floor(peak$minPos))
    upper <- ceiling(peak$maxPos)
    if (length(chrLen) == 1L && is.finite(chrLen)) upper <- min(upper, as.integer(chrLen))

    loessX <- as.numeric(loessObj$x)
    loessPeak <- .peakApexPosition(loessX, loessObj$fitted)

    # [LIMIT] This data model still represents one refined locus per chromosome.
    # Proper polygenic / multiple-same-chromosome support requires changing the
    # peaks structure and downstream output, not simply choosing a second maximum.
    list(seqname = seqname,
         start = lower,
         end = upper,
         densityFunction = densityFunction,
         peakPosition = peak$peakPos,       # resampling-density apex; preserves the old field meaning
         densityPeakPosition = peak$peakPos,
         loessPeakPosition = loessPeak,
         densityData = densityData,
         resampleSuccessRate = successRate,
         resampleSeed = seed,
         cutoff = inputList$cutoff,
         cutoffCenter = inputList$cutoffCenter,
         cutoffSpread = inputList$cutoffSpread,
         cutoffMethod = inputList$cutoffMethod,
         intervalMethod = peakIntervalMethod(mmapprData@param))
}


#' @title Identify chromosomes containing peaks
#'
#' @name prePeak
#'
#' Follows the \code{\link{loessFit}} step and precedes
#' \code{\link{peakRefinement}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \linkS4class{MmapprData} object with the \code{peaks}
#'   slot initialized.
#' @export
NULL


# [CHANGE — EXPLICIT INITIAL-PEAK CUTOFF SEMANTICS]
# The old implementation hard-coded a multiplier of 3 and used the mean of chromosome means plus
# 3 * sqrt(sum(var/n)). The new implementation exposes the multiplier and method: `legacy_current`
# applies that old formula to finite fitted values, while `global_sd` is an explicit alternative.
.calculatePeakCutoff <- function(snpDistance, method = c("legacy_current", "global_sd"), k = 3) {
    method <- match.arg(method)
    valid <- lapply(snpDistance, function(chr) {
        if (!is.list(chr) || !inherits(chr$loess, "loess")) return(numeric())
        chr$loess$fitted[is.finite(chr$loess$fitted)]
    })
    valid <- valid[vapply(valid, length, integer(1)) > 0L]
    if (length(valid) == 0L) stop("No finite LOESS values are available for peak thresholding")

    if (method == "legacy_current") {
        # [NOTE] This path uses the old cutoff formula: an unweighted mean of
        # chromosome means plus k * sqrt(sum(var/n)). The spread term is
        # standard-error-like rather than a literal genome-wide SD; keeping it
        # here preserves the old threshold semantics when this method is selected.
        centers <- vapply(valid, mean, numeric(1))
        varianceTerms <- vapply(valid, function(x) if (length(x) > 1L) stats::var(x) / length(x) else 0, numeric(1))
        center <- mean(centers)
        spread <- sqrt(sum(varianceTerms))
    } else {
        # [IMPROVE] Optional statistically clearer alternative: a robust global
        # center with the ordinary SD of all fitted marker values. This is still
        # a heuristic threshold, not a calibrated false-positive probability.
        allFitted <- unlist(valid, use.names = FALSE)
        if (length(allFitted) < 2L) stop("Insufficient finite LOESS values for global-SD thresholding")
        center <- stats::median(allFitted)
        spread <- stats::sd(allFitted)
    }
    list(cutoff = center + k * spread, center = center, spread = spread, method = method)
}


# [CHANGE — NA-SAFE PEAK DETECTION WITH PROVENANCE]
# The old implementation tested `any(fitted > cutoff)` directly and stored only the sequence name for
# chromosomes that passed. The new implementation ignores non-finite fitted values safely and stores
# the cutoff center, spread, value, and method with each initial peak for downstream diagnostics.
prePeak <- function(mmapprData) {
    mmapprData@peaks <- list()
    cutoffInfo <- .calculatePeakCutoff(
        mmapprData@snpDistance,
        method = peakCutoffMethod(mmapprData@param),
        k = peakCutoffSd(mmapprData@param)
    )

    for (i in seq_along(mmapprData@snpDistance)) {
        chr <- mmapprData@snpDistance[[i]]
        if (!is.list(chr) || !inherits(chr$loess, "loess")) next
        loessForChr <- chr$loess
        if (length(loessForChr$x) < 50L) next
        # [FIX] NA fitted values no longer make any() return NA and break if().
        containsPeak <- any(loessForChr$fitted > cutoffInfo$cutoff, na.rm = TRUE)
        chrName <- names(mmapprData@snpDistance)[[i]]
        if (containsPeak) {
            mmapprData@peaks[[chrName]] <- list(
                seqname = chrName,
                cutoff = cutoffInfo$cutoff,
                cutoffCenter = cutoffInfo$center,
                cutoffSpread = cutoffInfo$spread,
                cutoffMethod = cutoffInfo$method
            )
        }
    }

    .messageAndLog(sprintf("Peak cutoff (%s): center %.6g + %.3g*spread %.6g = %.6g",
                           cutoffInfo$method, cutoffInfo$center,
                           peakCutoffSd(mmapprData@param), cutoffInfo$spread,
                           cutoffInfo$cutoff),
                   outputFolder(mmapprData@param))
    mmapprData
}
