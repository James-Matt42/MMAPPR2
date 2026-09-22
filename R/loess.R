#' Perform optimized Loess regression for each chromosome
#'
#' Called after the \code{\link{calculateDistance}} step and before
#' \code{\link{prePeak}}.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{$loess} 
#'   element of the \code{snpDistance} slot list filled.
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
#'
#' postLoessMD <- loessFit(postCalcDistMD)
#' }
#' 
#' @export

# [CHANGE — NAMESPACE-STABLE PARALLEL LOESS]
# The old implementation called unqualified `bplapply`, relying on the symbol being attached/imported.
# The new implementation calls BiocParallel::bplapply explicitly and routes each chromosome through the
# hardened optimizer below.
loessFit <- function(mmapprData) {
    loessOptResolution <- mmapprData@param@loessOptResolution
    loessOptCutFactor <- mmapprData@param@loessOptCutFactor

    # each item (chr) of distance list has mutCounts, wtCounts,
    # distanceDf going in
    mmapprData@snpDistance <-
        BiocParallel::bplapply(mmapprData@snpDistance,
                 FUN = .loessFitForChr,   
                 loessOptResolution = loessOptResolution,
                 loessOptCutFactor = loessOptCutFactor)

    return(mmapprData)
}


# [CHANGE — LOESS FIT FAILURES BECOME OPTIMIZER DATA]
# The old implementation let many invalid/degenerate LOESS fits propagate as errors or
# non-finite AICc values. The new implementation wraps fitting so a bad span can be scored as failed
# and skipped without aborting the whole chromosome search.
.getLoess <- function(s, pos, eucDist, ...){
    suppressWarnings(try(stats::loess(eucDist ~ pos,
                                      span = s,
                                      degree = 1,
                                      family = "symmetric"),
                         silent = TRUE))
}


# [CHANGE — BOUNDARY-SAFE LOCAL SEARCH RESOLUTION]
# The old helper indexed diff() on both sides of a span and could step outside the vector at an
# endpoint. The new helper checks available neighbors explicitly and returns a safe spacing value.
.localResolution <- function(spans, span) {
    spans <- unique(sort(spans[is.finite(spans)]))
    index <- match(span, spans)
    if (is.na(index)) return(Inf)
    neighbors <- numeric()
    if (index > 1L) neighbors <- c(neighbors, span - spans[index - 1L])
    if (index < length(spans)) neighbors <- c(neighbors, spans[index + 1L] - span)
    if (length(neighbors) == 0L) Inf else max(neighbors)
}


# [CHANGE — AICC OPTIMIZER REWRITE]
# The old implementation's recursive search could evaluate the same span repeatedly and behaved
# poorly when fits were all invalid. The new implementation caches evaluated spans, performs an
# iterative coarse-to-fine search, handles an all-failed search explicitly, and
# adds an iteration cap. The objective remains AICc minimization.
# [IMPROVE] Coarse-to-fine AICc search with cached span evaluations. The old
# recursive implementation could refit the same span multiple times and handled
# all-NA fits poorly. This retains the same basic search idea while making failure
# explicit and avoiding duplicate LOESS work.
.aiccOpt <- function(distanceDf, spans, resolution, cutFactor) {
    evaluated <- data.frame(spans = numeric(), aiccValues = numeric())
    frontier <- unique(round(spans[spans > 0 & spans <= 1],
                             digits = .numDecimals(resolution)))

    converged <- FALSE
    for (iter in seq_len(50L)) {
        already <- evaluated$spans
        frontier <- setdiff(frontier, already)
        if (length(frontier) > 0L) {
            vals <- vapply(frontier, .aicc, FUN.VALUE = numeric(1),
                           eucDist = distanceDf$DISTANCE,
                           pos = distanceDf$POS)
            evaluated <- rbind(evaluated,
                               data.frame(spans = frontier, aiccValues = vals))
            evaluated <- evaluated[!duplicated(evaluated$spans), , drop = FALSE]
        }

        finite <- evaluated[is.finite(evaluated$aiccValues), , drop = FALSE]
        if (nrow(finite) == 0L)
            stop("All candidate LOESS/AICc fits failed for this chromosome")
        finite <- finite[order(finite$spans), , drop = FALSE]

        minimaIdx <- .localMinIndices(finite$aiccValues)
        if (length(minimaIdx) == 0L) minimaIdx <- which.min(finite$aiccValues)
        minimaIdx <- minimaIdx[order(finite$aiccValues[minimaIdx])]
        minimaIdx <- head(minimaIdx, 2L)
        minSpans <- finite$spans[minimaIdx]

        newSpans <- numeric()
        allSpans <- sort(unique(evaluated$spans))
        for (minSpan in minSpans) {
            localRes <- .localResolution(allSpans, minSpan)
            if (is.finite(localRes) && localRes > resolution) {
                step <- localRes * cutFactor
                addVector <- step * seq_len(9L)
                proposed <- c(minSpan - addVector, minSpan + addVector)
                proposed <- proposed[proposed > 0 & proposed <= 1]
                proposed <- round(proposed, digits = .numDecimals(resolution))
                newSpans <- c(newSpans, proposed)
            }
        }
        newSpans <- setdiff(unique(newSpans), evaluated$spans)
        if (length(newSpans) == 0L) {
            converged <- TRUE
            break
        }
        frontier <- newSpans
    }

    # [IMPROVE] Do not silently pretend the coarse-to-fine search converged if it
    # hit its defensive iteration cap. Returning the evaluated table is still useful,
    # but the warning makes an unexpectedly difficult objective visible to the user.
    if (!converged)
        warning("LOESS span optimization reached its 50-iteration safety cap before convergence")

    evaluated[order(evaluated$spans), , drop = FALSE]
}


# [CHANGE — FINITE/BOUNDARY-SAFE LOCAL MINIMA]
# The old helper returned minimum values from a diff-based expression, making repeated values, endpoints,
# and non-finite AICc entries fragile to match back to spans. The new helper returns indices directly,
# ignores non-finite values, and handles endpoint/repeated-extrema cases explicitly.
.localMinIndices <- function(x) {
    if (!length(x)) return(integer())
    out <- integer()
    for (i in seq_along(x)) {
        if (!is.finite(x[i])) next
        left <- if (i == 1L) Inf else x[i - 1L]
        right <- if (i == length(x)) Inf else x[i + 1L]
        if ((!is.finite(left) || x[i] <= left) &&
            (!is.finite(right) || x[i] <= right)) out <- c(out, i)
    }
    out
}


# [CHANGE — ROBUST NUMERIC RESOLUTION PARSING]
# The old helper inferred decimal places by applying string operations directly to the numeric value,
# which is fragile for scientific notation. The new helper formats a stable non-scientific string first.
.numDecimals <- function(x) {
    stopifnot(is.numeric(x), length(x) == 1L, is.finite(x), x > 0)
    txt <- format(x, scientific = FALSE, trim = TRUE, digits = 15)
    txt <- sub("0+$", "", txt)
    if (!grepl("\\.", txt)) 0L else nchar(sub("^.*\\.", "", txt))
}


# [CHANGE — DEGENERATE AICC GUARDS]
# The old AICc helper could return Inf/NaN for failed fits, zero residual variance, or a non-positive
# denominator. The new helper returns NA for those spans so the optimizer can exclude them explicitly.
.aicc <- function(s, eucDist, pos) {
    x <- .getLoess(s, pos, eucDist)
    if (inherits(x, "try-error")) return(NA_real_)
    n <- x$n
    traceL <- x$trace.hat
    denom <- n - traceL - 2
    if (!is.finite(denom) || denom <= 0 || n <= 1L) return(NA_real_)
    sigma2 <- sum(x$residuals^2, na.rm = TRUE) / (n - 1)
    if (!is.finite(sigma2) || sigma2 <= 0) return(NA_real_)
    value <- log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / denom
    if (is.finite(value)) value else NA_real_
}


# [CHANGE — DETERMINISTIC TIED-SPAN SELECTION]
# The old implementation could average tied optimum spans and return a span whose AICc had
# never actually been evaluated. The new implementation resolves floating-point ties to an evaluated
# span deterministically (choosing the smaller symmetric optimum).
.chooseBestAiccSpan <- function(aiccTable, resolution) {
    finite <- aiccTable[is.finite(aiccTable$aiccValues), , drop = FALSE]
    if (nrow(finite) == 0L) stop("No finite AICc values were produced")
    bestRows <- finite$aiccValues == min(finite$aiccValues)
    tiedSpans <- sort(unique(finite$spans[bestRows]))
    tiedMedian <- stats::median(tiedSpans)
    # [FIX] Return an ACTUALLY EVALUATED optimum. Averaging tied spans can create
    # an unevaluated value whose AICc is unknown. Ties equidistant from the median
    # resolve to the smaller evaluated span for deterministic behavior.
    distances <- abs(tiedSpans - tiedMedian)
    minDistance <- min(distances)
    # [FIX/R-VALIDATED] Binary floating-point can make mathematically symmetric
    # ties (for example 0.1 and 0.3 around 0.2) differ by ~1e-17. which.min()
    # could therefore choose the larger span despite the documented smaller-span
    # tie rule. Treat machine-scale differences as equal, then choose min().
    tol <- sqrt(.Machine$double.eps) * max(1, abs(tiedMedian), minDistance)
    tiedClosest <- tiedSpans[abs(distances - minDistance) <= tol]
    bestSpan <- min(tiedClosest)
    round(bestSpan, digits = .numDecimals(resolution))
}


# [CHANGE — CHROMOSOME-LEVEL LOESS HARDENING]
# The old chromosome wrapper attempted optimization without an explicit minimum-row check, did not record
# the chosen span separately, and did not verify that the final fit was a valid loess object. The new
# wrapper adds those checks/metadata and returns an explicit chromosome-level diagnostic on failure.
.loessFitForChr <- function(resultList, loessOptResolution, loessOptCutFactor){
    startTime <- proc.time()
    tryCatch({
        if (is(resultList, "character")) stop("-- Loess fit failed: ", resultList)
        if (is.null(resultList$distanceDf) || nrow(resultList$distanceDf) < 5L)
            stop("Too few informative positions for LOESS fitting")

        startSpans <- c(seq(.01, .16, .01), seq(.21, .91, .10))
        resultList$aicc <- .aiccOpt(distanceDf = resultList$distanceDf,
                                    spans = startSpans,
                                    resolution = loessOptResolution,
                                    cutFactor = loessOptCutFactor)

        bestSpan <- .chooseBestAiccSpan(resultList$aicc, loessOptResolution)
        resultList$bestSpan <- bestSpan

        resultList$loess <- .getLoess(bestSpan,
                                      resultList$distanceDf$POS,
                                      resultList$distanceDf$DISTANCE)
        if (inherits(resultList$loess, "try-error") ||
            !inherits(resultList$loess, "loess"))
            stop("Final LOESS fit failed at optimized span ", bestSpan)

        resultList$distanceDf <- NULL
        resultList$seqname <- NULL
        resultList$loessTime <- proc.time() - startTime
        resultList
    }, error = function(e) {
        if (is(resultList, "character")) paste0(resultList, ": ", e$message) else e$message
    })
}
