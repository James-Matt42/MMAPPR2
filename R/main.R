#' @title Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
#'
#' @name mmappr
#'
#' MMAPPR2 is designed to map the causative mutation in a forward genetics
#' screen. It analyzes aligned sequence files, calculates the per-base
#' Euclidean distance between the mutant and wild-type pools, performs
#' a Loess regression on that distance, and generates candidate variants
#' in regions of peak distance.
#'
#' @param mmapprParam A \code{\linkS4class{MmapprParam}} object containing
#'   desired parameters.
#'
#' @return A \code{\linkS4class{MmapprData}} object containing results
#'   and/or intermediate data.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly = TRUE)) {
#'
#'     # Specify parameters:
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     # Run pipeline:
#'     mmapprData <- mmappr(mmappr_param)
#' }
#' 
#' \dontrun{
#' ### Alternately, you can navigate the pipeline step by step.
#' ### This may be helpful for debugging.
#' md <- mmapprData(mmapprParam)
#' postCalcDistMD <- calculateDistance(md)
#' postLoessMD <- loessFit(postCalcDistMD)
#' postPrePeakMD <- prePeak(postLoessMD)
#' postPeakRefMD <- peakRefinement(postPrePeakMD)
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' outputMmapprData(postCandidatesMD)
#' }
#'
#' @seealso \code{\link{calculateDistance}}, \code{\link{loessFit}},
#'   \code{\link{prePeak}}, \code{\link{peakRefinement}},
#'   \code{\link{generateCandidates}}, \code{\link{outputMmapprData}}
#'   
NULL

mmappr <- function(mmapprParam) {
    startTime <- Sys.time()
    message("------------------------------------")
    message("-------- Welcome to MMAPPR2 --------")
    message("------------------------------------\n")

    # [FIX] The current implementation uses Rsamtools for BAM pileup and no
    # longer requires an external samtools executable at startup. The old check
    # was stale and could reject an otherwise valid Bioconductor installation.
    md <- mmapprData(mmapprParam)
    oF <- outputFolder(md@param)
    .messageAndLog(paste('Start time:', Sys.time()), oF)
    .messageAndLog(paste('Output folder:', oF, '\n'), oF)
    .log('Parameters:', oF)
    .log(mmapprParam, oF)
    .log('', oF)

    md <- tryCatch({
        .messageAndLog("Reading BAM files and generating Euclidean distance data...", oF)
        md <- calculateDistance(md)
        .messageAndLog("Generating optimal Loess curves for each chromosome...", oF)
        md <- loessFit(md)
        .messageAndLog("Identifying chromosome(s) harboring linkage region...", oF)
        md <- prePeak(md)
        if (length(peaks(md)) > 0)
            .messageAndLog("Peak regions successfully identified", oF)
        else stop("No peak regions identified")
        .messageAndLog("Refining peak characterization using SNP resampling...", oF)
        md <- peakRefinement(md)
        .messageAndLog('Generating, analyzing, and ranking candidate variants...', oF)
        md <- generateCandidates(md)
        .messageAndLog("Writing output plots and tables...", oF)
        outputMmapprData(md)

        md  # return for use after block
    },
    error = function(e) {
        .messageAndLog(paste('ERROR:', e$message), oF)
        .messageAndLog("MmapprData object up to the failing step is returned.", oF)
        .messageAndLog(paste0("You can also recover this object ",
            "from 'mmappr_data.RDS' in the '", outputFolder(param(md)),
            "' output folder"), oF)
        md
    })

    endTime <- Sys.time()
    .messageAndLog(paste('\nEnd time:', endTime), oF)
    .messageAndLog(paste("MMAPPR2 runtime:", format(endTime - startTime)), oF)
    .atomicSaveRDS(md, file.path(md@param@outputFolder, "mmappr_data.RDS"))
    .log('\nsessionInfo()', oF)
    .log(sessionInfo(), oF)

    return(md)
}


.messageAndLog <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file = logFile, sep = '\n', append = TRUE)
    msg <- paste(msg, collapse = '\n')
    message(msg)
}


.log <- function(msg, outputFolder) {
    logFile <- file.path(outputFolder, 'mmappr2.log')
    if (!is.character(msg)) msg <- capture.output(msg)
    cat(msg, file = logFile, sep = '\n', append = TRUE)
}


# [FIX] The former .checkDep('samtools') helper was removed. The revised current
# implementation performs BAM pileups through Rsamtools and no longer shells out
# to samtools, so requiring an unrelated executable was a stale startup failure.
