% Manually synchronized with R/main.R in the revised archive.
\name{mmappr}
\alias{mmappr}
\title{Mutation Mapping Analysis Pipeline for Pooled RNA-Seq}
\usage{mmappr(mmapprParam)}
\arguments{\item{mmapprParam}{A validated \code{MmapprParam} object.}}
\value{A \code{MmapprData} object containing results or, if a stage fails, the intermediate state reached before the failure.}
\description{Runs distance calculation, optimized LOESS smoothing, initial peak detection, resampling-based peak refinement, candidate generation/annotation, and output generation.}
\seealso{\code{\link{calculateDistance}}, \code{\link{loessFit}}, \code{\link{prePeak}}, \code{\link{peakRefinement}}, \code{\link{generateCandidates}}, \code{\link{outputMmapprData}}}
