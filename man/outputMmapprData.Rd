% Manually synchronized with R/output.R.
\name{outputMmapprData}
\alias{outputMmapprData}
\title{Generate plots and tables from MMAPPR2 data}
\usage{outputMmapprData(mmapprData)}
\arguments{\item{mmapprData}{A \code{MmapprData} object containing results to write.}}
\value{The same \code{MmapprData} object, invisibly.}
\description{Writes genome and peak plots plus candidate-variant, coding-effect, and descriptive expression-summary tables to the configured output folder. When \code{exportAiccPlots(param(mmapprData))} is TRUE and finite AICc evaluations are available, it also writes \code{aicc_plots.pdf} without rerunning LOESS fitting.}
