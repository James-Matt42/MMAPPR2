## ----setup, eval=FALSE--------------------------------------------------------
# param <- mmapprParam(
#   wtFiles = c("wt_pool.bam"),
#   mutFiles = c("mutant_pool.bam"),
#   refFasta = "reference.fa",
#   gtf = "annotation.gtf",
#   outputFolder = "mmappr2_results"
# )

## ----parameters, eval=FALSE---------------------------------------------------
# param <- mmapprParam(
#   wtFiles = "wt_pool.bam",
#   mutFiles = "mutant_pool.bam",
#   refFasta = "reference.fa",
#   gtf = "annotation.gtf",
#   minDepth = 20,
#   minBaseQuality = 20,
#   minMapQuality = 20,
#   maxPileupDepth = 1000,
#   candidateMinDepth = 1,
#   candidateMinAltDepth = 2,
#   candidateMinAltFreq = 0.80,
#   # These WT-specific hard filters are optional by default. WT evidence is
#   # still reported for every candidate but does not change frozen default ranking.
#   candidateMaxWtAltFreq = 1.00,
#   candidateMinDeltaAF = 0.00,
#   peakCutoffSd = 3,
#   # Compatibility-safe default: preserve the supplied main-branch cutoff.
#   peakCutoffMethod = "legacy_current",
#   peakIntervalWidth = 0.80,
#   # Conservative default: span all selected high-density modes.
#   peakIntervalMethod = "hpd_span",
#   peakResampleIterations = 1000,
#   randomSeed = 1,
#   pairedEnd = FALSE,
#   ignoreStrand = FALSE
# )

## ----run, eval=FALSE----------------------------------------------------------
# md <- mmappr(param)

## ----stages, eval=FALSE-------------------------------------------------------
# md <- mmapprData(param)
# md <- calculateDistance(md)  # pooled A/C/G/T frequencies and ED^p
# md <- loessFit(md)           # chromosome-wise robust LOESS
# md <- prePeak(md)            # genome-wide peak threshold
# md <- peakRefinement(md)     # reproducible half-marker resampling
# md <- generateCandidates(md) # mutant SNVs + WT evidence + coding effects
# outputMmapprData(md)

