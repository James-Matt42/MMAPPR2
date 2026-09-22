#' @title Read BAM files and generate Euclidean distance data
#'
#' @name calculateDistance
#'
#' First step in the MMAPPR2 pipeline. Precedes the \code{\link{loessFit}}
#' step.
#'
#' @param mmapprData The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{snpDistance}
#'   slot filled.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
#'     mmappr_param <- mmapprParam(wtFiles = MMAPPR2data::exampleWTbam(),
#'                                 mutFiles = MMAPPR2data::exampleMutBam(),
#'                                 refFasta = MMAPPR2data::goldenFasta(),
#'                                 gtf = MMAPPR2data::gtf(),
#'                                 outputFolder = tempOutputFolder())
#'
#'     md <- mmapprData(mmappr_param)
#'     postCalcDistMD <- calculateDistance(md)
#' }
#' 
NULL

# [CHANGE — DISTANCE ENTRY-POINT HARDENING]
# The new implementation fails explicitly when annotation processing yields no pileup ranges and
# qualifies BiocParallel::bplapply rather than depending on an attached symbol.
# The old implementation proceeded directly to bplapply with no empty-range diagnostic.
calculateDistance <- function(mmapprData) {
  chrList <- suppressWarnings(.getFileReadChrList(param(mmapprData)))

  if (length(chrList) == 0L) stop("No annotated genomic ranges were available for pileup")

  mmapprData@snpDistance <-
    BiocParallel::bplapply(chrList, .calcDistForChr, param = mmapprData@param)

  return(mmapprData)
}


# [CHANGE — ANNOTATION QUERY-RANGE REWRITE]
# The old implementation shell-decompressed the GTF and parsed `gene` rows with fread.
# The new implementation keeps the same preference for explicit annotation gene spans, but obtains
# them through .annotationGeneRanges(), adds GFF/rtracklayer/TxDb fallbacks, aligns
# annotation naming with the FASTA, checks shared sequence names, handles custom
# assemblies more defensively, and recognizes additional mitochondrial aliases.
.getFileReadChrList <- function(param) {
  # Preserve the old implementation's choice of explicit annotation `gene` spans.
  # Standard GTFs stay on the lightweight parser; .annotationGeneRanges() builds
  # and closes a TxDb only as a fallback when explicit gene rows are unavailable.
  genes <- .annotationGeneRanges(param)

  if (length(genes) == 0L) stop("No gene ranges could be derived from the annotation")

  # Harmonize annotation naming to the FASTA where possible.
  fa_si <- .faSeqinfo(param@refGenome)
  targetStyle <- .choose_target_style(fa_si)
  if (!is.na(targetStyle)) {
    # [CHANGE — CENTRALIZED SAFE STYLE TRANSLATION]
    # The old linkage-range code assumed annotation sequence names were already usable as BAM/FASTA
    # query names. The new implementation uses one guarded style-conversion helper and falls back to
    # the annotation's exact names when a safe conversion cannot be established.
    genes <- tryCatch(.setSeqlevelsStyleFrozen(genes, targetStyle),
                      error = function(e) genes)
  }

  common <- intersect(GenomeInfoDb::seqlevels(genes), GenomeInfoDb::seqlevels(fa_si))
  if (length(common) == 0L)
    stop("Annotation and FASTA do not share any sequence names after style harmonization")
  genes <- GenomeInfoDb::keepSeqlevels(genes, common, pruning.mode = "coarse")

  if (!isTRUE(param@includeScaffolds)) {
    # [NOTE] Standard-chromosome heuristics are useful for common model organisms
    # but can remove valid contigs from custom assemblies. Users of such assemblies
    # should set includeScaffolds=TRUE.
    genes <- GenomeInfoDb::keepStandardChromosomes(genes, pruning.mode = "coarse")
    if (length(genes) == 0L)
      stop("No standard chromosomes remained after filtering. For a custom/non-model assembly, try includeScaffolds=TRUE.")
  }

  # Preserve the old chrM/MT exclusions and recognize a few equivalent
  # mitochondrial labels without changing autosomal/nuclear behavior.
  mitoNames <- intersect(GenomeInfoDb::seqlevels(genes),
                         c("chrM", "MT", "M", "MtDNA", "mitochondrion_genome"))
  if (length(mitoNames) > 0L)
    genes <- GenomeInfoDb::dropSeqlevels(genes, mitoNames, pruning.mode = "coarse")
  if (length(genes) == 0L)
    stop("No non-mitochondrial annotated gene ranges remain for mapping")

  # [CHANGE — OVERLAPPING-GENE DOUBLE-COUNT FIX]
  # The old implementation passed overlapping gene ranges directly to Rsamtools; overlapping
  # `which` ranges can return the same alignment more than once. The new implementation reduces the
  # ranges first so overlapping annotations cannot inflate depth or allele counts.
  # [FIX] Query intervals must be disjoint. Rsamtools treats overlapping `which`
  # ranges independently; without reduce(), the same read can be returned twice
  # and its depth can be double-counted in overlapping genes.
  genes <- GenomicRanges::reduce(genes, ignore.strand = TRUE)
  split(genes, as.character(GenomicRanges::seqnames(genes)))
}


# [CHANGE — DISTANCE FILTER/AGGREGATION CORRECTIONS]
# The new implementation retains per-file identity while aggregating replicates, averages coverage
# across all supplied files, changes the minimum-depth boundary from old `>`
# to inclusive `>=`, checks emptiness with nrow(), joins WT/mutant tables explicitly
# on CHROM+POS, fixes the WT-homozygosity column selection, and removes non-finite
# distance rows. The Euclidean-distance formula itself is unchanged.
.calcDistForChr <- function(chrRange, param){
  tryCatch({
    stopifnot(length(unique(GenomicRanges::seqnames(chrRange))) == 1L)
    stopifnot(is(param, "MmapprParam"))
    stopifnot(length(chrRange) > 0L)

    # Keep sample identity long enough to make replicate aggregation semantics
    # explicit. In particular, coverage is averaged across *all* input files so
    # a replicate with zero coverage cannot silently inflate AVE.CVG.
    pileupWT <- lapply(seq_along(wtFiles(param)), function(i) {
      x <- .getPileup(wtFiles(param)[[i]], param, chrRange)
      if (nrow(x)) x[, FILE_ID := i]
      x
    })
    pileupWT <- data.table::rbindlist(pileupWT, use.names = TRUE, fill = TRUE)
    if (nrow(pileupWT) == 0L) stop("Insufficient data in wild-type file(s)")
    pileupWT <- .avgFiles(pileupWT,
                          fileAggregation = fileAggregation(param),
                          nFiles = length(wtFiles(param)))
    # [FIX] "minimum depth" now means >= threshold rather than strictly >.
    pileupWT <- pileupWT[AVE.CVG >= minDepth(param)]
    if (nrow(pileupWT) == 0L) stop("Insufficient wild-type depth after filtering")

    pileupMut <- lapply(seq_along(mutFiles(param)), function(i) {
      x <- .getPileup(mutFiles(param)[[i]], param, chrRange)
      if (nrow(x)) x[, FILE_ID := i]
      x
    })
    pileupMut <- data.table::rbindlist(pileupMut, use.names = TRUE, fill = TRUE)
    if (nrow(pileupMut) == 0L) stop("Insufficient data in mutant file(s)")
    pileupMut <- .avgFiles(pileupMut,
                           fileAggregation = fileAggregation(param),
                           nFiles = length(mutFiles(param)))
    pileupMut <- pileupMut[AVE.CVG >= minDepth(param)]
    if (nrow(pileupMut) == 0L) stop("Insufficient mutant depth after filtering")

    data.table::setkey(pileupWT, CHROM, POS)
    data.table::setkey(pileupMut, CHROM, POS)
    # [IMPROVE/R-VALIDATED] Make the scientific join key explicit. data.table's
    # current merge() defaults would infer the shared keyed columns, but naming
    # CHROM/POS here prevents a future extra shared column from silently changing
    # which WT and mutant observations are paired.
    # [CHANGE — EXPLICIT WT/MUTANT GENOMIC JOIN KEY]
    # The old implementation relied on merge() inferring all shared columns. The new implementation declares
    # CHROM/POS so future shared metadata cannot silently alter which observations
    # are paired.
    distanceDf <- merge(pileupWT, pileupMut, by = c("CHROM", "POS"),
                        suffixes = c(".WT", ".MT"))

    # [FIX] length(data.frame) counts columns, not rows.
    if (nrow(distanceDf) == 0L)
      stop("Empty dataframe after joining WT and mutant count tables")

    distanceDf[, DISTANCE := sqrt((AVE.A.FREQ.WT - AVE.A.FREQ.MT)^2 +
                                  (AVE.C.FREQ.WT - AVE.C.FREQ.MT)^2 +
                                  (AVE.G.FREQ.WT - AVE.G.FREQ.MT)^2 +
                                  (AVE.T.FREQ.WT - AVE.T.FREQ.MT)^2) ^
                                  distancePower(param)]

    # Retain markers that are polymorphic in the phenotypically WT F2 pool.
    # [FIX] Select the four WT frequency columns explicitly. In data.table,
    # `AVE.A.FREQ.WT:AVE.T.FREQ.WT` evaluates `:` on the column vectors; it is
    # not a safe column-range selector and can produce nonsensical input to apply().
    wtFreqCols <- c("AVE.A.FREQ.WT", "AVE.C.FREQ.WT",
                    "AVE.G.FREQ.WT", "AVE.T.FREQ.WT")
    homozygousWT <- .wtHomozygousRows(distanceDf, wtFreqCols,
                                      cutoff = homozygoteCutoff(param))
    distanceDf <- distanceDf[!homozygousWT & is.finite(DISTANCE)]
    if (nrow(distanceDf) == 0L) stop("No informative markers remain after filtering")

    resultList <- list(wtCounts = pileupWT,
                       mutCounts = pileupMut,
                       distanceDf = distanceDf)

    .messageAndLog(paste("Finished", as.character(GenomicRanges::seqnames(chrRange)[1])),
                   outputFolder = outputFolder(param))
    resultList
  }, error = function(e) {
    paste(toString(unique(GenomicRanges::seqnames(chrRange))), e$message, sep = ": ")
  })
}


# Keep the WT-homozygosity rule in a pure helper so the old column-selection bug
# can be regression-tested without constructing BAM files.
# [CHANGE — WT HOMOZYGOSITY FILTER BUG FIX]
# The old expression `AVE.A.FREQ.WT:AVE.T.FREQ.WT` applied `:` to vectors
# rather than selecting four columns. The new implementation passes the WT A/C/G/T columns explicitly,
# uses only finite evidence, and keeps the homozygosity rule independently testable.
.wtHomozygousRows <- function(distanceDf, wtFreqCols, cutoff) {
  stopifnot(length(wtFreqCols) == 4L, all(wtFreqCols %in% names(distanceDf)))
  mat <- as.matrix(distanceDf[, ..wtFreqCols])
  # Frequencies should already be finite, but make malformed rows conservative:
  # a row with no finite frequency evidence is not useful as a mapping marker.
  major <- apply(mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) Inf else max(x)
  })
  major > cutoff
}


# [FIX] Use one explicit alignment-filter policy anywhere reads contribute to
# allele counts or expression summaries. The old code did not explicitly reject
# secondary, supplementary, failed-QC, or unmapped records. The new filter rejects
# those records while leaving duplicate status unspecified, so duplicate-marked
# primary reads remain eligible just as they were before.
# [CHANGE — EXPLICIT PRIMARY-READ FILTER POLICY]
# The old pileup did not explicitly exclude secondary, supplementary,
# failed-QC, or unmapped records. The new implementation excludes those records anywhere reads
# contribute to linkage/candidate/expression evidence. Duplicate-marked primary
# reads remain eligible, preserving the old duplicate-read behavior.
.primaryMappedScanFlag <- function() {
  Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE
  )
}


# Imports per-position A/C/G/T data from one BAM file.
# [CHANGE — RNA-SEQ PILEUP AND BASE-TABLE HARDENING]
# The old implementation used simpleCigar=TRUE, which discards ordinary spliced reads with
# CIGAR `N`, fixed max_depth at 1000, and assumed every pileup had A/C/G/T columns.
# The new implementation accepts complex RNA-seq CIGARs, applies the shared primary-read filter,
# exposes linkage maxPileupDepth, ignores query Ns, handles empty pileups, and fills
# absent nucleotide categories explicitly before calculating frequencies.
.getPileup <- function(file, param, chrRange) {
  stopifnot(length(file) == 1L)

  # [FIX] simpleCigar=TRUE discarded ordinary spliced RNA-seq reads containing
  # CIGAR N operations, as well as clipped/indel-containing reads. pileup itself
  # understands reference skips, so do not throw away those reads wholesale.
  scanParam <- Rsamtools::ScanBamParam(
    flag = .primaryMappedScanFlag(),
    simpleCigar = FALSE,
    which = chrRange,
    mapqFilter = as.integer(minMapQuality(param))
  )

  pParam <- Rsamtools::PileupParam(
    max_depth = as.integer(maxPileupDepth(param)),
    min_mapq = as.integer(minMapQuality(param)),
    min_base_quality = as.integer(minBaseQuality(param)),
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,
    # [NOTE] The linkage statistic is intentionally A/C/G/T based. Indel support
    # belongs in candidate calling and requires exact allele reconstruction; merely
    # counting '+'/'-' pileup symbols would not provide a valid normalized indel.
    include_deletions = FALSE,
    include_insertions = FALSE
  )

  raw <- Rsamtools::pileup(file, scanBamParam = scanParam, pileupParam = pParam)
  if (NROW(raw) == 0L) return(.emptyPileupTable())

  pData <- data.table::as.data.table(raw)
  pData <- pData[nucleotide %in% c("A", "C", "G", "T")]
  if (nrow(pData) == 0L) return(.emptyPileupTable())

  pData <- data.table::dcast(pData, seqnames + pos ~ nucleotide,
                             value.var = "count", fun.aggregate = sum, fill = 0)

  data.table::setnames(pData, c("seqnames", "pos"), c("CHROM", "POS"))
  # [FIX] Small regions do not necessarily contain all four nucleotide categories.
  # Create absent columns explicitly rather than renaming by column position.
  for (base in c("A", "C", "G", "T")) {
    if (!base %in% names(pData)) pData[, (base) := 0]
  }
  data.table::setnames(pData, c("A", "C", "G", "T"),
                       c("A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ"))
  pData <- pData[, .(CHROM, POS, A.FREQ, C.FREQ, G.FREQ, T.FREQ)]

  pData[, CVG := A.FREQ + C.FREQ + G.FREQ + T.FREQ]
  pData <- pData[CVG > 0]
  pData[, c("A.FREQ", "C.FREQ", "G.FREQ", "T.FREQ") :=
          list(A.FREQ/CVG, C.FREQ/CVG, G.FREQ/CVG, T.FREQ/CVG)]
  pData
}

.emptyPileupTable <- function() {
  data.table::data.table(CHROM = character(), POS = integer(),
                         A.FREQ = numeric(), C.FREQ = numeric(),
                         G.FREQ = numeric(), T.FREQ = numeric(),
                         CVG = numeric())
}


# Aggregate replicate pileups.
# [CHANGE — DEFINED REPLICATE AGGREGATION SEMANTICS]
# The old implementation lost sample identity before averaging, so zero-coverage replicates
# disappeared from the coverage denominator. The new implementation defines `simple` as equal
# frequency weighting among replicates with evidence and `weighted` as pooled-read
# frequency weighting, while AVE.CVG always divides by the total number of supplied
# files so missing coverage contributes zero rather than vanishing.
.avgFiles <- function(chrDf, fileAggregation, nFiles) {
  stopifnot(fileAggregation %in% c("simple", "weighted"))
  stopifnot(nFiles >= 1L)
  data.table::setkey(chrDf, CHROM, POS)

  if (fileAggregation == "weighted") {
    # Count-weighted allele frequencies are equivalent to pooling the reads.
    chrDf <- chrDf[, .(
      AVE.CVG = sum(CVG) / nFiles,
      AVE.A.FREQ = sum(A.FREQ * CVG) / sum(CVG),
      AVE.C.FREQ = sum(C.FREQ * CVG) / sum(CVG),
      AVE.G.FREQ = sum(G.FREQ * CVG) / sum(CVG),
      AVE.T.FREQ = sum(T.FREQ * CVG) / sum(CVG)
    ), .(CHROM, POS)]
  } else {
    # Equal-replicate weighting: frequency is averaged over replicates with data,
    # while coverage is averaged over all supplied files (missing coverage = 0).
    chrDf <- chrDf[, .(
      AVE.CVG = sum(CVG) / nFiles,
      AVE.A.FREQ = mean(A.FREQ),
      AVE.C.FREQ = mean(C.FREQ),
      AVE.G.FREQ = mean(G.FREQ),
      AVE.T.FREQ = mean(T.FREQ)
    ), .(CHROM, POS)]
  }
  chrDf
}
