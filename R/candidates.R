#' @title Generate candidate mutations and consequences in peak regions
#'
#' @name generateCandidates
#'
#' Follows the \code{\link{peakRefinement}} step and produces a
#' \code{\linkS4class{MmapprData}} object ready for
#' \code{\link{outputMmapprData}}.
#'
#' @param md The \code{\linkS4class{MmapprData}} object to be analyzed.
#'
#' @return A \code{\linkS4class{MmapprData}} object with the \code{candidates}
#'   slot filled with a \code{\link[GenomicRanges]{GRanges}} object for each
#'   peak chromosome containing variants and predicted consequences.
#' @export
#'
#' @examples
#' if (requireNamespace('MMAPPR2data', quietly=TRUE)) {
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
#'
#' postCandidatesMD <- generateCandidates(postPeakRefMD)
#' }
#'
NULL

# ---- main API ------------------------------------------------------------------

generateCandidates <- function(md) {
  if (length(md@peaks) == 0L) stop("No refined peaks are available for candidate generation")
  .messageAndLog("Getting variants in refined peak interval(s)", outputFolder(param(md)))
  peakGRanges <- lapply(md@peaks, .getPeakRange)

  .messageAndLog("Peak Summary:", outputFolder(param(md)))
  for (seqname in names(md@peaks)) {
    peak <- md@peaks[[seqname]]
    log_text <- paste0(
      "seqname: ", seqname, "\n",
      "refined interval: ", as.integer(peak$start), "-", as.integer(peak$end), "\n",
      "full-data LOESS apex: ", as.integer(round(peak$loessPeakPosition)), "\n",
      "resampling-density apex: ", as.integer(round(peak$densityPeakPosition)), "\n",
      "resampling success rate: ", round(peak$resampleSuccessRate, 4), "\n",
      "resampling seed: ", peak$resampleSeed, "\n"
    )
    .messageAndLog(log_text, outputFolder(param(md)))
  }

  # Candidate pileup does not need a TxDb. Do that comparatively cheap work first
  # and construct transcript annotation only when coding-effect prediction actually
  # has variants to annotate. Ordinary GTF expression counting uses explicit gene
  # rows and therefore also avoids unnecessary SQLite construction on empty peaks.
  md@candidates$snps <- lapply(peakGRanges, .getVariantsForRange, param = md@param)
  haveVariants <- any(vapply(md@candidates$snps, length, integer(1)) > 0L)

  txdb <- NULL
  if (haveVariants) {
    txdb <- .buildTxDb(md@param)
    on.exit(.disconnectTxDb(txdb), add = TRUE)
  }
  annotationGenes <- .annotationGeneRanges(md@param, txdb)

  .messageAndLog("Predicting coding variant effects", outputFolder(param(md)))
  if (haveVariants) {
    md@candidates$effects <- lapply(md@candidates$snps, .predictEffects,
                                    param = md@param, txdb = txdb)
  } else {
    md@candidates$effects <- lapply(md@candidates$snps,
                                    function(x) GenomicRanges::GRanges())
  }

  .messageAndLog("Summarizing expression changes in peak interval(s)", outputFolder(param(md)))
  md@candidates$diff <- lapply(peakGRanges, .addDiff,
                               param = md@param, txdb = txdb,
                               annotationGenes = annotationGenes)

  .messageAndLog("Scoring candidates by peak position", outputFolder(param(md)))
  md@candidates <- .scoreVariants(md@candidates, md@peaks)
  md
}


# ---- internals -----------------------------------------------------------------

.getPeakRange <- function(peakList) {
  GenomicRanges::GRanges(seqnames = peakList$seqname,
                         ranges = IRanges::IRanges(start = as.integer(peakList$start),
                                                  end = as.integer(peakList$end)))
}


.getVariantsForRange <- function(inputRange, param) {
  fa <- param@refGenome
  fa_si <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)
  if (!is.na(tstyle))
    inputRange <- tryCatch(.setSeqlevelsStyleFrozen(inputRange, tstyle),
                           error = function(e) inputRange)

  # [IMPROVE] Pool replicate counts in memory instead of creating a fixed-name
  # merged.tmp.bam. This removes a concurrency hazard, avoids temporary BAM I/O,
  # and guarantees that samples with zero ALT reads still contribute to depth.
  mutPile <- .pooledBasePileup(
    bams = mutFiles(param), genome = fa, which = inputRange,
    minBaseQuality = minBaseQuality(param),
    minMapQuality = minMapQuality(param),
    # Original candidate calling left Rsamtools at its historical 250-read
    # pileup cap. Linkage maxPileupDepth is intentionally not reused here.
    maxDepth = .FROZEN_CANDIDATE_MAX_DEPTH
  )
  if (nrow(mutPile) == 0L) return(VariantAnnotation::VRanges())

  pile <- .basePileupToVRanges(mutPile, sampleName = "pooled_mutant")
  if (length(pile) == 0L) return(VariantAnnotation::VRanges())

  # [FIX/IMPROVE] The prior revision still sent these already-aggregated tallies
  # through VariantTools::callVariants(), which adds a hidden binomial/read-count
  # filter before MMAPPR2's much stricter ALT-frequency rule. Make every hard
  # candidate criterion explicit and user-visible instead. With the original
  # >80% ALT rule, VariantTools' ~4% likelihood filter is redundant; its default
  # two-ALT-read requirement is represented directly by candidateMinAltDepth.
  mutAF <- VariantAnnotation::altDepth(pile) / VariantAnnotation::totalDepth(pile)
  keep <- .candidateMutantKeep(
    totalDepth = VariantAnnotation::totalDepth(pile),
    altDepth = VariantAnnotation::altDepth(pile),
    minDepth = candidateMinDepth(param),
    minAltDepth = candidateMinAltDepth(param),
    minAltFreq = candidateMinAltFreq(param)
  )
  resultVr <- pile[keep]
  mutAF <- mutAF[keep]
  if (length(resultVr) == 0L) return(resultVr)

  # [FIX/IMPROVE] Candidate discovery used to compare mutant reads only with the
  # reference. Add the phenotypically WT F2 pool as evidence at every candidate.
  # We report WT depth/AF for all candidates and expose optional WT/delta-AF
  # filters. Defaults are deliberately permissive because a recessive F2 WT pool
  # legitimately contains heterozygotes (~1/3 mutant-allele frequency at the locus).
  resultVr <- .addWtAlleleSupport(resultVr, param)
  S4Vectors::mcols(resultVr)$mutAltFreq <- mutAF
  S4Vectors::mcols(resultVr)$mutRefDepth <- VariantAnnotation::refDepth(resultVr)
  S4Vectors::mcols(resultVr)$mutAltDepth <- VariantAnnotation::altDepth(resultVr)
  S4Vectors::mcols(resultVr)$mutTotalDepth <- VariantAnnotation::totalDepth(resultVr)
  S4Vectors::mcols(resultVr)$deltaAltFreq <- mutAF - S4Vectors::mcols(resultVr)$wtAltFreq
  # [IMPROVE] Preserve the WT evidence even when the optional hard filters are
  # disabled. A high-frequency ALT allele in both pools is a strong hint that the
  # site is simply a background difference from the reference genome.
  wtEnough <- S4Vectors::mcols(resultVr)$wtTotalDepth >= minDepth(param)
  S4Vectors::mcols(resultVr)$wtEvidenceSufficient <- wtEnough
  S4Vectors::mcols(resultVr)$sharedBackgroundLike <-
    wtEnough & is.finite(S4Vectors::mcols(resultVr)$wtAltFreq) &
    S4Vectors::mcols(resultVr)$wtAltFreq > candidateMinAltFreq(param)

  wtAF <- S4Vectors::mcols(resultVr)$wtAltFreq
  deltaAF <- S4Vectors::mcols(resultVr)$deltaAltFreq
  wtFilterRequested <- candidateMaxWtAltFreq(param) < 1
  deltaFilterRequested <- candidateMinDeltaAF(param) > 0

  keep <- .candidateWtKeep(
    wtAF = wtAF, deltaAF = deltaAF,
    maxWtAltFreq = candidateMaxWtAltFreq(param),
    minDeltaAF = candidateMinDeltaAF(param),
    wtFilterRequested = wtFilterRequested,
    deltaFilterRequested = deltaFilterRequested
  )
  resultVr[keep]
}


# Pure threshold helpers make the candidate rules explicit and independently
# testable. Keeping them free of BAM/annotation I/O also prevents a future caller
# refactor from silently changing the biological thresholds.
.candidateMutantKeep <- function(totalDepth, altDepth, minDepth, minAltDepth, minAltFreq) {
  af <- ifelse(totalDepth > 0, altDepth / totalDepth, NA_real_)
  is.finite(totalDepth) & is.finite(altDepth) & is.finite(af) &
    totalDepth >= 0 & altDepth >= 0 & altDepth <= totalDepth &
    totalDepth >= minDepth & altDepth >= minAltDepth & af > minAltFreq
}

.candidateWtKeep <- function(wtAF, deltaAF, maxWtAltFreq, minDeltaAF,
                             wtFilterRequested = maxWtAltFreq < 1,
                             deltaFilterRequested = minDeltaAF > 0) {
  keep <- rep(TRUE, length(wtAF))
  if (wtFilterRequested) keep <- keep & is.finite(wtAF) & wtAF <= maxWtAltFreq
  if (deltaFilterRequested) keep <- keep & is.finite(deltaAF) & deltaAF >= minDeltaAF
  keep
}


# Pooled A/C/G/T counts across one or more BAMs. Counts are summed before
# frequencies are calculated, which is the natural representation of a pooled
# sequencing library and avoids losing zero-ALT samples during replicate merge.
.pooledBasePileup <- function(bams, genome, which = NULL,
                              minBaseQuality = 0L,
                              minMapQuality = 0L,
                              maxDepth = 1000L) {
  # [FIX] BamFileList() is dots-style in Rsamtools. Construct it explicitly for
  # character vectors so multiple BAM paths cannot be interpreted as one list-like
  # argument on versions with stricter S4 coercion.
  if (is.character(bams)) {
    bams <- do.call(Rsamtools::BamFileList, lapply(bams, Rsamtools::BamFile))
  }
  if (is(bams, "BamFile")) bams <- Rsamtools::BamFileList(bams)
  if (!is(bams, "BamFileList") || length(bams) == 0L)
    stop("bams must contain at least one BAM file")

  pupar <- Rsamtools::PileupParam(
    max_depth = as.integer(maxDepth),
    distinguish_nucleotides = TRUE,
    distinguish_strands = FALSE,
    min_base_quality = as.integer(minBaseQuality),
    min_mapq = as.integer(minMapQuality),
    ignore_query_Ns = TRUE,
    include_insertions = FALSE,
    include_deletions = FALSE
  )
  # [FIX] Apply the same primary-alignment policy used by linkage pileup so a
  # chimeric/split read cannot contribute both its primary and supplementary
  # alignments to candidate allele depth. Secondary alignments are excluded too.
  sbpar <- if (!is.null(which))
    Rsamtools::ScanBamParam(flag = .primaryMappedScanFlag(),
                            which = which, simpleCigar = FALSE,
                            mapqFilter = as.integer(minMapQuality))
  else Rsamtools::ScanBamParam(flag = .primaryMappedScanFlag(),
                               simpleCigar = FALSE,
                               mapqFilter = as.integer(minMapQuality))

  pieces <- lapply(seq_along(bams), function(i) {
    x <- Rsamtools::pileup(bams[[i]], scanBamParam = sbpar, pileupParam = pupar)
    if (NROW(x) == 0L) return(NULL)
    x <- data.table::as.data.table(x)
    x[nucleotide %in% c("A", "C", "G", "T"), .(seqnames, pos, nucleotide, count)]
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0L) return(.emptyPooledBasePileup())

  long <- data.table::rbindlist(pieces)
  long <- long[, .(count = sum(count)), .(seqnames, pos, nucleotide)]
  wide <- data.table::dcast(long, seqnames + pos ~ nucleotide,
                            value.var = "count", fun.aggregate = sum, fill = 0)
  for (base in c("A", "C", "G", "T")) if (!base %in% names(wide)) wide[, (base) := 0]
  wide <- wide[, .(seqnames, pos, A, C, G, T)]
  wide[, totalDepth := A + C + G + T]
  wide <- wide[totalDepth > 0]
  if (nrow(wide) == 0L) return(.emptyPooledBasePileup())

  gr <- GenomicRanges::GRanges(seqnames = wide$seqnames,
                               ranges = IRanges::IRanges(wide$pos, width = 1L))
  wide[, ref := as.character(Rsamtools::getSeq(genome, gr))]
  wide
}

.emptyPooledBasePileup <- function() {
  data.table::data.table(seqnames = character(), pos = integer(),
                         A = integer(), C = integer(), G = integer(), T = integer(),
                         totalDepth = integer(), ref = character())
}


.basePileupToVRanges <- function(pile, sampleName = "pooled") {
  if (nrow(pile) == 0L) return(VariantAnnotation::VRanges())
  bases <- c("A", "C", "G", "T")
  countMatrix <- as.matrix(pile[, ..bases])
  refIndex <- match(pile$ref, bases)
  refDepthVec <- rep(0L, nrow(pile))
  validRef <- !is.na(refIndex)
  if (any(validRef))
    refDepthVec[validRef] <- countMatrix[cbind(which(validRef), refIndex[validRef])]

  # [FIX] data.table modifies by reference. Work on a copy so converting a pileup
  # to VRanges does not unexpectedly add an internal .row_id column to the caller's
  # object (important for tests and for future reuse of the pileup table).
  pile <- data.table::copy(pile)
  pile[, .row_id := .I]
  long <- data.table::melt(pile,
                           id.vars = c(".row_id", "seqnames", "pos", "totalDepth", "ref"),
                           measure.vars = bases,
                           variable.name = "alt", value.name = "altDepth")
  long[, alt := as.character(alt)]
  long <- long[alt != ref & altDepth > 0 & ref %in% bases]
  if (nrow(long) == 0L) return(VariantAnnotation::VRanges())
  # [FIX] Keep the source vector name distinct from the data.table column name.
  # This avoids data.table scope ambiguity during := evaluation.
  long[, refDepth := refDepthVec[.row_id]]

  VariantAnnotation::VRanges(
    seqnames = S4Vectors::Rle(long$seqnames),
    ranges = IRanges::IRanges(start = long$pos, width = 1L),
    ref = long$ref,
    alt = as.character(long$alt),
    refDepth = S4Vectors::Rle(as.integer(long$refDepth)),
    altDepth = S4Vectors::Rle(as.integer(long$altDepth)),
    totalDepth = S4Vectors::Rle(as.integer(long$totalDepth)),
    sampleNames = S4Vectors::Rle(sampleName, nrow(long))
  )
}


.matchVariantToPileup <- function(variants, pile) {
  if (length(variants) == 0L) return(integer())
  if (nrow(pile) == 0L) return(rep(NA_integer_, length(variants)))
  varGr <- GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(variants)),
    ranges = IRanges::IRanges(BiocGenerics::start(variants), width = 1L)
  )
  pileGr <- GenomicRanges::GRanges(
    seqnames = as.character(pile$seqnames),
    ranges = IRanges::IRanges(as.integer(pile$pos), width = 1L)
  )
  # [IMPROVE] Match by genomic ranges instead of concatenated "chr:position"
  # strings. This uses genomic semantics directly and avoids delimiter/formatting
  # assumptions in sequence names. Pooled pileup rows are unique by position.
  as.integer(GenomicRanges::findOverlaps(varGr, pileGr, type = "equal",
                                         select = "first", ignore.strand = TRUE))
}


.addWtAlleleSupport <- function(variants, param, which = NULL) {
  if (length(variants) == 0L) return(variants)
  if (is.null(which)) {
    which <- unique(GenomicRanges::GRanges(
      seqnames = as.character(GenomicRanges::seqnames(variants)),
      ranges = IRanges::IRanges(BiocGenerics::start(variants), width = 1L)
    ))
    which <- GenomicRanges::reduce(which, ignore.strand = TRUE)
  }
  # Query only the candidate loci by default. This metadata is additive to the
  # historical candidate model, so avoid rescanning the entire refined peak merely
  # to annotate WT support at a comparatively small set of positions.
  wtPile <- .pooledBasePileup(
    bams = wtFiles(param), genome = param@refGenome, which = which,
    minBaseQuality = minBaseQuality(param),
    minMapQuality = minMapQuality(param),
    maxDepth = .FROZEN_CANDIDATE_MAX_DEPTH
  )

  n <- length(variants)
  wtAltDepth <- integer(n)
  wtRefDepth <- integer(n)
  wtTotalDepth <- integer(n)
  if (nrow(wtPile) > 0L) {
    idx <- .matchVariantToPileup(variants, wtPile)
    found <- !is.na(idx)
    if (any(found)) {
      bases <- c("A", "C", "G", "T")
      counts <- as.matrix(wtPile[idx[found], ..bases])
      altIdx <- match(as.character(VariantAnnotation::alt(variants))[found], bases)
      refIdx <- match(as.character(VariantAnnotation::ref(variants))[found], bases)
      rr <- seq_len(nrow(counts))
      okAlt <- !is.na(altIdx)
      okRef <- !is.na(refIdx)
      tempAlt <- integer(nrow(counts)); tempRef <- integer(nrow(counts))
      if (any(okAlt)) tempAlt[okAlt] <- counts[cbind(rr[okAlt], altIdx[okAlt])]
      if (any(okRef)) tempRef[okRef] <- counts[cbind(rr[okRef], refIdx[okRef])]
      wtAltDepth[found] <- tempAlt
      wtRefDepth[found] <- tempRef
      wtTotalDepth[found] <- wtPile$totalDepth[idx[found]]
    }
  }

  wtAltFreq <- ifelse(wtTotalDepth > 0, wtAltDepth / wtTotalDepth, NA_real_)
  S4Vectors::mcols(variants)$wtRefDepth <- wtRefDepth
  S4Vectors::mcols(variants)$wtAltDepth <- wtAltDepth
  S4Vectors::mcols(variants)$wtTotalDepth <- wtTotalDepth
  S4Vectors::mcols(variants)$wtAltFreq <- wtAltFreq
  variants
}


.isKnownInternalRangeWarning <- function(w) {
  msg <- conditionMessage(w)
  call <- conditionCall(w)
  grepl("^GRanges object contains [0-9]+ out-of-bound ranges located on sequence ", msg) &&
    !is.null(call) && identical(as.character(call[[1L]]), "valid.GenomicRanges.seqinfo")
}

.validateRangesWithinReference <- function(x, refSeqinfo, label) {
  if (!length(x)) return(invisible(TRUE))
  seqn <- as.character(GenomeInfoDb::seqnames(x))
  refNames <- GenomeInfoDb::seqlevels(refSeqinfo)
  lens <- GenomeInfoDb::seqlengths(refSeqinfo)
  idx <- match(seqn, refNames)
  lim <- as.numeric(lens[idx])
  bad <- is.na(idx) | !is.finite(lim) |
    BiocGenerics::start(x) < 1L | BiocGenerics::end(x) > lim
  if (any(bad)) {
    example <- unique(paste0(seqn[bad], ":", BiocGenerics::start(x)[bad],
                             "-", BiocGenerics::end(x)[bad]))
    stop(label, " contains range(s) outside the indexed reference FASTA bounds: ",
         paste(utils::head(example, 5L), collapse = ", "))
  }
  invisible(TRUE)
}

.predictEffects <- function(inputVariants, param, txdb) {
  if (length(inputVariants) == 0L) return(GenomicRanges::GRanges())
  fa <- param@refGenome
  fa_si <- .faSeqinfo(fa)
  tstyle <- .choose_target_style(fa_si)

  # [FIX] Harmonize sequence naming *before* intersecting seqlevels. Intersecting
  # first would incorrectly discard everything when, for example, the GTF uses
  # "1" while the FASTA/candidates use "chr1".
  txdbUse <- txdb
  vars <- inputVariants
  if (!is.na(tstyle)) {
    txdbUse <- tryCatch(.setSeqlevelsStyleFrozen(txdbUse, tstyle),
                        error = function(e) txdbUse)
    vars <- tryCatch(.setSeqlevelsStyleFrozen(vars, tstyle),
                     error = function(e) vars)
  }
  # [FIX] Restrict the common set to all THREE objects. A FASTA/TxDb-wide list
  # can contain chromosomes that are not seqlevels of this peak's VRanges, and
  # passing those extra names to keepSeqlevels(vars, ...) can itself error.
  keep <- Reduce(intersect, list(GenomeInfoDb::seqlevels(txdbUse),
                                 GenomeInfoDb::seqlevels(vars),
                                 GenomeInfoDb::seqlevels(fa_si)))
  if (length(keep) == 0L)
    stop("No common sequence names between the annotation, reference FASTA, and candidate variants")
  txdbUse <- GenomeInfoDb::keepSeqlevels(txdbUse, keep, pruning.mode = "coarse")
  vars <- GenomeInfoDb::keepSeqlevels(vars, keep, pruning.mode = "coarse")
  if (length(vars) == 0L) return(GenomicRanges::GRanges())

  # Validate external genomic coordinates before predictCoding(). This turns
  # wrong-reference annotation problems into clear errors rather than opaque
  # internal transcript-mapping warnings.
  .validateRangesWithinReference(vars, fa_si, "Candidate variants")
  cds <- GenomicFeatures::cds(txdbUse)
  .validateRangesWithinReference(cds, fa_si, "Annotation CDS")

  # [LIMIT] predictCoding() only describes coding consequences. The full SNV
  # candidate table is therefore still written separately for noncoding variants.
  effects <- withCallingHandlers(
    VariantAnnotation::predictCoding(
      query = vars,
      subject = txdbUse,
      seqSource = fa,
      varAllele = Biostrings::DNAStringSet(VariantAnnotation::alt(vars))
    ),
    warning = function(w) {
      if (.isKnownInternalRangeWarning(w)) invokeRestart("muffleWarning")
    }
  )

  # [IMPROVE] Carry mutant/WT evidence into the coding-effect table. predictCoding()
  # identifies each source query with QUERYID but does not guarantee that custom
  # metadata added during candidate calling survives in the result.
  if (length(effects) > 0L && "QUERYID" %in% colnames(S4Vectors::mcols(effects))) {
    qid <- as.integer(S4Vectors::mcols(effects)$QUERYID)
    sourceMeta <- S4Vectors::mcols(vars)
    evidenceCols <- intersect(
      c("mutAltFreq", "mutRefDepth", "mutAltDepth", "mutTotalDepth",
        "wtRefDepth", "wtAltDepth", "wtTotalDepth", "wtAltFreq",
        "deltaAltFreq", "wtEvidenceSufficient", "sharedBackgroundLike"),
      colnames(sourceMeta)
    )
    validQid <- !is.na(qid) & qid >= 1L & qid <= length(vars)
    for (nm in evidenceCols) {
      out <- rep(NA, length(effects))
      out[validQid] <- sourceMeta[[nm]][qid[validQid]]
      S4Vectors::mcols(effects)[[nm]] <- out
    }
  }
  effects
}


.expressionLog2FC <- function(ave_mt, ave_wt, pseudocount = 0) {
  if (pseudocount > 0)
    log2((ave_mt + pseudocount) / (ave_wt + pseudocount))
  else
    log2(ave_mt / ave_wt)
}


.poolMeanCounts <- function(countMat, num_wt, num_mut) {
  if (!is.matrix(countMat)) countMat <- as.matrix(countMat)
  if (num_wt < 1L || num_mut < 1L || ncol(countMat) != num_wt + num_mut)
    stop("countMat columns must equal num_wt + num_mut, with both groups non-empty")
  wt_idx <- seq_len(num_wt)
  mut_idx <- num_wt + seq_len(num_mut)
  # [FIX] Keep group means isolated from one another. This helper exists partly as
  # a regression target for the original off-by-one/chained-mutation expression bug.
  list(wt = rowMeans(countMat[, wt_idx, drop = FALSE]),
       mut = rowMeans(countMat[, mut_idx, drop = FALSE]))
}


.extractQuotedAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]+\"([^\"]+)\"")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
         FUN.VALUE = character(1))
}

.extractEqualsAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]*=[[:space:]]*([^;]+)")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  out <- vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
                FUN.VALUE = character(1))
  trimws(out)
}

.annotationGeneIdentity <- function(attributes) {
  # GTF convention first; fill missing identifiers from standard GFF3 ID/Name.
  gene_id <- .extractQuotedAnnotationAttribute(attributes, "gene_id")
  gff_id <- .extractEqualsAnnotationAttribute(attributes, "ID")
  gene_id[is.na(gene_id) | !nzchar(gene_id)] <- gff_id[is.na(gene_id) | !nzchar(gene_id)]

  gene_name <- .extractQuotedAnnotationAttribute(attributes, "gene_name")
  gff_name <- .extractEqualsAnnotationAttribute(attributes, "Name")
  missing_name <- is.na(gene_name) | !nzchar(gene_name)
  gene_name[missing_name] <- gff_name[missing_name]
  missing_name <- is.na(gene_name) | !nzchar(gene_name)
  gene_name[missing_name] <- gene_id[missing_name]
  list(gene_id = gene_id, gene_name = gene_name)
}

.annotationGeneRanges <- function(param, txdb = NULL) {
  # Fast/fidelity path: the original MMAPPR2 expression summary used explicit
  # GTF `gene` records. Read only the fields it used, preserving those exact
  # coordinates and strand rather than deriving replacement gene models.
  tab <- tryCatch(
    data.table::fread(gtf(param), sep = "\t", header = FALSE,
                      select = c(1L, 3L, 4L, 5L, 7L, 9L), fill = TRUE,
                      quote = "", showProgress = FALSE, data.table = TRUE),
    error = function(e) NULL
  )
  if (!is.null(tab) && ncol(tab) == 6L) {
    data.table::setnames(tab,
                         c("seqnames", "type", "start", "end", "strand", "attributes"))
    tab[, seqnames := as.character(seqnames)]
    tab[, type := as.character(type)]
    tab[, start := suppressWarnings(as.numeric(start))]
    tab[, end := suppressWarnings(as.numeric(end))]
    tab[, strand := as.character(strand)]
    # GTF/GFF commonly uses "." for an unknown/unstranded feature. GRanges
    # represents that state as "*"; normalize it explicitly rather than allowing
    # an atypical annotation to fail during GRanges construction.
    tab[is.na(strand) | !(strand %in% c("+", "-", "*")), strand := "*"]
    tab[, attributes := as.character(attributes)]
    tab <- tab[!is.na(seqnames) & nzchar(trimws(seqnames)) &
               !grepl("^[[:space:]]*#", seqnames) & !is.na(type) & type == "gene" &
               is.finite(start) & is.finite(end) & start >= 1 & end >= start]
    if (nrow(tab)) {
      identity <- .annotationGeneIdentity(tab$attributes)
      df <- data.frame(seqnames = tab$seqnames,
                       start = as.integer(tab$start),
                       end = as.integer(tab$end),
                       strand = tab$strand,
                       gene_id = identity$gene_id,
                       gene_name = identity$gene_name,
                       stringsAsFactors = FALSE)
      return(GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE))
    }
  }

  # Compatibility fallback for GFF3, compressed/atypical annotations, or files
  # that fread cannot parse directly. Prefer explicit gene rows here as well.
  ann <- tryCatch(rtracklayer::import(gtf(param)), error = function(e) NULL)
  if (!is.null(ann) && length(ann)) {
    type <- S4Vectors::mcols(ann)$type
    if (!is.null(type)) {
      type <- as.character(type)
      genes <- ann[!is.na(type) & type == "gene"]
      if (length(genes)) {
        mc <- S4Vectors::mcols(genes)
        ids <- if ("gene_id" %in% colnames(mc)) as.character(mc$gene_id) else rep(NA_character_, length(genes))
        fallbackIds <- if ("ID" %in% colnames(mc)) as.character(mc$ID) else names(genes)
        if (is.null(fallbackIds)) fallbackIds <- rep(NA_character_, length(genes))
        missingId <- is.na(ids) | !nzchar(ids)
        ids[missingId] <- fallbackIds[missingId]

        geneNames <- if ("gene_name" %in% colnames(mc)) as.character(mc$gene_name) else rep(NA_character_, length(genes))
        fallbackNames <- if ("Name" %in% colnames(mc)) as.character(mc$Name) else rep(NA_character_, length(genes))
        missingName <- is.na(geneNames) | !nzchar(geneNames)
        geneNames[missingName] <- fallbackNames[missingName]
        missingName <- is.na(geneNames) | !nzchar(geneNames)
        geneNames[missingName] <- ids[missingName]

        S4Vectors::mcols(genes)$gene_id <- ids
        S4Vectors::mcols(genes)$gene_name <- geneNames
        return(genes)
      }
    }
  }

  # Last-resort hardening for annotations without explicit gene rows. Build a
  # TxDb only if the caller does not already have one; ordinary GTFs therefore
  # avoid this comparatively expensive SQLite construction in the linkage stage.
  ownTxdb <- is.null(txdb)
  if (ownTxdb) {
    txdb <- .buildTxDb(param)
    on.exit(.disconnectTxDb(txdb), add = TRUE)
  }
  genes <- GenomicFeatures::genes(txdb)
  if (!length(genes)) return(genes)
  S4Vectors::mcols(genes)$gene_id <- names(genes)
  S4Vectors::mcols(genes)$gene_name <- names(genes)
  genes
}

.addDiff <- function(peakGRange, param, txdb, annotationGenes = NULL) {
  genes <- if (is.null(annotationGenes)) .annotationGeneRanges(param, txdb) else annotationGenes
  if (!length(genes)) return(GenomicRanges::GRanges())

  # Match annotation/peak naming before overlap tests; original gene-span features
  # are retained instead of replacing them with exon-group features.
  peakUse <- peakGRange
  geneUse <- genes
  peakStyle <- .choose_target_style(.faSeqinfo(param@refGenome))
  if (!is.na(peakStyle)) {
    geneUse <- tryCatch(.setSeqlevelsStyleFrozen(geneUse, peakStyle), error = function(e) geneUse)
    peakUse <- tryCatch(.setSeqlevelsStyleFrozen(peakUse, peakStyle), error = function(e) peakUse)
  }
  common <- intersect(GenomeInfoDb::seqlevels(geneUse), GenomeInfoDb::seqlevels(peakUse))
  if (!length(common)) return(GenomicRanges::GRanges())
  geneUse <- GenomeInfoDb::keepSeqlevels(geneUse, common, pruning.mode = "coarse")
  peakUse <- GenomeInfoDb::keepSeqlevels(peakUse, common, pruning.mode = "coarse")
  geneUse <- IRanges::subsetByOverlaps(geneUse, peakUse, ignore.strand = TRUE)
  if (!length(geneUse)) return(geneUse)

  bamObjects <- c(
    lapply(seq_along(param@wtFiles), function(i) param@wtFiles[[i]]),
    lapply(seq_along(param@mutFiles), function(i) param@mutFiles[[i]])
  )
  readfiles <- do.call(Rsamtools::BamFileList, bamObjects)

  # Preserve original gene-span/peak-window counting semantics. The only counting
  # hardening here is excluding secondary/supplementary/unmapped records so the
  # same physical alignment cannot inflate expression evidence.
  counts <- GenomicAlignments::summarizeOverlaps(
    features = geneUse,
    reads = readfiles,
    mode = "Union",
    ignore.strand = ignoreStrand(param),
    singleEnd = !pairedEnd(param),
    fragments = pairedEnd(param),
    param = Rsamtools::ScanBamParam(which = peakUse,
                                    flag = .primaryMappedScanFlag())
  )
  countMat <- SummarizedExperiment::assays(counts)$counts
  num_wt <- length(param@wtFiles)
  num_mut <- length(param@mutFiles)
  if (ncol(countMat) != num_wt + num_mut)
    stop("Unexpected number of columns returned by summarizeOverlaps")

  # Retain the 9007 correction for the original off-by-one/chained-mutation bug.
  means <- .poolMeanCounts(countMat, num_wt = num_wt, num_mut = num_mut)
  ave_wt <- means$wt
  ave_mt <- means$mut
  # Frozen default is the original raw ratio. Keep 9007's pseudocount as an
  # explicit opt-in only, so a non-default public parameter is never silently ignored.
  pc <- expressionPseudocount(param)
  log2FC <- .expressionLog2FC(ave_mt, ave_wt, pseudocount = pc)

  S4Vectors::mcols(geneUse)$ave_wt <- ave_wt
  S4Vectors::mcols(geneUse)$ave_mt <- ave_mt
  S4Vectors::mcols(geneUse)$log2FC <- round(log2FC, 3)

  geneUse[(abs(S4Vectors::mcols(geneUse)$log2FC) > 1 |
             is.na(S4Vectors::mcols(geneUse)$log2FC)) &
            (S4Vectors::mcols(geneUse)$ave_wt > 10 |
             S4Vectors::mcols(geneUse)$ave_mt > 10)]
}


.scoreVariants <- function(candList, peaks) {
  for (groupName in names(candList)) {
    group <- candList[[groupName]]
    if (!is.list(group)) next
    for (seqname in names(group)) {
      if (is.null(peaks[[seqname]]) || length(group[[seqname]]) == 0L) next
      densityFunc <- peaks[[seqname]]$densityFunction
      positions <- BiocGenerics::start(group[[seqname]]) +
        ((BiocGenerics::width(group[[seqname]]) - 1) / 2)
      densityCol <- vapply(positions, densityFunc, FUN.VALUE = numeric(1))
      S4Vectors::mcols(candList[[groupName]][[seqname]])$peakDensity <- densityCol
      candList[[groupName]][[seqname]] <- .orderVariants(candList[[groupName]][[seqname]])
    }
  }
  candList
}


.orderVariants <- function(candidateGRanges) {
  if (!inherits(candidateGRanges, c("VRanges", "GRanges"))) {
    warning("Invalid data type for candidate sorting: ", class(candidateGRanges)[1])
    return(candidateGRanges)
  }
  if (length(candidateGRanges) == 0L) return(candidateGRanges)

  consequence <- S4Vectors::mcols(candidateGRanges)$CONSEQUENCE
  density <- S4Vectors::mcols(candidateGRanges)$peakDensity
  if (!is.null(consequence)) {
    # [NOTE] This remains a deliberately simple coding-consequence ordering; a
    # full ontology/functional-impact model is a larger scientific redesign.
    impactLevels <- c("synonymous", "nonsynonymous", "frameshift", "nonsense")
    severity <- match(as.character(consequence), impactLevels)
    severity[is.na(severity)] <- 0L
    # Frozen ordering: consequence severity first, then peak density. WT-pool
    # evidence is retained as metadata but does not redefine original ranking.
    orderVec <- order(severity, density, na.last = TRUE, decreasing = TRUE)
  } else {
    orderVec <- order(density, na.last = TRUE, decreasing = TRUE)
  }
  candidateGRanges[orderVec]
}


# [LIMIT] Exact indel reconstruction is intentionally NOT faked here. Rsamtools
# pileup can count insertion/deletion events, but insertion sequence is truncated
# to '+', which is insufficient to build normalized VCF-style alleles or predict
# coding effects correctly. A production indel upgrade should use a mature
# haplotype-aware caller or a thoroughly tested CIGAR-aware implementation.
