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

# [CHANGE — CANDIDATE PIPELINE REWRITE]
# The old implementation merged multiple mutant BAMs into a fixed temporary BAM, passed pooled
# variants through VariantTools::callVariants(), had no WT-pool evidence in final candidate calls,
# emitted only coding-effect rows in candidate tables, and sliced WT/mutant expression columns
# incorrectly. The new implementation pools base counts in memory, applies explicit SNV thresholds,
# records WT/delta-AF evidence, adds a full-SNV candidate table, and fixes expression grouping while
# retaining the old >0.80/two-ALT-read defaults.
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
  # [CHANGE — LAZY TXDB CONSTRUCTION/LIFECYCLE]
  # The old coding-effect path built a TxDb inside `.predictEffects()` for each peak that was processed.
  # The new implementation discovers SNVs first, builds at most one TxDb only when coding annotation is
  # needed, reuses it across peaks, and closes its transient SQLite connection deterministically.
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


# [CHANGE — EXPLICIT, USER-VISIBLE SNV CALLING]
# The old implementation used VariantTools::callVariants() after pileup, adding a hidden
# likelihood/read-count filter before MMAPPR2's strict >80% ALT rule. The new implementation calls
# SNVs directly from pooled A/C/G/T depths so every hard candidate criterion is an
# explicit MmapprParam setting and the evidence used is retained in output metadata.
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
  # [CHANGE — CANDIDATE PILEUP CAP MADE EXPLICIT]
  # The old candidate path used Rsamtools' implicit max_depth=250 after merging multiple mutant BAMs.
  # The new path names 250 explicitly and applies it to each mutant BAM before in-memory pooling, while
  # keeping this candidate-stage limit independent from the linkage-stage maxPileupDepth setting.
  mutPile <- .pooledBasePileup(
    bams = mutFiles(param), genome = fa, which = inputRange,
    minBaseQuality = minBaseQuality(param),
    minMapQuality = minMapQuality(param),
    # The old candidate path relied on Rsamtools' implicit 250-read pileup cap.
    # The new candidate path keeps an explicit 250-read cap per BAM and does not
    # reuse the separately configurable linkage maxPileupDepth setting.
    maxDepth = .FROZEN_CANDIDATE_MAX_DEPTH
  )
  if (nrow(mutPile) == 0L) return(VariantAnnotation::VRanges())

  pile <- .basePileupToVRanges(mutPile, sampleName = "pooled_mutant")
  if (length(pile) == 0L) return(VariantAnnotation::VRanges())

  # [FIX/IMPROVE] The old implementation passed pooled tallies through
  # VariantTools::callVariants(), which adds a binomial calling filter and a
  # two-ALT-read minimum before the separate >80% ALT rule. The new implementation
  # makes the effective hard criteria explicit: candidateMinAltDepth represents the
  # two-read minimum and candidateMinAltFreq retains the strict >80% boundary.
  mutAF <- VariantAnnotation::altDepth(pile) / VariantAnnotation::totalDepth(pile)
  # [CHANGE — MUTANT CALL RULES EXPOSED AS PARAMETERS]
  # The new implementation makes the old strict ALT-frequency rule (`> 0.80`) and two-ALT-read
  # requirement explicit. candidateMinDepth is new but defaults to 1, so it is
  # effectively inert under the default two-ALT-read rule and does not add a new
  # default scientific filter.
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
# [CHANGE — PURE/TESTABLE MUTANT CANDIDATE CRITERIA]
# The old candidate thresholds were split between VariantTools behavior and an
# inline >0.8 expression. The new implementation centralizes the explicit depth/ALT-depth/ALT-frequency
# checks in a pure helper; the old strict `>` frequency boundary is retained.
.candidateMutantKeep <- function(totalDepth, altDepth, minDepth, minAltDepth, minAltFreq) {
  af <- ifelse(totalDepth > 0, altDepth / totalDepth, NA_real_)
  is.finite(totalDepth) & is.finite(altDepth) & is.finite(af) &
    totalDepth >= 0 & altDepth >= 0 & altDepth <= totalDepth &
    totalDepth >= minDepth & altDepth >= minAltDepth & af > minAltFreq
}

# [CHANGE — OPTIONAL WT AND DELTA-AF FILTERS]
# The new implementation adds optional hard filters on WT ALT frequency and mutant-minus-WT ALT
# frequency. Their defaults (WT maximum 1; delta minimum 0) disable filtering, so
# old candidate inclusion is not narrowed unless the user opts in.
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
# [CHANGE — IN-MEMORY REPLICATE POOLING]
# The old implementation merged multiple mutant BAMs on disk and then piled the merged file. The new
# implementation piles each BAM separately and sums A/C/G/T counts in memory, eliminating the fixed
# temporary BAM. Counts are pooled before ALT selection, so reads with the reference base still
# contribute to total depth; for multiple BAMs, the explicit 250-read cap now applies per BAM.
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


# [CHANGE — EXPLICIT REF/ALT DEPTH CONSTRUCTION]
# The old path created per-ALT VRanges and then passed them through VariantTools::callVariants(). The new
# path creates one VRanges row for each observed non-reference A/C/G/T base directly from the pooled
# counts, retains ref/ALT/total depth explicitly, and applies the candidate filters separately.
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


# [CHANGE — GENOMIC-OVERLAP EVIDENCE MATCHING]
# The old candidate model had no WT-pileup evidence to attach to variants. The new WT-evidence path
# matches pooled base-count rows to candidate loci with exact GRanges overlaps, using genomic coordinates
# directly rather than inventing a text key for the newly added join.
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
  # [IMPROVE] The old candidate model had no WT-evidence join. The new WT-evidence
  # join matches exact genomic ranges directly, avoiding text-key formatting assumptions
  # in the newly added path. Pooled pileup rows are unique by position.
  as.integer(GenomicRanges::findOverlaps(varGr, pileGr, type = "equal",
                                         select = "first", ignore.strand = TRUE))
}


# [CHANGE — WT-POOL EVIDENCE ADDED TO CANDIDATES]
# The old final candidate call considered mutant reads against the reference
# only. The new implementation measures WT ref/ALT/total depth, WT ALT frequency, mutant ALT frequency,
# and delta AF for every candidate. By default this evidence is reported but does
# not exclude candidates; optional WT/delta filters can be enabled explicitly.
.addWtAlleleSupport <- function(variants, param, which = NULL) {
  if (length(variants) == 0L) return(variants)
  if (is.null(which)) {
    which <- unique(GenomicRanges::GRanges(
      seqnames = as.character(GenomicRanges::seqnames(variants)),
      ranges = IRanges::IRanges(BiocGenerics::start(variants), width = 1L)
    ))
    which <- GenomicRanges::reduce(which, ignore.strand = TRUE)
  }
  # Query only the candidate loci by default. WT support is new metadata layered onto
  # the old candidate model, so there is no reason to rescan the entire refined peak
  # merely to annotate a comparatively small set of candidate positions.
  # [CHANGE — WT PILEUP LIMITED TO CANDIDATE LOCI]
  # The old implementation did not pile up the WT pool during final candidate calling. The new WT
  # evidence pass queries only the candidate coordinates—not every base in the refined interval—and
  # uses the explicit candidate-stage 250-read-per-BAM cap.
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


# [CHANGE — NARROW WARNING HANDLING]
# The old predictCoding() call allowed every warning to propagate. The new implementation first
# validates reference bounds, then muffles only one known internal VariantAnnotation range warning;
# all other warnings still propagate.
.isKnownInternalRangeWarning <- function(w) {
  msg <- conditionMessage(w)
  call <- conditionCall(w)
  grepl("^GRanges object contains [0-9]+ out-of-bound ranges located on sequence ", msg) &&
    !is.null(call) && identical(as.character(call[[1L]]), "valid.GenomicRanges.seqinfo")
}

# [CHANGE — PRE-predictCoding REFERENCE-BOUNDS CHECK]
# The old implementation passed candidate/CDS ranges to predictCoding() without an explicit FASTA-bounds
# check. The new implementation validates them against indexed FASTA lengths first, turning wrong-build
# or out-of-bounds coordinates into a direct actionable error.
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

# [CHANGE — CODING-ANNOTATION HARDENING]
# The old path built a TxDb inside each effect call, used direct style assignment, and returned only
# predictCoding() consequences for this stage. The new path reuses a managed TxDb, harmonizes
# variant/TxDb/FASTA naming before intersecting seqlevels, validates FASTA/CDS bounds, keeps the full
# candidate-SNV set separately, and carries mutant/WT evidence into coding-effect rows.
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


# [CHANGE — EXPRESSION RATIO MADE EXPLICIT/CONFIGURABLE]
# The old implementation used raw log2(mean mutant / mean WT), including Inf/NaN at zero
# counts. The new implementation preserves that behavior by default (pseudocount=0) but factors it
# into a testable helper and permits an explicit positive pseudocount when desired.
.expressionLog2FC <- function(ave_mt, ave_wt, pseudocount = 0) {
  if (pseudocount > 0)
    log2((ave_mt + pseudocount) / (ave_wt + pseudocount))
  else
    log2(ave_mt / ave_wt)
}


# [CHANGE — WT/MUTANT EXPRESSION INDEXING BUG FIX]
# The old WT slice included the first mutant column and the mutant slice could
# overlap the WT boundary. The new implementation partitions columns exactly as WT 1..num_wt and
# mutant (num_wt+1)..(num_wt+num_mut) before calculating group means.
.poolMeanCounts <- function(countMat, num_wt, num_mut) {
  if (!is.matrix(countMat)) countMat <- as.matrix(countMat)
  if (num_wt < 1L || num_mut < 1L || ncol(countMat) != num_wt + num_mut)
    stop("countMat columns must equal num_wt + num_mut, with both groups non-empty")
  wt_idx <- seq_len(num_wt)
  mut_idx <- num_wt + seq_len(num_mut)
  # [FIX] Keep group means isolated from one another. This helper exists partly as
  # a regression target for the old implementation's off-by-one/chained-mutation expression bug.
  list(wt = rowMeans(countMat[, wt_idx, drop = FALSE]),
       mut = rowMeans(countMat[, mut_idx, drop = FALSE]))
}


# [CHANGE — ROBUST GTF ATTRIBUTE EXTRACTION]
# The old GTF path extracted gene_id/gene_name with rigid greedy substitutions over the full attribute
# string. The new implementation parses quoted attributes by key, preserving the intended values without
# depending on attribute order or unrelated trailing fields.
.extractQuotedAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]+\"([^\"]+)\"")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
         FUN.VALUE = character(1))
}

# [CHANGE — GFF3 ATTRIBUTE SUPPORT]
# The old expression path assumed GTF-style quoted attributes. The new implementation also
# understands GFF3 `key=value` attributes so ID/Name can be preserved.
.extractEqualsAnnotationAttribute <- function(x, key) {
  pattern <- paste0("(?:^|;)[[:space:]]*", key, "[[:space:]]*=[[:space:]]*([^;]+)")
  hit <- regexec(pattern, x, perl = TRUE)
  parts <- regmatches(x, hit)
  out <- vapply(parts, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_,
                FUN.VALUE = character(1))
  trimws(out)
}

# [CHANGE — PRESERVE ANNOTATION GENE IDENTIFIERS]
# The old direct GTF parser preserved gene_id/gene_name but had no equivalent alternate-format path. The
# new implementation centralizes identity extraction so GTF gene_id/gene_name and GFF3 ID/Name survive
# whichever supported annotation reader supplies the gene ranges.
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

# [CHANGE — GENE-RANGE PARSING WITH FALLBACKS]
# The old expression/linkage code directly parsed explicit GTF `gene` rows.
# The new implementation preserves those exact spans/strand/identifiers as the preferred path, adds
# normalization of unknown strand to `*`, supports GFF3/compressed/atypical input
# through rtracklayer, and uses TxDb-derived genes only as a last resort.
.annotationGeneRanges <- function(param, txdb = NULL) {
  # Fast/fidelity path: the old expression summary used explicit
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

  # New fallback for GFF3, compressed/atypical annotations, or files that fread
  # cannot parse directly. Explicit gene rows remain preferred here as well.
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

# [CHANGE — DESCRIPTIVE EXPRESSION SUMMARY REPAIR]
# The old expression summary counted annotation gene spans within the refined peak but sliced WT/mutant
# count columns incorrectly and left counting assumptions implicit. The new implementation preserves the
# gene-span/peak-window design, harmonizes seqnames, applies the primary-read filter, exposes paired-end
# and strand controls, fixes group indexing, and keeps raw log2(mutant/WT) with zero pseudocount by
# default. This remains a descriptive raw-count summary, not formal differential-expression inference.
.addDiff <- function(peakGRange, param, txdb, annotationGenes = NULL) {
  genes <- if (is.null(annotationGenes)) .annotationGeneRanges(param, txdb) else annotationGenes
  if (!length(genes)) return(GenomicRanges::GRanges())

  # Match annotation/peak naming before overlap tests; old gene-span features
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

  # Preserve old gene-span/peak-window counting semantics. The only counting
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

  # Correct the old implementation's off-by-one/chained-mutation group-indexing bug.
  # [CHANGE — CORRECT GROUP MEANS USED HERE]
  # The old code used overlapping/off-by-one column slices at this point. The new code uses the exact WT
  # and mutant partitions returned by .poolMeanCounts().
  means <- .poolMeanCounts(countMat, num_wt = num_wt, num_mut = num_mut)
  ave_wt <- means$wt
  ave_mt <- means$mut
  # The default remains the old implementation's raw ratio. A positive pseudocount is an
  # explicit opt-in only, so a non-default public parameter is never silently ignored.
  # [CHANGE — OPTIONAL PSEUDOCOUNT WITHOUT CHANGING DEFAULT]
  # The old implementation always used the raw mutant/WT ratio with no pseudocount. The new implementation
  # keeps zero as the default but exposes a positive pseudocount as an explicit recorded opt-in.
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


# [CHANGE — DEFENSIVE/EXPLICIT CANDIDATE ORDERING]
# The old invalid-type branch could fail while trying to log through an undefined object, and empty/NA
# consequence cases were not handled defensively. The new implementation returns safely for invalid or
# empty inputs, handles missing consequence values, and preserves the old ranking dimensions:
# coding-consequence severity first, then peak density; WT-background evidence remains metadata.
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
    # Preserve the old ordering: consequence severity first, then peak density.
    # WT-pool evidence is retained as metadata and does not redefine the ranking.
    orderVec <- order(severity, density, na.last = TRUE, decreasing = TRUE)
  } else {
    orderVec <- order(density, na.last = TRUE, decreasing = TRUE)
  }
  candidateGRanges[orderVec]
}


# [CHANGE — SNV-ONLY SCOPE MADE EXPLICIT]
# The old candidate path also operated on base-only substitution tallies. The new direct caller keeps
# that scope explicit by constructing exact A/C/G/T SNVs only; it does not turn pileup insertion/deletion
# symbols into incomplete VCF alleles that cannot be normalized or annotated reliably.
# [LIMIT] Exact indel reconstruction is intentionally NOT faked here. Rsamtools
# pileup can count insertion/deletion events, but insertion sequence is truncated
# to '+', which is insufficient to build normalized VCF-style alleles or predict
# coding effects correctly. A production indel upgrade should use a mature
# haplotype-aware caller or a thoroughly tested CIGAR-aware implementation.
