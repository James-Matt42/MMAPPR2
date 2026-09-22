# Cross-resource input concordance checks ----------------------------------------
#
# These checks intentionally stay small: no execution state, checkpoints, cached
# provenance, or alternate analysis modes. They run during mmapprParam()
# construction before an existing output directory can be cleared.

# [CHANGE — NEW CROSS-RESOURCE PREFLIGHT LAYER]
# The old implementation validated resources mostly in isolation and could prepare/clear output before
# proving that BAM headers, FASTA sequence dictionaries, and annotation coordinates belonged to the same
# reference build. The new preflight layer establishes that cross-resource contract first.
.bamSeqinfoFrozen <- function(bam) {
    header <- Rsamtools::scanBamHeader(bam)
    entry <- if (length(header) == 1L && is.list(header[[1]]) &&
                 !is.null(header[[1]]$targets)) header[[1]] else header
    targets <- entry$targets
    if (is.null(targets) || !length(targets))
        stop("BAM header contains no sequence dictionary: ", BiocGenerics::path(bam))
    GenomeInfoDb::Seqinfo(seqnames = names(targets),
                          seqlengths = as.numeric(targets))
}

# Apply a naming-style translation while muffling only GenomeInfoDb's harmless
# "more than one best map" warning. Exact-name comparison remains the fallback.
# [CHANGE — CONTROLLED SEQNAME-STYLE CONVERSION]
# The old code performed direct seqlevelsStyle assignments at individual call sites. The new helper
# centralizes those conversions and muffles only the known ambiguous-renaming warning; unrelated
# warnings remain visible.
.setSeqlevelsStyleFrozen <- function(x, style) {
    if (length(style) != 1L || is.na(style) || !nzchar(style)) return(x)
    withCallingHandlers({
        GenomeInfoDb::seqlevelsStyle(x) <- style
        x
    }, warning = function(w) {
        if (grepl("^found more than one best sequence renaming map compatible with seqname style ",
                  conditionMessage(w))) invokeRestart("muffleWarning")
    })
}

# [CHANGE — ANNOTATION-LIKE SEQINFO HARMONIZATION]
# The old implementation had ad-hoc style conversion at individual call sites. The new implementation may
# harmonize annotation-like sequence metadata to the FASTA style when a reliable
# style hint exists, while retaining exact names when style inference is unavailable.
.harmonizeSeqinfoFrozen <- function(si, refSi) {
    targetStyle <- .choose_target_style(refSi)
    if (!is.na(targetStyle)) {
        candidate <- tryCatch(.setSeqlevelsStyleFrozen(si, targetStyle),
                              error = function(e) NULL)
        if (!is.null(candidate)) si <- candidate
    }
    si
}

# [CHANGE — EXACT SHARED BAM/FASTA CONTIG MATCHING]
# The old implementation did not compare BAM headers with the reference FASTA before analysis. Because
# BAM contig names cannot be renamed for Rsamtools queries, the new preflight requires at least one exact
# shared contig name and matching lengths for shared contigs, and rejects style-only matches such as
# `1` versus `chr1`.
.compareBamSeqinfoFrozen <- function(label, si, refSi) {
    # BAM queries are issued with reference/annotation GRanges. Unlike annotation
    # ranges, a BAM header cannot be renamed in memory before scanBam()/pileup().
    # Therefore a chr1-vs-1 mismatch must *not* be accepted merely because
    # GenomeInfoDb can recognize the two naming styles as equivalent.
    refNames <- GenomeInfoDb::seqlevels(refSi)
    bamNames <- GenomeInfoDb::seqlevels(si)
    common <- intersect(refNames, bamNames)

    styled <- .harmonizeSeqinfoFrozen(si, refSi)
    styledCommon <- intersect(refNames, GenomeInfoDb::seqlevels(styled))
    styleOnly <- setdiff(styledCommon, common)
    if (length(styleOnly)) {
        stop(label,
             " uses sequence naming that differs from the reference FASTA ",
             "(for example ", paste(utils::head(styleOnly, 4L), collapse = ", "),
             "). BAM and FASTA sequence names must match exactly for frozen MMAPPR2 queries.")
    }
    if (!length(common))
        stop(label, " shares no exact sequence names with the reference FASTA")

    refLen <- GenomeInfoDb::seqlengths(refSi)[common]
    bamLen <- GenomeInfoDb::seqlengths(si)[common]
    comparable <- is.finite(refLen) & is.finite(bamLen)
    conflicts <- common[comparable & as.numeric(refLen[comparable]) != as.numeric(bamLen[comparable])]
    if (length(conflicts)) {
        stop(label, " has sequence-length conflict(s) with the reference FASTA: ",
             paste(utils::head(conflicts, 8L), collapse = ", "))
    }
    common
}

# Read only the three annotation fields needed for reference-concordance checks.
# fread is much lighter than constructing a TxDb. If a compressed/atypical GFF
# cannot be read this way, rtracklayer is used as the compatibility fallback.
# [CHANGE — LIGHTWEIGHT ANNOTATION COORDINATE SCAN]
# The old implementation had no annotation-coordinate preflight. The new code reads only seqname/start/end
# for validation, falling back to rtracklayer when needed, so malformed or out-of-reference coordinates
# are detected without constructing a TxDb solely for input checking.
.annotationObservedRanges <- function(path) {
    tab <- tryCatch(
        data.table::fread(path, sep = "\t", header = FALSE,
                          select = c(1L, 4L, 5L), fill = TRUE, quote = "",
                          comment.char = "#", showProgress = FALSE, data.table = TRUE),
        error = function(e) NULL
    )
    if (!is.null(tab) && ncol(tab) == 3L) {
        data.table::setnames(tab, c("seqname", "start", "end"))
        # fread may infer character columns for malformed/comment-heavy inputs.
        # Coerce explicitly before is.finite() so preflight reports concordance
        # problems instead of failing with a type error.
        tab[, seqname := as.character(seqname)]
        tab[, start := suppressWarnings(as.numeric(start))]
        tab[, end := suppressWarnings(as.numeric(end))]
        dataRows <- !is.na(tab$seqname) & nzchar(trimws(tab$seqname)) & !grepl("^[[:space:]]*#", tab$seqname)
        if (any(dataRows & (!is.finite(tab$start) | !is.finite(tab$end))))
            stop("Annotation contains non-numeric start/end coordinates: ", path)
        if (any(dataRows & (tab$start < 1 | tab$end < tab$start)))
            stop("Annotation contains invalid genomic coordinate bounds: ", path)
        tab <- tab[dataRows & is.finite(start) & is.finite(end)]
        if (nrow(tab)) {
            maxEnd <- tab[, .(max_end = max(end)), by = seqname]
            return(maxEnd)
        }
    }

    gr <- rtracklayer::import(path)
    if (!length(gr)) stop("Annotation contains no genomic features: ", path)
    data.table::data.table(
        seqname = as.character(GenomeInfoDb::seqnames(gr)),
        end = as.numeric(BiocGenerics::end(gr))
    )[, .(max_end = max(end)), by = seqname]
}

# [CHANGE — ANNOTATION/FASTA BOUNDS VALIDATION]
# After safe seqname harmonization, the new implementation verifies that the maximum observed
# annotation coordinate on every shared sequence fits inside the indexed FASTA.
# The old implementation did not perform this build-level bounds check before analysis.
.checkAnnotationAgainstReferenceFrozen <- function(path, refSi) {
    obs <- .annotationObservedRanges(path)
    annSi <- GenomeInfoDb::Seqinfo(seqnames = obs$seqname,
                                   seqlengths = as.numeric(obs$max_end))
    annSi <- .harmonizeSeqinfoFrozen(annSi, refSi)

    # Harmonization may rename the Seqinfo, so reconstruct observed maxima by the
    # translated names in the same order.
    annNames <- GenomeInfoDb::seqlevels(annSi)
    annMax <- GenomeInfoDb::seqlengths(annSi)
    refNames <- GenomeInfoDb::seqlevels(refSi)
    common <- intersect(annNames, refNames)
    if (!length(common))
        stop("Annotation shares no sequence names with the reference FASTA after harmonization")

    refLen <- GenomeInfoDb::seqlengths(refSi)[common]
    maxEnd <- annMax[common]
    bad <- common[is.finite(refLen) & is.finite(maxEnd) &
                  as.numeric(maxEnd) > as.numeric(refLen)]
    if (length(bad)) {
        detail <- paste0(bad, " (annotation end ", as.numeric(maxEnd[bad]),
                         " > FASTA length ", as.numeric(refLen[bad]), ")")
        stop("Annotation contains coordinates outside the indexed reference FASTA: ",
             paste(utils::head(detail, 8L), collapse = ", "))
    }
    common
}

# [CHANGE — ONE BUILD-CONCORDANCE GATE]
# The old constructor/setters had no single build-concordance gate. The new implementation routes
# construction and resource replacement through one check requiring every BAM to agree with the FASTA,
# WT and mutant pools to share reference sequence(s), and annotation coordinates to fit that FASTA.
.preflightInputResourcesFrozen <- function(wt, mut, ref, annotationPath) {
    refSi <- .faSeqinfo(ref)
    if (!length(refSi)) stop("Reference FASTA contains no indexed sequences")

    wtSets <- lapply(seq_along(wt), function(i) {
        .compareBamSeqinfoFrozen(paste0("WT BAM ", i), .bamSeqinfoFrozen(wt[[i]]), refSi)
    })
    mutSets <- lapply(seq_along(mut), function(i) {
        .compareBamSeqinfoFrozen(paste0("Mutant BAM ", i), .bamSeqinfoFrozen(mut[[i]]), refSi)
    })
    wtUnion <- Reduce(union, wtSets)
    mutUnion <- Reduce(union, mutSets)
    if (!length(intersect(wtUnion, mutUnion)))
        stop("WT and mutant BAM pools have no common reference sequence after harmonization")

    .checkAnnotationAgainstReferenceFrozen(annotationPath, refSi)
    invisible(TRUE)
}
