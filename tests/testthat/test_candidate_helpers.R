test_that("mutant candidate threshold preserves strict original >80 percent boundary", {
  keep <- MMAPPR2:::.candidateMutantKeep(
    totalDepth = c(20, 20, 19, 2, 20, 0),
    altDepth = c(16, 17, 16, 2, 1, 0),
    minDepth = 1, minAltDepth = 2, minAltFreq = .80
  )
  # 2/2 is intentionally valid: the historical caller's meaningful depth rule
  # was two ALT-supporting reads, not the linkage-stage 20-read depth floor.
  expect_identical(keep, c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("WT candidate filters can be disabled or applied independently", {
  wt <- c(.2, .9, NA)
  delta <- c(.7, .05, NA)
  expect_true(all(MMAPPR2:::.candidateWtKeep(wt, delta, 1, 0)))
  expect_identical(
    MMAPPR2:::.candidateWtKeep(wt, delta, .5, 0),
    c(TRUE, FALSE, FALSE)
  )
  expect_identical(
    MMAPPR2:::.candidateWtKeep(wt, delta, 1, .2),
    c(TRUE, FALSE, FALSE)
  )
})

test_that("base pileup conversion does not mutate input and preserves depths", {
  pile <- data.table::data.table(
    seqnames = c("chr1", "chr1"), pos = c(10L, 20L),
    A = c(2L, 0L), C = c(18L, 0L), G = c(0L, 20L), T = c(0L, 0L),
    totalDepth = c(20L, 20L), ref = c("A", "G")
  )
  before <- data.table::copy(pile)
  vr <- MMAPPR2:::.basePileupToVRanges(pile)
  expect_identical(pile, before)
  expect_false(".row_id" %in% names(pile))
  expect_length(vr, 1)
  expect_equal(as.character(VariantAnnotation::ref(vr)), "A")
  expect_equal(as.character(VariantAnnotation::alt(vr)), "C")
  expect_equal(as.integer(VariantAnnotation::refDepth(vr)), 2L)
  expect_equal(as.integer(VariantAnnotation::altDepth(vr)), 18L)
  expect_equal(as.integer(VariantAnnotation::totalDepth(vr)), 20L)
})

test_that("variant-to-pileup matching uses exact genomic coordinates", {
  vars <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 30), width = 1))
  pile <- data.table::data.table(seqnames = c("chr1", "chr1"), pos = c(30L, 10L))
  idx <- MMAPPR2:::.matchVariantToPileup(vars, pile)
  expect_identical(idx, c(2L, 1L))
})

test_that("expression group means cannot leak samples across groups", {
  mat <- matrix(c(10, 20, 100, 200,
                  30, 50, 300, 500), nrow = 2, byrow = TRUE)
  means <- MMAPPR2:::.poolMeanCounts(mat, num_wt = 2L, num_mut = 2L)
  expect_equal(means$wt, c(15, 40))
  expect_equal(means$mut, c(150, 400))
})

test_that("WT background metadata does not change frozen candidate ordering", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 20), width = 1))
  S4Vectors::mcols(gr)$peakDensity <- c(0.9, 0.5)
  S4Vectors::mcols(gr)$sharedBackgroundLike <- c(TRUE, FALSE)
  ordered <- MMAPPR2:::.orderVariants(gr)
  expect_equal(BiocGenerics::start(ordered)[1], 10)
})

test_that("reference range validation rejects coordinates beyond FASTA", {
  si <- GenomeInfoDb::Seqinfo("chr1", seqlengths = 100)
  ok <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 100))
  bad <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 101))
  expect_invisible(MMAPPR2:::.validateRangesWithinReference(ok, si, "test"))
  expect_error(MMAPPR2:::.validateRangesWithinReference(bad, si, "test"),
               "outside the indexed reference FASTA bounds")
})

test_that("expression fold-change helper preserves raw default and pseudocount opt-in", {
  raw <- MMAPPR2:::.expressionLog2FC(c(20, 10, 0, 0), c(10, 0, 10, 0), 0)
  expect_equal(raw[1], 1)
  expect_true(is.infinite(raw[2]) && raw[2] > 0)
  expect_true(is.infinite(raw[3]) && raw[3] < 0)
  expect_true(is.nan(raw[4]))

  opt <- MMAPPR2:::.expressionLog2FC(c(20, 10, 0, 0), c(10, 0, 10, 0), 0.5)
  expect_true(all(is.finite(opt)))
  expect_equal(opt[4], 0)
})
