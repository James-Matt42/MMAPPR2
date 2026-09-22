# MMAPPR2 0.99.24.9008 — frozen/fixed maintenance branch

This tree starts from the fully validated 0.99.24.9007 source and reconciles it against the supplied 0.99.24 implementation. The goal is deliberately narrow: keep MMAPPR2 recognizable and scientifically faithful while repairing identifiable bugs and adding low-risk input/reference hardening.

Read these first:

- `FROZEN_BEHAVIOR.md` — scientific/default behavior that is intentionally frozen.
- `FIXES_FROM_ORIGINAL.md` — accepted bug fixes and hardening.
- `9007_RECONCILIATION.md` — 9007 changes retained, reverted, or neutralized.
- `9018_BACKPORT_DECISIONS.md` — selected 9018 ideas accepted or rejected.
- `VALIDATION_STATUS_9008.md` — what has and has not been executed on this modified tree.
- `FROZEN_AUDIT_MATRIX.md` — line-item original/9007/9008 decisions and regression gates.
- `SOURCE_BASELINES.md` — supplied source/archive provenance and hashes.
- `tests/testthat/` — active regression and synthetic integration tests.
- `validation/run_frozen_validation.sh` — one-command R-enabled validation gate.

## Runtime baseline

The package targets R 4.4 or newer / Bioconductor 3.19 or newer and uses `txdbmaker`, the maintained home of `makeTxDbFromGFF()`.

## Important compatibility choices

The default path remains one sequential MMAPPR2 analysis. There is no 9018-style compatibility/modern mode, checkpoint runner, provenance state machine, or alternate statistical model.

Candidate calling preserves the original strict `>0.80` mutant ALT-frequency rule and two-ALT-read requirement. The 9007-added total candidate depth floor is neutralized by default (`candidateMinDepth=1`), and candidate pileup uses the historical 250-read Rsamtools cap while retaining MMAPPR2's configured base/map-quality thresholds.

The descriptive expression screen again uses explicit annotation gene spans, the refined peak as the BAM query window, strand-aware counting by default, and raw `log2(mutant/wildtype)`. A positive `expressionPseudocount` remains an explicit 9007-style opt-in, while the frozen default is zero. The serious 9007-fixed WT/mutant indexing bug remains fixed.

Duplicate status is intentionally left unspecified in the shared read filter, matching original MMAPPR2 and 9007: duplicate-marked primary reads remain eligible in linkage, candidate, and expression counting unless removed upstream.

Input preflight requires BAM sequence names to match the reference FASTA exactly, because regional BAM queries are issued with reference-coordinate names. Annotation sequence naming may still be harmonized in memory. Missing, empty, or stale BAM/FASTA sidecar indexes are rebuilt before output-folder preparation.

## Validation status

The inherited 9007 tree was previously fully validated under R 4.5.0. **Those logs are historical and do not constitute validation of 9008.** This execution environment does not contain R/Rscript, so 9008 currently has source-level/static validation and an expanded executable regression suite, but still requires the R-enabled gate in `validation/run_frozen_validation.sh` before release.
