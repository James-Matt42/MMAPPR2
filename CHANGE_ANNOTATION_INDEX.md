# MMAPPR2 — old-to-new change annotation index

## Comparison and annotation rule

This review describes the differences between the **old implementation** and the **new implementation**. The old implementation is the clean source-repository baseline used for comparison; the new implementation is the annotated package in this archive.

Every local explanatory annotation uses the form:

`[CHANGE — ...]`

The comments explain what the old functionality did and what the new functionality does. They intentionally avoid package revision numbers and intermediate-version history. Runtime R expressions are unchanged; only comments/documentation are added or reworded.

## Package/API surface

### `R/MMAPPR2.R`

- **R-check/NSE declarations:** the new implementation adds `utils::globalVariables()` declarations for the data.table symbols used by the expanded linkage, candidate, annotation, and preflight code. This is static-check metadata only.
- **Dependency/API modernization:** the old implementation VariantTools-oriented import surface is replaced by the current APIs actually used by the new implementation, including `txdbmaker` and additional `VariantAnnotation` accessors/types.

### `R/all_generics.R`

- **Public analysis-setting getters:** the new implementation exposes formerly hard-coded linkage, candidate, peak, RNG, and expression choices as S4 generics.
- **Replacement generics:** matching setters are added so those settings can be changed through validated methods rather than direct slot editing.

### `NAMESPACE` / `DESCRIPTION`

- The namespace exports the added parameter accessors/replacements and imports the current dependency set used by the runtime implementation.
- the minimum R version and Bioconductor-era dependencies are updated, the obsolete VariantTools dependency is removed, and `preflight.R` is added to collation.
- Package metadata/documentation now describes local coding annotation and WT evidence rather than the older external-caller/VEP-oriented description.

## Parameter construction, validation, and resource consistency — `R/param.R`

The new implementation makes parameter/resource behavior much more explicit than the old implementation.

- **Candidate pileup cap is named explicitly:** the old implementation candidate path's effective 250-read/BAM cap is preserved as `.FROZEN_CANDIDATE_MAX_DEPTH`, separate from the new linkage `maxPileupDepth` control.
- **New recorded settings:** `maxPileupDepth`, `candidateMinDepth`, `candidateMinAltDepth`, `candidateMinAltFreq`, `candidateMaxWtAltFreq`, `candidateMinDeltaAF`, `peakCutoffSd`, `peakCutoffMethod`, `peakIntervalMethod`, `peakResampleIterations`, `randomSeed`, `pairedEnd`, `ignoreStrand`, and `expressionPseudocount` are stored in `MmapprParam`.
- **Compatibility-preserving defaults:** the default candidate rule still has the historical 250-read/BAM candidate cap, two ALT reads, and strict `>0.80` mutant ALT frequency; expression remains strand-aware by default and uses raw `log2(mutant/wildtype)` with zero pseudocount.
- **Constructor rewrite:** documented BAM input types are normalized; scalar/cross-parameter validation happens before filesystem work; BAM/FASTA indexes are prepared and validated; annotation handling no longer depends on destructive shell decompression; overwrite is explicit; and build preflight occurs before output deletion.
- **Input-before-output ordering:** corrupt sequence inputs, stale indexes, incompatible builds, or out-of-bounds annotations fail before an existing output directory can be cleared.
- **Central scalar validation:** invalid types/ranges and impossible settings are rejected consistently, including candidate depth capacity and open upper bound for the strict ALT-frequency threshold.
- **Whole-object validity:** resource shapes, scalar settings, and `refFasta`/cached `FaFile` consistency are checked in one S4 validity path.
- **FASTA empty-path fix:** zero-length/NA/non-scalar FASTA paths are explicitly rejected.
- **Deep BAM/index validation:** construction and resource replacement can open BAM headers/indexes rather than checking filenames alone.
- **BAM type normalization:** character paths, `BamFile`, and `BamFileList` are normalized before path/index handling.
- **Index freshness:** an index must exist, be nonempty, and not be older than its BAM/FASTA; otherwise the new implementation rebuilds it.
- **Safer BAM index discovery:** sibling-index names are checked explicitly and `BamFileList` is built from real `BamFile` objects.
- **Custom assemblies:** failure to infer a canonical seqlevel style becomes a nonfatal “no hint”; exact-name matching remains possible.
- **FASTA Seqinfo:** names/lengths come directly from the indexed `FaFile`.
- **TxDb lifecycle:** transient SQLite-backed TxDb connections are closed explicitly.
- **TxDb warnings:** only the specifically recognized non-actionable genome-version warning is muffled; other annotation warnings propagate.
- **Annotation parsing:** TxDb construction uses Bioconductor GTF/GFF parsing instead of shell text processing.
- **New getters/setters:** all recorded settings are available through the S4 API.
- **Setter consistency:** replacement methods revalidate the object instead of leaving invalid combinations latent.
- **Resource setter preflight:** changing BAMs, FASTA, or annotation re-runs build-concordance checks.
- **Mutant setter bug:** `mutFiles<-` now validates the mutant replacement itself instead of accidentally re-checking WT files.
- **FASTA replacement synchronization:** `refFasta<-` updates both the visible path and cached `FaFile`, refreshes `.fai` when necessary, and preflights the new reference.
- **Annotation replacement:** `gtf<-` now checks the replacement against the existing reference contract.
- **Shared scalar setter path:** scalar replacement methods use one validation helper so rules cannot diverge between setters.

## New cross-resource preflight — `R/preflight.R`

The old implementation had no dedicated cross-resource preflight module; the new implementation adds one.

- **BAM sequence dictionaries** are read from headers for reference-build checks.
- **Controlled seqname-style conversion** centralizes GenomeInfoDb translation and narrowly handles only the known ambiguous-map warning.
- **Annotation-like Seqinfo harmonization** is attempted only when a reliable style hint exists.
- **Exact BAM/FASTA contract:** BAM contig names and shared lengths must agree exactly because BAM query names cannot be renamed in memory.
- **Lightweight annotation scan:** seqname/start/end are read without building a TxDb solely for preflight.
- **Annotation bounds:** maximum observed coordinates must fit the indexed FASTA after safe name harmonization.
- **One build-level gate:** both BAM pools, FASTA, and annotation must be mutually compatible before construction or resource replacement succeeds.

## Linkage distance stage — `R/distance.R`

- **Entry-point guard/namespacing:** the new implementation reports a direct error when annotation processing produces no ranges and calls `BiocParallel::bplapply` explicitly.
- **Annotation query-range rewrite:** explicit annotation `gene` spans remain the preferred mapping regions, but the new implementation adds robust GFF/rtracklayer/TxDb fallbacks, FASTA-name harmonization, shared-sequence checks, safer custom-assembly handling, and additional mitochondrial aliases.
- **Controlled seqname translation:** style conversion uses the shared helper and falls back to exact names on conversion failure.
- **Overlap double-count fix:** gene ranges are reduced before pileup so reads in overlapping genes cannot be returned twice through overlapping `which` ranges.
- **Replicate/depth corrections:** per-file identity is retained until aggregation; average coverage includes zero-coverage replicates; the minimum-depth boundary changes from `>` to inclusive `>=`; empty joins use `nrow()`; non-finite distances are removed.
- **Explicit genomic join:** WT and mutant tables merge only on `CHROM` + `POS`.
- **WT homozygosity bug fix:** the four WT frequency columns are selected explicitly instead of using a vector `:` expression.
- **Primary-read policy:** unmapped, QC-fail, secondary, and supplementary alignments are excluded; duplicate-marked primary reads remain eligible as in the old implementation.
- **RNA-seq pileup repair:** `simpleCigar=FALSE` restores spliced/complex RNA-seq alignments; linkage pileup depth is configurable; empty pileups and missing A/C/G/T categories are handled explicitly.
- **Defined replicate aggregation:** simple and depth-weighted frequency aggregation are formalized, while coverage always divides by the number of supplied files.

## LOESS fitting — `R/loess.R`

- **Explicit parallel namespace:** chromosome fitting uses `BiocParallel::bplapply` directly.
- **Fit failures are recoverable:** invalid spans/degenerate fits become failed evaluations rather than aborting the chromosome loop.
- **Boundary-safe local resolution:** sparse and edge span sets do not index outside `diff()` output.
- **AICc optimizer rewrite:** recursive duplicate-prone searching is replaced by cached iterative coarse-to-fine optimization with explicit all-failed handling and an iteration cap.
- **Finite local-minimum handling:** non-finite AICc values are ignored and minima are represented by stable indices.
- **Resolution parsing:** decimal precision is derived robustly from a non-scientific representation.
- **Degenerate AICc guards:** invalid residual variance/denominators return `NA` rather than poisoning optimization with Inf/NaN.
- **Tie handling:** tied minima resolve to an actually evaluated span instead of averaging to a never-evaluated span.
- **Chromosome fit validation:** the new implementation records `bestSpan`, checks sufficient informative markers and the final fit, and turns chromosome-specific failures into diagnostics.

## Peak detection/refinement — `R/peaks.R`

- **Parameterization/reproducibility:** cutoff, interval strategy, resampling count, and RNG seed are recorded settings.
- **Flat-apex fix:** a tied LOESS plateau is centered by median genomic position rather than taking the first maximum.
- **Resample failure tolerance:** failed half-marker fits return `NA` and contribute to success-rate metadata instead of aborting all refinement.
- **KDE probability mass:** density is integrated with local grid widths before selecting interval mass.
- **Explicit interval methods:** `hpd_span` spans all selected high-density modes; `shortest_contiguous` is an opt-in shortest single interval.
- **Deterministic per-chromosome RNG:** seeds derive from the recorded base seed and chromosome name and the caller RNG state is restored.
- **Refinement hardening:** density support is clipped to chromosome bounds; coordinates are valid integers; LOESS/density apex, resampling success rate, seed, and interval metadata are retained.
- **Initial cutoff control:** the old implementation cutoff calculation remains the default through `legacy_current`, with a configurable multiplier; `global_sd` is an explicit alternative.
- **NA-safe prePeak:** NA fitted values do not break linked-chromosome detection, and cutoff provenance is stored for plotting/output.

## Candidate calling, coding annotation, and expression — `R/candidates.R`

- **Candidate pipeline rewrite:** the fixed merged-BAM + VariantTools calling path is replaced by explicit in-memory pooled SNV evidence; WT evidence and delta AF are added; full SNVs are retained separately from coding effects; expression indexing is repaired.
- **Lazy TxDb lifecycle:** SNV discovery occurs before TxDb construction; a TxDb is built only when needed and is closed deterministically.
- **Explicit SNV criteria:** the hidden VariantTools calling filter is removed; every hard final candidate threshold is visible in `MmapprParam`.
- **Previous candidate defaults preserved:** candidate pileup remains capped at 250 reads/BAM; the mutant ALT boundary remains strict `>0.80`; the two-ALT-read default remains; `candidateMinDepth=1` adds no effective default gate.
- **Pure threshold helpers:** mutant and optional WT/delta-AF criteria are independently testable.
- **Optional WT filters:** WT ALT frequency and mutant-minus-WT delta AF can be used as explicit hard filters, but defaults leave them disabled.
- **In-memory replicate pooling:** mutant BAMs are piled independently and counts summed instead of creating a fixed `merged.tmp.bam`.
- **Explicit VRanges construction:** reference and strongest non-reference depths are retained directly.
- **Genomic evidence matching:** WT/base-count evidence is associated with variants via overlaps instead of string coordinates.
- **WT evidence reporting:** every mutant candidate receives WT ref/ALT/total depth, WT ALT frequency, mutant ALT frequency, and delta AF; sufficient-WT/shared-background flags are metadata by default.
- **WT query narrowing:** WT pileup runs only at candidate loci and uses the same candidate-stage cap.
- **Narrow warning handling:** only the recognized internal range warning is muffled.
- **Reference-bounds validation:** candidate/CDS ranges are checked against FASTA lengths before `predictCoding()`.
- **Coding annotation hardening:** sequence naming is harmonized across variants/TxDb/FASTA; only common levels are used; coding effects preserve mutant/WT evidence; full noncoding SNVs remain separately available.
- **Expression ratio helper:** raw `log2(mutant/wildtype)` remains the default; a positive pseudocount is explicit opt-in behavior.
- **Expression indexing bug fix:** WT and mutant columns are partitioned exactly rather than overlapping at the boundary.
- **GTF/GFF attribute support:** quoted GTF and key/value GFF3 attributes are parsed to preserve gene IDs/names.
- **Gene-range fidelity:** explicit annotation gene spans are preferred for linkage/expression, with rtracklayer and TxDb fallbacks; unknown strand symbols are normalized safely.
- **Expression summary repair:** counts are restricted to peak-overlapping supplied gene spans, use the shared primary-read filter, expose paired/strand controls, and use corrected group means. The output remains a descriptive raw-count screen rather than formal differential-expression inference.
- **Candidate ordering hardening:** empty/malformed inputs and missing consequence values are handled safely; coding severity and peak density remain the default ranking dimensions.
- **SNV-only scope is explicit:** the code refuses to invent normalized indel alleles from incomplete pileup insertion/deletion symbols.

## Analysis-object consistency — `R/mmapprData.R`

- **Parameter replacement invalidation:** replacing `param` after analysis results exist now warns and clears distances, peaks, and candidates.
- **Duplicate pileup implementation removed:** only `distance.R` owns `.getPileup()`, preventing read-filter/CIGAR fixes from diverging across copies.
- Documentation/printing is corrected to use the actual `snpDistance` terminology and current constructor names.

## Pipeline driver — `R/main.R`

- **External samtools requirement removed:** the R/Bioconductor implementation no longer aborts merely because an unrelated shell executable is absent.
- **Canonical output path logging:** the already-normalized output path is logged directly rather than prepending the current working directory.
- **Atomic final checkpoint:** `mmappr_data.RDS` is written through the new durable replacement helper.
- **Dead dependency helper removed:** the old executable-check helper is no longer part of the runtime code.

## Output, plotting, and serialization — `R/output.R`

- **Recovery/API fixes:** a missing output directory is recreated correctly; only devices opened by MMAPPR2 are closed; the documented `MmapprData` is returned invisibly.
- **Atomic RDS helper:** completed temporary serialization is verified before installation; the previous checkpoint is preserved/restored if replacement fails.
- **Unique output paths:** default/temp path generation gains collision-resistant suffixes.
- **Explicit safe overwrite:** no interactive prompt; overwrite requires `TRUE`; protected directories cannot be cleared; failures to empty/create the folder are detected.
- **Safe plot limits:** empty/constant/non-finite data receive valid finite axes.
- **Genome plot hardening:** numeric LOESS coordinates, explicit cutoff/peak overlays, safe device handling, and namespaced sequence ordering are used.
- **Peak plot diagnostics:** shading uses the actual baseline, cutoff is visible, and density overlays handle one-point/degenerate support.
- **TSV safety:** list/Bioconductor metadata are flattened only when necessary and embedded record-breaking whitespace is sanitized.
- **Legacy mutation-table compatibility:** the historical DetectedMutations schema is reconstructed while richer evidence remains in the additive full-candidate table.
- **Full candidate output:** `AllCandidateVariantsFor*.tsv` retains noncoding candidates and WT/delta-AF evidence while previous mutation/expression filenames continue to be written.

## Tests and user-facing documentation

- The old implementation tests are replaced/expanded by helper-level regression tests plus a synthetic integration fixture covering parameter validation, linkage helpers, LOESS selection, peak helpers, candidate rules/evidence, output helpers, preflight, and compatibility-preserving defaults.
- The vignette is rewritten around the current R-only BAM/FASTA/GTF workflow and the explicit new parameters rather than obsolete external-tool interfaces.
- Generated `.Rd` pages are updated to match the current class names, arguments, defaults, and output semantics.

## Intentional documentation differences

- The README intentionally uses the compact, original-style project introduction while updating the installation note to describe the current Rsamtools-based implementation.
- Generated help and vignette wording is kept version-neutral and describes old functionality versus new functionality where comparison context is useful.
- This change-annotation index is an added review aid and is not part of the executable package logic.

## Verification of this annotated copy

- **107** localized `[CHANGE — ...]` comments are present across the 11 runtime R files plus `NAMESPACE`.
- Removing comment-only/blank lines from every runtime R file yields the **same executable lines, in the same order**, as the supplied new implementation.
- Non-comment `NAMESPACE` directives are also identical to the supplied new implementation.
- The explanatory annotations use only old-functionality/new-functionality language and do not identify intermediate package revisions.
