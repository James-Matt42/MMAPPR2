# MMAPPR2
## Mutation Mapping Analysis Pipeline for Pooled RNA-Seq
### Authors
Kyle Johnsen, Nathaniel Jenkins, and Jonathon Hill

### Introduction
MMAPPR2 maps mutations resulting from pooled RNA-seq data from the F2
cross of forward genetic screens. Its predecessor is described in a paper published
in Genome Research (Hill et al. 2013). MMAPPR2 accepts aligned BAM files as well as
a reference genome as input, identifies loci of high sequence disparity between the
control and mutant RNA sequences, predicts variant effects, 
and outputs a ranked list of candidate mutations.

[See vignette for instructions](vignettes/MMAPPR2.Rmd)

Publication for the [original MMAPPR](http://genome.cshlp.org/content/23/4/687.full.pdf)

## Installation Notes
MMAPPR2 performs BAM and FASTA operations through the Bioconductor package Rsamtools. An external Samtools executable is not required.

## Candidate-pileup memory and diagnostics

Candidate SNV calling uses an adaptive memory strategy by default. On machines with
sufficient available RAM, ordinary refined intervals retain the fast one-shot pileup.
Large intervals or constrained environments are processed in bounded genomic chunks;
within those chunks BAMs are pooled incrementally so all per-BAM pileup tables do not
need to coexist in memory. WT candidate support uses the same adaptive strategy. If a
one-shot mutant pileup still encounters an allocation failure, MMAPPR2 retries with
smaller chunks and temporary RDS staging under R's `tempdir()`; staged files are
removed automatically. Set `candidatePoolMode` and `candidateChunkSize` in
`mmapprParam()` only when explicit control is needed.

Set `exportAiccPlots=TRUE` in `mmapprParam()` to have `outputMmapprData()` write
`aicc_plots.pdf` using the span/AICc evaluations already retained by `loessFit()`.

