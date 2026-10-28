# MMAPPR2

*Mutation Mapping Analysis Pipeline for Pooled RNA-Seq*

MMAPPR2 helps find DNA changes linked to a trait in pooled RNA-sequencing
(RNA-seq) data from second-generation (F2) genetic crosses. It compares groups
with and without the trait, highlights promising regions of the genome, and
lists possible single-base DNA changes.

## Installation

You need **R 4.4.0 or newer**. Unzip this source package (version
`0.99.24.9011`), set R's working directory to the folder containing
`DESCRIPTION`, and run:

```r
install.packages("remotes")
remotes::install_local(".", dependencies = NA, upgrade = "never")
```

This installs the package and its required dependencies. MMAPPR2 was
[removed from current Bioconductor releases](https://bioconductor.org/about/removed-packages/),
so `BiocManager::install("MMAPPR2")` will not install this version.

## Quick start

You need coordinate-sorted BAM files (aligned RNA-seq data) for the **mutant**
group (with the trait) and **wild-type** group (without it). You also need the
reference genome (FASTA file) used for alignment and a matching gene annotation
(GTF or GFF3 file).

Replace the example filenames with your own:

```r
library(MMAPPR2)

param <- mmapprParam(
  wtFiles = "wt_pool.bam",
  mutFiles = "mutant_pool.bam",
  refFasta = "reference.fa",
  gtf = "annotation.gtf",
  outputFolder = "mmappr2_results"
)

results <- mmappr(param)
```

Use a new or empty results folder. Missing BAM and FASTA indexes are created
automatically. See the [full workflow guide](vignettes/MMAPPR2.Rmd) for more options.

## Output files

Files are saved to `mmappr2_results` (or your chosen folder). Some appear only
when enough data are available.

| File | What it contains |
| --- | --- |
| `AllCandidateVariantsFor<seq>.tsv` | Possible single-base changes, with evidence from both groups. |
| `DetectedMutationsFor<seq>.tsv` | Predicted effects of changes in genes that make proteins. |
| `DifferentiallyExpressedGenesFor<seq>.tsv` | Gene-expression differences; **not** a statistical test, despite the filename. |
| `genome_plots.pdf` and `peak_plots.pdf` | Figures showing where the strongest evidence appears. |
| `aicc_plots.pdf` | Optional figures showing how the smoothing setting was chosen. |
| `mmappr_data.RDS` | Saved results for creating new figures later. |
| `mmappr2.log` | Run messages and warnings. |

`<seq>` means a chromosome or other genome sequence. Files ending in `.tsv` are
tables you can open in a spreadsheet. Start with the PDF figures, then review
the candidate tables. To save `aicc_plots.pdf`, add `exportAiccPlots = TRUE` to
`mmapprParam()`.

## Figure colors and replotting

To change the color scheme for a new analysis, replace `mmappr(param)` above
with `mmappr(param, palette = "colorblind")`. The options are `default`,
`colorblind`, and `grayscale`.

To make new PDF figures later **without running the analysis again**, use the
saved results file:

```r
plotMmapprData(
  "mmappr2_results/mmappr_data.RDS",
  plots = c("genome", "peaks"),
  palette = "grayscale"
)
```

New figures go into a `replots` folder by default. The saved results and
candidate tables stay unchanged. `overwrite = TRUE` in `plotMmapprData()`
replaces only the selected figures; **in `mmapprParam()`, it can clear the
existing results folder**. For individual color changes and other options,
see `?plotMmapprData` in R.

## Limitations

MMAPPR2 looks for single-base DNA changes, not insertions, deletions, or larger
changes. RNA-seq covers the genome unevenly, so some changes may be missed;
gene expression and how reads align can also affect results. The analysis
reports one main linked region per chromosome. The gene-expression summary is
not a statistical test.

## Documentation, citation, and support

Read the [full workflow guide](vignettes/MMAPPR2.Rmd), or use `?mmapprParam`,
`?mmappr`, and `?plotMmapprData` in R for detailed help. The predecessor of
MMAPPR2 was described in the
[original MMAPPR paper (Hill et al., 2013)](https://genome.cshlp.org/content/23/4/687.full.pdf).
To cite this package, run `citation("MMAPPR2")` in R.

Authors: Kyle Johnsen, Nathaniel Jenkins, and Jonathon Hill.
License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.html).
Report problems through [GitHub issues](https://github.com/jonathonthill/MMAPPR2/issues).
