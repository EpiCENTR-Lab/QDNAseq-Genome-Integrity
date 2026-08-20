# Genome integrity using low-pass whole-genome sequencing and QDNAseq

## Overview

This pipeline performs a genome integrity screen using low-pass whole-genome sequencing (WGS) BAM files.

It is designed for cultured cell quality control, including iPSCs, ESCs, fibroblasts, and engineered cell lines.

The script produces genome-wide copy number plots and summary tables that help detect:

- Whole chromosome gains or losses
- Chromosome arm-level changes
- Large structural abnormalities

This tool is intended for genome integrity QC, not full CNV discovery.

## Pipeline Steps

1. **Read binning** – reads are counted in genomic bins using QDNAseq
2. **Filtering** – problematic regions are removed (low mappability, blacklist regions)
3. **Bias correction** – GC and mappability bias correction
4. **Normalization** – counts converted to relative copy number
5. **Segmentation** – contiguous regions with similar copy number are identified
6. **Visualization** – genome-wide plots are generated

## Software Requirements

R version 4.1 or newer is recommended.

**Install packages in R:**

```bash
Rscript install_packages.R
```

### QDNAseq Genome Annotation Packages

Install the annotation matching your genome build.

For hg38:

```r
BiocManager::install("QDNAseq.hg38")
```

For hg19:

```r
BiocManager::install("QDNAseq.hg19")
```

### Reproducing the Exact Environment

`install_packages.R` pins dependencies to the exact versions used during development of the pipeline:

- **R**: 4.5.3 — the script checks your R version and warns (without stopping) if it differs
- **CRAN packages** (dplyr, tidyr, tibble, ggplot2, optparse): pinned to exact versions via `remotes::install_version()`
- **Bioconductor packages** (QDNAseq, Biobase, DNAcopy, Rsamtools, GenomeInfoDb): pinned via Bioconductor release `3.22`. Bioconductor doesn't support pinning individual packages to arbitrary versions — the release as a whole determines which package versions you get, and release 3.22 is what produced the versions used during development (QDNAseq 1.46.0, Biobase 2.70.0, DNAcopy 1.84.0, Rsamtools 2.26.0)
- **QDNAseq.hg38**: pinned to the exact GitHub commit used for the analysis (`cf7c07e39de0ac64a9c38cb030cba4626e2aae83`), since annotation packages aren't captured by `sessionInfo()`

Running `Rscript install_packages.R` on a clean R installation should reproduce this environment exactly. The full session details from the original analysis run — including the OS, BLAS/LAPACK libraries, and every loaded package with its version — are recorded in `SESSION_INFO.txt`, written automatically at the end of each pipeline run.

## Input Data

The script expects BAM files and BAM index files in a folder within your working directory called `bam`. BAM files must be indexed (`.bai` files), and the BAM filename becomes the sample name.

Example:

```
bam/
    sample1.bam
    sample1.bam.bai
    sample2.bam
    sample2.bam.bai
```

### Sex Annotation File

If your samples include both XX and XY genomes, modify the file called `sample_sex.csv` to replace the placeholder example text.

Example:

| sample      | sex |
| ----------- | --- |
| iPSC_A      | XX  |
| iPSC_B      | XY  |
| Fibroblast1 | XX  |

Sex must be `XX` or `XY`.

Sample names must match BAM filenames without `.bam`.

## Running the Pipeline

Run from RStudio's Terminal tab with:

```bash
Rscript run_qdnaseq_genome_integrity.R \
  --bam-dir bam \
  --out-dir qdnaseq_output \
  --bin-size 100 \
  --genome hg38 \
  --sex-map sample_sex.csv \
  --map-cutoff 80 \
  --bases-cutoff 99.9 \
  --seed 170826
```

This command reads BAM files from `bam/`, writes results to `qdnaseq_output/`, uses 100 kb genomic bins, and assumes the hg38 reference genome.

## Important Parameters

### Bin size

```
--bin-size 100
```

Recommended:

- 50 kb – higher resolution but noisier
- 100 kb – recommended default
- 500 kb – smoother but lower resolution

### QC thresholds

```
--map-cutoff 80
--bases-cutoff 99.9
```

### Sex-chromosome ploidy correction

```
--no-sex-ploidy-correction
```

XY samples naturally carry one copy of chrX and one copy of chrY, versus two copies of each autosome. Left uncorrected, this ~50% dosage difference reads on the genome-wide plot as an apparent whole-chromosome loss on X and Y, even in a chromosomally normal sample.

By default, the pipeline rescales chrX/chrY copy number by 2x for XY samples before computing the log2 ratio, so a normal XY sample sits at ~0 on X and Y, just like the autosomes. A genuine loss or gain on X/Y still shows up as a deviation from 0 after this correction. Correction is applied automatically to any sample classified as XY (via `--sex-map`, which must be provided as there is no default). XX samples are unaffected.

To see the raw, uncorrected QDNAseq output instead, pass `--no-sex-ploidy-correction`.

Whether the correction was applied is recorded per sample in the `sex_ploidy_corrected` column of `sample_qc_summary.tsv`.

### Plot range

```
--ymin -1.2
--ymax 1.2
```

### Reproducibility

DNAcopy's segmentation algorithm (used internally by `segmentBins()`) relies on permutation-based significance testing, so exact segment breakpoints can vary slightly between runs unless a seed is fixed.

The pipeline sets a random seed automatically via `--seed` (default: `170826`) before segmentation begins, and records the seed used in `SESSION_INFO.txt` for every run.

To reproduce the figures from the published analysis exactly, use `--seed 170826` (the default), or the seed value recorded in the `SESSION_INFO.txt` from that run.

## Output Files

Example output directory:

```
qdnaseq_output/
    genomewide_sample1.jpg
    sample1.qdnaseq.rds
    sample1_bins.tsv.gz
    sample1_segments.tsv
    genomewide_sample2.jpg
    sample2.qdnaseq.rds
    sample2_bins.tsv.gz
    sample2_segments.tsv
    sample_qc_summary.tsv
    SESSION_INFO.txt
```

## Output Descriptions

### Genome plot — `genomewide_SAMPLE.jpg`

- Black points = high-quality bins
- Purple lines = segmented copy number
- Dashed line = diploid baseline

Interpretation:

- Flat near 0 → diploid genome
- Whole chromosome above 0 → gain
- Whole chromosome below 0 → loss

### Bin-level data — `sample_bins.tsv.gz`

Contains one row per genomic bin including genomic coordinates, QC status, raw copy number, and segmented copy number.

### Segment table — `sample_segments.tsv`

Contains segmented copy number regions.

### QDNAseq object — `sample.qdnaseq.rds`

Stores the full QDNAseq analysis object for future analysis.

### QC summary — `sample_qc_summary.tsv`

One row per sample with QC statistics including number of bins, median log2 copy number, and fraction of bins failing QC, plus `sex_mode` (XX or XY, as used for that sample) and `sex_ploidy_corrected` (whether the chrX/chrY ploidy correction described above was applied).

### Session information — `SESSION_INFO.txt`

Records software versions and command used to run the pipeline.

## Recommended Sequencing Parameters

- Coverage: 0.2×–5×
- Read length: 75–150 bp
- Paired-end sequencing preferred
- Reference genome: hg38 recommended

## Limitations

This pipeline cannot reliably detect small CNVs.

Detection reliability:

| Size    | Reliability          |
| ------- | --------------------- |
| >10 Mb  | Very reliable          |
| 2–10 Mb | Usually detectable     |
| <1 Mb   | Unreliable with low-pass WGS |

Given 100 kb bins and observed genome-wide noise, fully clonal copy-number changes are statistically detectable below 1 Mb; however, we adopt 1 Mb as a conservative calling threshold to guard against locus-specific mapping artifacts and treat sub-Mb findings as requiring orthogonal investigation.

## Troubleshooting

**Missing BAM index**

Create indexes with:

```bash
samtools index sample.bam
```

**Cannot load QDNAseq annotations**

Install the annotation package matching your genome build:

```r
BiocManager::install("QDNAseq.hg38")
```

**X/Y chromosome issues**

A sample showing a clean whole-chromosome-scale offset (~±1.0 log2) confined to chrX/chrY, with autosomes flat, most likely indicates a sex-map/sample mismatch rather than a true copy-number event. Verify `sample_sex.csv` against the sample's known genotype.
