# ReplicationCollision

This repository contains code used in figure generation and analysis of the replication collision events assoociatied with the manuscript titled "Large tandem duplications in cancer result from transcription and DNA replication collision".

Author: Yang Yang yangyang0110@uchicago.edu, yangyangclover@gmail.com

Contact: Lixing Yang lixingyang@uchicago.edu

# Scripts Description
## PlotSignature
### plot-sigprofile.R
Using SV signature tables called by SigProfiler to plot simple SV signatures. This code generates Fig1.b-c, Fig,s1-s3

## GenerateRandomSV
### generate-random-sv.R
Generateing 4 (default) random SVs for each real/observed SVs according to the SV type, chromosome and SV size.
### run.sh
Running generate-random-sv.R in Linux system sample by sample to save time.

## ReplicationOrientation
### 01-replication-timing-direction-and-plot.R
Processing the raw cell line replication timing data, then determin the replication orientation of gene windows. This code generates Fig2.a,and Fig.s6
### 02-coding-gene-replication-direction.R
Annotating the replication orientation of cell lines for all coding genes.
### 03-add-coding-gene-replication-orientation-to-bkpt.R
Annotatiing the replication orientation of the coding gene where the breakpoint is locates.

## AnnodateBreakpoints
### add-signature-category-to-bkpt.R
Annotating all breakpoints with signature category.

## StrandBias
### 01-test-strand-bias.R
Testing strand bias for all signatures.
### 02-plot-strand-bias.R
Plotting strand bias for all signatures. This code generates Fig.2b-c, Fig.3, Fig.4a-b and Fig.s8

## SVEnrichmentGene
### 01-test-sv-enrichment.R
Testing the coding gene enrichment for all signatures.
### 02-plot-sv-enrichment.R
Plotting the coding gene enrichment for all signatures. This coed generates Fig.s4

## MultiVariateRegression
### test-and-plot-logistic-all-features-and-sig.R
Testing and Plotting logistic regression module, strand bias ~ SV + Expression + Replication timing + ... at all signatures. This code generates Fig.s9.

## PlotLargeTD
### 01-plot-large-td-freq.R
This code generates Fig.2d-e.
### 02-test-and-plot-large-td-and-gene-mutation.R
This code generates Fig.4c-f.
### 03-plot-large-td-gistic-hotspot.R
This code generates Fig.s10.
### 04-plot-large-td-and-survival.R
This code generates Fig.4g.
### 05-test-and-plot-large-td-and-brca-subtype.R
This code generates Fig.s7.

## DrugScreen
### ccle-cellline-drug-screen.R
This code generates Fig.6b.

# Additional packages
[dplyr 1.0.10](https://cran.r-project.org/web/packages/dplyr/index.html)
[ggplot2 3.4.0](https://rdocumentation.org/packages/ggplot2/versions/3.4.0)
[ggrepel 0.9.2](https://www.rdocumentation.org/packages/ggrepel/versions/0.9.2)
[data.table 1.14.6](https://www.rdocumentation.org/packages/data.table/versions/1.14.6)
[doParallel 1.0.17](https://www.rdocumentation.org/packages/doParallel/versions/1.0.17)
[foreach 1.5.2](https://www.rdocumentation.org/packages/foreach/versions/1.5.2)
[ggpubr 0.5.0](https://www.rdocumentation.org/packages/ggpubr/versions/0.5.0)
[stringr 1.5.0](https://www.rdocumentation.org/packages/stringr/versions/1.5.0)
[BSgenome.Hsapiens.UCSC.hg19 1.4.3](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)
[BSgenome.Hsapiens.UCSC.hg19.masked 1.3.993](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.masked.html)
[jtools 2.2.1](https://www.rdocumentation.org/packages/jtools/versions/2.2.1)
[ggstance 0.3.6](https://www.rdocumentation.org/packages/ggstance/versions/0.3.6)
[scales 1.2.1](https://www.rdocumentation.org/packages/scales/versions/1.2.1)
[facetscales 0.1.0.9000](https://github.com/zeehio/facetscales/tree/archived)
[ggh4x 0.2.3](https://cran.rstudio.com/web/packages/ggh4x/index.html)
[reshape 0.8.9](https://www.rdocumentation.org/packages/reshape/versions/0.8.9)
[tidyr 1.2.1](https://www.rdocumentation.org/packages/tidyr/versions/1.2.1)
[openxlsx 4.2.5.1](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.5.1)
[cowplot 1.1.1](https://www.rdocumentation.org/packages/cowplot/versions/1.1.1)
[GenomicRanges 1.44.0](http://bioconductor.riken.jp/packages/3.13/bioc/html/GenomicRanges.html)
[metR 0.13.0](https://www.rdocumentation.org/packages/metR/versions/0.13.0)
[patchwork 1.1.2](https://www.rdocumentation.org/packages/patchwork/versions/1.1.2)
[grid 4.1.1](https://www.rdocumentation.org/packages/grid/versions/3.6.2)
[egg 0.4.5](https://www.rdocumentation.org/packages/egg/versions/0.4.5)
[survival 3.4-0](https://www.rdocumentation.org/packages/survival/versions/3.4-0)
[survminer 0.4.9](https://rdocumentation.org/packages/survminer/versions/0.4.9)
[magrittr 2.0.3](https://www.rdocumentation.org/packages/magrittr/versions/2.0.3)
[preprocessCore 1.54.0](https://bioc.ism.ac.jp/packages/3.13/bioc/html/preprocessCore.html)
[gdata 2.18.0.1](https://www.rdocumentation.org/packages/gdata/versions/2.18.0)
[parallel 4.1.1](https://www.rdocumentation.org/packages/parallel/versions/3.6.2)
