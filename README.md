# The sum of two halves may be different from the whole. Effects of splitting sequencing samples across lanes.
Over the past two decades, the advances in high throughput sequencing (HTS) enabled the characterisation of biological processes at an unprecedented level of detail; as a result the vast majority of hypotheses in molecular biology rely on analyses of HTS data. However, achieving increased robustness and reproducibility of results remains one of the main challenges across analyses. Although variability in results may be introduced at various stages, such as alignment, summarisation or detection of differences in expression, one source of variability has been systematically omitted: the consequences of choices that influence the sequencing design which propagate through analyses and introduce an additional layer of technical variation.

In this study, we illustrate qualitative and quantitative differences in results arising from the splitting of samples across lanes, on bulk and single cell sequencing outputs. For bulk mRNAseq data, we focus on differential expression and enrichment analyses; for bulk ChIPseq data, we investigate the effect on peak calling, and the peaks' properties. At single cell level, we concentrate on the identification of cell subpopulations (cells clustered based on their expression profiles). We rely on the identity of markers used for assigning cell identities; both smartSeq and 10x data are presented.

We conclude that the observed reduction in the number of unique sequenced fragments reduces the level of detail on which the different prediction approaches depend. Further, the sequencing stochasticity adds in a weighting bias corroborated with variable sequencing depths.

Preprint: https://www.biorxiv.org/content/10.1101/2021.05.10.443429v1
