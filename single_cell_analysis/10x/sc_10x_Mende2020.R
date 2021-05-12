# The paper presenting this data can be accessed at https://doi.org/10.1101/2020.01.26.919753
# The data is available at https://www.ebi.ac.uk/biostudies/studies/S-SUBS4
#
# Subsampling for simulations was performed with seqtk:
#
# seqtk sample -s1234 [input] 0.5 | gzip > [lane1]
# seqtk sample -s1235 [input] 0.5 | gzip > [lane2]
# cat [lane1] [lane2] > [concat]
#
# and choosing a unique seed for each S sample. Alignment and protein-coding
# gene quantification were performed with 10x Cellranger v3.1.0 using the 10x
# GRCh38 v3 reference transcriptome.

# start analysis
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ClustAssess)
theme_set(theme_bw())
set.seed(777)

min.cells = 0
min.features = 0

# the processing is performed identically for each GT and S sample and each donor
# load in the data
cts = readRDS(...)
meta = readRDS(...)
so = CreateSeuratObject(counts=cts, min.cells=min.cells, min.features=min.features, meta.data=meta)

mt.genes=grep("^MT-", rownames(so), value=FALSE)
rp.genes=grep("^RP[SL]", rownames(so), value=FALSE)

so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[rp.genes])

Idents(so) = so@meta.data$organ
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size=0)

Idents(so) = so@meta.data$library
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size=0)

plot(so@meta.data$nCount_RNA, so@meta.data$nFeature_RNA, xlab = 'n total counts', ylab='n unique features', col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), main='Scatterplot of features-seq depth')
plot(so@meta.data$percent.mt, so@meta.data$percent.rp, xlab = 'MT %', ylab='RP %', col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), pch = 16, main='Scatterplot of MT-RP')

# filter away some cells
so = subset(so, nFeature_RNA>1e3 & percent.mt<10 & percent.rp>20)

so = so[-c(mt.genes, rp.genes),]
so[['raw.seq.depth.no.MTRP']] = Matrix::colSums(GetAssayData(so, assay='RNA', slot='counts'))

so = SCTransform(so, return.only.var.genes=FALSE, verbose=FALSE)

summary(as.factor(so@meta.data$organ))
barplot(table(so@meta.data$organ))

summary(as.factor(so@meta.data$library))
barplot(table(so@meta.data$library))

table(so@meta.data$organ, so@meta.data$library)

n.abundant = 2000
matrix.so = GetAssayData(so, assay='SCT', slot='counts')
abundant.genes=rownames(matrix.so)[order(Matrix::rowSums(matrix.so), decreasing=TRUE)[1:n.abundant]]
length(abundant.genes)
abundant.genes[1:200]
rm(matrix.so)

Idents(so) = so@meta.data$library
so[['percent.abundant']] = PercentageFeatureSet(so, features=abundant.genes)
VlnPlot(so, features = 'percent.abundant', pt.size=0) + theme(legend.position='none')

so <- RunPCA(so, features = abundant.genes, verbose = FALSE, npcs=30)
so <- RunUMAP(so, dims = 1:30, verbose = FALSE)

FeaturePlot(so, features=c('nCount_RNA', 'nCount_SCT', 'percent.mt', 'percent.rp', ))

so = FindNeighbors(so, reduction='pca', dims=1:30, k.param=20, verbose=FALSE)
so = FindClusters(so, algorithm=3, verbose=FALSE)
DimPlot(so, group.by='seurat_clusters')


# now, we use separate seurat objects for GT and S
so.gt = readRDS(...)
so.s = readRDS(...)
common.barcodes = intersect(colnames(so.gt), colnames(so.s))
so.gt@meta.data$ecs = NA
gt.clustering = create_clustering(so.gt@meta.data$seurat_clusters)
s.clustering = create_clustering(so.s@meta.data$seurat_clusters)
so.gt@meta.data[common.barcodes, 'ecs'] = element_sim_elscore(gt.clustering, s.clustering)
FeaturePlot(so.gt, 'ecs')

lfc.thresh = log(2^1)
Idents(so.gt) = so.gt@meta.data$seurat_clusters
gt.markers = FindAllMarkers(so.gt, logfc.threshold=lfc.thresh, min.pct=0.0, test.use='roc', verbose=FALSE)
Idents(so.s) = so.s@meta.data$seurat_clusters
s.markers = FindAllMarkers(so.s, logfc.threshold=lfc.thresh, min.pct=0.0, test.use='roc', verbose=FALSE)
so.gt@meta.data$jsi = NA
so.gt@meta.data[common.barcodes, 'jsi'] = marker_overlap(gt.markers, s.markers, gt.clustering, s.clustering, rank_by='power')
FeaturePlot(so.gt, 'jsi')
