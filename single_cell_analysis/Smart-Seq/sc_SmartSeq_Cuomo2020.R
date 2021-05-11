## This document summarises the data analysis performed in the Smart-Seq section of the manuscript:
## The sum of two halves may be different from the whole.Effects of splitting sequencing samples across lanes 
## Eleanor C. Williams, Ruben Chazarra-Gil, Arash Shahsavari, and Irina Mohorianu

### Data Generation ### 

## 1. Six experimental study cases were designed to illustrate the effect of the different covariates of the data: the donor, the specific time-point and clusters as depicted in Table 1 of the manuscript. 
# For each of these cases, the data generation steps comprised: 

## 2. Initial Trimming of the FASTQ files with Trim Galore (0.4.1) : trim 10bp from the 5’ end, and 40 bp from the 3' end. 
# ${trim_galore_path} ${file_ID}_1.fastq.gz ${file_ID}_2.fastq.gz --path_to_cutadapt $cutadapt_path \
# --gzip --paired \
# --clip_R1 10 --clip_R2 10 \
# --three_prime_clip_R1 40 --three_prime_clip_R2  40 \
# --adapter "X" # to not trim based on adapter sequence

## 3. Generation of simulations by: i) Subsampling the FASTQ files using the seqtk toolkit (https://github.com/lh3/seqtk) and ii) concatenation of the subsamples
## Calculate number of entries to subsample as 0.5: 
# n_entries_sub=$(($(zgrep -c "^@" ${fastq_path}${file_ID}_1.fastq.gz)/2))
## Subsample:
# seqtk sample -s100 ${fastq_path}${file_ID}_1.fastq.gz  $n_entries_sub > ${file_ID}_1.fastq.gz
## Concatenate subsampled FASTQs from $sub_dir_1 and $sub_dir_2 into $comb_r1_r2 dir:
# cat ${sub_dir_1}${file_ID}_1.fastq.gz ${sub_dir_2}${file_ID}_1.fastq.gz | gzip -c > ${comb_r1_r2}${file_ID}_1.fastq.gz

## 4. Alignment to the reference genome using STAR (2.7.0a)
# STAR --genomeDir $genome_dir --runThreadN 12 --readFilesIn ${fastq_path}${file_ID}_1.fastq.gz,${fastq_path}${file_ID}_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${out_dir}${file_ID}.

## 5. Generation of a count matrix with featureCounts (2.0.0)
# featureCounts -p  -a ${gtf_file} -t exon -g gene_id -o ${out_dir}counts.txt ${star_path}*.bam

## 6. Obtaining of read quality metrics with fastQC
# fastqc ${fastq_path}${file_ID}_1.fastq.gz -o ${out_dir}

## 7. Obtaining of alignment and quantification metrics with multiQC
# multiqc ${fastQC_path} ${star_path} ${featureCouts_path} -o ${multiQC_path}


### Data Analysis ### 

# The following functions were used for the generation of results of the different individual Study cases of SmartSeq section

### 1. Collection of fastQC data ---------------------------------------------------------------------------

get_multiqc_data <-function (multiqc_path, which_data){
  # Retrieve custom multiQC data (either fastQC [per read end], or STAR and featureCounts [per sample]) from 'multiqc_general_stats.txt' report
  # This report must include fastQC, STAR and featureCounts for the function to retrieve the desired info.
  data <- read.csv(paste0(multiqc_path, "multiqc_data/multiqc_general_stats.txt"), sep = "\t", header = T)
  if (! which_data %in% c("fastqc", "star")) stop(print("which_data arg has to be either 'fastqc' or 'satr'."))
  # fastQC data
  if(which_data == "fastqc"){ 
    inv <- FALSE
    extract_cols <- c("Sample","FastQC_mqc.generalstats.fastqc.percent_duplicates", "FastQC_mqc.generalstats.fastqc.percent_gc")
    rename_cols <-  c("sample", "dup.pc", "gc.pc")
  }else{
    # STAR data
    inv <- TRUE
    extract_cols <-c("Sample", "featureCounts_mqc.generalstats.featurecounts.percent_assigned",
                     "featureCounts_mqc.generalstats.featurecounts.Assigned", 
                     "STAR_mqc.generalstats.star.uniquely_mapped_percent", 
                     "STAR_mqc.generalstats.star.uniquely_mapped")
    rename_cols <-  c("sample", "fC.pc", "fC", "star.pc", "star")
  }
  
  # get fastq samples
  fastq_samples <- sort(grep(pattern = "\\_1$|\\_2$", data$Sample, invert = inv))
  # extract % of duplicates and % GC
  df <- data[fastq_samples, extract_cols]
  # edit colnames
  colnames(df) <- rename_cols
  df
}

# Fetch fastQC data
case_path = ""
gt_df <- get_multiqc_data(multiqc_path = paste0(case_path, "Ground_truth/5.multiQC/"), which_data = "fastqc")
s1_df <- get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r2/"), which_data = "fastqc")
s2_df <- get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r3/"), which_data = "fastqc")
s3_df <-get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r2-r3/"), which_data = "fastqc")

# Add N of Unique and Duplicated Reads from file: `multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt`
file.name = "multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt"
gt_reads <- read.csv(paste0(case_path, "Ground_truth/5.multiQC/", file.name), sep = "\t")
s1_reads <- read.csv(paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r2/", file.name), sep = "\t")
s2_reads <- read.csv(paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r3/", file.name), sep = "\t")
s3_reads <- read.csv(paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r2-r3/", file.name), sep = "\t")
reads.list = list("gt" = gt_reads, "s1" = s1_reads, "s2" = s2_reads, "s3" = s3_reads)

# Group results
df <- data.frame("id" = rep(c("gt", "s1", "s2", "s3"), each = nrow(gt_df)), 
                      "sample" = c(as.character(gt_df$sample), as.character(s1_df$sample), 
                                   as.character(s2_df$sample), as.character(s3_df$sample)), 
                      "gc.pc"= c(gt_df$gc.pc, s1_df$gc.pc, s2_df$gc.pc, s3_df$gc.pc),
                      "dup.pc"= c(gt_df$dup.pc, s1_df$dup.pc, s2_df$dup.pc, s3_df$dup.pc), 
                      "Unique.reads" = unlist(lapply(reads.list, function(x) x$Unique.Reads)), 
                      "Duplicate.reads" = unlist(lapply(reads.list, function(x) x$Duplicate.Reads))
)

# Save
save_path = ""
write.csv(df, paste0(save_path, "1.multiQC.fastQC.csv"))

### 2. Alignment and Quantification QC data ---------------------------------------------------------------------------

# Fetch alignment and quantification QC data from multiQC report

gt_star <- get_multiqc_data(multiqc_path = paste0(case_path, "Ground_truth/5.multiQC/"), which_data = "star")
s1_star <- get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r2/"), which_data = "star")
s2_star <- get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r1-r3/"), which_data = "star")
s3_star <-get_multiqc_data(multiqc_path = paste0(case_path, "Simulation/5.multiQC/multiQC-comb_r2-r3/"), which_data = "star")

# Group alignment and quantification QC data
map.df <- data.frame("id" = rep(c("gt", "s1", "s2", "s3"), each = nrow(gt_star)), 
                     "sample" = c(as.character(gt_star$sample), as.character(gt_star$sample), 
                                  as.character(gt_star$sample), as.character(gt_star$sample)), 
                     "fC.pc"= c(gt_star$fC.pc, s1_star$fC.pc, s2_star$fC.pc, s3_star$fC.pc),
                     "fC"= c(gt_star$fC, s1_star$fC, s2_star$fC, s3_star$fC), 
                     "star.pc"= c(gt_star$star.pc, s1_star$star.pc, s2_star$star.pc, s3_star$star.pc),
                     "star"= c(gt_star$star, s1_star$star, s2_star$star, s3_star$star) 
)

# Save
write.csv(map.df, "2.multiQC.STAR.fC.csv")

### 3. Generation of Seurat Objects ------------------------------------------------------

create_seu_obj <- function(fC_path, min_cells, min_features){
  # Create Seurat Object from counts matrix (from featureCounts)
  suppressPackageStartupMessages(require(Seurat))
  counts_mat <- read.table(paste0(fC_path, "counts.txt"), header = T, row.names = 1)[, -c(1:5)]
  seu_obj <- CreateSeuratObject(counts = counts_mat, assay = "RNA", min.cells = min_cells, min.features = min_features)
  seu_obj
}

# Create Seurat objects
gt_seu <- create_seu_obj(fC_path = paste0(case_path, "Ground_truth/4.featureCounts/"),  min_cells = 3, min_features = 50)
s1_seu <- create_seu_obj(paste0(case_path, "Simulation/4.featureCounts/fC-comb_r1-r2/"), min_cells = 3, min_features = 50)
s2_seu <- create_seu_obj(paste0(case_path, "Simulation/4.featureCounts/fC-comb_r1-r3/"), min_cells = 3, min_features = 50)
s3_seu <- create_seu_obj(paste0(case_path, "Simulation/4.featureCounts/fC-comb_r2-r3/"), min_cells = 3, min_features = 50)

# Group nCount and nFeature data
count.df <- data.frame("id" = rep(c("gt", "s1", "s2", "s3"), each = ncol(gt_seu)), 
                       "sample" = c(colnames(gt_seu), colnames(s1_seu), colnames(s2_seu), colnames(s3_seu)), 
                       "nCount" = c(gt_seu$nCount_RNA, s1_seu$nCount_RNA, s2_seu$nCount_RNA, s3_seu$nCount_RNA), 
                       "nFeatures" = c(gt_seu$nFeature_RNA, s1_seu$nFeature_RNA, s2_seu$nFeature_RNA, s3_seu$nFeature_RNA)
)
# Save nCount and nFeature data
write.csv(count.df, "3.Seurat.counts.csv")

### 4. Pre-Processing of Seurat Objects ------------------------------------------------------
pre_process <- function(seu_obj){
  ## Pre-procesing Seurat object
  seu_obj <- SCTransform(seu_obj, new.assay.name = "logcounts", do.scale = T, variable.features.n = 3000, verbose = FALSE)
  seu_obj <- RunPCA(seu_obj, assay = "logcounts", npcs = 10, verbose = F)
  seu_obj <- RunUMAP(seu_obj,  dims = 1:10, verbose = F)
  seu_obj
}

gt_seu <- pre_process(gt_seu)
s1_seu <- pre_process(s1_seu)
s2_seu <- pre_process(s2_seu)
s3_seu <- pre_process(s3_seu)

### 5. Clustering of Seurat Objects ------------------------------------------------------

seurat_clusters <- function(seu, clust_alg) {
  ## Cluster Seurat objects
  seu <- FindNeighbors(seu, reduction = "pca", k.param = 30, compute.SNN = F, verbose = F)
  # clustering
  if(!(clust_alg %in% c(1:4))) stop ("clust_alg must be one of: 1, 2, 3, 4")
  seu <- FindClusters(object = seu, graph.name = tail(names(seu@graphs), 1), algorithm = clust_alg, resolution = 0.8 , verbose = F)
  seu
}

## Run Clustering
seu_list <- list("gt_seu"=gt_seu, "s1_seu"=s1_seu, "s2_seu"=s2_seu, "s3_seu"=s3_seu)
seu_list <- lapply(seu_list, function(seu) seurat_clusters(seu = seu, clust_alg = 1))

### 6. Differential Expression Analysis ------------------------------------------------------
seurat_markers <- function(seu){
  markers <- FindAllMarkers(seu, assay = "logcounts", logfc.threshold = 0.25,  min.pct = 0.1, test.use = "wilcox")
  seu@misc[["markers"]] <- markers
  return(seu)
}

## Run
seu.list.markers <- lapply(seu_list, function(seu) seurat_markers(seu = seu))
# Save processed seurat objects
saveRDS(seu.list.markers, "seurat.processed.list.rds")

### 7. Differential Expression Data Processing ------------------------------------------------------

# 1. Load Seurat object list with computed Markers 
markers.list <- lapply(seu.list.markers, function(seu) seu@misc$markers)

# 2. Filter markers 
filt_markers <- function(markers_df, logFC.thres, p.val.thres, pct.clust.thres){
  ## Filter markers based on logFC, adj p.val and pct of cells in cluster that have the marker
  markers_df[markers_df$avg_logFC > logFC.thres & markers_df$p_val_adj < p.val.thres & markers_df$pct.1 > pct.clust.thres, ]
}
## Run
markers.list.filt <- lapply(markers.list, function(df) filt_markers(markers_df = df, logFC.thres = 0.5, p.val.thres = 0.05, pct.clust.thres = 0.25))

# 3. Split markers by cluster
split_df <- function(markers_df, cluster_col){
  split_list = split(markers_df, markers_df[[cluster_col]])
  # edit names
  names(split_list) = paste0("cluster_", names(split_list))
  split_list
}
## Run
markers.list.split <- lapply(markers.list.filt, function(df) split_df(df, cluster_col = "cluster"))

# 4. Order dataframes by feature (`logFC`)
order_dfs <- function(markers_df_list, order_col){
  # Order list of dfs by order_col
  lapply(markers_df_list, function(cl) cl[ order(cl[[order_col]], decreasing = T), ])
}
## Run
markers.list.ord <- lapply(markers.list.split, function(df)order_dfs(markers_df_list = df, order_col  = "avg_logFC"))

# 5. Get gene names
get_feat_names <- function(by_cl_markers, feat_col ){
  # Get gene names of list of dfs
  lapply(by_cl_markers, function(cl) factor(cl[[feat_col]]))
}
## Run
markers.list.names <- lapply(markers.list.split, function(df) get_feat_names(by_cl_markers = df, feat_col = "gene"))
## Flatten list
markers.list.flat = unlist(markers.list.names, recursive = F)

## 6. Compute JSI 
split.vec <- function(vec, split.char.1, pos.1, split.char.2, pos.2, recursive){
  # Split a vector by split.char and get the item in position pos
  split.1 = unlist(lapply(strsplit(vec, split = split.char.1), function(x) x[pos.1]))
  if(recursive == T){
    split.2 = split.vec(split.1, split.char.1 = split.char.2, pos.1 = pos.2, split.char.2 = "X", pos.2 = "X", recursive = FALSE )
    return(split.2)
  }else{
    return(split.1)
  }
}
# Compute JSI list
jsi.list = list()
for (i in names(markers.list.flat)){
  set_i = markers.list.flat[[i]]
  for (j in names(markers.list.flat)){
    set_j = markers.list.flat[[j]]
    jsi.list[[paste0(i, "-", j)]] = compute_jsi(set_a = set_i, set_b = set_j, min_n = 1000)
  }
}

## Group JSI results data.frame
df = data.frame(comparison = names(jsi.list), 
                comp_1 = split.vec(vec = names(jsi.list), split.char = "-", pos = 1, split.char.2 = "X", pos.2 = "X", recursive = F), 
                obj_1 = split.vec(vec = names(jsi.list), split.char = "-", pos = 1, split.char.2 = "[.]", pos.2 =  1, recursive = T),
                cl_1 = split.vec(vec = names(jsi.list), split.char = "-", pos = 1, split.char.2 = "[.]", pos.2 =  2, recursive = T),
                
                comp_2 = split.vec(vec = names(jsi.list), split.char = "-", pos = 2, split.char.2 = "X", pos.2 = "X", recursive = F),
                obj_2 = split.vec(vec = names(jsi.list), split.char = "-", pos = 2, split.char.2 = "[.]", pos.2 =  1, recursive = T),
                cl_2 = split.vec(vec = names(jsi.list), split.char = "-", pos = 2, split.char.2 = "[.]", pos.2 =  2, recursive = T),
                
                JSI = unlist(jsi.list), row.names = names(jsi.list)
) 
## Save data
export.path = ""
write.csv(df, paste0(export.path, "1.Case_1_1.Markers.JSI.csv"))