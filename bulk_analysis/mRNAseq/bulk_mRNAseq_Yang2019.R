# Bulk mRNAseq data

#The paper presenting this data can be accessed at https://pubmed.ncbi.nlm.nih.gov/31078527/
#The 0h and 12h samples were used and there are 2 biological replicates for each (SRR7624365,SRR7624366 and SRR7624371,SRR7624372 respectively)

#Preprocessing ------------------------------------------------------------------

#Subsampling (using the randomReadSubSample script from the metameta package: https://github.com/dsquintana/metameta)
#For k=2
#python3 randomReadSubSample.py -f1 [forward sample] -f2 [reverse sample] -s 0.5 -r 0 -g 1 -o [output prefix]

#Concatenating
#For k=2
#cat [subsample 1] [subsample 2] > [output filename]

#Alignment using STAR (paired-end)
#STAR --genomeDir [Mus musculus genome] --readFilesIn [forward sample] [reverse sample] --outSAMtype BAM SortedByCoordinate \
#--readFilesCommand zcat --runMode alignReads --outFileNamePrefix [output prefix]

#Expression quantification using featureCounts
#featureCounts -p -t exon -g gene_id -a [Mus musculus gtf] -o [output file] [BAM file]


#Creating count matrix ----------------------------------------------------------

library(preprocessCore)
creating_count_matrix <-function(samplenames,samplelocations){
  meta = data.frame(index=samplenames,n_replicate=c(1,2,1,2),type=c('0h','0h','12h','12h'))
  meta$type = as.factor(meta$type)
  meta$type=relevel(meta$type, '0h')
  gene_ids = read.csv(samplelocations[1],sep='\t',skip=1)
  cts = data.frame(gene_id = gene_ids$Geneid)
  rownames(cts) = cts$gene_id
  for (i in (1:4)){
    currentfile = read.csv(samplelocations[i],sep='\t',skip=1)
    columns = colnames(currentfile)
    cts[,paste('sam_',samplenames[i],sep='')]=currentfile[,columns[7]]
  }
  cts = subset(cts,select=-c(gene_id))
  cts.filtered = cts
  cts.qnorm=data.frame(normalize.quantiles(as.matrix(cts.filtered)),row.names=rownames(cts.filtered))
  colnames(cts.qnorm)=colnames(cts.filtered)
  return_list <- list("cts" = cts,"cts.qnorm"=cts.qnorm, "meta" = meta)
  return(return_list)
}
#gt=creating.count.matrix(c('0h_1','0h_2','12h_1','12h_2'),
#                         c('SRR7624365_counts.txt','SRR7624366_counts.txt','SRR7624371_counts.txt','SRR7624372_counts.txt'))

#MA -----------------------------------------------------------------------------

ma.plot <-function(sample1,sample2){
  counts = cbind(sample1,sample2)
  zero.mask = !(counts[,1] == 0 | counts[,2] == 0)
  l1 = log2(counts[zero.mask, 1])
  l2 = log2(counts[zero.mask, 2])
  m = l1 - l2
  a = 0.5 * (l1 + l2)
  data = data.frame(A = a, M = m)
  a.binned = cut(a, breaks = seq(0,20))
  data.binned = data.frame(A = a.binned, M = m)
  return(data.binned)
}
#ma.plot(gt$0h_1$cts,k2$0h_1$cts)

#DE -----------------------------------------------------------------------------

#find differentially expressed genes for a set of samples
library(edgeR)
de.set <-function(samples){
  edg = DGEList(samples,group=c(1,1,2,2))
  edg = estimateDisp(edg)
  exact =exactTest(edg, dispersion=0.2)
  exact$table$adjustedp = p.adjust(exact$table$PValue,method='BH')
  de<-rownames(exact$table[(abs(exact$table$logFC) > 0.5 & exact$table$adjustedp < 0.05),])
  return(list('de'=de,'exact'=exact$table))
}
#de.set(gt$cts.qnorm)

#Cross plots --------------------------------------------------------------------

crossplots <- function(cts1,cts2){
  dgeobj1 = DGEList(cts1,group=c(1,1,2,2))
  exact1 =exactTest(dgeobj1, dispersion=0.2)
  exact1$table$adjustedp = p.adjust(exact1$table$PValue,method='BH')
  ind1=which((abs(exact1$table$logFC) > 0.5 & exact1$table$adjustedp<0.05), arr.ind=TRUE)
  
  dgeobj2 = DGEList(cts2,group=c(1,1,2,2))
  exact2 =exactTest(dgeobj2, dispersion=0.2)
  exact2$table$adjustedp = p.adjust(exact2$table$PValue,method='BH')
  ind2=which((abs(exact2$table$logFC) > 0.5 & exact2$table$adjustedp<0.05), arr.ind=TRUE)
  
  names1=rownames(cts1)
  names2=rownames(cts2)
  intersectnames=intersect(names1,names2)
  exact1$table = exact1$table[intersectnames,]
  exact2$table = exact2$table[intersectnames,]
  
  logFC_frame = data.frame(sample1=exact1$table$logFC,
                           sample2=exact2$table$logFC,
                           abn=(rowMeans(cts1[intersectnames,])+rowMeans(cts2[intersectnames,]))/2,
                           gene_id=intersectnames)
  rownames(logFC_frame)=intersectnames
  
  return(logFC_frame)
}

#Enrichment analysis ------------------------------------------------------------

#finding enriched terms for a set of samples
library(gprofiler2)
enrichment.set <-function(samples){
  edg = DGEList(samples,group=c(1,1,2,2))
  edg = estimateDisp(edg)
  exact =exactTest(edg, dispersion=0.2)
  exact$table$adjustedp = p.adjust(exact$table$PValue,method='BH')
  de<-rownames(exact$table[(abs(exact$table$logFC) > 0.5 & exact$table$adjustedp < 0.05),])
  gost = gost(de,organism='mmusculus')
  gostterms = gost$result$term_id
  return(gostterms)
}
#enrichment.set(gt$cts.qnorm)

#Noise analysis -----------------------------------------------------------------

library(noisyr)
#Suppose the current working directory contains the 4 BAMs (sorted using samtools sort and indexed using samtools index) 
#for a single set (e.g. 0h rep1,2, 12h rep1,2 for GT or S) and the Mus musculus gtf (e.g. from Ensembl).
noise = calculate_expression_similarity_transcript(path.bams='.',path.gtf = '.')
plot_expression_similarity(noise)
