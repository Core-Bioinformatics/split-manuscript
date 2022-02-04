## This document summarises the data analysis performed in the ChIPseq section of the manuscript:
## The sum of two halves may be different from the whole. Effects of splitting sequencing samples across lanes 
## Eleanor C. Williams, Ruben Chazarra-Gil, Arash Shahsavari, and Irina Mohorianu

#The paper presenting this data can be accessed at https://pubmed.ncbi.nlm.nih.gov/31078527/
#The 0h and 12h samples were used for H3K4me3 and H3K27ac (SRR7624381,SRR7624384 and SRR7624389,SRR7624392 respectively)

#Preprocessing ------------------------------------------------------------------

#Subsampling (using the randomReadSubSample script from the metameta package: https://github.com/dsquintana/metameta)
#For k=2
#python3 randomReadSubSample.py -f [fastq file] -s 0.5 -r 0 -g 1 -o [output_prefix]

#Concatenating
#For k=2
#cat [subsample 1] [subsample 2] > [output_filename]

#Alignment using bowtie2
#bowtie2 -p 5 -q --local -x [Mus musculus genome] -U [sample] -S [output SAM name]
#samtools view -h -S -b -q 5 -o [output BAM name] [output SAM name]

#Peak calling using macs2
#macs2 callpeak -t [BAM file] -f BAM --outdir [output directory] -n [output name prefix]

#Matching peaks between 2 time-points -------------------------------------------

library(preprocessCore)
library(DESeq2)
collecting_peaks<-function(file1,file2){
  file1$midpoint = (file1$stop + file1$start)/2
  file2$midpoint = (file2$stop + file2$start)/2
  file1$appears = FALSE
  file2$appears = FALSE
  i=1
  j=1
  
  listofvals1 = unfactor(unique(file1$Chr))
  file2 = file2[file2$Chr %in% listofvals1,]
  listofvals2 = unfactor(unique(file2$Chr))
  file1 = file1[file1$Chr %in% listofvals2,]
  listofchrvals = unfactor(unique(file1$Chr))
  
  while (i<=nrow(file1)&&j<=nrow(file2)){
    if (file1[i,'Chr']!=file2[j,'Chr']){
      itemi = which(listofchrvals==file1[i,'Chr'])
      itemj = which(listofchrvals==file2[j,'Chr'])
      if (itemi<itemj){
        while (file1[i,'Chr']!=file2[j,'Chr']){
          file1[i,'appears']<-FALSE
          i=i+1
        }
      }
      else {
        while (file1[i,'Chr']!=file2[j,'Chr']){
          file2[j,'appears']<-FALSE
          j=j+1
        }
      }
    }
    else {
      if (file1[i,'midpoint']<file2[j,'start']){
        file1[i,'appears']<-FALSE
        i=i+1
      }
      else if (file2[j,'midpoint']<file1[i,'start']){
        file2[j,'appears']<-FALSE
        j=j+1
        
      }
      else if (file1[i,'midpoint']>=file2[j,'start']&&file1[i,'midpoint']<=file2[j,'stop']){
        file1[i,'appears']<-TRUE
        file2[j,'appears']<-TRUE
        file1[i,'match']<-file2[j,'id']
        file2[j,'match']<-file1[i,'id']
        i=i+1
        j=j+1
      }
      else if (file2[j,'stop']-file2[j,'start']<0.5*(file1[i,'stop']-file1[i,'start'])){
        file2[j,'appears']<-NA
        j=j+1
      }
      else if (file1[i,'stop']-file1[i,'start']<0.5*(file2[j,'stop']-file2[j,'start'])){
        file1[i,'appears']<-NA
        i=i+1
      }
      else if (file1[i,'midpoint']>=file2[j,'stop']){
        file2[j,'appears']<-NA
        j=j+1
      }
      else if (file2[j,'midpoint']>=file1[i,'stop']){
        file1[i,'appears']<-NA
        i=i+1
      }
    }
  }
  
  file1_match = na.omit(file1[file1$appears == TRUE,])
  file2_match = na.omit(file2[file2$appears == TRUE,])

  colnames(file1_match)=c('sample0','Chr','start','stop','midpoint','appears','sample1')
  colnames(file2_match)=c('sample1','Chr','start1','stop1','midpoint1','appears','sample0')
  merge_file1_file2=merge(file1_match,file2_match,by=c('sample0','sample1','appears','Chr'))
  
  merged_df = data.frame(id1=merge_file1_file2$sample0,
                         id2=merge_file1_file2$sample1,
                         Chr=merge_file1_file2$Chr,
                         overallstart=pmin(merge_file1_file2$start,merge_file1_file2$start1),
                         overallstop=pmax(merge_file1_file2$stop,merge_file1_file2$stop1),
                         start1=merge_file1_file2$start,
                         stop1=merge_file1_file2$stop,
                         start2=merge_file1_file2$start1,
                         stop2=merge_file1_file2$stop1)
  return(merged_df)
}

#Here sample1file and sample2file are the locations of the macs2 narrowPeak output
sample1 = read.csv(sample1file,sep='\t',header=FALSE)
sample1 = subset(sample1,select=c(V4,V1,V2,V3))
colnames(sample1)=c('id','Chr','start','stop')
sample2 = read.csv(sample2file,sep='\t',header=FALSE)
sample2 = subset(sample2,select=c(V4,V1,V2,V3))
colnames(sample2)=c('id','Chr','start','stop')
matchedpeaks = collecting_peaks(sample1,sample2)
matchedpeaks$peak_length_ratio = (matchedpeaks$stop1-matchedpeaks$start1)/(matchedpeaks$stop2-matchedpeaks$start2)

#Creating expression matrix -----------------------------------------------------
sample1 = read.csv(sample1file,sep='\t',header=FALSE)
rownames(sample1)=sample1$V4
sample1.matched = sample1[sample1$V4 %in% matchedpeaks$id1,]
sample2 = read.csv(sample2file,sep='\t',header=FALSE)
rownames(sample2)=sample2$V4
sample2.matched = sample2[sample2$V4 %in% matchedpeaks$id2,]
expression = data.frame('sam_0h'=sample1.matched$V5,'sam_12h'=sample2.matched$V5)

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
#ma.plot(expression$sam_0h,expression$sam_12h)
