##Supplementary figure 5A,B,C====================================================
library(ggplot2)

sample='2RA3'
#using matrix which is output after using perl/R2NR.pl and createMatrix_fromFasta.pl
expressionMatrix=read.csv(paste0(sample,'_matrix.txt'),skip=1,header=FALSE)
colnames(expressionMatrix)=c('fragment','split1','split2','whole')
expressionMatrix[is.na(expressionMatrix)]<-0

seqdepth_ratio=(sum(expressionMatrix$split1)+sum(expressionMatrix$split2))/sum(expressionMatrix$whole)

#take fragments which are non-zero in the whole and at least one split
nonzero.ratio=expressionMatrix[(expressionMatrix$split1+expressionMatrix$split2)!=0 & expressionMatrix$whole!=0,]
nonzero.ratio$abundances=ceiling(log2(nonzero.ratio$whole))

#compute ratios
nonzero.ratio$combined_vs_whole=(nonzero.ratio$split1+nonzero.ratio$split2)/nonzero.ratio$whole
nonzero.ratio$split1_vs_whole=(nonzero.ratio$split1)/nonzero.ratio$whole
nonzero.ratio$split2_vs_whole=(nonzero.ratio$split2)/nonzero.ratio$whole

#melt
nonzero.ratio.melt=reshape2::melt(nonzero.ratio,id=c('split1','split2','whole','abundances'))

#plot boxplot
A=ggplot(nonzero.ratio.melt,aes(y=log2(value),x=as.factor(abundances),color=variable))+
  geom_boxplot(outlier.shape=NA)+
  geom_hline(aes(yintercept = log2(seqdepth_ratio)),color='red')+
  theme_bw()+xlab('log2(abundance in whole sample)')+
  ylab('log2(ratio between abundances in two samples)')+
  ggtitle(sample)+ylim(-4,4)+ theme(legend.position = c(0.9, 0.8))

##Supplementary figure 5D,E,F======================================================
sample='2RA3'
#using matrix which is output after using perl/R2NR.pl and createMatrix_fromFasta.pl
expressionMatrix=read.csv(paste0(sample,'_matrix.txt'),skip=1,header=FALSE)
colnames(expressionMatrix)=c('fragment','split1','split2','whole')
expressionMatrix[is.na(expressionMatrix)]<-0

#take fragments which are non-zero in the whole and at least one split
whole.specific=expressionMatrix[(expressionMatrix$split1+expressionMatrix$split2)==0 & expressionMatrix$whole!=0,]
whole.specific$binnedabn=factor(cut(whole.specific$whole,breaks=c(seq(0,100,5),150,250),right=FALSE),
                                levels=unique(cut(1:200,breaks=c(seq(0,100,5),150,250),right=FALSE)))

wholespecific.grouped=data.frame(table(whole.specific$binnedabn))
D=ggplot(na.omit(wholespecific.grouped),aes(x=Var1,y=log2(Freq+1)))+
  geom_bar(stat='identity')+
  ylab('log2(frequency+1)')+
  xlab('Binned abundance')+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##Supplementary figure 5G,H,I=====================================================
sample='2RA3'
#using output from python/transcript_coverage_pyBigWig.py
diff=read.csv(paste0(sample,'_average_difference.csv'))
diff$scaleddifference=(seqdepth_ratio*diff$avgabn_whole - diff$avgabn_split)/(myratio*diff$avgabn_whole)
diff$sign=ifelse(diff$scaleddifference<0,'-','+')
diff$absolutedifference=abs(diff$scaleddifference)
diff$logabn=cut(log2(diff$avgabn_whole),breaks=0:20,right = FALSE)

diff=na.omit(diff)
diffsig=diff[diff$absolutedifference>=0.2,]
diffsig$logabn=factor(diffsig$logabn,
                      levels = levels(diff$logabn))
diff$significant='all'
diffsig$significant='>0.2'
difftotal=rbind(diff,diffsig)

label.n <- function(x){
  return(c(y = ifelse(sum(x)<0,1.5,2), label = length(x)))
}
G=ggplot(Gdatatotal,aes(x=logabn,y=scaleddiff,color=sign))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun.data = label.n, geom = "text", fun.y = median,size=3,position = position_dodge(width = 0.9))+
  ylab('Difference')+
  xlab('Binned log2 abundance in whole sample')+
  coord_cartesian(ylim=c(-3,2.1))+ 
  theme(legend.position = "none")+
  facet_grid(rows=vars(significant))

##Supplementary figure 5J,K,L
#See python/transcript_coverage_pyBigWig.py