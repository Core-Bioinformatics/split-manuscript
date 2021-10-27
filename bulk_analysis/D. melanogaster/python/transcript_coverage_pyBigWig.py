#calculate average in whole and combined split in each window across transcripts

import pandas as pd
import numpy as np
import pyBigWig

bwwhole = pyBigWig.open(bigwig_name_whole)
bwsplit = pyBigWig.open(bigwig_name_split)
#genetable is a table of transcripts with chromosome, start and end for each, generated using noisyr::cast_gtf_to_genes.R on the reference gtf

abntable=pd.DataFrame({'id':['AA','BB'],'overallstart':[0,0],'start':[0,0],'end':[0,0],'numwindows':[0,0],'avgabn_whole':[0,0],'avgabn_split':[0,0]})
for i in range(len(genetable)):
  #get information for current transcript
  chr=genetable.iloc[i]['seqid']
  start=genetable.iloc[i]['start']
  end=genetable.iloc[i]['end']
  id=genetable.iloc[i]['gene_id']
  #identify windows
  sequence=list(range(start,end,step//2))
  listofstarts=sequence[0:-2]
  for window in listofstarts:
    whole_dist=bwwhole.values(chr,int(window),int(window+100))
    split_dist=bwsplit.values(chr,int(window),int(window+100))
    vector={'id':id,'overallstart':start,'start':window,'end':window+100,'numwindows':len(listofstarts),'avgabn_whole':np.mean(whole_dist),'avgabn_split':np.mean(whole_dist)}
    abntable=abntable.append(vector,ignore_index=True)

abntable=abntable[2:]
#output is abntable with columns id (gene id), overallstart (transcript start), start (window start), end (window end), numwindows (number of windows in transcript), avgabn_whole (average abundance of whole on window), abnabn_split (average abundance of combined split on window)

#distribution of an transcript
chr='2L'
start=9904000
end=9908000
id='FBgn0032156'
bwwhole = pyBigWig.open(bigwig_name_whole)
bwcombined = pyBigWig.open(bigwig_name_combined)
bwsplit1 = pyBigWig.open(bigwig_name_split1)
bwsplit2 = pyBigWig.open(bigwig_name_split2)
whole_dist=bwwhole.values(chr,int(start),int(end))
combined_dist=bwcombined.values(chr,int(start),int(end))
split1_dist=bwsplit1.values(chr,int(start),int(end))
split2_dist=bwsplit2.values(chr,int(start),int(end))
