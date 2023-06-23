#!/usr/bin/env python3

import pandas as pd
import sys

def getData(f):
    names=['chrom','position','depth']
    df=pd.read_csv(f,sep='\t',names=names)
    return df

def coverageStats(df):
    cov1=df[df.depth > 1].groupby(['chrom']).count()
    cov10=df[df.depth > 10].groupby(['chrom']).count()

    bases=df.groupby(['chrom'])['depth'].sum()
    chromLens=df.groupby(['chrom'])['position'].count()
    avDepth=bases/cov1['depth']

    cov1.reset_index(inplace=True)
    cov10.reset_index(inplace=True)

    df2=cov1.merge(cov10,on='chrom',suffixes=[' cov1',' cov10'],how='outer')
    df2['length']=df2.chrom.map(chromLens)
    df2['covBreadth1x']=df2['depth cov1']/df2['length']
    df2['covBreadth10x']=df2['depth cov10']/df2['length']
    df2['avDepth']=df2.chrom.map(avDepth)
    df2['id']=sys.argv[2]
    df2['bases']=df2.chrom.map(bases)
    df2=df2[['id','chrom','length','bases','avDepth','position cov1','position cov10','covBreadth1x','covBreadth10x']]
    df2.to_csv('coverage_stats.csv',index=False)
    print(df2)
    

df=getData(sys.argv[1])
coverageStats(df)
