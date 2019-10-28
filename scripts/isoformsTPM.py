#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hier
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def getarg():
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-i", "--inFile", help = "whole path", required = True)
    parser.add_argument("-o", "--outFile", help = "whole path", required = True)
    parser.add_argument("-g", "--geneList", help = "path to a file with ensemble ids; one per line", required = True)
    arg = parser.parse_args()
    return arg.inFile, arg.outFile, arg.geneList

fin,fout,glist=getarg()
df=pd.read_csv(fin)
df=df.fillna('0.0')
df["gene_number"]=df.apply(lambda x: x['gene_id'].split('.')[0],1)
df['gene_id']=df["gene_number"]
genes = [line.rstrip('\n') for line in open(glist)]
plt.ioff()

c=["#25224c","#c7890f","#b4242f","#444036"]
   
df=df[df['gene_id'].isin(genes)]
df=df.sort_values(by='gene_name')
with PdfPages(fout) as pdf:
    data=df.groupby('gene_name').sum().filter(like='day')
    data = data.sort_index(1)
    linkageMatrix = hier.linkage(data.values,method='average',metric='correlation')
    plt.figure()
    dendro = hier.dendrogram(linkageMatrix, labels=list(data.index))#no_labels=True)#,truncate_mode='level')#,labels=list(data.index))
    plt.close()
    Z=data.transform(lambda x: (x - x.mean()) / x.std(),axis=1)
    plt.figure(figsize=(6,data.shape[0]*0.2+2))
    ax=sns.heatmap(Z.iloc[dendro['leaves']], cmap=plt.get_cmap('coolwarm'),cbar=True,yticklabels=True)#,vmax=np.percentile(data.values,95),vmin=np.percentile(data.values,5))
    plt.ylabel('')
    plt.subplots_adjust(left=0.25, right=0.9, top=0.95, bottom=0.1)
    pdf.savefig()
    data=data.sort_index()
    for i in range(len(genes)):
        data=df[df['gene_id']==genes[i]]
        if not data.shape[0]:
            continue
        fig, ax = plt.subplots(figsize=(5,5))
        cur_axes = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.plot([0,3,5],data.filter(like='day').sum(),color='black',linewidth=2)
        plt.xticks([0,3,5],['0','3','5'])
        plt.xlabel('Day')
        plt.ylabel('Normalized read count')
        plt.title(data['gene_name'].iloc[0])
        a=plt.axis()
        yl=plt.yticks()
        plt.xlim([-0.3,5.3])
        plt.ylim([max(-5,-0.01*a[3]),a[3]])  
        yl=plt.yticks()
        ax.spines['left'].set_bounds(max(-5,-0.01*a[3]),yl[0][yl[0]<a[3]][-1])
        ax.spines['bottom'].set_bounds(-0.3,5)
        plt.tight_layout()
        pdf.savefig()
            
        data['maximum']=data.filter(like="day").max(axis=1)
        data=data.sort_values('maximum',ascending=False)
        means=data.filter(like="day")
        labels=data['transcript_id']
        fig, ax = plt.subplots(figsize=(5,5)) 
        cur_axes = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        z=min(data.shape[0],4)
        for j in range(z):
            plt.plot(np.array([0,3,5])+np.random.normal(0, 0.03, 3),means.iloc[j,:],label=labels.iloc[j],color=c[j],linewidth=2)
        plt.xticks(np.array([0,3,5]),['0','3','5'])
        plt.xlabel('Day')
        plt.ylabel('TPM')
        plt.title(data['gene_name'].iloc[0])      
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.98, box.height])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),frameon=False)
        a=plt.axis()
        yl=plt.yticks()
        plt.xlim([-0.3,5.3])
        plt.ylim([max(-5,-0.01*a[3]),a[3]])  
        yl=plt.yticks()
        ax.spines['left'].set_bounds(max(-5,-0.01*a[3]),yl[0][yl[0]<a[3]][-1])
        ax.spines['bottom'].set_bounds(-0.3,5)
        plt.subplots_adjust(left=0.15, right=0.99, top=0.95, bottom=0.4)
        pdf.savefig()