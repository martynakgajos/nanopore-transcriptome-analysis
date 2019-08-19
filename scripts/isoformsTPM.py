#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

dr="/project/owlmayerTemporary/Martyna/Test/Results/"
df=pd.read_csv(dr+'Quantification/all_counts.txt')

plt.ioff()
ag=[]
c=["#25224c","#c7890f","#b4242f","#444036"]
dr='/project/owlmayer/Olga/spliceosome_components_lists/ENSIDed_lists.txt/20190802_2nd_round_lists/'
for fn in os.listdir(dr)+['all.txt']:
    dr='/project/Neurodifferentiation_System/'
    fn='BAFgenes.txt'
    if fn.endswith('.txt'):
        if fn!='all.txt':
            genes = [line.rstrip('\n').split('.')[0] for line in open(dr+fn)]
            ag=ag+genes
        else:
            genes=list(set(ag))
        with PdfPages('/home/gajos/Presentation/ForOlga/TCEA.pdf') as pdf:
            for i in range(len(genes)):
                #data=df[df['id']==genes[i]]
                data=df[df['gene']==genes[i]]
                if not data.shape[0]:
                    continue
                #data=data[data["type"]=='protein_coding']
                data['maximum']=data.filter(like="day").max(axis=1)
                data=data.sort_values('maximum',ascending=False)
                means=data.filter(like="day")
                labels=data.index
                fig, ax = plt.subplots(figsize=(5,5)) 
                cur_axes = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                z=min(data.shape[0],4)
                for j in range(z):
                    plt.plot(np.array([0,3,5])+np.random.normal(0, 0.03, 3),means.iloc[j,:],label=labels[j],color=c[j],linewidth=2)
                plt.xticks(np.array([0,3,5]),['0','3','5'])
                plt.xlabel('Day')
                plt.ylabel('TPM')
                plt.title(data['gene'].iloc[0])      
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
                plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.4)
                pdf.savefig()