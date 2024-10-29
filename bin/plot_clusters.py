#!/usr/bin/env python
import sys
import math
import random
import numpy as np
import scipy.stats
from scipy.spatial import distance
from collections import OrderedDict
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns
## usage
# 	$./plot_clusters.py outdir
#		outdir: the path to the result output directory of TimesVector
def printProgressBar (iteration, total, prefix='', suffix='', decimals=1, bar_length=50):
	str_format = "{0:." + str(decimals) + "f}"
	percents = str_format.format(100 * (iteration / float(total)))
	filled_length = int(round(bar_length * iteration / float(total)))
	bar = '█' * filled_length + '-' * (bar_length - filled_length)
	sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
	if iteration == total:
		sys.stdout.write('\n')
	sys.stdout.flush()
# result output directory of TimesVector
outdir=sys.argv[1]
# color map
cmaps = ['crimson', 'royalblue', 'limegreen', 'gold', 'sienna', 'teal', 'navy', 'magenta', 'kahki', 'coral', 'sienna']
# cmaps=["darkcyan", "blue", "crimson", "orange", "yellow", "saddlebrown", "turquoise", "teal", "navy", "kahki", "coral", "sienna"]
### 
### Plotting cluster centroids
###
print("Plotting cluster centroids...")
df_DEP=pd.read_csv("%s/plots/DEP_clusters_plot.dat"%(outdir), sep='\t')
df_SEP=pd.read_csv("%s/plots/SEP_clusters_plot.dat"%(outdir), sep='\t')
df=df_DEP.append(df_SEP)
dat=df.loc[df['Clid'] == 1]
pheno=pd.Series(df['Samples']).unique()
timepoints=pd.Series(df['TimePoint']).unique()
colors=dict(zip(pheno, cmaps[:len(pheno)]))
clust_pdf = matplotlib.backends.backend_pdf.PdfPages("%s/Clusters.pdf"%(outdir))
fs=10	# fontsize
row_n=3
col_n=3
pos=0
fig, ax=plt.subplots(figsize=(12, 10), nrows=row_n, ncols=col_n)
clusters=pd.unique(df['Clid'])
printProgressBar(0, len(clusters), prefix = 'Progress:', suffix = 'Complete')
for cidx, cid in enumerate(clusters):
	if pos//col_n==col_n:
		fig.tight_layout()
		clust_pdf.savefig(fig)
		fig.clf()
		fig, ax=plt.subplots(figsize=(12, 10), nrows=row_n, ncols=col_n)
		pos=0
	row=pos//col_n
	col=pos%row_n
	dat=df.loc[df['Clid'] == cid]
	ctype=np.unique(dat['Type'])[0]
	for pdx, p in enumerate(pheno):
		p_dat=dat.loc[dat['Samples']==p]
		ax[row][col].plot(p_dat['TimePoint'], p_dat['Value'], color=colors[p], linewidth=3, label=p)
		
	# plot annotations
	ax[row][col].set_title("C%d %s"%(cid, ctype), size=10)
	if pos==0:
		ax[row][col].set_ylabel("Normalized expression", fontsize=fs)
		ax[row][col].set_xlabel("Time points", fontsize=fs)
		handles, labels = ax[row][col].get_legend_handles_labels()
		by_label = OrderedDict(zip(labels, handles))
		leg=ax[row][col].legend(by_label.values(), by_label.keys(), loc="upper right", fancybox=False, fontsize=fs, edgecolor='black')
	pos+=1
	# ylim=min(1, list(p_dat['Max'])[0]+.01)
	# ax[row][col].set_ylim(0, ylim)
	printProgressBar(cidx + 1, len(clusters), prefix = 'Progress:', suffix = 'Complete')
fig.tight_layout()
clust_pdf.savefig(fig)
clust_pdf.close()
### 
### Plotting cluster genes with normalized expression
###
print("Plotting cluster genes (normalized)...")
df_DEP=pd.read_csv("%s/plots/DEP_genes_norm_plot.dat"%(outdir), sep='\t')
df_SEP=pd.read_csv("%s/plots/SEP_genes_norm_plot.dat"%(outdir), sep='\t')
df=df_DEP.append(df_SEP)
pheno=pd.Series(df['Samples']).unique()
timepoints=pd.Series(df['TimePoint']).unique()
colors=dict(zip(pheno, cmaps[:len(pheno)]))
geneclust_pdf = matplotlib.backends.backend_pdf.PdfPages("%s/Gene_clusters.pdf"%(outdir))
row_n=3
col_n=3
pos=0
fig, ax=plt.subplots(figsize=(12, 10), nrows=row_n, ncols=col_n)
clusters=pd.unique(df['Clid'])
printProgressBar(0, len(clusters), prefix = 'Progress:', suffix = 'Complete')
for cidx, cid in enumerate(clusters):
	if pos//col_n==col_n:
		fig.tight_layout()
		geneclust_pdf.savefig(fig)
		fig.clf()
		fig, ax=plt.subplots(figsize=(12, 10), nrows=row_n, ncols=col_n)
		pos=0
	row=pos//col_n
	col=pos%row_n
	dat=df.loc[df['Clid'] == cid]
	ctype=np.unique(dat['Type'])[0]
	genes_n=np.unique(dat['GeneID']).size//len(pheno)
	genes=np.unique(dat['GeneID'])
	for p in pheno:
		plot_dat=[]
		for g in genes:
			g_dat=dat.loc[(dat['GeneID']==g) & (dat['Samples']==p)]
			plot_dat.append(tuple(g_dat['TimePoint']))
			plot_dat.append(tuple(g_dat['Value']))			
		ax[row][col].plot(*plot_dat, color=colors[p], linewidth=0.5, alpha=.4, label=p)
		
	# plot annotations
	ax[row][col].set_title("C%d %s (%d genes)"%(cid, ctype, genes_n), size=10)
	if pos==0:
		ax[row][col].set_ylabel("Normalized expression", fontsize=fs)
		ax[row][col].set_xlabel("Time points", fontsize=fs)
		handles, labels = ax[row][col].get_legend_handles_labels()
		by_label = OrderedDict(zip(labels, handles))
		leg= ax[row][col].legend(by_label.values(), by_label.keys(), loc="upper right", fancybox=False, fontsize=fs, edgecolor='black')
		# leg = ax[row][col].legend()
		for line in leg.get_lines():
			line.set_linewidth(2.0)
			line.set_alpha(1)
	pos+=1
	printProgressBar(cidx + 1, len(clusters), prefix = 'Progress:', suffix = 'Complete')
fig.tight_layout()
geneclust_pdf.savefig(fig)
geneclust_pdf.close()
