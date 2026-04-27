#code in new env
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

sc.settings.verbosity = 3
sc.logging.print_header()

###################################################################################################
# Load combined dataset 
adata = sc.read_h5ad('WT_KO_DG_combined.h5ad')
print(adata)

###################################################################################################
# Apply Blank gene based filter to remove False discoveries
blank_genes = [gene for gene in adata.var_names if gene.startswith('Blank')]
print(len(blank_genes))
print("\n")
print('----------------------------------------------------------------------------')

total_blank_transcripts_per_cell = np.sum(adata[:,blank_genes].X,axis=1)
average_blank_transcripts_per_cell = np.mean(total_blank_transcripts_per_cell)
print('Average_blank_transcripts_per_cell = ',average_blank_transcripts_per_cell)
min_filter = int(12*average_blank_transcripts_per_cell)
print('Min filter = ',min_filter)
print("\n")
print('----------------------------------------------------------------------------')

###################################################################################################
# Calculate QC metrics after blank gene filtration
sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

print('Max volume', max(adata.obs['volume']))
fig, axs = plt.subplots(1, 3, figsize=(15, 4))

###################################################################################################
# Visualize the QC metrics
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
tx_lower_cutoff = min_filter
tx_upper_cutoff = 4000
axs[0].axvline(tx_lower_cutoff, color='red', linestyle='-')
axs[0].axvline(tx_upper_cutoff, color='red', linestyle='-')

axs[1].set_title("Genes per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
    bins=75
)
genes_cutoff = 5
axs[1].axvline(genes_cutoff, color='red', linestyle='-')

axs[2].set_title("Volume per cell")
sns.histplot(
    adata.obs["volume"],
    kde=False,
    ax=axs[2],
)
vol_lower_cutoff = 100
median_vol = np.median(adata.obs['volume'])
vol_upper_cutoff = 3*median_vol

axs[2].axvline(vol_lower_cutoff, color='red', linestyle='-')
axs[2].axvline(vol_upper_cutoff, color='red', linestyle='-')

plt.savefig('WT_KO_DG_combined_QC_Before_Filter.png',bbox_inches='tight')
print("\n")
print('----------------------------------------------------------------------------')

#sc.pl.violin(adata, keys='volume', ylabel='Volume per cell', groupby=None,stripplot=False)
# Adding horizontal cutoff lines
#plt.axhline(vol_lower_cutoff, color='red', linestyle='-')
#plt.axhline(vol_upper_cutoff, color='red', linestyle='-')

#plt.savefig('vol_dist.png')

###################################################################################################
# Apply blank gene, counts, genes and volume based filtering
print('----------------------------------------------------------------------------------')
print(adata)
print('----------------------------------------------------------------------------------')
adata = adata[:, ~adata.var_names.isin(blank_genes)]
print(adata)
print('----------------------------------------------------------------------------------')
adata = adata[(adata.obs['total_counts']>tx_lower_cutoff) & (adata.obs['total_counts']<tx_upper_cutoff)]
print(adata)
print('----------------------------------------------------------------------------------')
sc.pp.filter_cells(adata, min_genes=genes_cutoff)
print(adata)
print('----------------------------------------------------------------------------------')
adata = adata[(adata.obs['volume']>vol_lower_cutoff) & (adata.obs['volume']<vol_upper_cutoff)]
print(adata)
print('----------------------------------------------------------------------------------')

adata.write_h5ad('filtered_WT_KO_DG_combined.h5ad')

###################################################################################################
# Scrublet - doublet filtering
import scrublet as scr
import scipy.io
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

input_dir = '../filtered_WT_KO_combined.h5ad'
adata = sc.read_h5ad(input_dir)
print(adata)

from scipy.sparse import csc_matrix

# Convert to CSC sparse matrix
counts_matrix = csc_matrix(adata.X)
genes = np.array(adata.var_names)

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

print('----------------------------------------------------------------------') 
print(doublet_scores)
adata.obs['Doublets_predict'] = predicted_doublets
print(adata[adata.obs['Doublets_predict']==True])
scrub.plot_embedding('UMAP', order_points=True)
#plt.savefig('umap_doublets_Combined_cortical.png')

# predicted doublets using Scrublet's automatic threshold
n_doublets = predicted_doublets.sum()
n_cells = len(predicted_doublets)