###########################################################
## ccAF:  scvelo_analysis.py                             ##
##  ______     ______     __  __                         ##
## /\  __ \   /\  ___\   /\ \/\ \                        ##
## \ \  __ \  \ \___  \  \ \ \_\ \                       ##
##  \ \_\ \_\  \/\_____\  \ \_____\                      ##
##   \/_/\/_/   \/_____/   \/_____/                      ##
## @Developed by: Plaisier Lab                           ##
##   (https://plaisierlab.engineering.asu.edu/)          ##
##   Arizona State University                            ##
##   242 ISTB1, 550 E Orange St                          ##
##   Tempe, AZ  85281                                    ##
## @Author:  Chris Plaisier, Samantha O'Connor           ##
## @License:  GNU GPLv3                                  ##
##                                                       ##
## If this program is used in your analysis please       ##
## mention who built it. Thanks. :-)                     ##
###########################################################

# docker run -it -v '/home/cplaisier/Dropbox (ASU):/files' cplaisier/scrna_seq_velocity

import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
from matplotlib.colors import ListedColormap

# Load U5 data from Seurat
adata_Seurat = scv.read_loom('data/U5_hNSC/WT/U5_all_v3.loom')
ccAF1_tsne = pd.read_csv('data/U5_hNSC/WT/tsne_cell_embeddings_Perplexity_26.csv', header = 0, index_col = 0)
ccAF1_tsne = ccAF1_tsne.loc[[True if not i.find('-1_1')==-1 else False for i in ccAF1_tsne.index]]
ccAF1_tsne.index = [i.rstrip('-1_1') for i in list(ccAF1_tsne.index)]
adata_Seurat.obsm['tsne_cell_embeddings'] = np.array(ccAF1_tsne.loc[adata_Seurat.obs_names])

# Load U5 data from velocyto
adata_vc = scv.read_loom('data/U5_hNSC/WT/U5_all.loom')

# Merge data
adata = scv.utils.merge(adata_vc, adata_Seurat)
adata.obs.clusts_WT = adata.obs.clusts_WT.astype('category')
adata.obs.clusts_WT = adata.obs.clusts_WT.cat.rename_categories({'G1/Differentiation':'Neural G0'})
adata.obs.clusts_WT = adata.obs.clusts_WT.cat.reorder_categories(['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'], ordered=True)

# Calcualte moments
scv.pp.moments(adata, n_pcs=19, n_neighbors=10, use_rep='pca_cell_embeddings')

# Compute velocity
ccAdata = adata[list(adata.obs[adata.obs['clusts_WT'].isin(['Neural G0','G1','Late G1','S','S/G2','G2/M','M/Early G1'])].index)]
scv.pp.moments(ccAdata, n_pcs=19, n_neighbors=15, use_rep='pca_cell_embeddings')
scv.tl.velocity(ccAdata, groupby='clusts_WT', groups=['Neural G0', 'G1','Late G1','S','S/G2','G2/M','M/Early G1'], groups_for_fit=['Neural G0', 'G1','Late G1','S','S/G2','G2/M','M/Early G1'])
scv.tl.velocity_graph(ccAdata)

# Get same colors for each cell cycle phase
import webcolors
from scipy.spatial import KDTree
hexnames = webcolors.CSS3_HEX_TO_NAMES
names = []
positions = []
for hex, name in hexnames.items():
    names.append(name)
    positions.append(webcolors.hex_to_rgb(hex))
spacedb = KDTree(positions)
result = []
for querycolor in [(123,175,65),(201,149,43),(243,118,110),(31,189,194),(166,129,186),(225,109,170),(43,181,103),(74,161,217)]:
    dist, index = spacedb.query(querycolor)
    result.append(names[index])
cmap1 = dict(zip(['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'],result))

# Plot tSNE with stream lines
scv.pl.velocity_embedding_stream(ccAdata, basis='tsne_cell_embeddings', color=['clusts_WT'], save='ccAdata_velocity_stream_tsne.png', dpi=300, palette = cmap1, title='', legend_loc='right margin', figsize=(3,3), size=70)

