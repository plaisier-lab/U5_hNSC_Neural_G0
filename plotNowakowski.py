##########################################################
## ccAF:  plotNowakowski.py                             ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier, Samantha O'Connor          ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.use('Agg')
matplotlib.style.use('ggplot')
import webcolors
from scipy.spatial import KDTree

def confusionMatrix(a1, b1):
    v1 = set(a1)
    v2 = set(b1)
    dfTmp = pd.DataFrame([pd.Series(list(a1)),pd.Series(list(b1))]).T
    df1 = pd.DataFrame(index=v1, columns = v2)
    for i in v1:
        for j in v2:
            df1.loc[i,j] = dfTmp.loc[(dfTmp[0]==i) & (dfTmp[1]==j)].shape[0]
    return df1

# lets populate some names into spatial name database
hexnames = webcolors.CSS3_HEX_TO_NAMES
names = []
positions = []

for hex, name in hexnames.items():
    names.append(name)
    positions.append(webcolors.hex_to_rgb(hex))

spacedb = KDTree(positions)

result = []
# query nearest point
for querycolor in [(123,175,65),(201,149,43),(243,118,110),(31,189,194),(166,129,186),(225,109,170),(43,181,103),(74,161,217)]:
    dist, index = spacedb.query(querycolor)
    result.append(names[index])

cmap1 = dict(zip(['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'],result))

# Load ccAF predicted cell cycle states
df1 = pd.read_csv('results/ccAF_results_Nowakowski_norm.csv', header=0, index_col=0)
df1 = df1.loc[~df1['WGCNAcluster_restricted'].isin(['U','none','Glyc'])]
df1 = df1.dropna()
tmp = confusionMatrix(df1['Predictions'],df1['WGCNAcluster_restricted']).loc[['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1']]
df = tmp.div(tmp.sum(),axis=1)
df = df.T
df = df.sort_values('Neural G0',ascending=False)

fig, ax = plt.subplots(figsize=(10,3))  

states = ['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1']
margin_bottom = np.zeros(df.shape[0])
for state1 in states:
    values = list(df[state1])
    
    df[state1].plot.bar(x='Cell Cycle State',y='Percent of Cells', ax=ax, stacked=True, 
                                    bottom = margin_bottom, color=cmap1[state1], label=state1)
    margin_bottom += values

fig.tight_layout()
plt.savefig('Nowakowski_ccAF_plot_Refined.pdf')
plt.clf()

