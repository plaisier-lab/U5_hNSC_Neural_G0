##########################################################
## ccAF:  Whitfield_classification_ACTINN_analysis.py   ##
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

# General
import os
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.preprocessing import scale, robust_scale, quantile_transform
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

# Custom classes for classification
import classifiersV3 as cl

# Classifier loading and writing
#import torch

# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Function to translate predictions
def compute_whitfield_experimental(whitfield_pred, conversion = {'M/Early G1':'M','G2/M':'M','S':'S','S/G2':'S','Late G1':'S','G1':'S', 'Unknown':'Unknown', 'Neural G0':'Unknown', 'G1/other':'Unknown'}):
    whit2 = [conversion[i] for i in whitfield_pred]
    return whit2

def compute_whitfield_translation(whitfield_pred, conversion = {'M/Early G1':'M/G1','G2/M':'G2/M','S':'S','S/G2':'G2','Late G1':'S','G1':'G1/S', 'Unknown':'Unknown', 'Neural G0':'G1/S', 'G1/other':'Unknown'}):
    ## Convert predictions to ['M','S']
    whit2 = [conversion[i] for i in whitfield_pred]
    #
    # Compute and return accuracy
    return whit2

# Load up classifier as ccAF1 -> Genes are Ensmbl gene IDs
with open('results/ACTINN/ccAF_1536_smaller.pkl','rb') as pklFile:
    ccAF1 = pickle.load(pklFile)

g2e = pd.read_csv('data/Whitfield/geneConversions/ensembl_entrez.csv', index_col=1, header=0)
g2e = g2e.loc[g2e.index.dropna()]
g2e = g2e.loc[~g2e.index.duplicated()]

# Make folder to store downstream results
if not os.path.exists('results/Whitfield'):
    os.makedirs('results/Whitfield')
    print ("Directory created")
else:
    print("Directory already exists")

# Load up whitfield dataset & incorporate ideal vector data from Figure 1 to meta
whit1 = ['whitfield_dataPlusScores_6_30_2020_', '_1334.csv', '1334']
mod1 = 'quantile'
experiments = {'TT1':[0,12], 'TT2':[12,38], 'TT3':[38,86], 'TN':[86,105], 'SHAKE':[105,114]}
whitfield = {}
print(whit1, mod1)
for exp1 in experiments:
    whitfield[exp1] = sc.read_csv('data/Whitfield/'+whit1[0]+exp1+whit1[1], first_column_names=True).T
    var_names = [str(g2e.loc[float(i),'Gene stable ID']) for i in whitfield[exp1].var_names if float(i) in g2e.index]
    whitfield[exp1] = whitfield[exp1][:,[True if float(i) in g2e.index else False for i in whitfield[exp1].var_names]]
    whitfield[exp1].var_names = pd.Index(var_names)
    if mod1 == 'quantile':
        tmp = whitfield[exp1].X
        whitfield[exp1].X = quantile_transform(tmp, axis = 1)
    whitfield[exp1].var_names = [i.rstrip('.0') for i in whitfield[exp1].var_names]
    whitfield[exp1].var_names_make_unique()

whitfield_mg = pd.read_csv('data/Whitfield/markergenes_ForPlotting.csv', header=0, index_col=0, skiprows=[1,2,3])
whitfield_phase = whitfield_mg['Phase']
whitfield_mg_1 = whitfield_mg.drop('Phase', axis=1)
whitfield_vectors = whitfield_mg_1.iloc[20:25]
maxWhitfield_vector = [i.replace(' phase Vector','').replace(' vector','').replace(' Vector','') for i in whitfield_vectors.idxmax()]
maxWhitfield_vector = [k.replace('_', '/') for k in maxWhitfield_vector]
maxWhitfield_vector_2 = whitfield_mg.iloc[0:20].groupby('Phase').median().idxmax()
maxWhitfield_vector_2 = [k.replace('_', '/') for k in maxWhitfield_vector_2]

# Experimental and exprs3
# Load up experimental meta data
meta = pd.read_csv('data/Whitfield/metaInformation.csv', header=0, index_col=0)
experimental = meta['Cycle']
maxWhitfield_vector_3 = [maxWhitfield_vector_2[i] if not isinstance(meta['Cycle'].iloc[i], str) else {'M':'G2/M','S':'S'}[str(meta['Cycle'].iloc[i])] for i in range(len(list(maxWhitfield_vector_2)))]

# Reconstitute output
testPredLbls = {'ACTINN':[]}

##############
### ACTINN ###
##############

# Classification with ground truth dataset
for exp1 in ['TT1','TT2','TT3','TN','SHAKE']:
    testPredLbls['ACTINN'] += list(ccAF1.predict_labels(whitfield[exp1]))

# Set true labels
truelab = maxWhitfield_vector_3

# Compute scores for each classifier method
pred = compute_whitfield_translation(testPredLbls['ACTINN'])
pred_exp = compute_whitfield_experimental(testPredLbls['ACTINN'])

# Dataframe of true labels, predictions, probabilities for all iterations
DF = pd.DataFrame({'Experimental': meta['Cycle'],'maxWhitfield_mean':maxWhitfield_vector,'maxWhitfield_median': maxWhitfield_vector_2, 'maxWhitfield_median_exp': maxWhitfield_vector_3, 'Translated_Predictions_experimental':pred_exp, 'Translated_Predictions':pred, 'Predictions':testPredLbls['ACTINN']})
DF.to_csv('results/Whitfield/ACTINN_results_Whitfield_'+whit1[2]+'_'+mod1+'.csv')

# Get classification report for each iteration
performanceResults = classification_report(truelab, pred, output_dict=True, zero_division=0)

# Convert into a dataframe
performDF = pd.DataFrame(performanceResults).T
states1 = ['G1/S','G2', 'G2/M', 'M/G1', 'S']
performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
performDF['Classifier'] = 'ACTINN'
performDF.to_csv('results/Whitfield/ACTINN_classification_report_Whitfield_'+whit1[2]+'_'+mod1+'.csv')

# Make folder to store downstream plots
if not os.path.exists('results/Whitfield/plots'):
    os.makedirs('results/Whitfield/plots')
    print ("Directory created")
else:
    print("Directory already exists")

# Plot f1-score boxplot
ACTINN_CR_Whitfield = pd.read_csv('results/Whitfield/ACTINN_classification_report_Whitfield_'+whit1[2]+'_'+mod1+'.csv')
ACTINN_CR_Whitfield.rename(columns={'Unnamed: 0':'Cell Cycle State'}, inplace=True)
# hue = cell cycle state
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.barplot(x="Classifier", y="f1-score", hue = "Cell Cycle State", data=ACTINN_CR_Whitfield)
plt.savefig('results/Whitfield/plots/F1_scores_for_each_state_Whitfield_'+whit1[2]+'_'+mod1+'.pdf')
plt.clf()
# hue = classifier
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.barplot(hue="Classifier", y="f1-score", x = "Cell Cycle State", data=ACTINN_CR_Whitfield)
plt.savefig('results/Whitfield/plots/F1_scores_for_each_state_Whitfield_'+whit1[2]+'_'+mod1+'_2.pdf')
plt.clf()

# Plot heatmap
# Load up ACTINN classifications
ACTINN_Results = pd.read_csv('results/Whitfield/ACTINN_results_Whitfield_'+whit1[2]+'_'+mod1+'.csv')
pred_ACTINN = list(ACTINN_Results['Translated_Predictions'])

# Format heatmap
color_map = {'G1_S': '#fee090', 'S': '#91bfdb', 'G2': '#4574b4', 'G2_M': '#d73027', 'M_G1': '#fc8d59', 'nan': '#ffffff', 'Neural_G0': '#ffffff', 'Unknown':'#ffffff', 'M':'#d73027'}
with PdfPages('results/Whitfield/plots/heatmapsOfWhitfield_'+whit1[2]+'_'+mod1+'.pdf') as pp:
    nn1 = [color_map[str(i).replace('/','_').replace(' ','_')] for i in pred_ACTINN]
    experimental = [color_map[{'M':'M', 'S':'S', 'nan':'nan'}[str(i).replace('/','_').replace(' ','_')]] for i in list(meta['Cycle'])]
    exprs3 = [color_map[str(maxWhitfield_vector_3[i]).replace('/','_').replace(' ','_')] for i in range(len(list(maxWhitfield_vector_3)))]
    col_colors = pd.DataFrame({'Exp.':dict(zip(list(whitfield_vectors.columns),experimental)),  'Exprs.3':dict(zip(list(whitfield_vectors.columns),exprs3)), 'NeuralNet':dict(zip(list(whitfield_vectors.columns),nn1))})
    fig = sns.clustermap(whitfield_mg_1, cmap = matplotlib.colors.LinearSegmentedColormap.from_list('rkg',colors=[(0,0.75,0),(0,0,0),(0.75,0,0)]), cbar_kws={'label':'Expression Ratio'}, row_cluster=False, col_cluster=False, standard_scale=1, figsize=(11,4), xticklabels=True, yticklabels=True, row_colors=[[color_map[i] for i in whitfield_phase]], col_colors=col_colors)
    fig.ax_col_colors.set_yticklabels(fig.ax_col_colors.get_yticklabels(), fontsize=6)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
    for label in ['G1_S','S','G2','G2_M','M_G1']:
        fig.ax_col_dendrogram.bar(0, 0, color=color_map[label], label=label, linewidth=0)
    fig.ax_col_dendrogram.legend(loc="center", ncol=6)
    fig.savefig(pp, format='pdf')
    plt.close()
