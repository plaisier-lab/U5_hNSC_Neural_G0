##########################################################
## ccAF:  sensitivityAnalysis.py                        ##
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
from os.path import exists
import numpy as np
import pandas as pd
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

# Custom classes for classification
import classifiersV3 as cl

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

# Load up loom file
print('\nLoading data...')
ccAF1_scanpy = sc.read_loom('data/cellcycle_int_integrated.loom')

# Parameters for CV
numSamples = 1000
nfolds = 100

# Initialize helper vars/indices for subsetting data (train/test)
nCells = ccAF1_scanpy.shape[0]
allInds = np.arange(0, nCells)

# Precompute testing and training data sets
trainInds = []
truelab = []
testInds = []
mapMe_ACTINN = []
errorACTINN = []
for k in range(nfolds):
    samp1 = sample_without_replacement(nCells, numSamples, random_state = 1234 + k)
    testInds.append(samp1)
    truelab.extend(ccAF1_scanpy.obs['clusts_named'][samp1])
    trainInds.append(np.setdiff1d(allInds, samp1))
    mapMe_ACTINN.append([k, np.setdiff1d(allInds, samp1), samp1, ccAF1_scanpy])

# Subset dataset with marker genes
mg1 = list(set(list(pd.read_csv('data/highlyVarGenes_WT_sgTAOK1_1584.csv', header = 0, index_col = 0).loc[:,'x'])).intersection(ccAF1_scanpy.var_names))
ccAF1_scanpy_mg1 = ccAF1_scanpy[:,mg1]

# Percentage of genes to test for cross-validation classification
percs = [1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

#################
### ACTINN CV ###
#################

# Cross validation within dataset
pred = dict(zip(percs,[[] for i in range(len(percs))]))
print('\nACTINN cross-validation...')
for k in range(nfolds):
    print('ACTINN round '+str(k)+'...')
    ACTINN = cl.Classifier_ACTINN(train = ccAF1_scanpy_mg1[trainInds[k],:], label = 'clusts_named')
    for perc1 in percs:
        samp2 = sample_without_replacement(len(ACTINN.genes), int(len(ACTINN.genes)*perc1), random_state = 1234 + k)
        samp3 = pd.Series(ACTINN.genes).iloc[samp2]
        testPredLbls = ACTINN.predict_labels(ccAF1_scanpy_mg1[testInds[k],samp3])
        pred[perc1].extend(testPredLbls)

# Make folder to store downstream results
if not os.path.exists('results/SensitivityAnalysis'):
    os.makedirs('results/SensitivityAnalysis')
    print ("Directory created")
else:
    print("Directory already exists")

# Dataframe of true labels, predictions, probabilities for all iterations
DF = pd.concat([pd.DataFrame({'True Labels':truelab}), pd.DataFrame(pred)], axis=1)[['True Labels']+percs]
DF.to_csv('results/SensitivityAnalysis/ccAF_CV_sensitivity_analysis.csv')

## If import from csv: run sensitivityAnalysis_plot

for i in range(0,len(percs)):
    comparison_column_i = np.where(DF['True Labels'] == DF[percs[i]], True, False)
    DF['Equal_'+str(percs[i])] = comparison_column_i
    errorACTINN.append(DF['Equal_'+str(percs[i])].value_counts(normalize=True))
    #errorACTINN[percs[i]] = 1-sum(DF['True Labels']==DF[percs[i]])/DF.shape[0]

aggErrorACTINN = []
for i in range(0,len(percs)):
    comparison_column_i = np.where(DF['True Labels'] == DF[percs[i]], True, False)
    DF['Equal_'+str(percs[i])] = comparison_column_i
    tmp = []
    for k in np.arange(0,100001, 1000)[1:]:
        tmp.append(DF.loc[(k-1000):k,'Equal_'+str(percs[i])].value_counts(normalize=True))
    tmpDF = pd.DataFrame(tmp)
    tmpDF['Percent'] = percs[i]
    aggErrorACTINN.append(tmpDF)

aggErrorDF = pd.concat(aggErrorACTINN, axis=0)
aggErrorDF.columns = aggErrorDF.columns.map(str)

# Plot boxplot
aggErrorDF.rename(columns={'False':'Error Rate', 'Percent': 'Percent Genes'}, inplace=True)
aggErrorDF = aggErrorDF.drop(index = 'Equal_0.95')
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
b = sns.boxplot(x="Percent Genes", y = "Error Rate", data=aggErrorDF, palette = sns.light_palette("red", 11))
b.axes.set_title("Cross Validation Sensitivity Analysis",fontsize=25)
b.set_xlabel("Percent Genes",fontsize=25)
b.set_ylabel("Error Rate",fontsize=25)
b.tick_params(labelsize=18)
plt.savefig('results/SensitivityAnalysis/ccAF_CV_sensitivity_analysis_boxplot.pdf')
plt.clf()
