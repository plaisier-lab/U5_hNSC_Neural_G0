##########################################################
## ccAF:  cvClassification_FullAnalysis.py              ##
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
from os.path import exists
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc
from scvi.dataset import LoomDataset, DownloadableAnnDataset

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from clusim.clustering import Clustering, print_clustering
import clusim.sim as sim

# Custom classes for classification
import classifiersV3 as cl

# Load up loom file
print('\nLoading data...')
ccAF1_scanpy = sc.read_loom('data/cellcycle_int_integrated_V3_6_18_2020.loom')

# Load up conversion file
symEnsmbl = pd.read_csv('data/U5/filtered_gene_bc_matrices/hg19/genes.tsv', header = None, index_col=1, sep='\t')
tmp1 = pd.Series(symEnsmbl.index)
tmp1.loc[symEnsmbl.index.duplicated()] = [i+'.1' for i in symEnsmbl.loc[symEnsmbl.index.duplicated()].index]
symEnsmbl.index = pd.Index(tmp1)

# Convert ccAF1
ccAF1_scanpy.var_names = pd.Index(symEnsmbl.loc[ccAF1_scanpy.var_names,0], name='Ensembl')
ccAF1_scanpy.var_names_make_unique()

# HGNC -> downlaoded from HGNC website (https://www.genenames.org/download/custom/)
hgncEnsmbl = pd.read_csv('data/Whitfield/geneConversions/hgnc_geneSymbols_ensmbl.txt', index_col=1, header=0, sep='\t')
hgncEnsmbl = hgncEnsmbl.loc[~hgncEnsmbl['Ensembl ID(supplied by Ensembl)'].isnull()]

ensmblHgnc = pd.Series(hgncEnsmbl.index)
ensmblHgnc.index = list(hgncEnsmbl['Ensembl ID(supplied by Ensembl)'])

hgncPrevEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Previous symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Previous symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncPrevEnsmbl[j] = ensmbl

hgncAliasEnsmbl = {}
for i in hgncEnsmbl.loc[~hgncEnsmbl['Alias symbols'].isnull()].index:
    splitUp = hgncEnsmbl.loc[i,'Alias symbols'].split(', ')
    ensmbl = hgncEnsmbl.loc[i,'Ensembl ID(supplied by Ensembl)']
    for j in splitUp:
        hgncAliasEnsmbl[j] = ensmbl

#ccAF1_scanpy = ccAF1_scanpy[:,[True if i in hgncEnsmbl.index or i in hgncPrevEnsmbl or i in hgncAliasEnsmbl else False for i in ccAF1_scanpy.var_names]]
#ccAF1_scanpy.var_names = pd.Index([str(int(i)) for i in converted])

tmp1 = set(list(pd.read_csv('data/markerGenes/highlyVarGenes_WT_sgTAOK1_1584.csv', header = 0, index_col = 0).loc[:,'x']))
tmp2 = []
missed = []
for j in tmp1:
    if j in hgncEnsmbl.index:
        tmp2.append(hgncEnsmbl.loc[j,'Ensembl ID(supplied by Ensembl)'])
    elif j in hgncPrevEnsmbl:
        tmp2.append(hgncPrevEnsmbl[j])
    elif j in hgncAliasEnsmbl:
        tmp2.append(hgncAliasEnsmbl[j])
    else:
        missed.append(j)

mg1 = tmp2 # Only 1547, only garbage genes are lost that cannot be mapped anymore

sum([1 for i in mg1 if i in ccAF1_scanpy.var_names]) # We can map 1536 geness between mg1 and ccAF1_scanpy

# Subset with classifier genes
mg1 = list(set(mg1).intersection(ccAF1_scanpy.var_names))
ccAF1_scanpy_mg1 = ccAF1_scanpy[:,mg1]

# Parameters for CV
numSamples = 1000
nfolds = 100
ncores = 10

# Initialize helper vars/indices for subsetting data (train/test)
nCells = ccAF1_scanpy.shape[0]
allInds = np.arange(0, nCells)

# Make folder to store downstream results
for meth1 in ['SVMrej', 'RFpy', 'KNN', 'ACTINN']:
    if not os.path.exists('results/'+meth1):
        os.makedirs('results/'+meth1)
        print ("Directory created")
    else:
        print("Directory already exists")

# Precompute testing and training data sets
trainInds = []
truelab = []
testInds = []
mapMe_SVMrejRF = []
mapMe_KNN = []
mapMe_ACTINN = []
errorACTINN = []
for k in range(nfolds):
    samp1 = sample_without_replacement(nCells, numSamples, random_state = 1234 + k)
    testInds.append(samp1)
    truelab.extend(ccAF1_scanpy.obs['clusts_named'][samp1])
    trainInds.append(np.setdiff1d(allInds, samp1))
    mapMe_SVMrejRF.append([k, ccAF1_scanpy[trainInds[k],:], ccAF1_scanpy.obs['clusts_named'][trainInds[k]], ccAF1_scanpy[testInds[k],:]])
    mapMe_KNN.append([k, np.setdiff1d(allInds, samp1), samp1, ccAF1_scanpy])
    mapMe_ACTINN.append([k, np.setdiff1d(allInds, samp1), samp1, ccAF1_scanpy])


#################
### SVMrej CV ###
#################

if not exists('results/SVMrej/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv'):
    # SVMrej multiprocessable function
    def runSVMrej(params):
        print('SVMrej round '+str(params[0])+'...')
        svmRej = cl.Classifier_SVMrej(params[1], params[2])
        testPredLbls = svmRej.predict_labels(params[3])
        return [params[0], testPredLbls]

    # Cross validation within dataset
    print('\nSVMrej cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(runSVMrej, mapMe_SVMrejRF)

    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])

    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('results/SVMrej/ccAF_CV_results_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['G1','G1/other','G2/M','Late G1','M/Early G1','Neural G0','S', 'S/G2']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'SVMrej'
    performDF.to_csv('results/SVMrej/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')

#############
### RF CV ###
#############

if not exists('results/RFpy/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv'):
    # RF multiprocessable function
    def runRF(params):
        print('RF round '+str(params[0])+'...')
        RF = cl.Classifier_RF(params[1], params[2])
        testPredLbls = RF.predict_labels(params[3])
        return [params[0], testPredLbls]

    # Cross validation within dataset
    print('\nRF cross-validation (k-fold = '+str(nfolds)+')...')
    with Pool(ncores) as p:
        res1 = p.map(runRF, mapMe_SVMrejRF)

    # Reassemble predictions
    res1 = dict(res1)
    pred = []
    for k in range(nfolds):
        pred.extend(res1[k])

    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('results/RFpy/ccAF_CV_results_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['G1','G1/other','G2/M','Late G1','M/Early G1','Neural G0','S', 'S/G2']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'RFpy'
    performDF.to_csv('results/RFpy/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')


##############
### KNN CV ###
##############

if not exists('results/KNN/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv'):
    # Cross validation within dataset
    pred = []
    print('\nKNN cross-validation...')
    for k in range(nfolds):
        print('KNN round '+str(k)+'...')
        KNN = cl.Classifier_KNN(ccAF1_scanpy_mg1[trainInds[k],:], 'clusts_named')
        testPredLbls = KNN.predict_labels(ccAF1_scanpy_mg1[testInds[k],:])
        pred.extend(testPredLbls)

    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('results/KNN/ccAF_CV_results_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')

    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))

    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['G1','G1/other','G2/M','Late G1','M/Early G1','Neural G0','S', 'S/G2']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'KNN'
    performDF.to_csv('results/KNN/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')

#################
### ACTINN CV ###
#################
if not exists('results/ACTINN/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv'):
    # Cross validation within dataset
    pred = []
    print('\nACTINN cross-validation...')
    for k in range(nfolds):
        print('ACTINN round '+str(k)+'...')
        ACTINN = cl.Classifier_ACTINN(train = ccAF1_scanpy_mg1[trainInds[k],:], label = 'clusts_named')
        with open('method/ACTINN/ccAF_'+ str(k)+ '.pkl', 'wb') as pklFile:
            pickle.dump(ACTINN, pklFile)
        testPredLbls = ACTINN.predict_labels(ccAF1_scanpy_mg1[testInds[k],:])
        pred.extend(testPredLbls)
    # Dataframe of true labels, predictions, probabilities for all iterations
    DF = pd.DataFrame({'True Labels':truelab, 'Predictions':pred})
    DF.to_csv('results/ACTINN/ccAF_CV_results_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')
    # Get classification report for each iteration
    performanceResults = []
    for k in range(nfolds):
        performanceResults.append(classification_report(truelab[slice((numSamples)*k, (numSamples)*(k+1), 1)], pred[slice((numSamples)*k, (numSamples)*(k+1), 1)], output_dict=True, zero_division=0))
    # Convert into a dataframe
    performDF = pd.concat([pd.DataFrame(i) for i in performanceResults], axis=1).T
    states1 = ['G1','G1/other','G2/M','Late G1','M/Early G1','Neural G0','S', 'S/G2']
    performDF = performDF.loc[[True if i in states1 else False for i in list(performDF.index)]]
    performDF['Classifier'] = 'ACTINN'
    performDF.to_csv('results/ACTINN/CV_classification_report_'+str(ccAF1_scanpy_mg1._n_vars)+'.csv')
    comparison_column = np.where(DF['True Labels'] == DF['Predictions'], True, False)
    DF["Equal"] = comparison_column
    errorACTINN.append(DF['Equal'].value_counts(normalize=True))

if not exists('results/ACTINN/ccAF_'+str(ccAF1_scanpy_mg1._n_vars)+'.pkl'):
    ACTINN = cl.Classifier_ACTINN(train = ccAF1_scanpy_mg1, label = 'clusts_named')
    with open('results/ACTINN/ccAF_'+ str(ccAF1_scanpy_mg1._n_vars)+ '.pkl', 'wb') as pklFile:
        pickle.dump(ACTINN, pklFile)

    with open('results/ACTINN/ccAF_'+ str(ccAF1_scanpy_mg1._n_vars)+ '_smaller.pkl', 'wb') as pklFile:
        ACTINN.train = []
        pickle.dump(ACTINN, pklFile)


print(errorACTINN)
