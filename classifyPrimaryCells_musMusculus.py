##########################################################
## ccAF:  classifyPrimaryCells_musMusculus.py           ##
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
from scipy.sparse import csr_matrix
from copy import deepcopy
from multiprocessing import Pool
import pickle

# Single Cell Packages
import scanpy as sc

# Cross-validation and metrics
from sklearn.utils.random import sample_without_replacement
from sklearn.metrics import classification_report
from sklearn.preprocessing import scale, robust_scale, quantile_transform

# Custom classes for classification
import classifiersV3 as cl

def confusionMatrix(a1, b1):
    v1 = set(a1)
    v2 = set(b1)
    dfTmp = pd.DataFrame([pd.Series(list(a1)),pd.Series(list(b1))]).T
    df1 = pd.DataFrame(index=v1, columns = v2)
    for i in v1:
        for j in v2:
            df1.loc[i,j] = dfTmp.loc[(dfTmp[0]==i) & (dfTmp[1]==j)].shape[0]
    return df1

# Load up classifier as ccAF1 -> Genes are human Ensmbl gene IDs
with open('../classifierTester/method/ACTINN/ccAF_1536_smaller.pkl','rb') as pklFile:
    ccAF1 = pickle.load(pklFile)

# Conversion of mouse Ensembl IDs to human Ensembl IDs
hsmmEnsmbl = pd.read_csv('data/geneConversions/mart_export.txt', sep='\t', index_col=1)

# Load up loom file
for set1 in ['GSE67833','PRJNA324289']:
    if not exists('results/table1_'+str(set1)+'.csv'): 
        print('\nLoading '+set1+' data...')
        set1_scanpy = sc.read_loom('data/forClassification/'+set1+'.loom')
        set1_scanpy = set1_scanpy[:,[True if i in hsmmEnsmbl.index else False for i in set1_scanpy.var_names]]
        # Convert to ensmbl
        convertMe_hgnc = pd.DataFrame(np.nan, index=set1_scanpy.var_names, columns=['ensembl'])
        numConverted_hgnc = 0
        missed = []
        for j in set1_scanpy.var_names:
            if j in hsmmEnsmbl.index:
                tmp = hsmmEnsmbl.loc[j,'Gene stable ID']
                if isinstance(tmp, pd.core.series.Series):
                    tmp = list(tmp)[0]
                convertMe_hgnc.loc[j,'ensembl'] = tmp
            else:
                missed.append(j)
        
        set1_scanpy.var_names = pd.Index(convertMe_hgnc['ensembl'])
        
        #############################
        ### ACTINN Classification ###
        #############################
        if not exists('results/CV_classification_report_'+set1+'.csv'):
            # Cross validation within dataset
            print('\nACTINN '+set1+'...')
            set1PredLbls = ccAF1.predict_labels(set1_scanpy)
            set1PredProbs = ccAF1.predict_probs(set1_scanpy)
            set1PredProbs.columns = set1_scanpy.obs_names
            set1_scanpy.obs['ccAF'] = set1PredLbls
            # Dataframe of true labels, predictions, probabilities for all iterations
            DF = pd.concat([set1_scanpy.obs, pd.DataFrame({'Predictions':set1PredLbls},index=set1_scanpy.obs_names),set1PredProbs.T], axis=1)
            DF.to_csv('results/ccAF_results_'+set1+'.csv')
            # Make Table 1 data
            t1 = pd.DataFrame(pd.Series(set1PredLbls).value_counts()).T
            t1.index = ['Counts']
            t2 = confusionMatrix(set1_scanpy.obs['cell_type'], set1PredLbls)
            table1 = pd.concat([t1, t2[list(t1.columns)]])
            table1.to_csv('results/table1_'+str(set1)+'.csv')

