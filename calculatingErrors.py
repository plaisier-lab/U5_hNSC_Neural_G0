##########################################################
## ccAF:  calculatingErrors.py                          ##
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
import pandas as pd
import numpy as np

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load CV results and concatenate for error comparisons
meths = ['SVMrej', 'RFpy', 'KNN', 'ACTINN']
genes = [228, 298, 342, 439, 571, 738, 1584]
errors = pd.DataFrame(index=genes, columns=meths)
for meth1 in meths:
    for gene1 in genes:
        if not gene1 == 1584:
            tmp1 = pd.read_csv('method/'+meth1+'/ccAF_CV_results_'+str(gene1)+'.csv', index_col=0, header=0)
        else:
            tmp1 = pd.read_csv('method/'+meth1+'/ccAF_CV_results_'+str(gene1)+'_2.csv', index_col=0, header=0)
        errors.loc[gene1,meth1] = 1-sum(tmp1['True Labels']==tmp1['Predictions'])/tmp1.shape[0]
errors.to_csv('errors_cross_validation.csv')


# Load CV results and concatenate for whitfield error comparisons
meths = ['SVMrej','RFpy','KNN','NN']
errors_noOther = pd.DataFrame(columns=meths)
data = ['no_manip', 'pow10', 'pow2', 'FC', '1334', 'pow10_1334', 'pow2_1334', 'FC_1334']
mods = ['none', 'scale', 'robust', 'quantile']
errors = pd.DataFrame(index=data, columns=mods)
errors_experimental = pd.DataFrame(index=data, columns=mods)
for mod1 in mods:
    for data1 in data:
        tmp1 = pd.read_csv('method/NN/whitfield/ccAF_CV_results_Whitfield_'+data1+'_pad_'+mod1+'.csv', index_col=0, header=0)
        errors.loc[data1,mod1] = 1-sum(tmp1['maxWhitfield_median_exp']==tmp1['Translated_Predictions'])/tmp1.shape[0]
        errors_experimental.loc[data1,mod1] = 1-sum(tmp1.dropna()['Experimental']==tmp1.dropna()['Translated_Predictions_experimental'])/tmp1.dropna().shape[0]

errors.to_csv('errors_Whitfield.csv')
errors_experimental.to_csv('errors_Whifield_experimental.csv')
