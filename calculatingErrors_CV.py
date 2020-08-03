##########################################################
## ccAF:  calculatingErrors_CV.py                       ##
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
genes = [1536]
errors = pd.DataFrame(index=genes, columns=meths)
for meth1 in meths:
    for gene1 in genes:
        tmp1 = pd.read_csv('results/'+meth1+'/ccAF_CV_results_'+str(gene1)+'.csv', index_col=0, header=0)
        errors.loc[gene1,meth1] = 1-sum(tmp1['True Labels']==tmp1['Predictions'])/tmp1.shape[0]
    
errors.to_csv('results/errors_cross_validation.csv')
