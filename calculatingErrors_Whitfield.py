##########################################################
## ccAF:  calculatingErrors_Whitfield.py                ##
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

# Load CV results and concatenate for whitfield error comparisons
data = ['1334']
mods = ['quantile']
errors = pd.DataFrame(index=data, columns=mods)
errors_experimental = pd.DataFrame(index=data, columns=mods)
for mod1 in mods:
    for data1 in data:
        tmp1 = pd.read_csv('results/Whitfield/ACTINN_results_Whitfield_'+data1+'_'+mod1+'.csv', index_col=0, header=0)
        errors.loc[data1,mod1] = 1-sum(tmp1['maxWhitfield_median_exp']==tmp1['Translated_Predictions'])/tmp1.shape[0]
        errors_experimental.loc[data1,mod1] = 1-sum(tmp1.dropna()['Experimental']==tmp1.dropna()['Translated_Predictions_experimental'])/tmp1.dropna().shape[0]

errors.to_csv('results/errors_Whitfield.csv')
errors_experimental.to_csv('results/errors_Whifield_experimental.csv')
