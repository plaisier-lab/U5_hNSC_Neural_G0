##########################################################
## ccAF:  sensitivityAnalysis_plot.py                   ##
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

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

# Import data
DF = pd.read_csv('results/SensitivityAnalysis/ccAF_CV_sensitivity_analysis.csv')
del DF['Unnamed: 0']

# Percentage of genes to test for cross-validation classification
percs = [1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

# Initialize
errorACTINN = []

for i in range(0,len(percs)):
    comparison_column_i = np.where(DF['True Labels'] == DF["%.1f" %(percs[i])], True, False)
    DF['Equal_'+str(percs[i])] = comparison_column_i
    errorACTINN.append(DF['Equal_'+str(percs[i])].value_counts(normalize=True))
    #errorACTINN[percs[i]] = 1-sum(DF['True Labels']==DF[percs[i]])/DF.shape[0]

aggErrorACTINN = []
for i in range(0,len(percs)):
    comparison_column_i = np.where(DF['True Labels'] == DF["%.1f" %(percs[i])], True, False)
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
