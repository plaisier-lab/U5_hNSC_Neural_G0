##########################################################
## ccAF:  plottingClassifiers.py                        ##
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

# Plotting
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

# Make folder to store downstream results
if not os.path.exists('results/CV_plots'):
    os.makedirs('results/CV_plots')
    print ("Directory created")
else:
    print("Directory already exists")

## Load datasets and concatenate
SVMrej_CR_1536 = pd.read_csv('results/SVMrej/CV_classification_report_1536.csv')
NN_CR_1536 = pd.read_csv('results/ACTINN/CV_classification_report_1536.csv')
RF_CR_1536 = pd.read_csv('results/RFpy/CV_classification_report_1536.csv')
KNN_CR_1536 = pd.read_csv('results/KNN/CV_classification_report_1536.csv')

Classifiers_CR_1536 = pd.concat([SVMrej_CR_1536, RF_CR_1536, KNN_CR_1536, NN_CR_1536], axis = 0)
Classifiers_CR_1536.rename(columns={'Unnamed: 0':'Cell Cycle State'}, inplace=True)

# Write out means
CR_1536 = Classifiers_CR_1536.groupby(['Cell Cycle State','Classifier']).mean()
CR_1536['Subset'] = 1536
CR_1536.to_csv('results/cell_cycle_states_CV_stats.csv')

## Plot f1-score boxplot
# hue = cell cycle state
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.boxplot(x="Classifier", y="f1-score", hue = "Cell Cycle State", data=Classifiers_CR_1536)
plt.savefig('results/CV_plots/F1_scores_for_each_state_per_classifier_1536.pdf')
plt.clf()
# hue = classifier
sns.set(style="darkgrid")
fig, ax = plt.subplots(figsize=(20,10))
sns.boxplot(hue="Classifier", y="f1-score", x = "Cell Cycle State", data=Classifiers_CR_1536)
plt.savefig('results/CV_plots/F1_scores_for_each_state_per_classifier_1536_2.pdf')
plt.clf()
