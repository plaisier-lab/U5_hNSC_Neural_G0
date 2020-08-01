##########################################################
## ccAF:  classifiersV3.py                              ##
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

##########################################
## Load Python packages for classifiers ##
##########################################

# General
import numpy as np
import pandas as pd
import os
from scipy.sparse import isspmatrix

# sklearn: Support Vector Machine (SVM) rejection
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import BaggingClassifier

# sklearn: Random Forest (RF)
from sklearn.ensemble import RandomForestClassifier

# ACTINN
import actinn

# scanpy: K-Nearest Neighbors (KNN)
import scanpy as sc
import anndata

#################
## Classifiers ##
#################

class Classifier_SVMrej:
    """A class designed to facilitate SVMrej classifier construction
    and use."""
    def __init__(self, data, labels, cutoff = 0.7):
        self.data = data
        self.labels = labels
        self.cutoff = cutoff
        self.genes = data.var_names
        self.classifier = self.__build_classifier(data, labels)

    # Build classifier
    def __build_classifier(self, data, labels):
        train = self.__prep_data(data)
        Classifier = LinearSVC(max_iter = 100000)
        clf = CalibratedClassifierCV(Classifier)
        clf.fit(train, labels)
        return clf

    # Predict probability
    def predict_prob(self, newData):
        return np.max(self.classifier.predict_proba(newData), axis=1)

    # Prepare data
    def __prep_data(self, data):
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names)
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if isspmatrix(test_data.X):
            data = pd.DataFrame(test_data.X.todense(), index = test_data.obs_names, columns = test_data.var_names)
        else:
            data = pd.DataFrame(test_data.X, index = test_data.obs_names, columns = test_data.var_names)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(0,index=data.index, columns = missing)], axis=1)
        return data[list(self.genes)]

    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff = None):
        pred_data = self.__prep_predict_data(new_data)
        labels = self.classifier.predict(pred_data)
        if cutoff == None and not self.cutoff == None:
            cutoff = self.cutoff
        if not cutoff == None:
            probs = self.predict_prob(pred_data)
            unlabeled = np.where(probs < cutoff)
            labels[unlabeled] = 'Unknown'
        return labels


class Classifier_KNN:
    """A class designed to facilitate scanpy ingetst based
    K-nearest neighbor classifier construction and use."""
    def __init__(self, data, label):
        self.data = data
        self.genes = data.var_names
        self.label = label
        # self.classifier = self.__build_classifier(data, labels)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if len(missing)>0:
            data = pd.concat([pd.DataFrame(test_data.X, index=test_data.obs_names, columns=test_data.var_names), pd.DataFrame(0,index=test_data.obs_names, columns = missing)], axis=1)
            data = data[list(self.genes)]
            data_sc = anndata.AnnData(X=data.to_numpy())
            data_sc.var_names = data.columns
            data_sc.obs_names = data.index
            return data_sc
        else:
            return test_data
    
    # Predict labels
    def predict_labels(self, new_data):
        # Subset based on common gene names
        adata_query = self.__prep_predict_data(new_data)
        var_names = self.data.var_names.intersection(adata_query.var_names)
        adata_ref = self.data[:,var_names]
        adata_query = adata_query[:,var_names]

        # Run embedding anlaysis using subset
        sc.pp.pca(adata_ref)
        sc.pp.neighbors(adata_ref)
        sc.tl.umap(adata_ref)

        # Map the identifiers from the reference dataset to the query dataset
        sc.tl.ingest(adata_query, adata_ref, obs=self.label)

        # Save the results in whitfield data object
        return adata_query.obs[self.label]


class Classifier_RF:
    """A class designed to facilitate RF classifier construction
    and use. Can also be used for """
    def __init__(self, data, labels, cutoff = None):
        self.data = data
        self.genes = data.var_names
        self.labels = labels
        self.cutoff = cutoff
        self.classifier = self.__build_classifier()

    # Build classifier
    def __build_classifier(self):
        train = self.__prep_data(self.data)
        clf = RandomForestClassifier(n_estimators=500, oob_score=True)
        clf.fit(train, self.labels)
        return clf

    # Predict probability
    def predict_prob(self, newData):
        return np.max(self.classifier.predict_proba(newData), axis=1)
    
    # Prepare data
    def __prep_data(self, data):
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names)
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names)

    # Prepare test data for predicting
    def __prep_predict_data(self, test_data):
        missing = set(self.genes).difference(test_data.var_names)
        if isspmatrix(test_data.X):
            data = pd.DataFrame(test_data.X.todense(), index = test_data.obs_names, columns = test_data.var_names)
        else:
            data = pd.DataFrame(test_data.X, index = test_data.obs_names, columns = test_data.var_names)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(0,index=data.index, columns = missing)], axis=1)
        return data[list(self.genes)]
    
    # Predict labels with rejection
    def predict_labels(self, new_data, cutoff = None):
        pred_data = self.__prep_predict_data(new_data)
        labels = self.classifier.predict(pred_data)
        if cutoff == None and not self.cutoff == None:
            cutoff = self.cutoff
        if not cutoff == None:
            probs = self.predict_prob(pred_data)
            unlabeled = np.where(probs < cutoff)
            labels[unlabeled] = 'Unknown'
        return labels


class Classifier_ACTINN:
    """A class designed to facilitate ACTINN classifier construction
    and use. Can also be used for """
    def __init__(self, train, label, learning_rate = 0.0001, num_epochs = 200, minibatch_size = 128, print_cost = True, output_probability = False):
        self.train = train
        self.label = label
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.minibatch_size = minibatch_size
        self.print_cost = print_cost
        self.output_probability = output_probability
        self.label = label
        self.classifier, self.label_to_type_dict, self.genes = self.__build_classifier()

    # Prepare data
    def __prep_data(self, data):
        # Make indicies unique for
        data.var_names_make_unique()
        # Remove all genes with zero counts
        sc.pp.filter_genes(data, min_cells=1)
        if isspmatrix(data.X):
            return pd.DataFrame(data.X.todense(), index = data.obs_names, columns = data.var_names).T
        else:
            return pd.DataFrame(data.X, index = data.obs_names, columns = data.var_names).T
    
    # Prepare test data for predicting
    def __prep_predict_data(self, data):
        missing = set(self.genes).difference(data.index)
        if len(missing)>0:
            data = pd.concat([data, pd.DataFrame(0,index=missing, columns = data.columns)])
        return data.loc[list(self.genes)]

    # Build classifier
    def __build_classifier(self):
        train = self.train
        # Convert into pandas DataFrame
        train_data = self.__prep_data(train)
        labels = self.train.obs[self.label]
        clf, label_to_type_dict, genes = actinn.train_model(train_data, labels, learning_rate = self.learning_rate, num_epochs = self.num_epochs, minibatch_size = self.minibatch_size, print_cost = self.print_cost)
        return clf, label_to_type_dict, genes

    # Predict labels with rejection
    def predict_labels(self, newData):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        labels = actinn.predict_new_data(pred_data, self.classifier, self.label_to_type_dict, self.genes)
        return list(labels['celltype'])

    # Predict labels with rejection
    def predict_probs(self, newData, axis = -1):
        test_data = self.__prep_data(newData)
        pred_data = self.__prep_predict_data(test_data)
        probs = actinn.predict_probabilities(pred_data, self.classifier, self.genes, axis=axis)
        probs.index = pd.Index([self.label_to_type_dict[i] for i in list(probs.index)])
        return probs

