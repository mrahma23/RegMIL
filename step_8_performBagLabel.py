import numpy as np
from BTrees.OOBTree import BTree
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score

####################################################################

import os
from os import walk
import gc
import shutil
import glob 
from shutil import copyfile
import errno
import sys
import itertools
import time
import math
from operator import itemgetter

####################################################################

metadata = os.path.join(os.getcwd(),"metadata.txt")
kmer_data_path = os.path.join(os.getcwd(), "metagenomic_data_kmer")
all_samples = []
for (dirpath, dirnames, filenames) in walk(kmer_data_path):
    all_samples.extend(filenames)
    break

train_source_data = os.path.join(os.getcwd(), "train_set_positivity_scores")
train_samples = []
for (dirpath, dirnames, filenames) in walk(train_source_data):
    train_samples.extend(filenames)
    break

test_source_data = os.path.join(os.getcwd(), "regression_output")
test_samples = []
for (dirpath, dirnames, filenames) in walk(test_source_data):
    test_samples.extend(filenames)
    break

sampleClassLabels = BTree()

####################################################################

def classLabelCreator():
    global all_samples
    global sampleClassLabels
    global metadata
    for s in all_samples:
        sampleClassLabels[s] = 0
    with open(metadata,'r') as mf:
        for line in mf:
            line  = line.strip().split()
            if ("CD" in line[0]) or ("UC" in line[0]):
                for key in sampleClassLabels.keys():
                    if line[0] in key:
                        sampleClassLabels[key] = 1

####################################################################

def bag_feature_maker(path, samples, max_features, step_size, start_points):
    features = []
    labels = []       
    for sample in samples:
        feature_dictionary = {}
        for i in start_points:
            feature_dictionary[i] = float(0)
        #sample_name = sample.split('.')[0]
        file_path = os.path.join(path,sample)
        lineCount = 0
        with open(file_path,'r') as f:
            for line in f:
                line = line.strip()
                if line!='' and line!=' ' and line!='\n':
                    lineCount = lineCount + 1                    
                    sequence_number, score = line.split()
                    sequence_number = int(sequence_number)
                    score = float(score)                    
                    for i in start_points:
                        if i<=score<i+step_size:
                            count = feature_dictionary[i]
                            feature_dictionary[i] = count + 1.00        
        feature_vector = []
        for i in start_points:
            feature_vector.append(feature_dictionary[i]/float(lineCount))            
        class_labels = -1 ## HEALTHY    
        if sampleClassLabels[sample]==1:
            class_labels = 1 ## DISEASE            
        labels.append(class_labels)
        features.append(feature_vector)    
    features = np.asarray(features)
    labels = np.asarray(labels)    
    return features, labels

####################################################################

def main():
    classLabelCreator() 
    max_features = 20
    step_size = float(1)/float(max_features)
    start_points = np.arange(float(0.0),float(1.0),float(step_size))
    
    train_bag_features, train_bag_labels = bag_feature_maker(train_source_data, train_samples, max_features, step_size, start_points)
    test_bag_features, test_bag_labels = bag_feature_maker(test_source_data, test_samples, max_features, step_size, start_points)
    
    clf = RandomForestClassifier(n_estimators = 100, max_features=None)
    clf.fit(train_bag_features, train_bag_labels)
    
    predicted_test_bag_labels = clf.predict(test_bag_features)
    rf_acc = accuracy_score(test_bag_labels, predicted_test_bag_labels) 
    rf_auc = roc_auc_score(test_bag_labels, predicted_test_bag_labels)

    with open("bag_prediction_results.txt", "w") as f:
        for smpl, smpl_pred in zip(test_samples, predicted_test_bag_labels):
            p = 0
            if smpl_pred==1:
                p = 1
            f.write(smpl + " ----> predicted: " + str(p) + " ----> actual: " + str(sampleClassLabels[smpl]) + "\n")
        f.write("accuracy: " + str(rf_acc) + "\n")
        f.write("AUC-ROC value: " + str(rf_auc) + "\n")

 
if __name__ == "__main__": 
    main()


            