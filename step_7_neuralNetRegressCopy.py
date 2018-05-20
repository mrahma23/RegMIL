import numpy as np
from BTrees.OOBTree import BTree
from sklearn.neural_network import MLPRegressor
from scipy.sparse import coo_matrix
from sklearn import preprocessing
import os
import gc
import shutil
import glob 
from shutil import copyfile
import errno
import sys
import itertools
import time
import math

################## global vars ############################################

train_data_path = os.path.join(os.getcwd(), "nn_train_file.txt")
test_data_path = os.path.join(os.getcwd(), "metagenomic_data_kmer")
regression_output_Path = os.path.join(os.getcwd(), "regression_output")

k = 6
total_cols = 4**k

##########################################################################

def train():
	row = []
	col = []
	data = [] 
	y_train = []
	row_counter = 0
	with open(train_data_path,'r') as f:
	    for line in f:
	        line = line.strip()
	        if line!=' ' and line!='' and line!='\n':
	            items = line.split()
	            sequence_id = items[0]
	            score = float(items[1])
	            y_train.append(score)
	            for i in items[2:]:
	                col_id, data_value = i.strip().split("-")
	                col_id = int(col_id)
	                data_value = float(data_value)
	                row.append(row_counter)
	                col.append(col_id)
	                data.append(data_value)
	            row_counter = row_counter + 1
	#########################################################
	row = np.asarray(row)
	col = np.asarray(col)
	data = np.asarray(data)
	x_train = coo_matrix((data, (row, col)), shape=(row_counter, total_cols))
	scaler = preprocessing.StandardScaler(with_mean=False).fit(x_train)
	x_train = scaler.transform(x_train)
	y_train = np.asarray(y_train)
	#########################################################
	total_layes = 4
	num_neurons_per_layer = 300
	hidden_layers_collection = (num_neurons_per_layer,)
	for ly in range(total_layes):
	    hidden_layers_collection = hidden_layers_collection + (num_neurons_per_layer,)
	#########################################################
	regressor = MLPRegressor(hidden_layer_sizes=hidden_layers_collection)
	regressor.fit(x_train,y_train)
	return regressor, scaler


def test(regressor, scaler):
	test_samples = []
	with open(os.path.join(os.getcwd(), "test_samples.txt"), "r") as f:
		line = f.readline()
		test_samples = line.strip().split()
	
	for sample in test_samples:
	    regressed_values = BTree()
	    test_row = []
	    test_col = []
	    test_data = []
	    path = os.path.join(test_data_path, sample)
	    row_counter = 0
	    with open(path,'r') as f:
	        for line in f:
	            line = line.strip()
	            if line!=' ' and line!='' and line!= '\n':
	                items = line.split()
	                sequence_number = int(items[0].strip())
	                regressed_values[row_counter] = sequence_number
	                for i in items[1:]:
	                    col_id, data_value = i.strip().split("-")
	                    col_id = int(col_id)
	                    data_value = float(data_value)
	                    test_row.append(row_counter)
	                    test_col.append(col_id)
	                    test_data.append(data_value)
	                row_counter = row_counter + 1
	    test_row = np.asarray(test_row)
	    test_col = np.asarray(test_col)
	    test_data = np.asarray(test_data)
	    x_test = coo_matrix((test_data, (test_row, test_col)), shape=(row_counter, total_cols))
	    x_test = scaler.transform(x_test)

	    y_test = regressor.predict(x_test)
	    y_test[y_test < 0] = 0.0
	    y_test[y_test > 1] = 1.0
	    

	    toBeSorted = []
	    for r in range(0,len(y_test)):
	        toBeSorted.append((regressed_values[r], y_test[r]))
	    toBeSorted = sorted(toBeSorted, key = lambda x: x[1], reverse=True)
	    gc.collect()
	    out = os.path.join(regression_output_Path,sample)
	    with open(out,'w') as f:
	        for value in toBeSorted:
	            f.write(str(value[0]) + ' ' + str(value[1]) + '\n')

 

def main():
	directory = os.path.join(os.getcwd(), "regression_output")
	if not os.path.exists(directory):
		os.makedirs(directory)
	regressor, scaler = train()
	test(regressor, scaler)


if __name__ == "__main__":
	main()