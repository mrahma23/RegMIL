import numpy as np
from BTrees.OOBTree import BTree

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

####################################################################

import multiprocessing
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import Queue
from Cython.Shadow import address

####################################################################

sampleClassLabels = BTree()
train_instances = BTree()
metadata = os.path.join(os.getcwd(),"metadata.txt")
kmer_data_path = os.path.join(os.getcwd(), "metagenomic_data_kmer")
canopy_output = os.path.join(os.getcwd(), "canopy_output.txt")
positivity_outputs = os.path.join(os.getcwd(), "positivity_outputs.txt")
train_set_positivity_scores = os.path.join(os.getcwd(), "train_set_positivity_scores")
samples = []
for (dirpath, dirnames, filenames) in walk(kmer_data_path):
	samples.extend(filenames)
	break

####################################################################

def classLabelCreator():
	global samples
	global sampleClassLabels
	global metadata
	for s in samples:
		sampleClassLabels[s] = 0
	with open(metadata,'r') as mf:
		for line in mf:
			line  = line.strip().split()
			if ("CD" in line[0]) or ("UC" in line[0]):
				for key in sampleClassLabels.keys():
					if line[0] in key:
						sampleClassLabels[key] = 1
####################################################################

def membertag_to_class(canopy_member):
	return sampleClassLabels[canopy_member.strip().split("_")[0]]
####################################################################

def positivity(address):
	################# read canopy output #################
	positivity_scores = BTree()
	with open(address,'r') as f:
	    for line in f:
	        canopy_members = line.strip().split()
	        total = float(len(canopy_members))
	        number_of_positive = 0
	        for member in canopy_members:
	            if membertag_to_class(member) == 1:
	                number_of_positive = number_of_positive + 1
	        p_score = float(number_of_positive)/total
	        for member in canopy_members:
	            try:
	                score_list = positivity_scores[member]
	                score_list.append(p_score)
	                positivity_scores[member] = score_list
	            except:
	                positivity_scores[member] = [p_score]
	############### take average score ###################
	for key in positivity_scores.keys():
	    score_list = positivity_scores[key]
	    sum = float(0)
	    for item in score_list:
	        sum = sum + item
	    final_score = sum/float(len(score_list))
	    positivity_scores[key] = final_score
	################ prepare output ####################
	output = ''
	for key,value in positivity_scores.items():
	    output = output + key + ' ' + str(value) + '\n'
	with open(positivity_outputs, 'w') as f:
	    f.write(output)


####################################################################
def train_set_positivity():
	data = BTree()
	with open("positivity_outputs.txt", "r") as f:
		for line in f:
			line = line.strip().split()
			sample = line[0].split("_")[0]
			seq_number = line[0].split("_")[1]
			score = line[1]
			try:
				curr = data[sample]
				curr.append((seq_number,score))
				data[sample] = curr
			except:
				data[sample] = [(seq_number,score)]
	for key, val in data.items():
		with open(os.path.join(train_set_positivity_scores, key), "w") as f:
			for item in val:
				f.write(item[0] + " " + item[1] + "\n")
		

 
####################################################################

def main():
	classLabelCreator()
	positivity(canopy_output)
	train_set_positivity()
            
if __name__ == "__main__": 
    main()
