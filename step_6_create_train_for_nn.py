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

####################################################################

positivity_score_output = os.path.join(os.getcwd(),"positivity_outputs.txt")
nn_train_file = os.path.join(os.getcwd(), "nn_train_file.txt")
kmer_file = os.path.join(os.getcwd(), "metagenomic_data_kmer")
scores = BTree()

def read_scores():
	with open(positivity_score_output, "r") as f:
		for line in f:
			line = line.strip().split()
			sample_sequence = line[0]
			sample = sample_sequence.split("_")[0]
			sequence = (int)(sample_sequence.split("_")[1])
			score = (float)(line[1])
			try:
				current = scores[sample]
				current.append((sequence,score))
				scores[sample] = current
			except:
				scores[sample] = [(sequence,score)]
			gc.collect()

def sort_sequences():
	for key, val in scores.items():
		val.sort(key=lambda x: x[0])
		scores[key] = val

def check(seq_number, sequence_list):
	for item in sequence_list:
		if item[0]>seq_number:
			return False
		if item[0]==seq_number:
			return True
	return False


def join_kmers():
	nn_file = open(nn_train_file, "w")	
	for sample, sequence_list in scores.items():		
		with open(os.path.join(kmer_file, sample), "r") as f:					
			for line in f:
				line = line.strip()
				seq_number = (int)(line.split()[0])
				kmers = line.split(" ", 1)[1]
				for item in sequence_list:
					if item[0]>seq_number:
						break
					elif item[0]==seq_number:
						output = sample + "_" + str(item[0]) + " " + str(item[1]) + " " + kmers + "\n"
						nn_file.write(output)	
	nn_file.close()


def main():
	read_scores()
	sort_sequences()
	join_kmers()



if __name__=="__main__":
	main()

