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

k = 6

####################################################################

kmer_index_file = os.path.join(os.getcwd(),'kmer_feat_index.txt')
f = open(kmer_index_file, 'w')
allowed = ['a','c','g','t']
features = BTree()
feature_count = -1
for ft_string in map(''.join, itertools.product('acgt', repeat=k)):	
	feature_count = feature_count + 1
	features[ft_string]=feature_count
	f.write(ft_string + ' ' + str(feature_count) + '\n')
f.close()

####################################################################
	
def reader(input_file, inq):
	sequence_number=0
	f = open(input_file, 'r')
	sequence=''
	for line in f:
		line=line.strip()
		if line!='' and line!='\n' and line!=' ':
			if line.startswith('>'):
				if sequence!='':
					item=(sequence_number,sequence)
					inq.put(item)
					sequence_number = sequence_number + 1
					sequence=''
			else:
				sequence=sequence+line    
	if sequence!='':
		item=(sequence_number,sequence)
		inq.put(item)  
	f.close()

####################################################################

def modifier(inq, outq):
	while True:
		item=inq.get()
		if item is None:
			break
		else:
			sequence_number=item[0]
			sequence=item[1].strip().lower()
			sequence = ''.join([i for i in sequence if i in allowed])
			local_features = BTree()
			for i in range(len(sequence)-k+1):
				kmer = sequence[i:i+k]
				feature_index = features[kmer]
				index = feature_index
				try:					
					current_count = local_features[index] 
					local_features[index] = current_count + 1
				except KeyError:
					local_features[index] = 1            
			#normalization
			sum = 0
			for value in local_features.values():
				sum = sum + (value*value)
			sum = math.sqrt(sum)
			for key in local_features.keys():
				value = local_features[key]
				local_features[key] = float(value)/sum			
			out_string = str(sequence_number)
			for key,val in local_features.iteritems():
				out_string = out_string + ' ' + str(key) + '-' + str(val)
			out_string = out_string + '\n'
			outq.put(out_string)
			del local_features
 
####################################################################

def writer(output_file, outq):
	f = open(output_file, 'w')
	while True:
		s = outq.get()
		if s is None:
			break
		else:
			f.write(s)
	f.close()
	
####################################################################

def main():
	data_path = os.path.join(os.getcwd(), "metagenomic_data")
	dest_data_path = os.path.join(os.getcwd(), "metagenomic_data_kmer")
	samples = []
	for (dirpath, dirnames, filenames) in walk(data_path):
		samples.extend(filenames)
		break	
	####################################################################
	for item in samples:			
		input_file = os.path.join(data_path, item)
		output_file = os.path.join(dest_data_path, item)			
		####################################################################
		num_workers = multiprocessing.cpu_count()
		workers=[]		
		inq = multiprocessing.Queue()
		outq = multiprocessing.Queue()
		####################################################################
		for i in range(num_workers):
			tmp = multiprocessing.Process(target=modifier, args=(inq, outq, ))
			tmp.daemon=True
			tmp.start()
			workers.append(tmp)	   
		####################################################################
		fileWriteProcess=multiprocessing.Process(target=writer, args=(output_file, outq,))
		fileWriteProcess.daemon=True
		fileWriteProcess.start()	
		####################################################################
		fileReadProcess=multiprocessing.Process(target=reader, args=(input_file, inq,))
		fileReadProcess.daemon=True
		fileReadProcess.start()	
		####################################################################
		fileReadProcess.join()
		####################################################################
		for i in range(num_workers):
			inq.put(None)			
		####################################################################
		for worker in workers:
			worker.join()	  
		####################################################################
		outq.put(None)		
		####################################################################
		fileWriteProcess.join()	 
		####################################################################			
		inq.close()
		outq.close()
		####################################################################	
		del num_workers
		del workers
		del inq
		del outq
		gc.collect()
		####################################################################

if __name__ == "__main__":
	##### create destination directory ##########
	directory = os.path.join(os.getcwd(), "metagenomic_data_kmer")
	if not os.path.exists(directory):
		os.makedirs(directory) 
	main()
