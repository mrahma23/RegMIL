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

d = 1000

####################################################################

hyper_planes = BTree()
ft_count = 0
kmer_feature_file = os.path.join(os.getcwd(), "kmer_feat_index.txt")
with open(kmer_feature_file, 'r') as f:
	for line in f:
		line = line.strip()
		if line!='' and line!=' ' and line!='\n':
			line = line.split()[1]
			ft_count = ft_count + 1
			try:
				index = int(line)
				hyper_planes[index] = np.random.randn(d)
			except ValueError:
				print('\nkmer feature index not valid\n')
				sys.exit(1)
if index!=ft_count-1:
	print('\nTotal number of kmers does not match!\n')
	sys.exit(1)

####################################################################
	
def reader(input_file, inq):
	sequence_number=0
	f = open(input_file, 'r')
	sequence=''
	for line in f:
		if inq.qsize()>100000:
			while(inq.qsize()>100000):
				pass 
		line=line.strip()
		if line!='' and line!='\n' and line!=' ':
			inq.put(line)  
	f.close()
 
####################################################################
 
def modifier(inq, outq):
	while True:
		item=inq.get()
		if item is None:
			break
		else:
			sum = np.zeros(d)
			item = item.split()
			for ft in item[1:]:
				ft = ft.split('-')
				try:
					ind = int(ft[0])
					value = float(ft[1])
				except ValueError:
					print('\nkmer feature (or value) invalid in sample data!\n')
					sys.exit(1)
				ft = np.empty(d)
				ft.fill(value)
				ft = np.multiply(ft, hyper_planes[ind])
				sum = np.add(ft,sum)
			new_data = 0
			for v in np.nditer(sum):
				new_data = new_data << 1
				if v >= 0:
					new_data |= 1
			out_string = item[0].strip() + ' ' + str(new_data) + '\n'
			outq.put(out_string)
  
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
	data_path = os.path.join(os.getcwd(), "metagenomic_data_kmer")
	dest_data_path = os.path.join(os.getcwd(), "metagenomic_data_lsh")
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
	directory = os.path.join(os.getcwd(), "metagenomic_data_lsh")
	if not os.path.exists(directory):
		os.makedirs(directory)  
	main()
