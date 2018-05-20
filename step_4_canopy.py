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
import random

####################################################################

import multiprocessing
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import Queue

####################################################################################################
t1 = 0.41
t2 = 0.38
d = float(1000.00)

def modifier(canopyCenterHash, inq, outq, deleteq):
	while True:
		s=inq.get()
		if s is None:
			break
		else:
			seq = s[0]
			hv = s[1]
			xor = canopyCenterHash^hv
			distance = (bin(xor)[2:].count('1'))/d 
			if distance <= t2:
				deleteq.put(seq)
				outq.put(seq)
			elif distance <= t1:
				outq.put(seq)

####################################################################################################

def writer(canopyIndex, outq, st_curr_canopy_q):
	counter = 0
	# f = open("temp_file.txt", "w")
	# f.write(canopyIndex)
	st_curr_canopy = canopyIndex
	while True:
		s = outq.get()
		if s is None:
			break
		else:
			#f.write(" " + s)
			st_curr_canopy = st_curr_canopy + " " + s
			counter += 1
	# f.write("\n")
	# f.close()
	st_curr_canopy = st_curr_canopy + "\n"
	st_curr_canopy_q.put(st_curr_canopy)
	print("\n\t# of members = " + str(counter) + "\n")

####################################################################################################

def data_updater(current_dataset, deleteq, current_dataset_q):
	counter = 0
	while True:
		s = deleteq.get()
		if s is None:
			break
		else:
			current_dataset.remove(s)
			counter += 1
	current_dataset_q.put(current_dataset)
	print("\n\t# of deleted members = " + str(counter) + "\n")

####################################################################################################

def main():
	hash_file = os.path.join("combined_train_lsh.txt")
	canopy_output_file = os.path.join("canopy_output.txt")

	hashed_sequence = BTree()
	current_dataset = []
	dataCounter = 0

	with open(hash_file,'r') as f:
		for line in f:
			line = line.strip().split()
			try:
				seq_id = line[0]
				hash_value = (int)(line[1])
				current_dataset.append(seq_id)
				dataCounter = dataCounter + 1
				hashed_sequence[seq_id] = hash_value
			except ValueError:
				print("\nhash value invalid\n")
				sys.exit(1)  
	print("total # of instances in training set: " + str(dataCounter))
	gc.collect()
	time.sleep(2)
	
	canopy_counter = 0
	with open(canopy_output_file,'w') as f:
		while dataCounter>0:
			canopyIndex = random.choice(current_dataset)
			canopyCenterHash = hashed_sequence[canopyIndex]
			current_dataset.remove(canopyIndex)
			dataCounter = dataCounter - 1
			remove_list = []
			####################################################################
			num_workers = (int)(multiprocessing.cpu_count()/2)
			workers = []      
			inq = multiprocessing.Queue()
			outq = multiprocessing.Queue()
			deleteq = multiprocessing.Queue()
			st_curr_canopy_q = multiprocessing.Queue()
			current_dataset_q = multiprocessing.Queue()
			####################################################################
			for i in range(num_workers):
				tmp = multiprocessing.Process(target=modifier, args=(canopyCenterHash, inq, outq, deleteq, ))
				tmp.daemon=True
				tmp.start()
				workers.append(tmp)    
			####################################################################
			fileWriteProcess=multiprocessing.Process(target=writer, args=(canopyIndex, outq, st_curr_canopy_q, ))
			fileWriteProcess.daemon=True
			fileWriteProcess.start()
			####################################################################
			for seq in current_dataset:
				hv = hashed_sequence[seq]
				inq.put((seq,hv))
			####################################################################
			datasetUpdater=multiprocessing.Process(target=data_updater, args=(current_dataset, deleteq, current_dataset_q, ))
			datasetUpdater.daemon=True
			datasetUpdater.start()
			####################################################################
			for i in range(num_workers):
				inq.put(None)
			for worker in workers:
				worker.join()     
			####################################################################
			outq.put(None)
			deleteq.put(None)

			while True:
				try:
					curr_canopy = st_curr_canopy_q.get()
					break
				except:
					pass
			f.write(curr_canopy)
			fileWriteProcess.join()
			
			while True:
				try:
					c_s = current_dataset_q.get()
					current_dataset = c_s
					break
				except:
					pass
			datasetUpdater.join()  
			####################################################################
			# with open("temp_file.txt", "r") as t:
			# 	s = t.readline()
			# f.write(s)
			####################################################################
			inq.close()
			outq.close()
			#final_output_string.close()
			deleteq.close()
			del inq
			del outq
			#del final_output_string
			del deleteq
			del workers
			gc.collect()			
			####################################################################
			dataCounter = len(current_dataset)
			canopy_counter = canopy_counter + 1
			print(str(canopy_counter) + ". Remaining # of instances: " + str(dataCounter))
			

		  

if __name__=="__main__":
	main()