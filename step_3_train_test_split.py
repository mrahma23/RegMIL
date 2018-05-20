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
from random import shuffle

####################################################################

import multiprocessing
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import Queue

####################################################################

def main():
	data_path = os.path.join(os.getcwd(), "metagenomic_data_lsh")
	dest_data_path = os.path.join(os.getcwd(), "combined_train_lsh.txt")
	samples = []
	for (dirpath, dirnames, filenames) in walk(data_path):
		samples.extend(filenames)
		break
	shuffle(samples)
	two_third = (int)(len(samples)*(2/3))
	one_third = len(samples) - two_third
	
	with open(dest_data_path,'w') as f:
		for smpl in samples[:two_third]:
			with open(os.path.join(data_path,smpl), 'r') as f2:
				for line in f2:
					line = line.strip()
					line = smpl + "_" + line + "\n" 
					f.write(line)

					
	with open("test_samples.txt", "w") as f:
		for smpl in samples[two_third:]:
			f.write(smpl + " ")


if __name__=="__main__":
	main()