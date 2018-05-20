import urllib.request
import numpy as np
import hashlib
####################################################################
import os
import gc
import shutil 
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
####################################################################

def downloader(inq):
    while True:
        s = inq.get()
        if s is None:
            break
        else:
            ######## information ##############
            sample_name = s[0]
            country = s[1]
            gender = s[2]
            age = s[3]
            bmi = s[4]
            is_ibd = s[5]
            url = s[6]
            destination_path = s[7]
            try:
                urllib.request.urlretrieve(url, destination_path)
                print("\nCompleted: " + url + "\n")
            except:
                print('\n\nError occurred for sample: ' + sample_name + '\n')

   
            
def main():    
    start = time.time()        
    num_workers=multiprocessing.cpu_count()*2
    workers=[]    
    inq = multiprocessing.Queue()
    for i in range(num_workers):
        tmp = multiprocessing.Process(target=downloader, args=(inq,))
        tmp.daemon=True
        tmp.start()
        workers.append(tmp)

    ##### create destination directory ##########
    directory = os.path.join(os.getcwd(), "metagenomic_data")
    if not os.path.exists(directory):
        os.makedirs(directory)

    ######### prepare URL's for samples using the metadata file and pass to multiprocess ##########
    metaData = os.path.join(os.getcwd(),'metadata.txt')
    url_prefix = "ftp://public.genomics.org.cn/BGI/gutmeta/Single_Sample_contig/"
    url_suffix = ".scaftig.more500.gz"
    ##### open the metadata file ############################
    with open(metaData,'r') as f:
        for line in f:
            line = line.strip().split()

            ######## information ##############
            sample_name = line[0]
            country = line[1]
            gender = line[2]
            age = (int)(line[3])
            bmi = (float)(line[4])
            is_ibd = line[5]

            ######## prepare URL ############
            # ftp://public.genomics.org.cn/BGI/gutmeta/Single_Sample_contig/MH0001.scaftig.more500.gz
            url = url_prefix + sample_name + url_suffix
            destination_path = os.path.join(directory, sample_name + url_suffix)
            it = (sample_name, country, gender, age, bmi, is_ibd, url, destination_path)
            inq.put(it)
    for i in range(num_workers):
        inq.put(None)  
    for worker in workers:
        worker.join()
    inq.close()
    duration = time.time() - start

    del num_workers
    del workers
    del inq
    gc.collect()
    print(str(duration) + ' in seconds ' + str(duration/60.00) + ' in minutes')
    
if __name__ == "__main__": 
    main()