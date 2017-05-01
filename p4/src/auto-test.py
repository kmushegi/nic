"""
Neural Networks for Digit Recognition - Project 4
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Digit Recognition, Project 4. This file contains the 
implementation of our automated testing framework used to derive optimal parameters for the
neural network.
"""

import os
import sys
from joblib import Parallel, delayed
import multiprocessing

stats_dir = '../stats/'

data_representations = [0,1] #downsampled(0), bitmap(1)
output_neurons = [1,10]
learning_rates = [0.1,0.5,1]

#stats file format: (data_rep)-(output_neurons)-(learning_rate).txt

def clean_stats_folder():
	print("Cleaning stats folder")
	cmd = "cd ../stats && rm *.txt"
	os.system(cmd)

def clean_pyc_files():
	print("Cleaning .pyc files")
	cmd = "rm *.pyc"
	os.system(cmd)

def work(data_format,n_outputs,learning_rate):
	msg = ("Bitmap" if data_format == 1 else "Downsampled") + "\t::\tN_Outputs: " + str(n_outputs) \
		+ "\t::\tLearning Rate: " + str(learning_rate)
	print msg

	cmd = "python NNRunner.py " + str(data_format) + " " + str(n_outputs) \
		+ " " + str(learning_rate) + ">> ../stats/" + str(data_format) + "-" + str(n_outputs) \
		+ "-" + str(learning_rate) + ".txt"
	os.system(cmd)

clean_stats_folder()
clean_pyc_files()

if not os.path.isdir(stats_dir):
	cmd = "cd .. && mkdir stats"
	os.system(cmd)

job_count = multiprocessing.cpu_count()
print("Spawning " +str(job_count) + " jobs")

Parallel(n_jobs=job_count)(delayed(work)(data_format=df,n_outputs=no,learning_rate=lr) \
	for df in data_representations for no in output_neurons for lr in learning_rates)