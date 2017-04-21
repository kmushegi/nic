import numpy as np
import sys

training_data_file = "digit-recognition-examples/32x32-bitmaps/optdigits-32x32.tra"

def read_data():
	f = open(training_data_file)

	for i in xrange(3): #skip first 3 lines
		next(f)

	inputs = []
	desired_outputs = []
	counter = 0

	current_input = ""
	for line in f:
		if(counter < 32):
			current_input += line
			counter += 1
		else:
			inputs.append(current_input)
			desired_outputs.append(int(line))

			current_input = ""
			counter = 0

read_data()