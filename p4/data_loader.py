import numpy as np
import sys

training_data_file = "digit-recognition-data/32x32-bitmaps/optdigits-32x32.tra"
testing_data_file = "digit-recognition-data/32x32-bitmaps/optdigits-32x32.tes"

def read_data(data_file):
	f = open(data_file)

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

	return (inputs, desired_outputs)

(tra_inputs,tra_desired_outputs) = read_data(training_data_file)
(tes_inputs,tes_desired_outputs) = read_data(testing_data_file)
