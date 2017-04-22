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
			current_input += line.strip()
			counter += 1
		else:
			inputs.append(current_input)
			desired_outputs.append(int(line))

			current_input = ""
			counter = 0

	return (inputs, desired_outputs)

def format_data_10_output_neurons():

	(tr_i, tr_o) = read_data(training_data_file)
	training_inputs = [np.reshape(list(x), (1024, 1)) for x in tr_i]
	training_outputs = [vectorize_output(x) for x in tr_o]
	training_data = zip(training_inputs,training_outputs)

	(te_i,te_o) = read_data(testing_data_file)
	testing_inputs = [np.reshape(list(x),(1024,1)) for x in te_i]
	#print(te_o[0])
	testing_outputs = [vectorize_output(x) for x in te_o]
	testing_data = zip(testing_inputs,testing_outputs)

	return (training_data, testing_data)

def vectorize_output(j):
	e = np.zeros((10,1))
	e[j] = 1.0
	return e

(tr_data, te_data) = format_data_10_output_neurons()

'''
tr_data is a tuple. First element is array of arrays of size 1024 containing
the binary image; the second element is the array of arrays containing the
information about the "desired" output formatted in a 10 output neuron way.
[1,0,0,0,0,0,0,0,0,0] ==> 0
[0,0,1,0,0,0,0,0,0,0] ==> 2
[0,0,0,0,0,0,0,0,0,1] ==> 9 and etc.
'''

assert(len(tr_data) == 3823)

#print(np.array_str(tr_data[2][1]))
