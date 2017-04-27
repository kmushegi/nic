import numpy as np

bitmap_training_data_file = "data/bitmaps/optdigits-32x32.tra"
bitmap_testing_data_file = "data/bitmaps/optdigits-32x32.tes"

downsampled_training_data_file = "data/downsampled/optdigits-8x8-int.tra"
downsampled_testing_data_file = "data/downsampled/optdigits-8x8-int.tes"

NUM_TRAINING_IMAGES = 3823

np.set_printoptions(threshold=np.nan)

def read_data_down_sampled(data_file):
	f = open(data_file)

	inputs = []
	desired_outputs = []

	for line in f:
		tokens = line.split(',')

		inputs.append(" ".join(tokens[:-1]))
		desired_outputs.append(int(tokens[-1]))

	return (inputs, desired_outputs)

def read_data_bit_map(data_file):
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

'''
*_data is a tuple. First element is array of arrays of size 1024 containing
the binary image; the second element is the array of arrays containing the
information about the "desired" output formatted in a 10 output neuron way.
[1,0,0,0,0,0,0,0,0,0] ==> 0
[0,0,1,0,0,0,0,0,0,0] ==> 2
[0,0,0,0,0,0,0,0,0,1] ==> 9 and etc.
'''

def format_data(bit_map, num_output_neurons):
	if bit_map:
		(tr_i, tr_o) = read_data_bit_map(bitmap_training_data_file)
		(te_i,te_o) = read_data_bit_map(bitmap_testing_data_file)

		training_inputs = [np.reshape(list(x), (1024, 1)) for x in tr_i]
		testing_inputs = [np.reshape(list(x), (1024, 1)) for x in te_i]


		if num_output_neurons == 10:
			training_outputs = [digit_to_vector_representation(x) for x in tr_o]
			testing_outputs = [digit_to_vector_representation(x) for x in te_o]
		else:
			training_outputs = [int(x) for x in tr_o]
			testing_outputs = [int(x) for x in te_o]

		training_data = zip(training_inputs,training_outputs)
		testing_data = zip(testing_inputs,testing_outputs)

		return (training_data, testing_data)
	else:
		(tr_i, tr_o) = read_data_down_sampled(downsampled_training_data_file)
		(te_i, te_o) = read_data_down_sampled(downsampled_testing_data_file)

		training_inputs = [np.reshape(x.split(' '), (64, 1)) for x in tr_i]
		testing_inputs = [np.reshape(x.split(' '), (64, 1)) for x in te_i]

		if num_output_neurons == 10:
			training_outputs = [digit_to_vector_representation(x) for x in tr_o]
			testing_outputs = [digit_to_vector_representation(x) for x in te_o]
		else:
			training_outputs = [int(x) for x in tr_o]
			training_outputs = [int(x) for x in te_o]

		training_data = zip(training_inputs, training_outputs)
		testing_data = zip(testing_inputs, testing_outputs)

		return (training_data, testing_data)

def digit_to_vector_representation(j):
	e = np.zeros((10,1))
	e[j] = 1.0
	return e


#(tr_data, test_data) = format_data(0,10)
#print(np.array_str(tr_data[0][0]))
#assert(len(tr_data) == NUM_TRAINING_IMAGES)
