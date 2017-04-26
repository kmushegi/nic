import data_loader as dl
import network
import sys

bit_map = 1 #1 - bitmap, 0 - downsampled
num_output_neurons = 10 # or 1

assert(num_output_neurons == 10 or num_output_neurons == 1)

(training_data, test_data) = dl.format_data(0,10)

net = network.Network(1024,10,1,0.01)
net.train(training_data)

