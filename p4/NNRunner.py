import data_loader as dl
import network
import numpy

(training_data, test_data) = dl.format_data_10_output_neurons()

net = network.Network(1024,10,1,0.01)
net.train(training_data)

