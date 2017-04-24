import data_loader as dl
import network
import numpy

(training_data, test_data) = dl.format_data_10_output_neurons()

net = network.Network([1024,10,10])
net.squaredGradientDescent(training_data, 100, 0, 0.5)

