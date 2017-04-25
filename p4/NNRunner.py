import data_loader as dl
import Network
import numpy

(training_data, test_data) = dl.format_data_10_output_neurons()

net = Network.Network(1024,10,1,1)
net.train(training_data)

