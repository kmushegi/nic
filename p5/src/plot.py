from __future__ import print_function
import os, sys

import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Conv2D, MaxPooling2D
from keras.layers import Dense, Dropout, Activation, Flatten
from keras import regularizers
from keras import backend as K
from keras.callbacks import History
from keras.regularizers import l2
from keras.datasets import cifar10
from keras.optimizers import Optimizer

import numpy as np
from matplotlib import pyplot
from scipy.misc import toimage

n_classes = 10
n_epochs = 20
steps_per_epoch = 32
batch_size = 32
data_augmentation = False

#read training & testing data
(x_train,y_train),(x_test,y_test) = cifar10.load_data()
print('Training Samples: ',x_train.shape[0])
print('Testing Samples: ',x_test.shape[0])

#Visualizing CIFAR 10
fig, axes1 = pyplot.subplots(5,5,figsize=(3,3))
for j in range(5):
    for k in range(5):
        i = np.random.choice(range(len(x_train)))
        axes1[j][k].set_axis_off()
        axes1[j][k].imshow(x_train[i:i+1][0])
pyplot.show()
 
'''
for i in range(0,9):
	pyplot.subplot(330+1+i)
	pyplot.imshow(toimage(x_train[i]))
pyplot.show()
'''
