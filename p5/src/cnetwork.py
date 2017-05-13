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
import time

outdir = '../stats/'

class CNetwork(object):

	def __init__(self,n_epochs,n_layers,dropout,batch_size,optimizer,data_augmentation,
					convActivation,denseActivation,x_train,y_train,x_test,y_test):
		self.n_epochs = n_epochs
		self.n_layers = n_layers
		self.data_augmentation = data_augmentation
		self.dropout = dropout
		self.batch_size = batch_size
		self.optimizer = optimizer
		self.convActivation = convActivation
		self.denseActivation = denseActivation

		self.x_train = x_train
		self.y_train = y_train
		self.x_test = x_test
		self.y_test = y_test

		self.model = Sequential()

		self.n_classes = 10
		self.steps_per_epoch = 32

	def build_network():
		self.model.add(Conv2D(32,(3,3), padding='same', \
			input_shape=x_train.shape[1:], activation=self.convActivation))
		self.model.add(Conv2D(32,(3,3), activation=self.convActivation))
		self.model.add(MaxPooling2D(pool_size=(2,2)))

		if(dropout):
			self.model.add(Dropout(0.25))

		self.model.add(Conv2D(64,(3,3), padding='same', activation=self.convActivation))
		self.model.add(Conv2D(64,(3,3), activation=self.convActivation))
		self.model.add(MaxPooling2D(pool_size=(2,2)))

		if(dropout):
			self.model.add(Dropout(0.5))

		self.model.add(Flatten())
		self.model.add(Dense(512, activation=self.convActivation))
		self.model.add(Dense(256, activation=self.convActivation))

		if(dropout):
			self.model.add(Dropout(0.5))

		self.model.add(Dense(n_classes, activation=self.denseActivation))

		#possibly add data augmentation here

		#opt = keras.optimizers.rmsprop(lr=0.0001,decay=1e-6)
		#use optimizer with default parameters for learning rate, decay and etc and find best optimizer
		self.model.compile(loss='categorical_crossentropy',optimizer=self.optimizer,metrics=['accuracy'])

	def train():
		hist = model.fit(x_train, y_train, 
			batch_size=self.batch_size,
			epochs=self.n_epochs,
			validation_data=(self.x_test,self.y_test),
			shuffle=True)

		write_out_training_history(hist)

	def write_out_training_history(h):
		fn = outdir+str(time.time())+".losses.txt"
		f = open(fn,'w')
		f.write('loss\tval_loss\tacc\tval_acc\ttop_k\tval_top_k\n')
		for i in np.arange(self.n_epochs):
			f.write('%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f'%(
				h.history['loss'][i],
				h.history['val_loss'][i],
				h.history['acc'][i],
				h.history['val_acc'][i],
				h.history['top_k_categorical_accuracy'][i],
				h.history['val_top_k_categorical_accuracy'][i]))



								
#this will be in format data called before the call to the network
#so that data is only read once
#read training & testing data
(x_train,y_train),(x_test,y_test) = cifar10.load_data()
print('Training Samples: ',x_train.shape[0])
print('Testing Samples: ',x_test.shape[0])

y_train = keras.utils.to_categorical(y_train, n_classes)
y_test = keras.utils.to_categorical(y_test, n_classes)

x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train /= 255.0
x_test /= 255.0

#########
