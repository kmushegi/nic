"""
Evolving Neural Networks Using Genetic Algorithms - Project 5
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Project 5. This file contains the implementation of a Convolutional Neural
Network based on the Keras API.
"""

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

	def __init__(self,n_epochs,dropout,batch_size,optimizer,data_augmentation,
					convActivation,denseActivation,x_train,y_train,x_test,y_test):
		self.n_epochs = n_epochs
		self.data_augmentation = data_augmentation
		self.dropout = dropout
		self.batch_size = batch_size
		self.optimizer = optimizer
		self.convActivation = convActivation
		self.denseActivation = denseActivation
		self.n_classes = 10

		self.parameters = [n_epochs, dropout, batch_size, optimizer, data_augmentation, convActivation, denseActivation]

		self.x_train = x_train
		self.y_train = y_train
		self.x_test = x_test
		self.y_test = y_test

		self.model = Sequential()

	#stack convolutional and pooling layers to create the CNN structure
	def build_network(self):
		self.model.add(Conv2D(32,(3,3), padding='same', \
			input_shape=self.x_train.shape[1:], activation=self.convActivation))
		self.model.add(Conv2D(32,(3,3), activation=self.convActivation))
		self.model.add(MaxPooling2D(pool_size=(2,2)))

		if(self.dropout):
			self.model.add(Dropout(0.25))

		self.model.add(Conv2D(64,(3,3), padding='same', activation=self.convActivation))
		self.model.add(Conv2D(64,(3,3), activation=self.convActivation))
		self.model.add(MaxPooling2D(pool_size=(2,2)))

		if(self.dropout):
			self.model.add(Dropout(0.5))

		self.model.add(Flatten())
		self.model.add(Dense(512, activation=self.convActivation))
		self.model.add(Dense(256, activation=self.convActivation))

		if(self.dropout):
			self.model.add(Dropout(0.5))

		self.model.add(Dense(self.n_classes, activation=self.denseActivation))

		#use optimizer with default parameters for learning rate, decay and etc. goal: find best optimizer
		self.model.compile(loss='categorical_crossentropy',optimizer=self.optimizer,metrics=['accuracy'])

	def train(self):
		hist = self.model.fit(self.x_train, self.y_train, 
			batch_size=self.batch_size,
			epochs=self.n_epochs,
			validation_data=(self.x_test,self.y_test),
			shuffle=True)

		return hist.history['acc'][-1] #return the accuracy on last epoch
