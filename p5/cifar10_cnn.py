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

import numpy as np
from matplotlib import pyplot
from scipy.misc import toimage

n_classes = 10
n_epochs = 200
steps_per_epoch = 32
batch_size = 32
data_augmentation = False

#read training & testing data
(x_train,y_train),(x_test,y_test) = cifar10.load_data()
print('Training Samples: ',x_train.shape[0])
print('Testing Samples: ',x_test.shape[0])
 
'''
for i in range(0,9):
	pyplot.subplot(330+1+i)
	pyplot.imshow(toimage(x_train[i]))
pyplot.show()
'''

y_train = keras.utils.to_categorical(y_train, n_classes)
y_test = keras.utils.to_categorical(y_test, n_classes)

model = Sequential()

model.add(Conv2D(32,(3,3),padding='same',input_shape=x_train.shape[1:]))
model.add(Activation('relu'))
model.add(Conv2D(32,(3,3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2,2)))
model.add(Dropout(0.25))

model.add(Conv2D(64,(3,3),padding='same'))
model.add(Activation('relu'))
model.add(Conv2D(64,(3,3)))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(2,2)))
model.add(Dropout(0.5))

model.add(Flatten())
model.add(Dense(512))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(n_classes))
model.add(Activation('softmax'))

opt = keras.optimizers.rmsprop(lr=0.0001,decay=1e-6)

model.compile(loss='categorical_crossentropy',optimizer=opt,metrics=['accuracy'])

x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train /= 255
x_test /= 255

hist = model.fit(x_train,y_train,
	batch_size=batch_size,
	epochs=n_epochs,
	validation_data=(x_test,y_test),
	shuffle=True)

