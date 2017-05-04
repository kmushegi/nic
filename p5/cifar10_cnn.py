from keras.datasets import cifar10
from matplotlib import pyplot
from scipy.misc import toimage

(x_train,y_train),(x_test,y_test) = cifar10.load_data()

print("data loaded")

for i in range(0,9):
	pyplot.subplot(330+1+i)
	pyplot.imshow(toimage(x_train[i]))
pyplot.show()