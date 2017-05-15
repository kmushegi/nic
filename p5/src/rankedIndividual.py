"""
Neural Networks for Image Classification - Final Project
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

This file is part of Neural Networks for Image Classification, Final Project. This file
contains the implementation of NEED TO WRITE THIS
"""

import network as nn
import cnetwork as cnn

class rankedIndividual(object):

	def __init__(self , network):
		self.individual = network
		self.fitness = 0
		self.probability = 0

	def setFitness(self, newFitness):
		self.fitness = newFitness

	def setProbability(self, newProbability):
		self.probability = newProbability

	def getIndividual(self):
		return self.individual

	def getFitness(self):
		return self.fitness

	def getProbability(self):
		return self.probability