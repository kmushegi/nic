"""
Evolving Neural Networks Using Genetic Algorithms - Project 5
Nature Inspired Computation
Spring 2017
Stephen Majercik

Ernesto Garcia, Marcus Christiansen, Konstantine Mushegian

A classic setup.py script typically included with Python projects to manage
dependencies and store author information.
"""

from setuptools import setup, find_packages

setup(
	name='Project 5',
	version='1.0',
	description='Evolving Neural Networks Using Genetic Algorithms',
	url='https://github.com/kmushegi/nic/tree/master/p5', #currently private
	author='Konstantine Mushegian',
	author_email='kmushegian@gmail.com',
	install_requires=['numpy','keras'],
	)
