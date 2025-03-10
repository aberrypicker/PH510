#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
from mpi4py import MPI
import math
import numpy as np


class Position:
	"""
	A class which uses the dimensions and number of points to determine where the point lies 
	and relate this to the position of a unit circle inside a square, and so forth.
	"""
	def __init__(self, n, d):
		"""
		Initialise the position using the number of samples and number of dimensions.
		"""
		self.n = n
		self.d = d
		j = 0 # dimension loop iterator
		samples = np.zeros((self.n, self.d))
		while j < self.n:
			i = 0 # sample loop iterator
			samples_randomisation = np.random.uniform(-1,1,self.d)
			samples_randomisation = samples_randomisation.reshape((1, self.d))
			while i < self.d:
				samples[j][i] = samples_randomisation[0][i]
				i = i + 1
			j = j + 1
		self.samples = samples

	def r(self):
		"""
		Perform the necessary calculations to determine the vector position depending on dimensionality
		"""
		
sample_1 = Position(20, 2)
print(sample_1.samples)

class MonteCarlo:
	def __init__(self, nd, ns):
		self.nd = nd
		self.ns = ns
		self.data_array = np.zeros((self.ns,self.nd))

		i = 0 # sample iterator
		j = 0 # dimension iterator 

		while i < self.ns:
			j = 0
			while j < self.nd:
				self.data_array[i][j] = np.random.uniform(-1,1,1)
				j+=1
			i += 1

#a =  MonteCarlo(2,20)
#print(a.data_array)







