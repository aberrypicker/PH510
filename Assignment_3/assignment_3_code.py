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
	def __init__(self, n, d, x, y, z, a, b):
		"""
		Initialise the position using the number of samples and number of dimensions.
		"""
		self.n = n
		self.d = d
		self.x, self.y, self.z, self.a, self.b = np.random.uniform(-1,1,5)
		if d == 1: r = x
		if d == 2: r = math.sqrt(x**2 + y**2)
		if d == 3: r = math.sqrt(x**2 + y**2 + z**2)
		if d == 4: r = math.sqrt(x**2 + y**2 + z**2 + a**2)
		if d == 5: r = math.sqrt(x**2 + y**2 + z**2 + a**2 + b**2)

Pos_1 = Position(2, 2, 0, 0, 0, 0, 0) 
print(Pos_1)


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

a =  MonteCarlo(2,10)
print(a.data_array)







