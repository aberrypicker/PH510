#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
from mpi4py import MPI


class Position:
	"""
	A class which uses the dimensions and number of points to determine where the point lies 
	and relate this to the position of a unit circle inside a square, and so forth.
	"""
	def __init__(self,n,d):
		"""
		Initialise the position using the number of samples and number of dimensions.
		"""
		self.n = n
		self.d = d

	


class MonteCarlo:
	"""
	A system for integration that uses the expectation value of a function, and 'random' points
	along the path to determine the area?
	"""
	def __init__(self,function, a, b):
		

