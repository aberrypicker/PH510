#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
from mpi4py import MPI

class MonteCarlo:
	"""
	A system for integration that uses the expectation value of a function, and 'random' points
	along the path to determine the area?
	"""
	def __init__(self,

