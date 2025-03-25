#!/usr/bin/env python3
"""
Version: Python 3.12.8

Licensed and Copyrghted 2025.

This file contains classes which allow Monte Carlo simulations to be ran using random
number sampling, and further classes for specific processes relating to assignment 3:
The Position class and the Gaussian class. Each class has documentation explaining
its' subclasses and processes.

"""
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
delta = 1/nproc

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
        self.samples = np.random.uniform(-1, 1, (self.n, self.d))

    def __str__(self):
        """
        Gives the printed information for the samples.
        """
        return f"Position:({self.samples})"

    def r(self):
        """
        Depending on the number of dimensions, turns the samples into vector positions
        in the form of r, for use in determining if inside or outside of unit shape.
        """
        r = (np.linalg.norm(self.samples, axis=1)).reshape(-1,1)
        check = (r <=  1).astype(int)
        return check



class MonteCarlo:
    """
    Approximates a given function by finding the expected value (average), then integrates this
    average across the limits of the function 'a' and 'b'. Finally, the function's variance (error)
    is found.
    """
    def __init__(self, function, classification, a, b):
        self.a = a
        self.b = b
        self.classification = classification
        self.function = function
        self.value = function()


    def __str__(self):
        """
        Confirms which function the Monte Carlo simulation is approximating
        """
        return f"Monte Carlo simulation of the function {self.function}"


    def average(self):
        """
        Calculates the average value of a given function.
        """
        value_2 = self.value**2
        average = 1/self.classification.n * np.sum(self.value)
        average_2 = 1/self.classification.n * np.sum(value_2)
        return average, average_2

    def parallel_version(self):
        """
        Combines integral calculation of function, variance of function, and other function
        calculations into single sub-class.
        """
        mean = comm.reduce(self.average()[0], op=MPI.SUM, root=0)
        mean_squared = comm.reduce(self.average()[1], op=MPI.SUM, root=0)

        if rank==0:
            integral = ((self.b - self.a)**self.classification.d * mean)/ nproc

            variance_1 = (delta * mean)**2
            variance_2 = delta * mean_squared
            variance = 1/self.classification.n * (variance_2 - variance_1)
            uncertainty = np.sqrt(variance) * ((self.b - self.a)**self.classification.d)

            print(f"Average = {mean * delta}, Integral = {integral},",
            f"Variance = {variance}")
            print(f"The {self.classification.d}D Integral is {integral:.4f} ± {uncertainty:.4f}",
            f"units**{self.classification.d}")
            print()
            return integral
        return None

class Gaussian:
    """
    Classifies the Gaussian distribution function to interpret infitinte integral stuff.
    """
    def __init__(self, x0, sigma, n, d):
        """
        Create conditions for the Gaussian distribution.
        """
        self.x0 = x0
        self.sigma = sigma
        self.n = n
        self.d = d
        self.t = np.random.uniform(-1, 1, (self.n, self.d))
        self.x = self.t/(1 - self.t**2)


    def __str__(self):
        """
        Gives the printed information for the Distribution.
        """
        return f"Gaussian:({self.x0}, {self.sigma})" #placeholder

    def transformation(self, t):
        """
        Integration by Substition multiplication term.
        """
        return (1 + t**2)/((1 - t**2)**2)

    def integral(self):
        """
        Performs integration by substitution on Gaussian Function.
        """
        term_1 = (1/(self.sigma * np.sqrt(2*np.pi)))
        term_2 = np.exp(np.sum(-(abs(self.x - self.x0))**2, axis = 1)/(2 * self.sigma**2))
        return term_1 * term_2 * np.prod(self.transformation(self.t), axis = 1)
