#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
#from mpi4py import MPI
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
        w = 0 # dimension loop iterator
        samples = np.zeros((self.n, self.d))
        while w < self.n:
            v = 0 # sample loop iterator
            samples_randomisation = np.random.uniform(-1,1,self.d)
            samples_randomisation = samples_randomisation.reshape((1, self.d))
            while v < self.d:
                samples[w][v] = samples_randomisation[0][v]
                v = v + 1
            w = w + 1
        self.samples = samples

    #def r(self):
        """
        Perform the necessary loop calculations to determine the vector position
        of each sample row depending on dimensionality.
        """
        #w = 0 
        #while w < self.n:


sample_1 = Position(20, 3)
print(sample_1.samples)

class MonteCarlo:
    """
    Approximates a given function by finding the expected value (average), then integrates this
    average across the limits of the function 'a' and 'b'. Finally, the function's variance (error)
    is found.
    """
    def __init__(self, function, a, b, n):
        self.a = a
        self.b = b
        self.n = n
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
        average = 1/self.n * np.sum(self.value)
        average_2 = 1/self.n * np.sum(value_2)
        return average, average_2


    def integral(self):
        """
        Calculates the integral of a given function  between limits 'a' and 'b'.  
        """
        integral = (self.b - self.a) * self.average()[0]
        return integral


    def variance(self):
        """
        Calculates the variance (error) of a given function.
        """
        variance = 1/self.n * (self.average()[1] - self.average()[0]**2)
        return variance


    def calculations(self):
        """
        Returns all three desired values at once.
        """
        return self.average()[0], self.integral(), self.variance()

print(sample_1.calculations)
