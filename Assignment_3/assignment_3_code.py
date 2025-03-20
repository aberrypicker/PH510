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
        self.samples = np.random.uniform(-1, 1, (self.n, self.d))

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


    def integral(self):
        """
        Calculates the integral of a given function  between limits 'a' and 'b'.  
        """
        integral = (self.b - self.a)**self.classification.d * self.average()[0]
        return integral


    def variance(self):
        """
        Calculates the variance (error) of a given function.
        """
        variance = 1/self.classification.n * (self.average()[1] - self.average()[0]**2)
        return variance


    def calculations(self):
        """
        Returns all three desired values at once.
        """
        return self.average()[0], self.integral(), self.variance()

class Gaussian:
    """
    Classifies the Gaussian distribution function to interpret infitinte integral stuff.
    """
    def __init__(self, x0, sigma, n):
        """
        Create conditions for distribution.
        """
        self.x0 = x0
        self.sigma = sigma
        self.n = n

sample_3D = Position(100000000, 3)
#print(sample_1.samples)
#print()
#print(Position.r(sample_1))

Monte_3D = MonteCarlo(sample_3D.r, sample_3D, -1, 1)
print("For a 3D Monte Carlo, the respective average, integral, and variance are:"f"{Monte_3D.calculations()[0]:.4f},",
      f"{Monte_3D.calculations()[1]:.4f},", f"& {Monte_3D.calculations()[2]:.4f}")
