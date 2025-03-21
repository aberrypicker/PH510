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
        Returns all four desired values at once.
        """
        return self.average()[0], self.integral(), self.variance(), np.sqrt(self.variance())

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

dist_1 = Gaussian(0, 1, 100000000, 1)

MonteGaus_1 = MonteCarlo(dist_1.integral, dist_1, -1, 1)
print()
print("For a 1D Gaussian Function, the respective average, integral, variance, & uncertainty are: ",
      f"{MonteGaus_1.calculations()[0]:.4f},", f"{MonteGaus_1.calculations()[1]:.4f},",
      f"{MonteGaus_1.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")
print()

dist_6 = Gaussian(0, 1, 100000000, 6)

MonteGaus_6 = MonteCarlo(dist_6.integral, dist_6, -1, 1)
print()
print("For a 6D Gaussian Function, the respective average, integral, variance, & uncertainty are: ",
      f"{MonteGaus_6.calculations()[0]:.4f},", f"{MonteGaus_6.calculations()[1]:.4f},",
      f"{MonteGaus_6.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")
print()

sample_2D = Position(1000, 2)
sample_3D = Position(1000, 3)
sample_4D = Position(1000, 4)
sample_5D = Position(1000, 5)

Monte_2D = MonteCarlo(sample_2D.r, sample_2D, -1, 1)
print("For a 2D Circle, the respective average, integral, variance, & uncertainty are: ",
      f"{Monte_2D.calculations()[0]:.4f},", f"{Monte_2D.calculations()[1]:.4f},",
      f"{Monte_2D.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")

Monte_3D = MonteCarlo(sample_3D.r, sample_3D, -1, 1)
print("For a 3D Sphere, the respective average, integral, variance, & uncertainty are: ",
      f"{Monte_3D.calculations()[0]:.4f},", f"{Monte_3D.calculations()[1]:.4f},",
      f"{Monte_3D.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")

Monte_4D = MonteCarlo(sample_4D.r, sample_4D, -1, 1)
print("For a 4D Sphere, the respective average, integral, variance, & uncertainty are: ",
      f"{Monte_4D.calculations()[0]:.4f},", f"{Monte_4D.calculations()[1]:.4f},",
      f"{Monte_4D.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")

Monte_5D = MonteCarlo(sample_5D.r, sample_5D, -1, 1)
print("For a 5D HyperSphere, the respective average, integral, variance, & uncertainty are: ",
      f"{Monte_5D.calculations()[0]:.4f},", f"{Monte_5D.calculations()[1]:.4f},",
      f"{Monte_5D.calculations()[2]:.4f},", f"& {MonteGaus_1.calculations()[3]:.4f},")
