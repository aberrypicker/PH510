#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

This file contains a class which allow Monte Carlo simulations to be ran using random
number sampling. Each subclass has documentation explaining its' subclasses and processes.

"""
import time
from mpi4py import MPI
import numpy as np
import gaussian_class as ga

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
delta = 1/nproc

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
        Combines integral calculation of the function, variance of function, and other function
        calculations into single sub-class, which is coded with MPI in mind such that it allows
        the code to be ran using multiple processors in parallel. It involves determinations
        for the number of processors involved to divide tasks fairly and efficiently. What it 
        returns will be the integral of the function, its variance, using the expectation value
        squared and the expectation of the squared mean value, to give the overall variance. It
        then uses this to determine the overall uncertainty for the function's integral. 
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
    def importance_sampling(self):
        """
        Makes use of the idea of importance sampling to reduce error in variance calculation versus
        regular Monte Carlo.
        """
        y = 1 + self.value
        mean = comm.reduce(self.average()[0], op=MPI.SUM, root=0)
        mean_squared = comm.reduce(self.average()[1], op=MPI.SUM, root=0)

        if rank==0:
            term = (self.b - self.a)**self.classification.d * mean
            weighted_integral1 = term * np.linalg.norm(y, 0)/np.linalg.norm(y, 0)
            weighted_integral2 = (term * np.linalg.norm(y, 0)/np.linalg.norm(y, 0))/nproc

            variance_1 = (delta * mean)**2
            variance_2 = delta * mean_squared
            variance = 1/self.classification.n * (variance_2 - variance_1)
            uncertainty = np.sqrt(variance) * term

            print(f"Average = {mean * delta}, Weighted_Integral = {weighted_integral2},",
            f"Variance = {variance}")
            print(f"The {self.classification.d}D Integral is {weighted_integral2:.4f} ± {uncertainty:.4f}",
            f"units**{self.classification.d}")
            print()
            return weighted_integral2
        return None


if rank==0:
    print(f"{nproc} Processor(s)")
    print()
    start_time = time.time()

# Define a number of samples to be more efficient rather than creating 100000000
# different sample arrays.
num_samples = np.int32(100000000/nproc)


dist_10 = ga.Gaussian(0, 1, num_samples, 1)

Monte_1D0 = MonteCarlo(dist_10.integral, dist_10, -1, 1)
if rank==0:
    print()
    print("Test Case 2")
    print()
    print()
    print(f"1D Normal Distribution, [sigma = {dist_10.sigma}], [x0 = {dist_10.x0}]")
Monte_1D0P = Monte_1D0.parallel_version()

Monte_1D0I = Monte_1D0.importance_sampling()
