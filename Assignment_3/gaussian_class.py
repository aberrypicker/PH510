#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

This file contains the code and documentation for the class which defines the normal distribution,
how it varies dimensionally, and how it can be substitutied to perform infinite integral calculations.

"""

import numpy as np

class Gaussian:
    """
    Classifies the Gaussian distribution function to interpret infitinte integral stuff.
    """
    def __init__(self, x0, sigma, n, d):
        """
        Create conditions for the Gaussian distribution. This includes the variables which
        make up the equation for the distribution's function such as its' initial location x0
        and distribution width sigma, the parameters which determine how many samples are taken
        across the function. Then n the number of random samples to be taken across the distribution
        and d the number of dimensions in which the distribution is situated. The transformation 
        variables are also defined in this initialisation, covering t and x.
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
