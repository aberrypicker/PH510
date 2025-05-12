#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 4: Using Green's Function and Random Walks to solve Poisson's Equation
			  for multiple examples.

"""
import time
from mpi4py import MPI
import numpy as np
import monte_carlo_class as mc


class Poisson_Grid:
    """
    Grid thing
    """
    def __init__(self, length, n_points):
        self.l = length
        self.n = n_points
        self.h = self.l / (self.n - 1)
        self.phi = np.random.uniform(0, 1000, (self.n, self.n))
        self.f = np.zeros((self.n, self.n))
        self.fixed_potential = set()
        self.fixed_charge = set()

    def grid_potential(self, x, y, potential):
        """
        Looks like volatge but isnt voltage, why is the unit a V fissics mfckers
        """
        if 0 <= x < self.n and 0 <= y < self.n:
            self.phi[x,y] = potential
            self.fixed_potential.add((x,y))
        else:
            raise ValueError(f"({x}, {y}) Coordinates not in Grid.")
        return self.phi

    def boundary_condition(self, bc_type):
        """
        Fucntion behaviviour at edge.
        """
        if bc_type == 'Q4-a':
            self.phi[-1,:] = 1
            self.phi[0,:] = 1
            self.phi[:,0] = 1
            self.phi[:,-1] = 1
        elif bc_type == 'Q4-b':
            self.phi[-1,:] = 1
            self.phi[0,:] = 1
            self.phi[:,0] = -1
            self.phi[:,-1] = -1
        elif bc_type == 'Q4-c':
            self.phi[-1,:] = 2
            self.phi[0,:] = 0
            self.phi[:,0] = 2
            self.phi[:,-1] = 4
        else:
            raise ValueError(f"Boundary Condition not recognised: {bc_type}")

        for i in range(self.n):
            self.fixed_potential.add((self.n - 1, i))
            self.fixed_potential.add((0, i))
            self.fixed_potential.add((i, 0))
            self.fixed_potential.add((i, self.n - 1))

        return self.phi

def relax/overrelax

boundary check

grid plot
