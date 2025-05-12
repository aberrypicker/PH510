#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 4: Using Green's Function and Random Walks to solve Poisson's Equation
			  for multiple examples.

"""
#import time
#from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
#import monte_carlo_class as mc


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

    def overrelaxation_method(self, max_iteration=10000, tolerance = 1e-10):
        """
        Meditation or smthn i have legit no scoobies.
        """
        omega = 2/(1 + np.sin(np.pi/self.n))
        for iteration in range(max_iteration):
            max_delta = 0
            for i in range(0, self.n):
                for j in range(0, self.n):
                    if (i, j) in self.fixed_potential:
                        continue

                    adjacent_charges = []

                    if i + 1 < self.n:
                        adjacent_charges.append(self.phi[i+1, j])

                    if i - 1 >= 0:
                        adjacent_charges.append(self.phi[i-1, j])

                    if j + 1 < self.n:
                        adjacent_charges.append(self.phi[i, j+1])

                    if j - 1 >= 0:
                        adjacent_charges.append(self.phi[i, j-1])

                    initial_phi = self.phi[i, j]
                    adjacent_calculation = -(self.h**2 * f[i, j]) + np.mean(adjacent_charges)
                    self.phi[i, j] = (omega * adjacent_calculation) + ((1 - omega) * initial_phi)
                    max_delta = max(max_delta, abs(self.phi[i, j] - initial_phi))

            if max_delta < tolerance:
                print(f"The Grid has converged in {iteration} iterations.")
                print(np.round(self.phi, 2))
                break

        else:
             print(f"Maximum iterations ({max_iteration}) reached without convergence.")
        return self.phi


    def grid_plot(self):
        """
        Allows the poisson grid to be visualised.
        """
        extent = [0, self.l * 100, 0, self.l * 100] #cm conversion
        plt.imshow(np.round(self.phi, 4), origin='lower', extent=extent, cmap='inferno')
        plt.colorbar(label='Potential (V)')
        plt.title("Potential Distribution")
        plt.xlabel("x (cm)")
        plt.ylabel("y (cm)")
        plt.grid(False)
        plt.show()

    def boundary_check(self, i, j):
        """
        performs a check to determine if the position of the walker is at the boundary of the grid.
        """
        return i == 0 or j == 0 or i == self.n - 1 or j == self.n - 1


    def random_walker(self,initial_i, initial_j, n_walks = 100000):
        """
        Simulates random chance for walk to travel in any cardinal direction when making its
        way to boundary.
        """
        values = []
        for k in range(n_walks):
            i, j = initial_i, initial_j
            while not self.boundary_check(i, j):
                direction = np.random.choice(['up', 'down', 'left', 'right'])
                if direction == 'up':
                    i += 1
                elif direction == 'down':
                    i -= 1
                elif direction == 'left':
                    j -= 1
                elif direction == 'right':
                    j += 1
            if self.boundary_check(i, j):
                values.append(self.phi[i, j])
        return np.mean(values), np.std(values)

example = Poisson_Grid(0.1, 50)
phi, f = example.phi, example.f
phi = example.boundary_condition('Q4-a')
phi = example.grid_potential(25, 25, 5)
print(phi)
print(example.fixed_potential)
phi = example.overrelaxation_method()
#example.grid_plot()
example.random_walker(25,25)
print(example.random_walker(20,10)) 





