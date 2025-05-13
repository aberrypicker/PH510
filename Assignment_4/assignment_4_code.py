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


class PoissonGrid:
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
        k = 0
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
                k = k + 1
        return np.mean(values), np.std(values)

    def potential_check(self, point_i, point_j, n_walks_per_point=500):
        """
        Code to perform random walks to determine the potential at a specific desired point.
        """
        total = 0.0
        for i in range(1, self.n - 1):
            for j in range(1, self.n - 1):
                if self.f[i, j] == 0:
                    continue
                green_function_value, _ = self.random_walker(point_i, point_j, n_walks_per_point)
                total += green_function_value * self.f[i, j] * self.h**2
        return total

example1= PoissonGrid(0.1, 50)
phi1, f = example1.phi, example1.f
phi1 = example1.boundary_condition('Q4-a')
phi1 = example1.grid_potential(25, 25, 5)
phi1 = example1.overrelaxation_method()
example1.grid_plot()

example2= PoissonGrid(0.1, 50)
phi2, f = example2.phi, example2.f
phi2 = example2.boundary_condition('Q4-b')
phi2 = example2.grid_potential(25, 25, 5)
phi2 = example2.overrelaxation_method()
example2.grid_plot()

example3= PoissonGrid(0.1, 50)
phi3, f = example3.phi, example3.f
phi3 = example3.boundary_condition('Q4-c')
phi3 = example3.grid_potential(25, 25, 5)
phi3 = example3.overrelaxation_method()
example3.grid_plot()
#print(phi)
#print(example.fixed_potential)



#example.random_walker(25,25)
#walk = example.random_walker(20,10, 10000)
#print(walk)
