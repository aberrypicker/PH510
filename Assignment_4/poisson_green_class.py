#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 4: Individual Code file for the Class created for the Poisson Grid,
and Green's function calculations.

"""
#import time
#from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt


class PoissonGrid:
    """
    Defines a 2D version of the Poisson Equation, in a grid form of n by n. Between each spot
    lies areas where potential can be determined, and charges can be placed. These all combine 
    to affect Green's function, and random walks performed from points to the boundary.
    """
    def __init__(self, length, n_points, n_samples):
        self.l = length
        self.n = n_points
        self.h = self.l / (self.n - 1)
        self.phi = np.random.uniform(0, 10, (self.n, self.n))
        self.f = np.zeros((self.n, self.n))
        self.fixed_potential = set()
        self.fixed_charge = set()
        self.d = 2
        self.n_samples = n_samples

    def grid_potential(self, x, y, potential):
        """
        Function sets a fixed potential at each coordinate called within the grid, or raises
        an error if the called coordinate does not fall within the grid.
        """
        if 0 <= x < self.n and 0 <= y < self.n:
            self.phi[x,y] = potential
            self.fixed_potential.add((x,y))
        else:
            raise ValueError(f"({x}, {y}) Coordinates not in Grid.")
        return self.phi

    def boundary_condition(self, bc_type):
        """
        Function sets behaviour at edge of grid, termed as boundary conditions for the function.
        Mainly focuses on the 3 boundary conditions highlighted within the assignment, with 'Q4-a'
        representing the case that all borders are at a potential of 1V, 'Q4-b' representing the
        case that the top & bottom sides of the grid have a potential of 1V, and the left & right
        hand sides have a potential of -1V. 'Q4-c' is the third case, where the top and left hand
        hand sides have a potential of 1V, the bottom has a 2V potential, and the right hand side
        has a potential of -4V. The function only contains these three cases but has the 
        functionality to carry as many boundary condition cases as desired, so long as they are 
        defined by a title. Any condition not in the function that is attempted to be called 
        will return an error and statement that it does not recognise it.
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
            self.phi[:,-1] = -4
        else:
            raise ValueError(f"Boundary Condition not recognised: {bc_type}")

        for i in range(self.n):
            self.fixed_potential.add((self.n - 1, i))
            self.fixed_potential.add((0, i))
            self.fixed_potential.add((i, 0))
            self.fixed_potential.add((i, self.n - 1))

        return self.phi


    def charge_distribution_scenario(self, distribution_type):

        """
        Function changes the charge throughout the grid depending on desired cicumstances, mainly 
        focusing on the desired cases for the latter part of task 4.

        """
        x = np.linspace(0, self.l, self.n)
        y = np.linspace(0, self.l, self.n)
        x, y = np.meshgrid(x, y, indexing='ij')

        if distribution_type == 'uniform_10C':
            self.f[:, :] = 10

        elif distribution_type == 'linear_gradient_top_to_bottom':
            for i in range(self.n):
                self.f[i, :] = 1 - (i/(self.n - 1))

        elif distribution_type == 'exp_decay':
            x0, y0 = self.l/2, self.l/2
            r = np.sqrt((x - x0)**2 + (y - y0)**2)
            self.f[:, :] = np.exp(-2000 * np.abs(r))
        return self.f


    def overrelaxation_method(self, max_iteration=10000, tolerance = 1e-10):
        """
        Function which iteratively changes the potential at each grid point, by taking an average
        from each adjacent point, in all 4 cardinal directions relative to each point. It also
        accounts for boundary points, which only take points within the grid into consideration
        when averaging. The iteration process will only stop when the changes in potential between
        iterative steps are deemed to be minute enough that the potential across the grid has
        converged, or if the maximum number of iterations has taken place.
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
                    adj_term = 0.25 * self.h**2 * self.f[i, j]
                    adjacent_calculation = adj_term + np.mean(adjacent_charges)
                    self.phi[i, j] = (omega * adjacent_calculation) + ((1 - omega) * initial_phi)
                    max_delta = max(max_delta, abs(self.phi[i, j] - initial_phi))

            if max_delta < tolerance:
                print(f"The Grid has converged in {iteration} iterations.")
                print(np.round(self.phi, 2))
                break

        else:
            print(f"Maximum iterations ({max_iteration}) reached without convergence.")
        return self.phi


    def boundary_check(self, i, j):
        """
        Performs a check to determine if the walk is at the boundary of the grid.
        """
        return i == 0 or j == 0 or i == self.n - 1 or j == self.n - 1


    def random_walker(self,initial_i, initial_j):
        """
        Simulates the random chance for walk to travel in any cardinal direction when making its
        way to a boundary site.
        """
        values = []
        for _ in range(self.n_samples):
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
        return np.mean(values)


    def random_walk_probabilities(self, initial_i, initial_j):
        """
        Simulates random walks starting at the given set of coordinates (given by initial_i &
        initial_j), and returns the empirical probabilities of reaching each boundary point.
        """
        prob_grid = np.zeros((self.n, self.n))
        self.site_visits = np.zeros((self.n, self.n))
        boundary_hits = {}

        # Initialize count for each boundary point
        for i in range(self.n):
            boundary_hits[(0, i)] = 0       # Bottom
            boundary_hits[(self.n - 1, i)] = 0  # Top
            boundary_hits[(i, 0)] = 0       # Left
            boundary_hits[(i, self.n - 1)] = 0  # Right

        for _ in range(self.n_samples):
            i, j = initial_i, initial_j
            while not self.boundary_check(i, j):
                self.site_visits[(i, j)] += 1
                direction = np.random.choice(['up', 'down', 'left', 'right'])
                if direction == 'up':
                    i += 1
                elif direction == 'down':
                    i -= 1
                elif direction == 'left':
                    j -= 1
                elif direction == 'right':
                    j += 1
            self.site_visits[(i, j)] += 1
            boundary_hits[(i, j)] += 1

# These 'hits' are then used to fill the 2D probability grid, using the
# number of hits as a proportion of the total walks to give a probability.

        for (i, j), count in boundary_hits.items():
            prob_grid[i, j] = count / self.n_samples
        return prob_grid


    def get_potential(self, x, y):

        """
        By providing x and y coordinates, this function returns the known potential
        at the specified point so long as it is within the grid boundaries, with an
        error statement in case an outside coordinate is given.
        """

        if 0 <= x < self.n and 0 <= y < self.n:
            return self.phi[x, y]
        else:
            raise ValueError(f"Invalid coordinates: ({x}, {y}) outside grid bounds.")


    def greens_charge(self, initial_i, initial_j):
        """
        This function determines the charge at a given point in the grid for the purpose of
        determining its' potential using Green's function, as the charge component makes up
        half of the final determination of the potential Phi.
        """
        green_charge = np.zeros((self.n, self.n))
        i, j = initial_i, initial_j
        site_visits = self.random_walk_probabilities(i, j)[1]
        for p in range(0, self.n):
            for q in range(0, self.n):
                green_charge[p, q] = self.h**2/self.n_samples * self.site_visits[p, q]
        return green_charge

    def greens_function(self, starting_point_i, starting_point_j):
        """
        Uses above charge component and probability component of the grid to deliver 
        the Green's function.
        """
        i, j = starting_point_i, starting_point_j
        return self.random_walk_probabilities(i, j) + self.greens_charge(i, j)

    def greens_potential(self, initial_i, initial_j):
        """
        Function to determine the other part of Green's function, using the boundary
        probabilities as a summation, and the potential.
        """
        i, j = initial_i, initial_j
        greens_laplace = self.random_walk_probabilities(i, j)[0]
        term1 = np.zeros((self.n, self.n))
        for x_b in range(0, self.n):
            for y_b in range(0, self.n):
                if self.boundary_check(x_b, y_b):
                    term1[x_b, y_b] = greens_laplace[x_b, y_b] * self.phi[x_b, y_b]

        term1_sum = np.sum(term1)
        term2 = np.sum(self.greens_charge(i, j) * self.f)
        phi_greens = term1_sum + term2
        return phi_greens

    def grid_plot(self, value, title, dp):
        """
        Function to allow the various grids to be visualised in a colour mapped figure, and takes
        a title from the input to allow dynamic changing based on the circumstance. Also contains
        an SI conversion to ensure the grid is in cm.
        """
        plt.figure()
        extent = [0, self.l * 100, 0, self.l * 100] #cm conversion
        plt.imshow(np.round(value, dp), origin='lower', extent=extent, cmap='inferno')
        plt.colorbar(label='Potential (V)')
        plt.title(title)
        plt.xlabel("x (cm)")
        plt.ylabel("y (cm)")
        plt.grid(False)
        plt.show()
