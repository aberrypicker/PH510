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
    Defines a 2D version of the Poisson Equation, in a grid form of n by n. Between each spot
    lies areas where potential can be determined, and charges can be placed. These all combine 
    to affect Green's function, and random walks performed from points to the boundary.
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
        Function sets the potential at each coordinate within the grid.
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


    def charge_distribution_scenario(self, distribution_type):

        """
        Function changes the charge throughput the grid depending on desired cicumstances.

        """
        x = np.linspace(0, self.l, self.n)
        y = np.linspace(0, self.l, self.n)
        x, y = np.meshgrid(x, y, indexing='ij')

        if distribution_type == 'uniform_10C':
            charge = 10
            self.f[:, :] = charge

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
                    adjacent_calculation = (0.25 * self.h**2 * f[i, j]) + np.mean(adjacent_charges)
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


    def random_walker(self,initial_i, initial_j, n_walks = 1000):
        """
        Simulates random chance for walk to travel in any cardinal direction when making its
        way to boundary.
        """
        values = []
        for _ in range(n_walks):
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


    def random_walk_probabilities(self, initial_i, initial_j, n_walks=100000):
        """
        Simulates random walks starting at the given set of coordinates (given by initial_i &
        initial_j), and returns the empirical probabilities of reaching each boundary point.
        """
        prob_grid = np.zeros((self.n, self.n))
        site_visits = np.zeros((self.n, self.n))
        boundary_hits = {}

        # Initialize count for each boundary point
        for i in range(self.n):
            boundary_hits[(0, i)] = 0       # Bottom
            boundary_hits[(self.n - 1, i)] = 0  # Top
            boundary_hits[(i, 0)] = 0       # Left
            boundary_hits[(i, self.n - 1)] = 0  # Right

        for _ in range(n_walks):
            i, j = initial_i, initial_j
            while not self.boundary_check(i, j):
                site_visits[(i, j)] += 1
                direction = np.random.choice(['up', 'down', 'left', 'right'])
                if direction == 'up':
                    i += 1
                elif direction == 'down':
                    i -= 1
                elif direction == 'left':
                    j -= 1
                elif direction == 'right':
                    j += 1
                site_visits[(i, j)] += 1
            boundary_hits[(i, j)] += 1

# These 'hits' are then used to fill the 2D probability grid, using the
# number of hits as a proportion of the total walks to give a probability.

        for (i, j), count in boundary_hits.items():
            prob_grid[i, j] = count / n_walks
        return prob_grid, site_visits


#    def potential_check(self, point_i, point_j, n_walks_per_point=500):
#        """
#        Code to perform random walks to determine the potential at a specific desired point.
#        """
#        total = 0.0
#        for i in range(1, self.n - 1):
#            for j in range(1, self.n - 1):
#                if self.f[i, j] == 0:
#                    continue
#                green_function_value, _ = self.random_walker(point_i, point_j, n_walks_per_point)
#                total += green_function_value * self.f[i, j] * self.h**2
#        return total

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


    def greens_function(self, initial_i, initial_j, n_walks=100000):
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
                green_charge[p, q] = self.h**2/n_walks * site_visits[p, q]
        return green_charge

    def greens_potential(self, initial_i, initial_j, n_walks=100000):
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
        term2 = np.sum(self.greens_function(i, j) * self.f)
        phi_greens = term1_sum + term2

        return phi_greens, greens_laplace, self.phi, term1, term2

    def grid_plot(self, value, title):
        """
        Function to allow the various grids to be visualised in a colour mapped figure.
        """
        plt.figure()
        extent = [0, self.l * 100, 0, self.l * 100] #cm conversion
        plt.imshow(np.round(value, 4), origin='lower', extent=extent, cmap='inferno')
        plt.colorbar(label='Potential (V)')
        plt.title(title)
        plt.xlabel("x (cm)")
        plt.ylabel("y (cm)")
        plt.grid(False)
        plt.show()


# Question 3

# Performs random walk and evaluates the Green's function via probability for each specified
# point from the assignment.

print("Exercise 3")
init_grid = PoissonGrid(0.10, 21)

phi, f = init_grid.phi, init_grid.f
print("Green's function evaluation for a square grid of side length 10cm:")

# (a)
a = init_grid.random_walk_probabilities(10, 10)
print(f"At centre point (5cm, 5cm):\n{a[0]}")

#init_grid.grid_plot(a[0], "Green's Function")
#init_grid.grid_plot(a[1], 'Number of Site Visits')

# Q3b
b = init_grid.random_walk_probabilities(5, 5)
print(f"At (2.5cm, 2.5cm):\n{b[0]}")

#init_grid.grid_plot(b[0], "Green's Function")
#init_grid.grid_plot(b[1], 'Number of Site Visits')


# Q3c
c = init_grid.random_walk_probabilities(1, 5)
print(f"At (0.1cm, 2.5cm):\n{c[0]}")

#init_grid.grid_plot(c[0], "Green's Function")
#init_grid.grid_plot(c[1], 'Number of Site Visits')


# Q3d
d = init_grid.random_walk_probabilities(1, 1)
print(f"At (0.1cm, 0.1cm):\n{d[0]}")

#init_grid.grid_plot(d[0], "Green's Function")
#init_grid.grid_plot(d[1], 'Number of Site Visits')
print()

# Question 4

# Now evaluates the same cases as Question 3 but now with additional Potentials added
# via boundary conditions, with the latter part of the question considering further
# permutations with charge placed within the grid in various configurations.

print("Exercise 4")
print("Potential calculation via Green's function for a square grid of side length 10cm:")

# Q4a
print("Q4a with boundary conditions: All edges uniformly at +1V")
init_grid = PoissonGrid(0.10, 21)

phi_1 = init_grid.phi
phi_1 = init_grid.boundary_condition('Q4-a')
phi_1 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]}V")
print()


# Q4b
print("Q4b) with boundary conditions: Top and bottom edges at +1V, left and right edges at -1V")
init_grid = PoissonGrid(0.10, 21)


phi_2 = init_grid.phi
phi_2 = init_grid.boundary_condition('Q4-b')
phi_2 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]}V")
print()

 
# Q4c
print("Q4c) with boundary conditions: Top and left edges at +2V, bottom edge at 0V and right edge at -4V")
init_grid = PoissonGrid(0.10, 21)


phi_3 = init_grid.phi
phi_3 = init_grid.boundary_condition('Q4-c')
phi_3 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]}V")
print()

# Question 4, Part 1a
#example1= PoissonGrid(0.1, 9)
#phi1, f = example1.phi, example1.f
#phi1 = example1.boundary_condition('Q4-a')
#phi1 = example1.grid_potential(4, 4, 0)
#phi1 = example1.overrelaxation_method()
#example1.grid_plot()

# Question 4, Part 1b
#example2= PoissonGrid(0.1, 9)
#phi2, f = example2.phi, example2.f
#phi2 = example2.boundary_condition('Q4-b')
#phi2 = example2.grid_potential(4, 4, 0)
#phi2 = example2.overrelaxation_method()
#example2.grid_plot()

# Question 4, Part 1c
#example3= PoissonGrid(0.1, 9)
#phi3, f = example3.phi, example3.f
#phi3 = example3.boundary_condition('Q4-c')
#phi3 = example3.grid_potential(25, 25, 0)
#phi3 = example3.overrelaxation_method()
#example3.grid_plot()
#print(phi)
#print(example.fixed_potential)

#example1.random_walker(4,4)
#walk = example1.random_walker(4,4, 100)
#print(walk)
