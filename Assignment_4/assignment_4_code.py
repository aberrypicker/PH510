#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 4: Using Green's Function and Random Walks to solve Poisson's Equation
			  for multiple examples.

"""
#import time
#from mpi4py import MPI
#import numpy as np
#import matplotlib.pyplot as plt
#import monte_carlo_class as mc
import poisson_green_class as pg

# Question 3

# Performs random walk and evaluates the Green's function via probability for each specified
# point from the assignment.

print("Exercise 3")
init_grid = pg.PoissonGrid(0.10, 21)

phi, f = init_grid.phi, init_grid.f
print("Green's function evaluation for a square grid of side length 10cm:")

# (a)
a = init_grid.random_walk_probabilities(10, 10)
print(f"At centre point (5cm, 5cm):\n{a[0]}")

#init_grid.grid_plot(a[0], "Green's Function: (5cm, 5cm)")
#init_grid.grid_plot(a[1], 'Number of Site Visits: (5cm, 5cm)')

# Q3b
b = init_grid.random_walk_probabilities(5, 5)
print(f"At (2.5cm, 2.5cm):\n{b[0]}")

#init_grid.grid_plot(b[0], "Green's Function: (2.5cm, 2.5cm)")
#init_grid.grid_plot(b[1], 'Number of Site Visits: (2.5cm, 2.5cm)')


# Q3c
c = init_grid.random_walk_probabilities(1, 5)
print(f"At (0.1cm, 2.5cm):\n{c[0]}")

#init_grid.grid_plot(c[0], "Green's Function: (0.1cm, 2.5cm)")
#init_grid.grid_plot(c[1], 'Number of Site Visits: (0.1cm, 2.5cm)')


# Q3d
d = init_grid.random_walk_probabilities(1, 1)
print(f"At (0.1cm, 0.1cm):\n{d[0]}")

#init_grid.grid_plot(d[0], "Green's Function: (0.1cm, 0.1cm)")
#init_grid.grid_plot(d[1], 'Number of Site Visits: (0.1cm, 0.1cm)')
print()

# Question 4

# Now evaluates the same cases as Question 3 but now with additional Potentials added
# via boundary conditions, with the latter part of the question considering further
# permutations with charge placed within the grid in various configurations.

print("Exercise 4")
print("Potential calculation via Green's function for a square grid of side length 10cm:")

# Q4a
print("Q4a with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_1 = init_grid.phi
phi_1 = init_grid.boundary_condition('Q4-a')
phi_1 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# Q4b
print("Q4b) with boundary conditions: Top and bottom edges at +1V, left and right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21)


phi_2 = init_grid.phi
phi_2 = init_grid.boundary_condition('Q4-b')
phi_2 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# Q4c
print("Q4c) with boundary conditions: Top/left edges at +2V, bottom edge at 0V, right edge at -4V")
init_grid = pg.PoissonGrid(0.10, 21)


phi_3 = init_grid.phi
phi_3 = init_grid.boundary_condition('Q4-c')
phi_3 = init_grid.overrelaxation_method()

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# Question 4 - Part 2: These are the same calclualtions but now with charge added
# into consideration. The cases are labelled the same as above, but with "d, e, f"
# used to distinguish from "a, b, c" which are used in the assignment instructions,
# and "i, ii, iii" to distinguish for the 3 boundary conditions

# Q4d
print("Q4d) Repeat, with uniform charge of 10C throughout the grid")
print()

# (i)
print("i) with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_4 = init_grid.phi
phi_4 = init_grid.boundary_condition('Q4-a')
phi_4 = init_grid.overrelaxation_method()

f_4 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()

# (ii)
print("ii) with boundary conditions: Top and bottom edges at +1V, left and right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21)


phi_5 = init_grid.phi
phi_5 = init_grid.boundary_condition('Q4-b')
phi_5 = init_grid.overrelaxation_method()

f_5 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# (iii)
print("iii) with boundary conditions: Top/left edges at +2V, bottom edge at 0V, right edge at -4V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_6 = init_grid.phi
phi_6 = init_grid.boundary_condition('Q4-c')
phi_6 = init_grid.overrelaxation_method()

f_6 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# Q4e
print("Q4e) Repeat, with uniform charge gradient from 1C at top to 0C at bottom")
print()

# (i)
print("i) with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21)


phi_7 = init_grid.phi
phi_7 = init_grid.boundary_condition('Q4-a')
phi_7 = init_grid.overrelaxation_method()

f_7 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# (ii)
print("ii) with boundary conditions: Top/bottom edges at +1V, left and right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21)


phi_8 = init_grid.phi
phi_8 = init_grid.boundary_condition('Q4-b')
phi_8 = init_grid.overrelaxation_method()

f_8 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()



# (iii)
print("iii) with boundary conditions: Top/left edges at +2V, bottom edge at 0V, right edge at -4V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_9 = init_grid.phi
phi_9 = init_grid.boundary_condition('Q4-c')
phi_9 = init_grid.overrelaxation_method()

f_9 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# Q4f
print("f) Repeat, with exponentially decaying charge exp(-2000|r|) placed at centre of grid")
print()

# (i)
print("i) with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_10 = init_grid.phi
phi_10 = init_grid.boundary_condition('Q4-a')
phi_10 = init_grid.overrelaxation_method()

f_10 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# (ii)
print("ii) with boundary conditions: Top/bottom edges at +1V, left and right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_11 = init_grid.phi
phi_11 = init_grid.boundary_condition('Q4-b')
phi_11 = init_grid.overrelaxation_method()

f_11 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()


# (iii)
print("iii) with boundary conditions: Top/left edges at +2V, bottom edge at 0V, right edge at -4V")
init_grid = pg.PoissonGrid(0.10, 21)

phi_12 = init_grid.phi
phi_12 = init_grid.boundary_condition('Q4-c')
phi_12 = init_grid.overrelaxation_method()

f_12 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50 = init_grid.greens_potential(10, 10)
potential_25_25 = init_grid.greens_potential(5, 5)
potential_1_25 = init_grid.greens_potential(1, 5)
potential_1_1 = init_grid.greens_potential(1, 1)

print(f"At (5cm, 5cm): {potential_50_50[0]:.4f}V")
print(f"At (2.5cm, 2.5cm): {potential_25_25[0]:.4f}V")
print(f"At (0.1cm, 2.5cm): {potential_1_25[0]:.4f}V")
print(f"At (0.1cm, 0.1cm): {potential_1_1[0]:.4f}V")
print()
