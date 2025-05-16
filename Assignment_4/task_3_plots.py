#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 4: Code to provide graphical plots for task 3 as the job script does not
like attempting to show plots.

"""

import numpy as np
import matplotlib.pyplot as plt
import poisson_green_class as pg


# Question 3

# Performs random walk and evaluates the Green's function via probability for each specified
# point from the assignment.
n_samples = np.int32(100000)

print("Exercise 3")
print("Green's function evaluation for a square grid of side length 10cm:")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

#phi, f = init_grid.phi, init_grid.f


# (a)
a = init_grid.greens_function(10,10)

init_grid.grid_plot(a, "Edge Potential Green's Function at (5.0cm, 5.0cm)", 10)
init_grid.grid_plot(init_grid.greens_charge(10,10), "Charge Potential at (5.0cm, 5.0cm)", 10)

# Q3b
b = init_grid.greens_function(5,5)

init_grid.grid_plot(b, "Edge Potential Green's Function at (2.5cm, 2.5m)", 10)
init_grid.grid_plot(init_grid.greens_charge(5,5), "Charge Potential at (2.5cm, 2.5cm)", 10)


# Q3c
c = init_grid.greens_function(1,5)

init_grid.grid_plot(c, "Edge Potential Green's Function at (0.1cm, 2.5cm)", 4)
init_grid.grid_plot(init_grid.greens_charge(1,5), "Charge Potential at (0.1cm, 2.5cm)", 10)


# Q3d
d = init_grid.greens_function(1,1)

init_grid.grid_plot(d, "Edge Potential Green's Function at (0.1cm, 0.1cm)", 10)
init_grid.grid_plot(init_grid.greens_charge(1,1), "Charge Potential at (0.1cm, 0.1cm)", 10)
#print()

print(f"At Central Point (5.0cm, 5.0cm): \n{a}")
print(f"At (2.5cm, 2.5cm): \n{b}")
print(f"At (0.1cm, 2.5cm): \n{c}")
print(f"At (0.1cm, 0.1cm): \n{d}")
