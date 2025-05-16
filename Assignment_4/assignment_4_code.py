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
#import matplotlib.pyplot as plt
import monte_carlo_class as mc
import poisson_green_class as pg


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
delta = 1/nproc

# Measures runtime of code
if rank==0:
    print(f"{nproc} Processor(s)")
    print()
    start_time = time.time()

# Define a number of samples to be more efficient rather than creating many
# different sample arrays.
n_samples = np.int32(100000/nproc)


# Question 3

# Performs random walk and evaluates the Green's function via probability for each specified
# point from the assignment.

if rank==0:
    print("Exercise 3")
    print("Green's function evaluation for a square grid of side length 10cm:")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

#phi, f = init_grid.phi, init_grid.f


# (a)
a_monte = mc.MonteCarlo(init_grid, init_grid.greens_function, -1, 1, 10, 10)
a_calc = a_monte.parallel_array()

#print(f"At centre point (5cm, 5cm):\n{a[0]}")
#init_grid.grid_plot(a[0], "Green's Function: (5cm, 5cm)")
#init_grid.grid_plot(a[1], 'Number of Site Visits: (5cm, 5cm)')

# Q3b
b_monte = mc.MonteCarlo(init_grid, init_grid.greens_function, -1, 1, 5, 5)
b_calc = b_monte.parallel_array()

#print(f"At (2.5cm, 2.5cm):\n{b[0]}")
#init_grid.grid_plot(b[0], "Green's Function: (2.5cm, 2.5cm)")
#init_grid.grid_plot(b[1], 'Number of Site Visits: (2.5cm, 2.5cm)')


# Q3c
c_monte = mc.MonteCarlo(init_grid, init_grid.greens_function, -1, 1, 1, 5)
c_calc = c_monte.parallel_array()

#print(f"At (0.1cm, 2.5cm):\n{c[0]}")
#init_grid.grid_plot(c[0], "Green's Function: (0.1cm, 2.5cm)")
#init_grid.grid_plot(c[1], 'Number of Site Visits: (0.1cm, 2.5cm)')


# Q3d
d_monte = mc.MonteCarlo(init_grid, init_grid.greens_function, -1, 1, 1, 1)
d_calc = d_monte.parallel_array()

#print(f"At (0.1cm, 0.1cm):\n{d[0]}")
#init_grid.grid_plot(d[0], "Green's Function: (0.1cm, 0.1cm)")
#init_grid.grid_plot(d[1], 'Number of Site Visits: (0.1cm, 0.1cm)')
#print()

if rank==0:
    print(f"At centre point (5cm, 5cm):\n{a_calc[1]}")
    print(f"At (2.5cm, 2.5cm):\n{b_calc[1]}")
    print(f"At (0.1cm, 2.5cm):\n{c_calc[1]}")
    print(f"At (0.1cm, 0.1cm):\n{d_calc[1]}")
    print()


# Question 4

# Now evaluates the same cases as Question 3 but now with additional Potentials added
# via boundary conditions, with the latter part of the question considering further
# permutations with charge placed within the grid in various configurations.

    print("Exercise 4")
    print("Potential calculation via Green's function for a square grid of side length 10cm:")

# Q4a
    print("Q4a with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_1 = init_grid.phi
phi_1 = init_grid.boundary_condition('Q4-a')
phi_1 = init_grid.overrelaxation_method()

potential_50_50_4a = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4a_calc = potential_50_50_4a.parallel_version()

potential_25_25_4a = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4a_calc = potential_25_25_4a.parallel_version()

potential_1_25_4a = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4a_calc = potential_1_25_4a.parallel_version()

potential_1_1_4a = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4a_calc = potential_1_1_4a.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4a_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4a_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4a_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4a_calc[0]:.4f}V")
    print()


# Q4b
    print("Q4b with boundary conditions: Top and bottom edges at +1V, left & right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)


phi_2 = init_grid.phi
phi_2 = init_grid.boundary_condition('Q4-b')
phi_2 = init_grid.overrelaxation_method()

potential_50_50_4b = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4b_calc = potential_50_50_4b.parallel_version()

potential_25_25_4b = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4b_calc = potential_25_25_4b.parallel_version()

potential_1_25_4b = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4b_calc = potential_1_25_4b.parallel_version()

potential_1_1_4b = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4b_calc = potential_1_1_4b.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4b_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4b_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4b_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4b_calc[0]:.4f}V")
    print()


# Q4c
    print("Q4c with boundary conditions: Top/left edges +2V, bottom edge 0V, right edge -4V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)


phi_3 = init_grid.phi
phi_3 = init_grid.boundary_condition('Q4-c')
phi_3 = init_grid.overrelaxation_method()

potential_50_50_4c = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4c_calc = potential_50_50_4b.parallel_version()

potential_25_25_4c = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4c_calc = potential_25_25_4b.parallel_version()

potential_1_25_4c = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4c_calc = potential_1_25_4b.parallel_version()

potential_1_1_4c = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4c_calc = potential_1_1_4b.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4c_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4c_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4c_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4c_calc[0]:.4f}V")
    print()


# Question 4 - Part 2: These are the same calclualtions but now with charge added
# into consideration. The cases are labelled the same as above, but with "d, e, f"
# used to distinguish from "a, b, c" which are used in the assignment instructions,
# and "i, ii, iii" to distinguish for the 3 boundary conditions

# Q4d
    print("Q4d) Repeat, with uniform charge of 10C throughout the grid")
    print()

# (i)
    print("i with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_4 = init_grid.phi
phi_4 = init_grid.boundary_condition('Q4-a')
phi_4 = init_grid.overrelaxation_method()

f_4 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50_4di = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4di_calc = potential_50_50_4di.parallel_version()

potential_25_25_4di = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4di_calc = potential_25_25_4di.parallel_version()

potential_1_25_4di = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4di_calc = potential_1_25_4di.parallel_version()

potential_1_1_4di = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4di_calc = potential_1_1_4di.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4di_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4di_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4di_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4di_calc[0]:.4f}V")
    print()

# (ii)
    print("ii with boundary conditions: Top and bottom edges at +1V, left and right edges at -1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)


phi_5 = init_grid.phi
phi_5 = init_grid.boundary_condition('Q4-b')
phi_5 = init_grid.overrelaxation_method()

f_5 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50_4dii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4dii_calc = potential_50_50_4dii.parallel_version()

potential_25_25_4dii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4dii_calc = potential_25_25_4dii.parallel_version()

potential_1_25_4dii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4dii_calc = potential_1_25_4dii.parallel_version()

potential_1_1_4dii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4dii_calc = potential_1_1_4dii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4dii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4dii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4dii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4dii_calc[0]:.4f}V")
    print()


# (iii)
    print("iii with boundary conditions: Top/left edges +2V, bottom edge 0V, right edge -4V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_6 = init_grid.phi
phi_6 = init_grid.boundary_condition('Q4-c')
phi_6 = init_grid.overrelaxation_method()

f_6 = init_grid.charge_distribution_scenario('uniform_10C')

potential_50_50_4diii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4diii_calc = potential_50_50_4diii.parallel_version()

potential_25_25_4diii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4diii_calc = potential_25_25_4diii.parallel_version()

potential_1_25_4diii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4diii_calc = potential_1_25_4diii.parallel_version()

potential_1_1_4diii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4diii_calc = potential_1_1_4diii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4diii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4diii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4diii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4diii_calc[0]:.4f}V")
    print()


# Q4e
    print("Q4e) Repeat, with uniform charge gradient from 1C at top to 0C at bottom")
    print()

# (i)
    print("i) with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)


phi_7 = init_grid.phi
phi_7 = init_grid.boundary_condition('Q4-a')
phi_7 = init_grid.overrelaxation_method()

f_7 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50_4ei = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4ei_calc = potential_50_50_4ei.parallel_version()

potential_25_25_4ei = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4ei_calc = potential_25_25_4ei.parallel_version()

potential_1_25_4ei = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4ei_calc = potential_1_25_4ei.parallel_version()

potential_1_1_4ei = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4ei_calc = potential_1_1_4ei.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4ei_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4ei_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4ei_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4ei_calc[0]:.4f}V")
    print()


# (ii)
    print("ii) with boundary conditions: Top/bottom edges +1V, left and right edges -1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)


phi_8 = init_grid.phi
phi_8 = init_grid.boundary_condition('Q4-b')
phi_8 = init_grid.overrelaxation_method()

f_8 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50_4eii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4eii_calc = potential_50_50_4eii.parallel_version()

potential_25_25_4eii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4eii_calc = potential_25_25_4eii.parallel_version()

potential_1_25_4eii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4eii_calc = potential_1_25_4eii.parallel_version()

potential_1_1_4eii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4eii_calc = potential_1_1_4eii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4eii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4eii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4eii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4eii_calc[0]:.4f}V")
    print()



# (iii)
    print("iii with boundary conditions: Top/left edges +2V, bottom edge 0V, right edge -4V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_9 = init_grid.phi
phi_9 = init_grid.boundary_condition('Q4-c')
phi_9 = init_grid.overrelaxation_method()

f_9 = init_grid.charge_distribution_scenario('linear_gradient_top_to_bottom')

potential_50_50_4eiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4eiii_calc = potential_50_50_4eiii.parallel_version()

potential_25_25_4eiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4eiii_calc = potential_25_25_4eiii.parallel_version()

potential_1_25_4eiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4eiii_calc = potential_1_25_4eiii.parallel_version()

potential_1_1_4eiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4eiii_calc = potential_1_1_4eiii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4eiii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4eiii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4eiii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4eiii_calc[0]:.4f}V")
    print()


# Q4f
    print("f) Repeat, with exponentially decaying charge exp(-2000|r|) placed at centre of grid")
    print()

# (i)
    print("i with boundary conditions: All edges uniformly at +1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_10 = init_grid.phi
phi_10 = init_grid.boundary_condition('Q4-a')
phi_10 = init_grid.overrelaxation_method()

f_10 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50_4fi = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4fi_calc = potential_50_50_4fi.parallel_version()

potential_25_25_4fi = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4fi_calc = potential_25_25_4fi.parallel_version()

potential_1_25_4fi = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4fi_calc = potential_1_25_4fi.parallel_version()

potential_1_1_4fi = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4fi_calc = potential_1_1_4fi.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4fi_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4fi_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4fi_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4fi_calc[0]:.4f}V")
    print()


# (ii)
    print("ii with boundary conditions: Top/bottom edges +1V, left and right edges -1V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_11 = init_grid.phi
phi_11 = init_grid.boundary_condition('Q4-b')
phi_11 = init_grid.overrelaxation_method()

f_11 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50_4fii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4fii_calc = potential_50_50_4fii.parallel_version()

potential_25_25_4fii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4fii_calc = potential_25_25_4fii.parallel_version()

potential_1_25_4fii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4fii_calc = potential_1_25_4fii.parallel_version()

potential_1_1_4fii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4fii_calc = potential_1_1_4fii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4fii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4fii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4fii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4fii_calc[0]:.4f}V")
    print()


# (iii)
    print("iii with boundary conditions: Top/left edges +2V, bottom edge 0V, right edge -4V")
init_grid = pg.PoissonGrid(0.10, 21, n_samples)

phi_12 = init_grid.phi
phi_12 = init_grid.boundary_condition('Q4-c')
phi_12 = init_grid.overrelaxation_method()

f_12 = init_grid.charge_distribution_scenario('exp_decay')

potential_50_50_4fiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 10, 10)
potential_50_50_4fiii_calc = potential_50_50_4fiii.parallel_version()

potential_25_25_4fiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 5, 5)
potential_25_25_4fiii_calc = potential_25_25_4fiii.parallel_version()

potential_1_25_4fiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 5)
potential_1_25_4fiii_calc = potential_1_25_4fiii.parallel_version()

potential_1_1_4fiii = mc.MonteCarlo(init_grid, init_grid.greens_potential, 0, 10, 1, 1)
potential_1_1_4fiii_calc = potential_1_1_4fiii.parallel_version()

if rank==0:
    print(f"At (5cm, 5cm): {potential_50_50_4fiii_calc[0]:.4f}V")
    print(f"At (2.5cm, 2.5cm): {potential_25_25_4fiii_calc[0]:.4f}V")
    print(f"At (0.1cm, 2.5cm): {potential_1_25_4fiii_calc[0]:.4f}V")
    print(f"At (0.1cm, 0.1cm): {potential_1_1_4fiii_calc[0]:.4f}V")
    print()

# Record ending time and determine total runtime based on no. of processors used.
if rank==0:
    end_time = time.time()
    print(f"The runtime was {(end_time - start_time):.4f} seconds for {nproc} processors")
    print()
    print()

MPI.Finalize()
