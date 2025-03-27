#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

Assignment 3: Code which utilises created classes to run Monte Carlo
simulations for two test cases, beginning with random points in a space
enclosed by a dimensioned square which contains a circular shape, and
then determines whether or not the points fall within the circular shape
or outside the shape, and return this and other determinations. The second
case 

"""
import time
from mpi4py import MPI
import numpy as np
import monte_carlo_class as mc
import shape_points_class as sp
import gaussian_class as ga


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
delta = 1/nproc

if rank==0:
    print(f"{nproc} Processor(s)")
    print()
    start_time = time.time()

# Define a number of samples to be more efficient rather than creating 100000000
# different sample arrays.
num_samples = np.int32(100000000/nproc)

# Use the Position class to initialise an object defined by a number of samples(rows)
# and a number of dimensions(columns) for the 4 necessary dimensional tests.
sample_2D = sp.Position(num_samples, 2)
sample_3D = sp.Position(num_samples, 3)
sample_4D = sp.Position(num_samples, 4)
sample_5D = sp.Position(num_samples, 5)

# Use the Monte Carlo class to process the data and return the required integral, variance,
# average & uncertainty.
Monte_2D = mc.MonteCarlo(sample_2D.r, sample_2D, -1, 1)
if rank==0:
    print("Test Case 1")
    print()
    print()
    print("2D Circle inside Square")
Monte_2D = Monte_2D.parallel_version()

Monte_3D = mc.MonteCarlo(sample_3D.r, sample_3D, -1, 1)
if rank==0:
    print("3D Sphere inside Cube")
Monte_3D = Monte_3D.parallel_version()

Monte_4D = mc.MonteCarlo(sample_4D.r, sample_4D, -1, 1)
if rank==0:
    print("4D Shape inside Tesseract")
Monte_4D = Monte_4D.parallel_version()

Monte_5D = mc.MonteCarlo(sample_5D.r, sample_5D, -1, 1)
if rank==0:
    print("5D Shape inside Hypersquare")
Monte_5D = Monte_5D.parallel_version()

# use the Gaussian class to initialise class objects, varying x0 and sigma to provide
# different scenarios to simulate the normal distribution.
dist_10 = ga.Gaussian(0, 1, num_samples, 1)
dist_11 = ga.Gaussian(1, 3, num_samples, 1)
dist_12 = ga.Gaussian(-1, 5, num_samples, 1)
dist_60 = ga.Gaussian(0, 1, num_samples, 6)
dist_61 = ga.Gaussian(1, 3, num_samples, 6)
dist_62 = ga.Gaussian(-1, 5, num_samples, 6)

# Again use the Monte Carlo class to process the data and return the required integral,
# variance, average & uncertainty.

Monte_1D0 = mc.MonteCarlo(dist_10.integral, dist_10, -1, 1)
if rank==0:
    print()
    print("Test Case 2")
    print()
    print()
    print(f"1D Normal Distribution, [sigma = {dist_10.sigma}], [x0 = {dist_10.x0}]")
Monte_1D0 = Monte_1D0.parallel_version()

Monte_1D1 = mc.MonteCarlo(dist_11.integral, dist_11, -1, 1)
if rank==0:
    print(f"1D Normal Distribution, [sigma = {dist_11.sigma}], [x0 = {dist_11.x0}]")
Monte_1D1 = Monte_1D1.parallel_version()

Monte_1D2 = mc.MonteCarlo(dist_12.integral, dist_12, -1, 1)
if rank==0:
    print(f"1D Normal Distribution, [sigma = {dist_12.sigma}], [x0 = {dist_12.x0}]")
Monte_1D2 = Monte_1D2.parallel_version()

Monte_6D0 = mc.MonteCarlo(dist_60.integral, dist_60, -1, 1)
if rank==0:
    print(f"6D Normal Distribution, [sigma = {dist_60.sigma}], [x0 = {dist_60.x0}]")
Monte_6D0 = Monte_6D0.parallel_version()

Monte_6D1 = mc.MonteCarlo(dist_61.integral, dist_61, -1, 1)
if rank==0:
    print(f"6D Normal Distribution, [sigma = {dist_61.sigma}], [x0 = {dist_61.x0}]")
Monte_6D1 = Monte_6D1.parallel_version()

Monte_6D2 = mc.MonteCarlo(dist_62.integral, dist_62, -1, 1)
if rank==0:
    print(f"6D Normal Distribution, [sigma = {dist_62.sigma}], [x0 = {dist_62.x0}]")
Monte_6D2 = Monte_6D2.parallel_version()

if rank==0:
    end_time = time.time()
    print(f"The runtime was {(end_time - start_time):.4f} seconds for {nproc} processors")
    print()
    print()

MPI.Finalize()
