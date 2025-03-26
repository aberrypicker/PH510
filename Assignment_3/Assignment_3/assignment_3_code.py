#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
import time
from mpi4py import MPI
import monte_carlo_class as mc
import shape_points_class as sp
import gaussian_class as ga


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
delta = 1/nproc

if rank==0:
    start_time = time.time()

# Use the Position class to initialise an object defined by a number of samples(rows)
# and a number of dimensions(columns) for the 4 necessary dimensional tests.
sample_2D = sp.Position(100000000, 2)
sample_3D = sp.Position(100000000, 3)
sample_4D = sp.Position(100000000, 4)
sample_5D = sp.Position(100000000, 5)

Monte_2D = mc.MonteCarlo(sample_2D.r, sample_2D, -1, 1)
if rank==0:
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

# use the Gaussian class to initialise class objects
dist_1 = ga.Gaussian(0, 1, 100000000, 1)
dist_6 = ga.Gaussian(0, 1, 100000000, 6)

MonteGaus_1 = mc.MonteCarlo(dist_1.integral, dist_1, -1, 1)
if rank==0:
    print("1D Normal Distribution")
MonteGaus_1 = MonteGaus_1.parallel_version()

MonteGaus_6 = mc.MonteCarlo(dist_6.integral, dist_6, -1, 1)
if rank==0:
    print("6D Normal Distribution")
MonteGaus_6 = MonteGaus_6.parallel_version()

if rank==0:
    end_time = time.time()
    print(f"The runtime was {end_time - start_time} seconds for {nproc} processors")
