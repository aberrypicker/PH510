#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
import monte_carlo_class as mc
import shape_points_class as sp
import gaussian_class as ga


# Use the Position class to initialise an object defined by a number of samples(rows)
# and a number of dimensions(columns) for the 4 necessary dimensional tests.
sample_2D = sp.Position(1000000, 2)
sample_3D = sp.Position(1000000, 3)
sample_4D = sp.Position(1000000, 4)
sample_5D = sp.Position(1000000, 5)

Monte_2D = mc.MonteCarlo(sample_2D.r, sample_2D, -1, 1)
Monte_2D = Monte_2D.parallel_version()

Monte_3D = mc.MonteCarlo(sample_3D.r, sample_3D, -1, 1)
Monte_3D = Monte_3D.parallel_version()

Monte_4D = mc.MonteCarlo(sample_4D.r, sample_4D, -1, 1)
Monte_4D = Monte_4D.parallel_version()

Monte_5D = mc.MonteCarlo(sample_5D.r, sample_5D, -1, 1)
Monte_5D = Monte_5D.parallel_version()

# use the Gaussian class to initialise class objects
dist_1 = ga.Gaussian(0, 1, 1000000, 1)
dist_6 = ga.Gaussian(0, 1, 1000000, 6)

MonteGaus_1 = mc.MonteCarlo(dist_1.integral, dist_1, -1, 1)
MonteGaus_1 = MonteGaus_1.parallel_version()

MonteGaus_6 = mc.MonteCarlo(dist_6.integral, dist_6, -1, 1)
MonteGaus_6 = MonteGaus_6.parallel_version()
