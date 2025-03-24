#!/usr/bin/env python3
"""
Version: Python 3.12.8

Assignment 3: Monte Carlo system creation.

"""
import assignment_3_classes as cl

sample_2D = cl.Position(1000000, 2)
sample_3D = cl.Position(1000000, 3)
sample_4D = cl.Position(1000000, 4)
sample_5D = cl.Position(1000000, 5)

Monte_2D = cl.MonteCarlo(sample_2D.r, sample_2D, -1, 1)
Monte_2D = Monte_2D.parallel_version()

Monte_3D = cl.MonteCarlo(sample_3D.r, sample_3D, -1, 1)
Monte_3D = Monte_3D.parallel_version()

Monte_4D = cl.MonteCarlo(sample_4D.r, sample_4D, -1, 1)
Monte_4D = Monte_4D.parallel_version()

Monte_5D = cl.MonteCarlo(sample_5D.r, sample_5D, -1, 1)
Monte_5D = Monte_5D.parallel_version()

dist_1 = cl.Gaussian(0, 1, 1000000, 1)
dist_6 = cl.Gaussian(0, 1, 1000000, 6)

MonteGaus_1 = cl.MonteCarlo(dist_1.integral, dist_1, -1, 1)
MonteGaus_1 = MonteGaus_1.parallel_version()

MonteGaus_6 = cl.MonteCarlo(dist_6.integral, dist_6, -1, 1)
MonteGaus_6 = MonteGaus_6.parallel_version()
