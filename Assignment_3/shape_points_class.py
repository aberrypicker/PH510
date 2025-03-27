#!/usr/bin/env python3
"""
Version: Python 3.12.8

Copyright 2025: Aaron Berryman. Licensed under MIT license.

This file contains the code and documentation for the class which defines points
inside shapes of varying dimensions.

"""

import numpy as np

class Position:
    """
    A class which uses the dimensions and number of points to determine where the point lies 
    and relate this to the position of a unit circle inside a square, and so forth.
    """
    def __init__(self, n, d):
        """
        Initialise the position using the number of samples and number of dimensions.
        """
        self.n = n
        self.d = d
        self.samples = np.random.uniform(-1, 1, (self.n, self.d))

    def __str__(self):
        """
        Gives the printed information for the samples.
        """
        return f"Position:({self.samples})"

    def r(self):
        """
        Depending on the number of dimensions, turns the samples into vector positions
        in the form of r, for use in determining if inside or outside of unit shape.
        """
        r = (np.linalg.norm(self.samples, axis=1)).reshape(-1,1)
        check = (r <=  1).astype(int)
        return check

