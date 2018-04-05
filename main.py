# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 15:12:12 2018

@author: Trace
"""

import numpy as np
import random
from helper import QuickHull, plot_hull

# This block tests the given numbers in the assignment
# =============================================================================
# x = np.array([1, 4, 7, 10, 11, 11, 11, 12, 15, 16, 18, 18, 19, 22])
# y = np.array([6, 15, 7, 13, 6, 18, 21, 10, 18, 6, 3, 12, 15, 19])
# h = QuickHull(x, y)
# plot_hull(x, y, h)
# =============================================================================

# This block tests a random set of 25 points within the range of -100:100
# Create x and y values
x1 = np.array([random.randint(-100, 100) for i in range(25)])
y1 = np.array([random.randint(-100, 100) for i in range(25)])
# Calculate the hull coordinates
h1 = QuickHull(x1, y1)
# Plot the points and the conves hull
plot_hull(x1, y1, h1)