# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 15:09:48 2018

@author: Trace
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def line_from_points(x1, y1, x2, y2):
    # Determine slope-intercept equation of line connecting two points
    # Input: (x1, y1), (x2, y2) - Pair of points 
    # Output: m - slope of line
    #         b - y-intercept
    m = (y2 - y1) / (x2 - x1)
    b = y1 - (m * x1)
    return m, b

def dist_line_point(a, b, c, x0, y0):
    # Calculate perpendicular distance from point to line
    # Input: a, b, c - Coefficients from equation of line as ax + by + c = 0
    #        (x0, y0) - x and y coordinates of point
    # Output: The distance from the point to the line
    return abs(a * x0 + b * y0 + c) / math.sqrt(a**2 + b**2)

def area_triangle(x1, y1, x2, y2, x3, y3):
    # Calculate area of triangle
    # Input: (x1, y1), (x2, y2), (x3, y3) - Coordinates of points making up
    #               vertices of triangle
    # Output: area of triangle
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)

def is_outside(x1, y1, x2, y2, x3, y3, x, y):
    # Determine whether a point lies within a triangle
    # Inputs: (x1, y1), (x2, y2), (x3, y3) - Coordinates of points making up
    #               vertices of triangle
    #         (x, y) - Coordinates of point to check
    # Output: true if outside triangle, otherwise false
    A = area_triangle(x1, y1, x2, y2, x3, y3)
    A1 = area_triangle(x, y, x2, y2, x3, y3)
    A2 = area_triangle(x1, y1, x, y, x3, y3)
    A3 = area_triangle(x1, y1, x2, y2, x, y)
    return A < (A1 + A2 + A3)

def is_above_line(y, m, x, b):
    # Determine whether or point is above or below a line
    # Input: y - y coordinates to be tested
    #        m - slope of line
    #        x - x coordinates used in equation of line
    #        b - y-intercept of line
    # Output: true if point is above line, otherwise false
    yy = m*np.array(x) + b
    return y - yy > 0

def QuickHullRecursive(x, y, x0, y0, x1, y1, m, b):
    # Recursive part of QuickHull function
    # Input: x - array of x coordinates of set of points
    #        y - array of y coordinates of set of points
    #        x0, y0, x1, y1 - x and y coordinates of line used to split points
    #        m, b - slope and y-intercept, respectively, of line used to split points
    # Output: Set of points making up the convex hull
    
    # Base Case
    # If x is empty, return an empty list
    if len(x) == 0:
        return []
    # If x contains only 1 element, it is part of the hull, return (x, y)
    elif len(x) == 1:
        return [[x[0], y[0]]]
    # Recursive case
    else:
        # Find point with furthest perpendicular distance from input line
        # This point is part of the convex hull
        max_d = 0
        index = 0
        # Loop through each point
        for i in range(len(x)):
            # Calculate distance from point to input line
            dist = dist_line_point(m, -1, b, x[i], y[i])
            # If the distance is greater than the max, update max 
            # and record the index of the max
            if dist > max_d:
                max_d = dist
                index = i
        
        # We only want to keep the points outside of the triangle 
        # made by x0, y0, x1, y1, and x[index], y[index]
        out_tri = []
        # Loop through each point and check to see if it is outside the triangle
        # Add it's boolean to the mask to be used to maintain only needed points
        for i in range(len(x)):
            out_tri.append(is_outside(x0, y0, x1, y1, x[index], y[index], x[i], y[i]))
        # Using out_tri as a boolean mask, we can keep only the points we need=
        new_x = np.array(x)[out_tri]
        new_y = np.array(y)[out_tri]

        # If there were no points remaning outside the triangle, we just need
        # to return x[index], y[index]
        if len(new_x) == 0:
            return [[x[index], y[index]]]
        # Similarly, if there is only one point outside the triangle, it is also
        # part of the hull, so return x[index], y[index], new_x[0], new_y[0]
        elif len(new_x) == 1:
            return [[x[index], y[index]]] +  [[new_x[0], new_y[0]]]
        # Otherwise, we must split the remaining points based on which side of triangle
        # they lie and then call the function recursively
        else:
            
            # This determines the equation of the line for the side of the triangle 
            # from (x0, y0) to (x[index], y[index]). We are considering this the left side
            ml, bl = line_from_points(x0, y0, x[index], y[index])
            # This determines the equation of the line for side of the triangle
            # from (x[index], y[index]) to (x1, y1). We are considering this the right side
            mr, br = line_from_points(x[index], y[index], x1, y1)

            # Now that we have the equations for the sides of the triangle, we split
            # the points by where they lie in relation to x[index] (left or right)
            left_x = new_x[new_x<x[index]]
            left_y = new_y[new_x<x[index]]
            right_x = new_x[new_x>x[index]]
            right_y = new_y[new_x>x[index]]
            
            # Finally, we return (x[index], y[index]) and make 2 recursive calls
            # The first one is passed the points outside the left side of the triangle,
            # the points (x0, y0) and (x[index], y[index]), along with the slope and y-intercept
            # of the line that forms the left side of the triangle
            # Similarly, the second one is passed the points outside the right side of the triangle,
            # the points (x1, y1) and (x[index], y[index]), along with the slope and y-intercept
            # of the line that forms the right side of the triangle

            return ([[x[index], y[index]]] + 
                        QuickHullRecursive(left_x, left_y, x0, y0, x[index], y[index], ml, bl) + 
                        QuickHullRecursive(right_x, right_y, x1, y1, x[index], y[index], mr, br))
            
def QuickHull(x, y):
    # Function to find the smallest convex polygon for a set of points that contains all the given points
    # Input: x - array of x coordinates of set of points
    #        y - array of y coordinates of set of points
    # Output: Array of points making up the convex hull
    
    # Start by finding the x-coordinates with smallest and largest values, as these are always included
    min_index = np.argmin(x)
    max_index = np.argmax(x)
    
    # Now find the equation of the line connecting the smallest and largest x-coordinates
    m, b = line_from_points(x[min_index], y[min_index], 
                            x[max_index], y[max_index])
    
    # Remove the min and max x-coordinates from the set of points
    if min_index < max_index:
        new_x = np.append(np.append(x[:min_index], x[min_index+1:max_index]), x[max_index+1:])
        new_y = np.append(np.append(y[:min_index], y[min_index+1:max_index]), y[max_index+1:])
    else:
        new_x = np.append(np.append(x[:max_index], x[max_index+1:min_index]), x[min_index+1:])
        new_y = np.append(np.append(y[:max_index], y[max_index+1:min_index]), y[min_index+1:])
        
    # Split the data by which points fall above and below the line
    above = is_above_line(new_y, m, new_x, b)
    top_x = new_x[above]
    top_y = new_y[above]
    bottom_x = new_x[np.logical_not(above)]
    bottom_y = new_y[np.logical_not(above)]
    
    # Create a 2-d array holding the x and y coordinates of the convex hull
    hull = np.array([[x[min_index], y[min_index]], [x[max_index], y[max_index]]])
    
    # Create a list of x and y coordinates of the convex hull to be populated
    # through recursive calls 
    pts = (QuickHullRecursive(top_x, top_y, 
                             x[min_index], y[min_index], 
                             x[max_index], y[max_index], 
                             m, b) + 
          QuickHullRecursive(bottom_x, bottom_y, 
                             x[min_index], y[min_index], 
                             x[max_index], y[max_index], 
                             m, b))

    # Return the array of hull values with the list of hull values appended to it
    return np.append(hull, pts, axis=0)

def order_points(points):
    # Function to find order points in a clockwise order
    # Input: points - array of (x, y) coordinates to be ordered
    # Output: points - original array sorted by angle between point and
    #       line made by average of all y-values
    
    # Compute average x and y to find center of points
    avg_x = np.mean(points[:, 0])
    avg_y = np.mean(points[:, 1])
    
    # Initialize theta to an empty list
    thetas = []
    
    # Loop through all points
    for i in range(len(points)):
        # Since angle is calculated using arctan, adjacent leg is x coordinate minus avg x
        # and opposite leg is y coordinate minus avg y
        a = points[i][0] - avg_x
        o = points[i][1] - avg_y
        theta = math.atan2(o, a) * 180 / math.pi
        # We want only positive angles, so we check
        if theta > 0:
            thetas.append(theta)
        # If theta is negative, add 360 to make positive
        else:
            thetas.append(360+theta)
            
    # Make an array from the list of thetas
    thetas = np.array([thetas])
    # Add the array to the points array
    points = np.append(points, thetas.T, 1)
    # Sort the points array by theta
    points = points[np.argsort(points[:, 2])]
    # Return the sorted points array without thetas
    return points[:, :-1]

def plot_hull(x, y, hull):
    # Function to plot the convex hull of a set of points
    # Input: x, y - arrays containing x and y coordinates of points
    #        hull - array containing x and y coordinates of convex hull

    # Call the order_points function to order points making up convex hull
    hull = order_points(hull)
    # Loop through the points and connect them with lines
    for i in range(len(hull)-1):
        plt.plot([hull[i][0], hull[i+1][0]], [hull[i][1], hull[i+1][1]], 'r')
    # Connect the last point to the first
    plt.plot([hull[-1][0], hull[0][0]], [hull[-1][1], hull[0][1]], 'r')
    # Create a scatter plot of all the points
    plt.scatter(x, y)