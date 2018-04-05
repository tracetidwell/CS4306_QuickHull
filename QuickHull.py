# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 15:09:48 2018

@author: Trace
"""
#comment
import numpy as np
import math
import matplotlib.pyplot as plt
import random

def line_from_points(x1, y1, x2, y2):
    m = (y2 - y1) / (x2 - x1)
    b = y1 - (m * x1)
    return m, b

def dist_line_point(a, b, c, x0, y0):
    ### Equation of line given by ax + by + c = 0
    ### Point defined as (x0, y0)
    return abs(a * x0 + b * y0 + c) / math.sqrt(a**2 + b**2)

def area_triangle(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)

def is_outside(x1, y1, x2, y2, x3, y3, x, y):
    A = area_triangle(x1, y1, x2, y2, x3, y3)
    A1 = area_triangle(x, y, x2, y2, x3, y3)
    A2 = area_triangle(x1, y1, x, y, x3, y3)
    A3 = area_triangle(x1, y1, x2, y2, x, y)
    return A != (A1 + A2 + A3)

def is_above_line(y, m, x, b):
    yy = m*np.array(x) + b
    return y - yy > 0

def intersection_of_lines(m0, b0, m1, b1):
    x = (b1 - b0) / (m0 - m1)
    y = m0 * x + b0
    return x, y

def dist_btw_pts(x1, y1, x2, y2):
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def QuickHullRecursive3(x, y, x0, y0, x1, y1, m, b):
    # Recursive part of QuickHull function
    # Input: x - array of x coordinates of set of points
    #        y - array of y coordinates of set of points
    #        x0, y0, x1, y1 - x and y coordinates of line used to split points
    #        m, b - slope and y-intercept, respectively, of line used to split points
    # Output: Set of points making up the convex hull
    

    # Find point with furthest perpendicular distance from line
    # This point is part of the convex hull
    if len(x) == 0:
        print('none')
        return []
    elif len(x) == 1:
        print('one')
        print([[x[0], y[0]]])
        return [[x[0], y[0]]]
    else:
        print('many')
        max_d = 0
        index = 0
        for i in range(len(x)):
            dist = dist_line_point(m, -1, b, x[i], y[i])
            if dist > max_d:
                max_d = dist
                index = i

        # All of the points inside the triangle made by x0, y0, x1, y1, and the point at the
        # furthest distance cannot be a part of the hull, so remove them
        out_tri = []
        for i in range(len(x)):
            out_tri.append(is_outside(x0, y0, x1, y1, x[index], y[index], x[i], y[i]))
        new_x = np.array(x)[out_tri]
        new_y = np.array(y)[out_tri]

        if len(new_x) == 0:
            print('many none')
            print(index)
            print([[x[index], y[index]]])
            return [[x[index], y[index]]]
        elif len(new_x) == 1:
            print('many one')
            print([[x[index], y[index]]] + [[new_x[0], new_y[0]]])
            return [[x[index], y[index]]] +  [[new_x[0], new_y[0]]]
        else:
            # We need the equation for the line going from the input line to the furthest point
            #m0 = -1/m
            #b0 = y[index] - (m0 * x[index])
            #x2, y2 = intersection_of_lines(m, b, m0, b0)
            
            ml, bl = line_from_points(x0, y0, x[index], y[index])
            
            mr, br = line_from_points(x[index], y[index], x1, y1)

            # Now that we have the new line, we split the data by which points fall above and below the line
            #above = is_above_line(new_y, m0, new_x, b0)
            #top_x = new_x[above]
            #top_y = new_y[above]
            #bottom_x = new_x[np.logical_not(above)]
            #bottom_y = new_y[np.logical_not(above)]
            left_x = new_x[new_x<x[index]]
            left_y = new_y[new_x<x[index]]
            right_x = new_x[new_x>x[index]]
            right_y = new_y[new_x>x[index]]

            return ([[x[index], y[index]]] + 
                        QuickHullRecursive3(left_x, left_y, x0, y0, x[index], y[index], ml, bl) + 
                        QuickHullRecursive3(right_x, right_y, x1, y1, x[index], y[index], mr, br))
            
def QuickHull(x, y):
    # Function to find the smallest convex polygon for a set of points that contains all the given points
    # Input: x - x coordinates of set of points
    #        y - y coordinates of set of points
    # Output: Set of points making up the convex hull
    
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
    
    #-y = -mx - b
    #mx - y + b = 0
    hull = np.array([[x[min_index], y[min_index]], [x[max_index], y[max_index]]])
    pts = (QuickHullRecursive3(top_x, top_y, 
                              x[min_index], y[min_index], 
                              x[max_index], y[max_index], 
                              m, b) + 
          QuickHullRecursive3(bottom_x, bottom_y, 
                              x[min_index], y[min_index], 
                              x[max_index], y[max_index], 
                              m, b))
#    hull = np.append(hull, 
#                     QuickHullRecursive3(top_x, top_y, 
#                                         x[min_index], y[min_index], 
#                                         x[max_index], y[max_index], 
#                                         m, b), axis=0)# + 
                     #QuickHullRecursive3(bottom_x, bottom_y, 
                      #                   x[min_index], y[min_index], 
                       #                  x[max_index], y[max_index], 
                        #                 m, b), axis=0)
#    hull = np.append(hull, QuickHullRecursive3(bottom_x, bottom_y, 
#                                              x[min_index], y[min_index], 
#                                              x[max_index], y[max_index], 
#                                              m, b), axis=0)
    return np.append(hull, pts, axis=0)