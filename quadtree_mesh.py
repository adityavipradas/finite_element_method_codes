# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:21:14 2013

@author: dell
"""

#spatial decomposition quadtree mesh generation method

#geometry creation

from __future__ import division
import matplotlib.pyplot as plt


class quad_elements(object):
    def __init__(self, node1, node2, node3, node4, number):
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.number = number
    def getNode1(self):
        return self.node1
    def getNode2(self):
        return self.node2
    def getNode3(self):
        return self.node3
    def getNode4(self):
        return self.node4
    def getnumber(self):
        return self.number
        
        
def generate_points(x, y, step):
    x.append(x[0])
    y.append(y[0])
    x_coord = []
    y_coord = []
    for i in range(len(x)-1):        
        x_coord.append(x[i])
        y_coord.append(y[i])
        x_new = x[i]
        y_new = y[i]
        y_stop = int(abs(y[i+1] - y[i])/step)
        x_stop = int(abs(x[i+1] - x[i])/step)
        if x[i+1] == x[i]:
            for j in range(1, y_stop):
                x_coord.append(x[i])
                if (y[i+1] > y[i]):
                    y_new = y_new + step
                    y_coord.append(y_new)
                else:
                    y_new = y_new - step
                    y_coord.append(y_new)
        else:
            for j in range(1, x_stop):
                if x[i] < x[i+1]:
                    x_new =  x_new + step
                else:
                    x_new = x_new - step
                y_new = (x_new - x[i]) * (y[i+1] - y[i])/(x[i+1] - x[i]) + y[i]
                x_coord.append(x_new)
                y_coord.append(y_new)
    return x_coord, y_coord, x, y

def generate_first_element(xl, yl):
    elements = []
    elements.append(quad_elements([xl[0], yl[0]], [xl[0], yl[1]], \
    [xl[1], yl[1]], [xl[1], yl[0]], 1))
    
    ###plot start##############################################################
    plt.scatter([xl[0], xl[0], xl[1], xl[1]], [yl[0], yl[1],\
    yl[1], yl[0]], color = 'red')

    plt.plot([xl[0], xl[1]], [yl[0], yl[0]], \
    [xl[1], xl[1]], [yl[0], yl[1]], \
    [xl[1], xl[0]], [yl[1], yl[1]], \
    [xl[0], xl[0]], [yl[1], yl[0]])
    print (xl, yl)
    ###plot end################################################################
    
    return elements
    
def quad_decompose(xl, yl, num, indx, elements):
    n = 4 * num - 2
    discretize = []
    add = 0
    for i in range(0, 2):
        for j in range(0, 2):
            discretize.append(quad_elements([xl[i], yl[j]], [xl[i], yl[j+1]], \
            [xl[i+1], yl[j+1]], [xl[i+1], yl[j]], n + add))
            add = add + 1
    elements[indx] = discretize
    
    
    
if __name__ == "__main__":
    x_coord, y_coord, x, y = generate_points([1, 1, 3, 4, 3, 0], [1, 2, 2, 3, \
    -6, 1], 0.1)

    avg_y = (max(y) + min(y))/2
    avg_x = (max(x) + min(x))/2
    mod_y = (max(y) - min(y) + 2)/2
    mod_x = (max(x) - min(x) + 2)/2

    if (max(x) - min(x)) > (max(y) - min(y)):
        p1, p2, p3, p4 = min(x) - 1, max(x) + 1, avg_y - mod_x, avg_y + mod_x
        elements = generate_first_element([p1, p2], [p3, p4])
    else:
        p1, p2, p3, p4 = avg_x - mod_y, avg_x + mod_y, min(y) - 1, max(y) + 1
        elements = generate_first_element([p1, p2], [p3, p4])
        
    for i in range(0, len(elements)):
        count = 0
        for j in range(0, len(x_coord)):
            if x_coord[j] > p1 and x_coord[j] < p2 and y_coord[j] > p3 and \
            y_coord[j] < p4:
                count = count + 1
            
            if count > 1:
                num = elements[i].getnumber()
                indx = elements.index(elements[i])
                px = (p2 - p1)/2
                py = (p4 - p3)/2
                quad_decompose([p1, px, p2], [p3, py, p4], num, indx, elements)
        
    ###plot start##############################################################    
    print(x_coord, y_coord)
    plt.plot(x, y)
    plt.scatter(x_coord, y_coord)
    plt.show()            
    ###plot end################################################################
