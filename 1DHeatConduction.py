"""Code for solving a case of one-dimensional heat conduction in a rod
with internal heat generation using 1D quadratic elements"""
from __future__ import division
import numpy as np
from pylab import*

"""Inputs"""
d = float(raw_input("Enter the diameter of the cross-sectional area of the rod(mm): "))
L = float(raw_input("Enter the length of the rod(mm): "))
k = float(raw_input("Enter the thermal conductivity of the material(W/mK): "))
n = int(raw_input("Enter the number of quadratic elements: "))
Q = float(raw_input("Enter the amount of heat generation(W/m3): "))
q = float(raw_input("Enter the heat flux(W/m2): "))
Temp = float(raw_input("Enter the temperature of the rod end(celsius): "))

"""Factors"""
fac = (k * ((pi/4)*((d/1000)**2)))/(3*((L/1000)/n))
gen = ((pi/4)*((d/1000)**2)*((L/1000)/n)*Q)

"""Matrices Declaration"""
G = np.zeros([(2*n+1), (2*n+1)])
mat = np.zeros([(2*n+1), (2*n+1)])
trim = np.zeros([(2*n), (2*n)])
fq = np.zeros([(2*n+1), (1)])
FQ = np.zeros([(2*n+1), (1)])
FG = np.zeros([(2*n), (1)])
trimFQ = np.zeros([(2*n), (1)])
TPlot = []
Local = []

"""Form the element conductance matrix"""
matrix = np.matrix([[(7*fac), (-8*fac), fac], [(-8*fac), (16*fac), (-8*fac)], [(fac), (-8*fac), (7*fac)]])

"""Form the heat generation vector"""
H = np.matrix([[(gen/6)], [(2*gen)/3], [gen/6]])

"""Direct Stiffness Method"""
"""Formation of global heat generation vector"""
l = 0
m = 0
p = 0
for k in range(0, 2*n-1):
    if(k%2 == 0):
        for i in range(k, k+3):
            fq[i,0] = H[p,0]
            p = p + 1
            for j in range(k, k+3):
                mat[i,j] = matrix[l,m]
                m = m + 1
            l = l + 1
            m = 0
        G = G + mat
        mat = np.zeros([(2*n+1), (2*n+1)])
        FQ = FQ + fq
        fq = np.zeros([(2*n+1), (1)])
        l = 0
        p = 0

"""Solve using elimination approach"""
for i in range(0, 2*n):
    trimFQ[i,0] = FQ[i,0]
    if(i == 0):
        FG[i,0] = ((pi/4)*((d/1000)**2))*q
    elif(i == (2*n-1) or i == (2*n-2)):
        FG[i,0] = -1 * G[i,2*n] * Temp
    else:
        FG[i,0] = 0
    for j in range(0, 2*n):
        trim[i,j] = G[i,j]
trim = np.matrix(trim)
T = (trim.I) * (trimFQ + FG)
finalT = np.vstack((T,Temp))
print "Nodal temperatures are: \n", finalT

"""Plot using the interpolation functions"""
val = 0
dist = 0
carry = 0
for k in range(0, 2*n-1):
    if (k%2 == 0):
        for p in range(0, 11, 1):
            s = p/10
            val = (2*finalT[k,0]*(s-0.5)*(s-1)) + (-4*finalT[(k+1),0]*s*(s-1)) + (2*finalT[(k+2),0]*s*(s-0.5))
            TPlot.append(val)
            dist = (carry*(L/n)) + (s*(L/n)) 
            Local.append(dist)
        carry = carry + 1 
plot(Local, TPlot)
xlabel("Rod length(mm)")
ylabel("Temperature(celsius)")
grid()
show()
