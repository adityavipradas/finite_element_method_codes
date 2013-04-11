"""
Created on Fri Jan  4 09:07:05 2013

@author: adityavipradas
"""
"""2D FEM Solver for a rectangular fin subjected to a certain inlet temperature"""
"""4-node quadrilateral elements are used"""
"""Leftmost face of the fin is isothermal""" 
import sympy as sy
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

#inputs
l = float(raw_input("Enter the length of the fin(mm): "))
w = float(raw_input("Enter the breadth of the fin(mm): "))
t = float(raw_input("Enter the thickness of the fin(mm): "))
T = float(raw_input("Enter the inlet temperature of the fin(deg): "))
Ta = float(raw_input("Enter the ambient temperature(deg): "))
kx = float(raw_input("Enter the thermal conductivity in x-direction(W/mK): "))
ky = float(raw_input("Enter the thermal conductivity in y-direction(W/mK): "))
h = float(raw_input("Enter the heat transfer coefficient(W/m2K): "))
Q = float(raw_input("Enter the internal heat generation(W/m^3): "))

#determine number of elements
if (l!=1 or w!=1):
    nrow = int(l/2.)
    ncol = int(w/2.)
else:
    nrow = 1
    ncol = 1
print"\n"
print nrow," x ",ncol," quadrilateral elements"
print"\n"


#declare local and global matrices
k = np.zeros((4,4))
K = np.zeros((((nrow + 1)*(ncol + 1)),((nrow + 1)*(ncol + 1))))
KQ = np.zeros(((nrow + 1)*(ncol + 1),1))
KH = np.zeros(((nrow + 1)*(ncol + 1),1))
khs = np.zeros(((nrow + 1)*(ncol + 1),1))
Temperature = np.zeros(((nrow + 1)*(ncol + 1),1))

#declare the interpolation functions in terms of normalized coordinates r,s
r = sy.Symbol('r')
s = sy.Symbol('s')
xr = sy.Matrix([[1-r, 1+r]])
ys = sy.Matrix([[1-s], [1+s]])

#inter array has all the interpolation functions
inter = 0.25*ys*xr
swap1 = inter[1,0]
swap2 = inter[1,1]
inter[1,0], inter[1,1] = swap2, swap1
inter.rows = 4
inter.cols = 1
print "Interpolation functions"
print"\n"
print inter
print"\n"

#differentiated interpolation functions
diffr = inter.diff(r)
diffs = inter.diff(s)
print"dN/dr"
print"\n"
print diffr
print"\n"
print"dN/ds"
print"\n"
print diffs
print"\n"

#declare a and b concerned with r and s respectively
a = l/(nrow * 2.)
b = w/(ncol * 2.)

#element conductance matrix
k = (t*((kx * diffr * diffr.T * (b/a)) + ((ky * diffs * diffs.T * (a/b))))).integrate((r, -1, 1),(s, -1, 1)) + (2 * h * a * b * inter * inter.T).integrate((r, -1, 1),(s, -1, 1))
print"Element Conductance Matrix"
print"\n"
print k
print"\n"

#edge convection matrices and edge convection force vectors
print"Edge convection conductance matrices and edge convection force vectors"""
print"\n"
#1.lower edge
s = -1
yslow = sy.Matrix([[1-s], [1+s]])
interlow = 0.25*yslow*xr
swap1 = interlow[1,0]
swap2 = interlow[1,1]
interlow[1,0], interlow[1,1] = swap2, swap1
interlow.rows = 4
interlow.cols = 1
klow = (h * t * a * interlow * interlow.T).integrate((r, -1, 1))
flow = (h * Ta * t * a * interlow).integrate((r, -1, 1))
print klow
print"\n"
print flow
print"\n"
print interlow
print"\n"
s = sy.Symbol('s')

#2. upper edge
s = 1
ysup = sy.Matrix([[1-s], [1+s]])
interup = 0.25*ysup*xr
swap1 = interup[1,0]
swap2 = interup[1,1]
interup[1,0], interup[1,1] = swap2, swap1
interup.rows = 4
interup.cols = 1
kup = (h * t * a * interup * interup.T).integrate((r, -1, 1))
fup = (h * Ta * t * a * interup).integrate((r, -1, 1))
print kup
print"\n"
print fup
print"\n"
s = sy.Symbol('s')

#3. right edge
r = 1
xrgh = sy.Matrix([[1-r, 1+r]])
intergh = 0.25*ys*xrgh
swap1 = intergh[1,0]
swap2 = intergh[1,1]
intergh[1,0], intergh[1,1] = swap2, swap1
intergh.rows = 4
intergh.cols = 1
krgh = (h * t * b * intergh * intergh.T).integrate((s, -1, 1))
frgh = (h * Ta * t * b * intergh).integrate((s, -1, 1))
print krgh
print"\n"
print frgh
print"\n"
r = sy.Symbol('r')

#internal heat generation vector
kq = (Q * inter * t).integrate((r, -1, 1),(s, -1, 1))
print"Internal heat generation vector"
print"\n"
print kq
print"\n"

#lateral convection force vector
kh = (2 * h * Ta * inter).integrate((r, -1, 1),(s, -1, 1))
print"Lateral convection force vector"
print"\n"
print kh
print"\n"

#global conductance and internal heat generation matrices
l = 0
m = 0 

for i in range(0,((nrow + 1)*(ncol + 1))):
#((ncol+1) * (nrow+1))-(ncol+1)) = (ncol+1)*nrow
    if(((i+1)%(ncol+1)!=0 or i==0) and i < (ncol+1)*nrow):
        for j in (i, i+(ncol+1), i+(ncol+2), i+1):
            KQ[j,0] = KQ[j,0] + kq[l,0]
            KH[j,0] = KH[j,0] + kh[l,0]
            if (i % (ncol+1) == 0):
                khs[j,0] = khs[j,0] + flow[l,0]
            elif ((i + 2) % (ncol + 1) == 0):
                khs[j,0] = khs[j,0] + fup[l,0]
            if (i >= (ncol+1)*(nrow-1)):
                khs[j,0] = khs[j,0] + frgh[l,0]
            for p in (i, i+(ncol+1), i+(ncol+2), i+1):
                K[j,p] = K[j,p] + k[l,m]
                if (i % (ncol+1) == 0):
                    K[j,p] = K[j,p] + klow[l,m]
                elif ((i + 2) % (ncol + 1) == 0):
                    K[j,p] = K[j,p] + kup[l,m]
                if (i >= (ncol+1)*(nrow-1)):
                    K[j,p] = K[j,p] + krgh[l,m]
                m = m + 1
            l = l + 1
            m = 0
    l = 0
    m = 0
print"Global conductance matrix"
print"\n"
print K
print"\n"
print K.shape
print"\n"

#elimination approach
Ktrim = K[ncol+1:,0:ncol+1]
Ttrim = np.ones((ncol+1,1))
force = (KQ[ncol+1:,] + KH[ncol+1:,] + khs[ncol+1:,])-np.dot(Ktrim,(T*Ttrim))

#nodal temperatures
Temp = np.dot(np.linalg.inv(K[ncol+1:,ncol+1:]),force)
for i in range(0,((nrow + 1)*(ncol + 1))):
    if (i < (ncol+1)):
        Temperature[i,0] = T
    else:
        Temperature[i,0] = Temp[i-(ncol+1)]
print"Nodal temperatures"
print"\n"
print Temperature
print"\n"

"""Plot"""
InterTemp = []
rSide = []
sSide = []
prevR = a
prevS = b

for i in range(0,((nrow + 1)*(ncol + 1))):
    if(((i+1)%(ncol+1))==0):
        prevR = prevR + (2*a)
        prevS = b
    if(((i+1)%(ncol+1)!=0 or i==0) and i < (ncol+1)*nrow):
        for s in range(-10, 11, 1):
            m = s/10.
            for r in range(-10, 11, 1):
                n = r/10.
                InterTemp.append((((-1*n + 1)*(-0.25*m + 0.25))*(Temperature[i,0]))+(((n + 1)*(-0.25*m + 0.25))*(Temperature[i+(ncol+1),0]))+(((n + 1)*(0.25*m + 0.25))*(Temperature[i+(ncol+2),0]))+(((-1*n + 1)*(0.25*m + 0.25))*(Temperature[i+1,0])))
                rSide.append(prevR + (n*a))
                sSide.append(prevS + (m*b))
        prevS = prevS + (2*b)
rSide = np.array(rSide)
sSide = np.array(sSide)
InterTemp = np.array(InterTemp)
xi, yi = np.linspace(rSide.min(), rSide.max(), 300), np.linspace(sSide.min(), sSide.max(), 300)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((rSide, sSide), InterTemp, (xi, yi), method = 'linear')
plt.imshow(zi, vmin=InterTemp.min(), vmax=InterTemp.max(), origin = 'lower', extent=[rSide.min(), rSide.max(), sSide.min(), sSide.max()])
plt.title("Temperature distribution")
plt.xlabel("length")
plt.ylabel("width")
plt.colorbar()
plt.show()

