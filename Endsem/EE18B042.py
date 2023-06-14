# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 11:23:03 2021

@author: Venkat Sai
"""

import os,sys
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np  


M = 30
N = 25                          #dimensions of boundary
i = 12
j = 10                          #capacitor location
Nit = 5000                      #given no. of iterations



x = linspace(int(0), int(N), N+1)
y = linspace(int(0), int(M), M+1)
X, Y = meshgrid(x, y)                   #creating a meshgrid

def getphi(i, Niter):
    phi = zeros([M+1, N+1])             #creating simulation array
    tt = where((Y == i) & (X >= j))     #points of top capacitor
    phi[tt] = 1                         #voltage = 1V
    mm = where(Y == M)                  #points of lower capacitor
    phi[mm] = 0                         #voltage = 0V
    errors = zeros(Niter)               #initialising errors
    x_scale = []
    for k in range(Niter):
        oldphi = phi.copy() 
        phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
        phi[1:-1,0] = phi[1:-1,1]       # boundary conditions(first column)
        phi[1:-1,-1] = phi[1:-1,-2]     #last column
        phi[0,1:-1] = phi[1,1:-1]       #top row
        phi[mm] = 0                     #reassigning potentials
        phi[tt] = 1
        errors[k] = (abs(phi-oldphi)).max()     #error calculation
    return phi, errors
   
mas, error = getphi(i, Nit)             #obtaining phi and error
x_scale = linspace(int(1), int(Nit), Nit)               
semilogy(x_scale, error, '-r', label='errors semilogy')    #semilog plot for errors
legend(loc='upper right')
xlabel(r'$N$',size=10)
ylabel(r'$errors$',size=10)
show()

s = []
errors1 = []
for k in range(Nit):
    s.append([1, x_scale[k]])
for k in range(Nit):
    errors1.append(log(error[k]))
errors1 = array(errors1)
s = array(s)
m = array(lstsq(s, errors1, rcond=None)[0])         # finding best fit1
A = exp(m[0])
B = m[1]
print(A,B)

Nnew = floor(np.log(-0.01*B/A)/B - 0.5)             #iteration no. for error<0.01
print(f"The new iteration number is {Nnew}")

mas, error = getphi(i, int(Nnew))       #new phi and error

x1 = linspace(int(0), int(N), N+1)      #creating a mesh grid of inverse dimensions
y1 = linspace(int(M), int(0), M+1)
X1,Y1 = meshgrid(x1,y1)

c = contourf(X1, Y1, mas)     #contour plot
colorbar(c)
xlabel(r'$x$',size=10)
ylabel(r'$30-y$',size=10)
show()

Jx = zeros([M+1, N+1])
Jy = zeros([M+1, N+1])
Jx[1:-1,1:-1] = 0.5*(mas[1:-1, 0:-2] - mas[1:-1, 2:])   #finding x-field vector
Jy[1:-1,1:-1] = 0.5*(mas[0:-2, 1:-1] - mas[2:, 1:-1])   #finding y-field vector
quiver(X1, Y1, Jx, -Jy)     #plotting field vectors
xlabel(r'$x$',size=10)
ylabel(r'$30-y$',size=10)
show()

e = 8.85e-12                #value of epsilon
Dy = e*Jy[i+1,j:]           #calculating Dy
Q = 0                       #initiating charge
for k in range(len(Dy)):    #summing Dy
    Q += Dy[k]/len(Dy)      #averaging
print(f"The value of Dy is {Q}")
Q = Q*(N-j+1)               #multipying by area
print(f"The value of Q is {Q}")
print(f"The vector Dy is {Dy}")

def Qt(l):                  #function to calculate Qtheory
    Qt = e*(N-j)/(l)
    return Qt


def cap(l):                 #function to calculate Q
    las, sec = getphi(M-l, int(Nnew))       #new phi and error
    Jy[1:-1,1:-1] = 0.5*(las[0:-2, 1:-1] - las[2:, 1:-1])   #correspondin Jy or Ey
    e = 8.85e-12
    Dy = e*Jy[M-l+1,j:]             #corresponding Dy
    Q = 0
    for k in range(len(Dy)):
        Q += Dy[k]
    Q = Q/len(Dy)             #averaging
    Q = Q*(N-j+1)             #multipying by area
    return Q

def dev(l):                 #function to calculate deviation
    d = cap(l) - Qt(l)
    return d

print(f"The value of Q, Qt, deviation for given M-i value are in order")
print(cap(5), Qt(5), dev(5))        #M-i=5
print(cap(10), Qt(10), dev(10))     #M-i=10
print(cap(20), Qt(20), dev(20))     #M-i=20
print(cap(30), Qt(30), dev(30))     #M-i=30

