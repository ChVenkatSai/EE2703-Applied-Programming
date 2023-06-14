# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 15:10:36 2020

@author: venkat sai
"""
import os,sys
from pylab import * 
import mpl_toolkits.mplot3d.axes3d as p3

if (len(sys.argv) == 5):        #importing parameters
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    radius = int(sys.argv[3])
    Niter = int(sys.argv[4])
else:                           #default values
    Nx=25; # size along x
    Ny=25; # size along y 
    radius=8;   # radius of central lead 
    Niter=1500;   # number of iterations to perform
    
if (radius >= int((Ny-1)/2)):       #checking for validness
    print("The radius is larger than (Ny-1)/2. Taking default values.")
    Nx=25; 
    Ny=25; 
    radius=8;
    Niter=1500;

phi = zeros([Ny, Nx])       #creating phi

x = linspace(-int((Nx-1)/2),int((Nx-1)/2),Nx)
y = linspace(int((Ny-1)/2),-int((Ny-1)/2),Ny)
X,Y = meshgrid(x,y)         # meshgrid
print(x)

ii = where(X*X + Y*Y <= radius*radius)      # finding electrode location
phi[ii] = 1.0           # setting potential to 1
errors = zeros(Niter)
x_scale = []

for k in range(Niter):
    oldphi = phi.copy()         # using laplace to solve
    phi[1:-1,1:-1] = 0.25*(phi[1:-1,0:-2]+phi[1:-1,2:]+phi[0:-2,1:-1]+phi[2:,1:-1])
    phi[1:-1,0] = phi[1:-1,1]       # boundary conditions
    phi[1:-1,-1] = phi[1:-1,-2]
    phi[0] = phi[1]
    phi[ii] = 1.0
    errors[k] = (abs(phi-oldphi)).max()     # finding max error
    if errors[k] == 0:
        print('The error has become 0')
        exit()
    x_scale.append(k)
x_scale = array(x_scale)
semilogy(x_scale[::50], errors[::50], '-r', label='errors semilogy')        # plotting error
s = []
errors1 = []
for k in range(Niter):
    s.append([1, x_scale[k]])
for k in range(Niter):
    errors1.append(log(errors[k]))
errors1 = array(errors1)
s = array(s)
m = array(lstsq(s, errors1, rcond=None)[0])         # finding best fit1
A = exp(m[0])
B = m[1]
semilogy(x_scale[::50], A*exp(B*x_scale[::50]),'ob', label='f1')    #plotting
m = array(lstsq(s[500:], errors1[500:], rcond=None)[0])     #finding best fit2
A = exp(m[0])
B = m[1]
semilogy(x_scale[::50], A*exp(B*x_scale[::50]),'oy',markersize = 4, label='f2')     #plotting
legend(loc='upper right')
show()
loglog(x_scale[::50], errors[::50], 'og', label='errors loglog')        #plotting error loglog
legend(loc='upper right')
show()
c = contourf(X, Y, phi)     #plotting contour
colorbar(c)
show()
fig4 = figure(4) # open a new figure 
ax = p3.Axes3D(fig4) # Axes3D is the means to do a surface plot
title(r'The 3-D surface plot of the potential')     # plotting surface 3-D
surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
show()
x1 = linspace(0, Nx, Nx)
y1 = linspace(Ny, 0, Ny)
print(x1)
X1,Y1 = meshgrid(x1,y1)
Jx = zeros([Ny, Nx])
Jy = zeros([Ny, Nx])
Jx[1:-1,1:-1] = 0.5*(phi[1:-1, 0:-2] - phi[1:-1, 2:])   #finding current vector
Jy[1:-1,1:-1] = 0.5*(phi[0:-2, 1:-1] - phi[2:, 1:-1])
figure(5)
plot(X1[ii],Y1[ii], 'r.')           #plotting electode points
quiver(X1, Y1, -Jx, -Jy, scale=0.85, headwidth=2, scale_units='inches')     #plotting vector current
show()











