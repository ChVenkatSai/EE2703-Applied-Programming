# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:19:54 2020

@author: venkat sai
"""
import math
from pylab import *
import scipy.integrate as teg

m = 0

def exponen(a):           # defining exponential function
    b = []
    if isinstance(a, float):        # for float values
        return math.exp(a)
    else:
        for i in range(len(a)):         # for vectors
            b.append(math.exp(a[i]))
        return b


def coscos(a):          # defining coscos
    b = []
    if isinstance(a, float):
        return math.cos(math.cos(a))        # for float vlued
    else:
        for i in range(len(a)):             # for vectors
            b.append(math.cos(math.cos(a[i])))
        return b
    

while(m!=2):            # 1 for exponen 2 for coscos
    V = []
    N = []
    
    x = linspace(0, 2*pi, 401)
    x = x[:-1]              # drop last term
    if m == 0:
        b = exponen(x) 
    elif m == 1:
        b = coscos(x)
    A = zeros((400, 51))            # creating A matrix
    A[:, 0] = 1                     # 1st column 1
    for k in range(1, 26):
        A[:, 2*k-1] = cos(k*x)
        A[:, 2*k] = sin(k*x)
    if m == 0:
        c1 = lstsq(A, b, rcond=None)[0]         # best fits
    if m == 1:
        c2 = lstsq(A, b, rcond=None)[0]
    if m == 0:
        u = lambda x, k : exponen(float(x))*cos(k*x)        # inteegral functions
        v = lambda x, k : exponen(float(x))*sin(k*x)
    elif m ==1:
        u = lambda x, k : coscos(float(x))*cos(k*x)
        v = lambda x, k : coscos(float(x))*sin(k*x)   

    for k in range(26):
        a = teg.quad(u, 0, 2*math.pi, args=(k))[0]      # inntegrating
        if k == 0:
            a = a/(2*math.pi)
        else:
            a = a/(math.pi)
        N += [[k]]
        V += [[a]]                  # appending correspong value
        if k != 0:
            b = teg.quad(v, 0, 2*math.pi, args=(k))[0]
            b = b/(math.pi)
            N += [[k]]
            V += [[b]]
    V_1 = array(V)
    N_1 = array(N)
    if m == 0:                  #       calculating deviation
        dev1 = max(abs(c_[c1]-V_1))
        print(f'The maximum deviation in the case of exponential is {dev1}')
    if m == 1:
        dev2 = max(abs(c_[c2]-V_1))
        print(f'The maximum deviation in the case of coscos is {dev2}')
    if m  == 0:                 # plotting bestfit vs coeff
        w = semilogy(N_1, abs(V_1), 'or', label='exponen coeff')
        semilogy(N_1, abs(c1), 'og',markersize=3, label='exponen bestfit')
    elif m == 1:
        semilogy(N_1, abs(V_1), 'or', label='coscos coeff')
        semilogy(N_1, abs(c2), 'og',markersize=3, label='coscos bestfit')
    legend(loc = 'upper right')
    grid(True)
    show()
    if m == 0:              # using loglog
        loglog(N_1, abs(V_1), 'or', label='exponen coeff')
        loglog(N_1, abs(c1), 'og',markersize=3, label='exponen bestfit')
    elif m == 1:
        loglog(N_1, abs(V_1), 'or', label='coscos coeff')
        loglog(N_1, abs(c2), 'og',markersize=3, label='coscos bestfit')
    legend(loc = 'upper right')
    grid(True)
    show()
    m += 1
    
b1 = dot(A, c1)             # fourier fn with bestfit coeff.
b2 = dot(A, c2)
x_scale = linspace(-2*math.pi, 4*math.pi, 100)
x_scale = array(x_scale)        # plotting actual
semilogy(x_scale, (exponen(x_scale)), '-r', label='Exponential')
x_scale_1 = linspace(0, 2*pi)
x_scale_1 = array(x_scale_1)            # plotting expected
semilogy(x_scale_1+2*pi, exponen(x_scale_1), '--b')
semilogy(x_scale_1-2*pi, exponen(x_scale_1), '--b')
semilogy(x_scale_1, exponen(x_scale_1), '--b', label='Exponen expect')
semilogy(x, b1, 'og', label='Exponen series')       # plotting fourier with bestfit
legend(loc = 'upper right')
grid(True)
show()
plot(x_scale, coscos(x_scale), '-g',label='CosCos')
plot(x_scale, coscos(x_scale), '--b', label='CosCos expect')
plot(x, b2, 'or', label='Coscos series')
plot()
legend(loc = 'upper right')
grid(True)
show()
    

    
    
    
    
