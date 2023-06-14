# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 00:12:11 2020

@author: venkat sai
"""
from pylab import *
import scipy.signal as sp

p = poly1d([1,1,2.5])       # defining polynomial
q = poly1d([1,0,2.25])
w = polymul(p,q)        # multiplying poly
X = sp.lti([1,0.5],w)       # transfer function
t, x = sp.impulse(X, None, linspace(0, 50, 2000))   # getting x(t)
plot(t, x, label='x with 0.5 decay')
legend(loc= 'upper right')
xlabel(r'$time->$', size=15)
ylabel(r'$x(t)->$', size=15)
show()
p_1 = poly1d([1,0.1,2.2525])        # with decay = 0.05
w_1 = polymul(p_1,q)
X_1 = sp.lti([1,0.5],w_1)
t, x_1 = sp.impulse(X_1, None, linspace(0, 50, 2000))
plot(t, x_1, label='x with 0.05 decay')
legend(loc= 'upper right')
xlabel(r'$time->$', size=15)
ylabel(r'$x(t)->$', size=15)
show()
H = sp.lti([1], [1,0,2.25])     # impulse response
def f(t, w):                # defining f
    f = cos(w*t)*exp(-0.05*t)
    return f

w = 1.4
for k in range(1, 5):
    t, y, svec = sp.lsim(H,f(t, w),linspace(0, 50, 2000))       # finding output
    w += 0.05                                   # varying frequency
    plot(t, y, label=f"w = {w}")
    legend(loc= 'upper right')
xlabel(r'$time->$', size=15)
ylabel(r'$x(t)->$', size=15)
show()

Y = sp.lti([2],[1,0,3,0])           # plotting y from spring prob
t, y = sp.impulse(Y, None, linspace(0, 20, 2000))
plot(t, y, label='y')
X = sp.lti([1,0,2],[1,0,3,0])       # plotting x from spring prob
t, x = sp.impulse(X, None, linspace(0, 20, 2000))
plot(t, x, label='x')
xlabel(r'$time->$', size=15)
ylabel(r'$x,y ->$', size=15)
legend(loc= 'upper right')
show()

H = sp.lti([1],[10**(-12), 10**(-4), 1])    #H for circuit
w, S, phi = H.bode()
subplot(2,1,1)
semilogx(w,S,label='magnitude')         # plotting bode
legend(loc= 'upper right')
subplot(2,1,2)
semilogx(w,phi,label='phase')       # plotting phase
xlabel(r'$time->$', size=15)
legend(loc= 'upper right')
show()

def v(t):               # defining input
    v = cos((10**3)*t) - cos((10**6)*t)
    return v

t = linspace(0, 10**(-2), 10**(5))
t1, y, svec = sp.lsim(H, v(t), t)   #plotting output
plot(t, y)
xlabel(r'$time->$', size=15)
ylabel(r'$output->$', size=15)
show()


