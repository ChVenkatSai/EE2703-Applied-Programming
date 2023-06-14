# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:59:09 2020

@author: venkat sai
"""
from sympy import *
import pylab as p
import scipy.signal as sp

s = symbols('s')
def lowpass(R1,R2,C1,C2,G,Vi):      
    A = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0], [0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b = Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A, b, V)


A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vout = V[3]
lp = lambdify(s,Vout,'numpy')
ww = p.logspace(0,8,801)
ss = 1j*ww
v = lp(ss)
h2 = simplify(Vout)
numer,denom = h2.as_numer_denom()
num = poly(numer, s)
den = poly(denom, s)
Hlp = sp.lti([float(i) for i in num.all_coeffs()],[float(i) for i in den.all_coeffs()])
p.loglog(ww,abs(v),lw=2)
p.title("Low pass filter")
p.grid(True)
p.show()


def highpass(R1,R3,C1,C2,G,Vi):
    A = Matrix([[0,0,1,-1/G],[-s*C2*R3/(1+s*C2*R3),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    B = Matrix([0,0,0,-s*C1*Vi])
    V = A.inv()*B
    return (A,B,V)


A,B,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vout = V[3]
hf = lambdify(s,Vout,'numpy')
ww = p.logspace(0,8,801)
ss = 1j*ww
v = hf(ss)
h3 = simplify(Vout)
numer,denom = h3.as_numer_denom()
num = poly(numer, s)
den = poly(denom, s)
Hhp = sp.lti([float(i) for i in num.all_coeffs()],[float(i) for i in den.all_coeffs()])
p.loglog(ww,abs(v),lw=2)
p.title("High pass filter")
p.grid(True)
p.show()

t = p.linspace(0,1.5e-3,50000)
v = p.cos(2*(10**6)*p.pi*t) + p.sin(2000*p.pi*t)
y = sp.lsim(Hlp, v, t)[1]   
p.plot(t,y)
p.title("Input response to low pass filter(t)")
p.show()
y1 = sp.lsim(Hhp, v, t)[1]   
p.plot(t,y1)
p.title("Input response to high pass filter(t)")
p.show()
v = p.array([1 for i in t])
y = sp.lsim(Hlp, v, t)[1]  
p.title("Step response to low pass filter(t)") 
p.plot(t,y)
p.show()
y1 = sp.lsim(Hhp, v, t)[1]   
p.title("Step response to high pass filter(t)") 
p.plot(t,y1)
p.show()

v = p.sin(500*t)*(p.e**(-0.05*t))
y = sp.lsim(Hhp, v, t)[1] 
p.title("Damped sin response to high pass filter(t)") 
p.plot(t,y)
p.show()

