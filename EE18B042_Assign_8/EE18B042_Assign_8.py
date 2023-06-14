# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:43:00 2020

@author: venkat sai
"""
from pylab import *
x=linspace(0,2*pi,129)
x=x[:-1]
y=sin(5*x)
Y=fftshift(fft(y))/128.0
w=linspace(-64,63,128)
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin(5t)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
show()

y=sin(x)*sin(x)*sin(x)
Y=fftshift(fft(y))/128
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
ylabel(r"$[Y]$",size=16)
title("Spectrum of sin**3(t)")
xlim([-10,10])
grid(True)
subplot(2,1,2)
plot(w,angle(w),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
ylabel(r"Phase of $Y$",size=16)
xlim([-10,10])
grid(True)
show()

y=cos(x)*cos(x)*cos(x)
Y=fftshift(fft(y))/128
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
ylabel(r"$[Y]$",size=16)
title("Spectrum of cos**3(t)")
xlim([-10,10])
grid(True)
subplot(2,1,2)
plot(w,angle(w),'ro',lw=2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
ylabel(r"Phase of $Y$",size=16)
xlim([-10,10])
grid(True)
show()

t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513)
w=w[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()


y=cos(20*t + cos(5*t))
Y=fftshift(fft(y))/512.0
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-20,20])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $cos(20t + cos(5t))$")
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-20,20])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()


y=exp(-t**2/2)
Y=fftshift(fft(y))/512
Y=Y*8*pi
figure()
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-20,20])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $e(-t**2/2)$")
grid(True)
subplot(2,1,2)
ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'ro',lw=2)
xlim([-20,20])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()