import numpy as np
import matplotlib.pyplot as plt
import numba as nb
import time

alpha = 1/137
D = 1./1.79**3 
dtda_part2 = 2*np.pi/3
f_pi = 131 
Gf = 1.166*10**-11 
me = .511        
mpi_neutral = 135  
mpi_charged = 139.569  
mPL = 1.124*10**22 
mu = 105.661  
eps_e = me/mu
Enumax = (mu/2)*(1-(eps_e)**2)
n = 10
x_values, w_values = np.polynomial.laguerre.laggauss(n)  
x_valuese, w_valuese = np.polynomial.legendre.leggauss(n)

@nb.jit(nopython=True)
def I1(eps,x): #Energy Density
    numerator = (np.e**eps)*(eps**2)*((eps**2+x**2)**.5)
    denominator = np.e**((eps**2+x**2)**.5)+1
    return numerator/denominator

@nb.jit(nopython=True)
def I2(eps,x): #Pressure
    numerator = (np.e**eps)*(eps**4)
    denominator = ((eps**2+x**2)**.5)*(np.e**((eps**2+x**2)**.5)+1)
    return numerator/denominator

@nb.jit(nopython=True)
def dI1(eps,x): #Derivative of Energy Density
    numerator = (np.e**eps)*((eps**2+x**2)**.5)
    denominator = np.e**((eps**2+x**2)**.5)+1
    return (-x)*numerator/denominator

@nb.jit(nopython=True)
def dI2(eps,x): #Derivative of Pressure
    numerator = (np.e**eps)*3*(eps**2)
    denominator = ((eps**2+x**2)**.5)*(np.e**((eps**2+x**2)**.5)+1)
    return (-x)*numerator/denominator

@nb.jit(nopython=True)
def calculate_integral(I,x): #I is the function to integrate over, x is me/temp 
    return np.sum(w_values*I(x_values,x))  

@nb.jit(nopython=True)
def trapezoid(y_array,x_array):
    total = np.sum((x_array[1:]-x_array[:-1])*(y_array[1:]+y_array[:-1])/2)
    return total

@nb.jit(nopython=True)
def derivatives(a,y): #only calculates dTda (y[0]) and dtda (y[1])
    d_array = np.zeros(len(y))
    Tcm = 1/a 

    dtda_part1 = mPL/(2*a)
    dtda_part3 = (y[0]**4*np.pi**2)/15
    dtda_part4 = 2*y[0]**4*calculate_integral(I1,me/y[0])/np.pi**2
    dtda_part7 = (7*np.pi**2)/40 * Tcm**4
    dtda = dtda_part1/(dtda_part2*(dtda_part3+dtda_part4+dtda_part7))**.5
    d_array[1] = dtda 

    dTda_constant1 = (4*np.pi**2/45)+(2/np.pi**2)*(calculate_integral(I1,me/y[0]) + (1/3)*calculate_integral(I2,me/y[0]))
    dTda_constant2 = 2*me*y[0]*a**3/(np.pi**2)
    dTda_numerator1 = -3*a**2*y[0]**3*dTda_constant1
    dTda_numerator2 = 0 #dQda/y[0]
    dTda_denominator = (3*y[0]**2*a**3*dTda_constant1) - (dTda_constant2*(calculate_integral(dI1,me/y[0]) + (1/3)*calculate_integral(dI2,me/y[0])))
    dTda = (dTda_numerator1 + dTda_numerator2)/dTda_denominator
    d_array[0] = dTda

    return d_array