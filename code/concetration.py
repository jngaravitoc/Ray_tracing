"""
Code to compute the concentration (c_200) of a halo given the mass
M200 and the redshift z derived by Prada et al.(2012).

function:

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Global constants
c_0 = 3.681
c_1 = 5.033
alpha = 6.948
beta = 7.386
x_0 = 0.424
x_1 = 0.526
i_sigma_0 = 1.047
i_sigma_1 = 1.646

A = 2.881
b = 1.257
c = 1.022
d = 0.060
h = 0.7

Omega_m0 = 0.3
Omega_L0 = 0.7

# Functions:

# scale factor
def a_f(z):
    return 1/(1+z)

# x
def x_f(A):
    xx = (Omega_L0 / Omega_m0)**(1.0/3.0) * A
    return xx

def integral_x(X):
    func = lambda xx: xx**(3.0/2.0)/(1.0+xx**3.0)**(3.0/2.0)
    I , err = quad(func, 0, X)
    return I

# D(a)  eq. 12
def D_f(x):
    D_a = 5.0/2.0*(Omega_m0/Omega_L0)**(1.0/3.0)*(1.0+x**3.0)**(0.5)/(x**(3.0/2.0))*integral_x(x)
    return D_a

def sigma(M, x):
    y = (1E12/M)*h
    Da = D_f(x)
    return Da*16.9*y**0.41/(1.0+1.102*y**0.2+6.22*y**0.333)

def c_min(x):
    cmin = c_0 + (c_1 - c_0)*(1/np.pi * np.arctan(alpha*(x-x_0)) + 0.5)
    return cmin

def s_min_1(x):
    smin_1 = i_sigma_0+(i_sigma_1 - i_sigma_0)*(1.0/np.pi*np.arctan(beta*(x-x_1)) + 0.5)
    return smin_1

def B_0(x):
    B0 =  c_min(x)/c_min(1.393)
    return B0

def B_1(x):
    B1 = s_min_1(x)/s_min_1(1.393)
    return B1

def sigma_prime(M, x):
    return B_1(x)*sigma(M, x)

def CC(M, x):
    s_p = sigma_prime(M, x)
    return A*((s_p/b)**c + 1.0)*np.exp(d/s_p**2.0)

def cc(M , z):
    a = a_f(z)
    x_x = x_f(a)
    return B_0(x_x) * CC(M, x_x)


M = np.logspace(11, 15, 40)
con = cc(M, 0.0)

#z = np.linspace(0, 6, 20)
#a_test = a_f(z)
#x_test = x_f(a_test)
#c_min_t = c_min(x_test)
plt.plot(M, np.log10(con))
plt.ylim(0.55, 1.0)
#plt.yscale('log')
plt.xscale('log')
plt.show()
