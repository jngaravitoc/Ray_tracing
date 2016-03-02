import numpy as np
from scipy.integrate import quad


def a_f(M, z):
    return 1/(1+z)

def x_f(M, z, Omega_m, Omega_L):
    A = a_f(M, z)
    xx = (Omega_L / Omega_m)**(1/3.0) * A
    return xx

def model(x):
    return x**(3.0/2.0)/(1.0+x**3.0)**(3.0/2.0)

def integral_x(X):
    x = np.linspace(0, 1, 100)
    I, err = quad(model, 0, X
    return I

def D_f(a):
    xx = x_f(M, z, Omega_m, Omega_L)
    D_a = 5.0/2.0*(Omega_m0/Omega_L0)**(1.0/3.0)*(1+x**3.0)**(0.5)/(x**(3.0/2.0))*integral_x(xx)
    return D_a

def sigma(M, a):
    y = 1/(M/1E12)
    Da = D_f(a)
    return Da*16.9*y**(0.41)/(1+1.102*y**(0.2)+6.22*y**(0.333))


