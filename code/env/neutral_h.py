import numpy as np

NSSH = 4.0E-3 # See table 2 of Rahmati et al 2013
GUVB = 4.5E-13 # See table 2 of Rahmati et al 2013
#Z = 6.0

#def NHtonh(NH):


def lamda_t(t): # cm3/s equation a6, rahmati et al 2013
    return 1.17E-10 * t**(0.5) * np.exp(-157809.0/t) / (1.0 + np.sqrt(t/10**5))

def alpha_a(t): # cm3 / s equation a3
    l = 315614/t
    return 1.269E-13 * l**(1.503) / (1 + (l/0.522)**(0.47))**1.923

def gammaphot(gammauvb, nh, nhssh):
    gammap =  gammauvb*(0.98*(1.0+ (nh /nhssh)**(1.64))**(-2.28) + 0.22 * (1 + nh  / nssh)**(-0.84))
    return gammap

def abc(t, g, nh, nssh):
    a_a = alpha_a(t)
    lt = lamda_t(t)
    gphot = gammaphot(g, nh, nssh)
    a = a_a + lt
    b = 2.0*a_a + gphot / nh + lt
    c = a_a
    return a, b, c

def eta(a, b, c):
    eta = (b - np.sqrt(b**2.0 - 4.0*a*c)) / (2.0*a)
    return eta

def tvir(m, z):
    t = 2554.0 * (m / 1E6)**(2/3.)*((1 + z) / 31.0)
    return t
