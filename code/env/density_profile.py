import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy.special import gamma
from scipy.integrate import quad
#Constants:
G = const.G
G = G.to(u.cm**3/u.g/u.s**2)
b = 0.7
f_gas = 1.0
Omega_m = 0.32
Omega_b = 0.0455102
H_0 = 2.19E-18*u.Hz #1/s
rho_c = 3.0*H_0**2.0 / (8.0*np.pi*G) # critical density
cc = 3.75 # Concentration of the halo a z=6
mp = const.m_p
mp = mp.to(u.g)

def concentration(M):
    #This values from c0 and M0 come from Klypin et al 2011. for z=5
    c0 = 2.3
    M0 = 6E11
    cc = c0*(M/1E12)**(-0.075)*(1.0 + (M/M0)**(0.26))
    return cc

def H(z, Omega0, H_0):
    Lambda0 = 1. - Omega0
    return H_0*(Omega0*(1.0+z)**3.0 - (Omega0+Lambda0-1.0)*(1.0+z)**2.0 + Lambda0)**0.5

def Omega_z(z,Omega0, H_0):
    return Omega0*(1.0+z)**3.0*(H_0/H(z,Omega0, H_0))**2.0

def model(x):
    return x**2.0*((1.0+x)**(27.0*b/(2.0*x)))

def integral(c):
    x = np.linspace(0, c, 100)
    integral_scipy,err = quad(model, 0, c)
    return I

def rho_g0(zform, c): #g/cm3
    d_c = 3000.0*Omega_m*(1.0+zform)**3
    return (((f_gas * d_c * rho_c * Omega_b) / Omega_m) * np.exp(27.*b/2.0)*(np.log(1.0+c) - c/(1.0+c)))/integral(c)

def rho(r_s, z):
    #r_s = r_vir/c
    total_rho = rho_g0(z) * np.exp**(-27.0*b/2.0) * (1.0+(r/r_s))**(27.0*b/(2.0*r/r_s))
    nh = rrho / mp
    return total_rho, nh

#B is the impact parameter.
#Column density computation see eq.5 of the document
def nh(rvir, B, M):
    c = concentration(M)
    rhog0 = rho_g0(6.0, c)#gm/cm3
    B = B*u.kpc
    B = B.to(u.cm)
    rvir = rvir*u.kpc
    rvir = rvir.to(u.cm)
    A = -0.178*b + 0.982
    beta = 0.9*b
    r_c =  0.22*(rvir / c)
    #if len(B)==0:
    #    NH = np.zeros(1)
    #else:
    NH = (np.sqrt(np.pi)*(1.0/r_c**2.0)**(-3.0*beta/2.0)*(B**2.0+r_c**2.0)**(1.0/2.- 3*beta/2.) * gamma(-0.5 + 3.0*beta/2.))/(2.0*gamma(3.0*beta/2.))
    return (A*NH*rhog0/mp).value
