"""

"""

import numpy as np
import matplotlib.pyplot as plt 
import sys
import warnings
from scipy.special import gamma
from scipy.integrate import quad
from sklearn.neighbors import KDTree

warnings.simplefilter('ignore')

if len(sys.argv)<2:
    sys.stderr.write("Usage: " + sys.argv[0] + " catalogue.txt"\
                     " 10000(period) 1(Number of random halos) 1(Number "\
                     "of random directions) 3(Number of nearest neighbors: env)\n")
    sys.exit('Error: Put all the input parameters')


data = sys.argv[1]
#P = float(sys.argv[2])
#N = int(sys.argv[3])
#K = int(sys.argv[4])
#env = int(sys.argv[5])


# P = ray length
# N = number of emmiters
# M = number of traced rays in a radom direction

halos = np.loadtxt("../data/" + data)
h = 0.7
ids = halos[:,0]
x = halos[:,1]
y = halos[:,2]
z = halos[:,3]
M = halos[:,4]
R = halos[:,5]

NSSH = 4.0E-3 # See table 2 of Rahmati et al 2013
GUVB = 4.5E-13 # See table 2 of Rahmati et al 2013
Z = 6.0




def H(z, Omega0, H_0):
    Lambda0 = 1. - Omega0
    return H_0*(Omega0*(1+z)**3 - (Omega0+Lambda0-1)*(1+z)**2 + Lambda0)**0.5

def Omega_z(z,Omega0, H_0):
    return Omega0 * (1+z)**3 * (H_0/H(z,Omega0, H_0))**2

def model(x):
    b = 0.7
    return x**2 * ((1+x)**(27*b/(2.*x)))

def integral(c):
    x = np.linspace(0, c, 100)
    integral_scipy,err = quad(model, 0, c)
    return integral_scipy

def rho_g0(zform): #g/cm3
    b = 0.7
    f_gas = 1.
    Omega_m = 0.32
    Omega_b = 0.022/(0.67**2)
    H_0 = 2.19E-18 #1/s
    G = 6.67E-8 #
    rho_c = 3*H_0**2 / (8.*np.pi*G)
    c = 3.75
    #zform = 6.0
    d_c = 3000.0 * Omega_m * (1 + zform)**3 # Check this d_c
    #d_c = 200/3. * (c**3 / (log(1 + c) - c/(1+c)))
    return (((f_gas * d_c * rho_c * Omega_b) / Omega_m) * np.exp(27. * b / 2.) * (np.log(1+c) - c/(1+c)))  /  integral(c)

# see eq 10 of the document
def rho(r, z):
    c = 3.75
    r_s = r_vir/c
    mp =  1.67262158E-24
    rrho = rho_g0(z) * e**(-27*b/2.) * (1 + (r/r_s))**(27*b/(2*r/r_s))
    nh = rrho / mp
    return rrho, nh

def nh(rvir, B):
    rho = rho_g0(6)#gm/cm3
    kpc = 3.085E21
    B = B*kpc
    rvir = rvir*kpc
    mp =  1.67262158E-24
    c = 3.75
    b = 0.7
    A = -0.178*b + 0.982
    beta = 0.9*b
    r_c =  rvir / c
    if len(B)==0:
        NH = np.zeros(1)
    else:
        NH = (np.sqrt(np.pi)*(1/r_c**2)**(-3*beta/2.) * (B**2+ r_c**2)**(1/2. - 3*beta/2.) * gamma(-0.5 + 3*beta/2.))/(2*gamma(3*beta/2.))
    return A*NH*rho/mp


# --------------------------------------------------------------------------------------------------------

# -----------------------------------   NHI computation --------------------------------------------------

def LamdaT(T): # cm3/s Equation A6, Rahmati et al 2013
    return 1.17E-10 * T**(0.5) * np.exp(-157809/T) / (1 + np.sqrt(T/10**5))

def alphaA(T): # cm3 / s Equation A3
    l = 315614/T
    return 1.269E-13 * l**(1.503) / (1 + (l/0.522)**(0.47))**1.923

def Gammaphot(GammaUVB, nh, NHSSH):
    Gammap =  GammaUVB * (0.98 *  (1 + (nh  / NHSSH )**(1.64))**(-2.28) + 0.22 * (1 + nh  / NSSH)**(-0.84))
    return Gammap

def ABC(T, G, nh, NSSH):
    a_A = alphaA(T)
    LT = LamdaT(T)
    GPhot = Gammaphot(G, nh, NSSH)
    A = a_A + LT
    B = 2*a_A + GPhot / nh + LT
    C = a_A
    return A, B, C

def Eta(A, B, C):
    eta = (B - np.sqrt(B**2 - 4.0*A*C)) / (2.0*A)
    return eta

def tvir(M, z):
    T = 2554 * (M / 1E6)**(2/3.) * ((1 + z) / 31.0)
    return T

# -----------------------------------------------------------------------------------------

idsh, x_hh, y_hh, z_hh, D3_mean = host_halos(M, x, y, z, ids)
#print D3_mean
#environment(x_hh[0], y_hh[0], z_hh[0], x_hh, y_hh, z_hh, D3_mean)

#rint "#emmiter ID, Delta3. Kpc, Total NH(1/cm2), NHI(1/cm2)"

x_in, y_in, z_in, id_in, D_env = random_halo(idsh, x_hh, y_hh, z_hh, D3_mean)
#print D_env
x_cube, y_cube, z_cube, R_cube, M_cube, ids_cube = selecting_halos(x_in, y_in, z_in, 10000, x, y, z, R, M, ids)
x_out, y_out, z_out = random_direction(x_in, y_in, z_in, 10000.0)

#plt.scatter(x, y, s=0.1, c='k')
#plt.scatter(x_hh, y_hh, s=5, c='k')
plt.scatter(x_cube, y_cube, s=1, c='g')
plt.scatter(x_in, y_in, s=40, c='b')
plt.scatter(x_out, y_out, s=60, c='r')
plt.show()

    #NHT = []
    #NHI = []
"""
        #print 'for one emitter'
        
	for j in range(K):
                x_out, y_out, z_out = random_direction(x_in, y_in, z_in, P)
		babs, xabs, yabs, zabs, R_abs, M_abs, ids_abs = impact_parameter(ids_cube, x_cube, y_cube, z_cube, P, x_in, y_in, z_in, R_cube, M_cube, x_out , y_out, z_out)
	        #for k in range(len(xabs)):
                #rint D_env, len(xabs), j 
                 
		NH = nh(R_abs, babs)
                if (len(NH==1) & (NH[0]==0)):
                   NHT.append(0)
                   NHI.append(0)
		else:
		    A = np.zeros(len(NH))
		    B = np.zeros(len(NH))
		    C = np.zeros(len(NH))
  		    E = np.zeros(len(NH))
		    T = tvir(M_abs*1e10, Z)
		    for k in range(len(A)):
			A[k], B[k], C[k] = ABC(T[k], GUVB, NH[k], NSSH )
    			E[k] = Eta(A[k], B[k], C[k])
			#rint E[k], A[k], B[k], C[k], NH[k], E[k]*NH[k]
			#print E[k]
		# Some times Lya photons dont find any halos
		#f (NH.sum() > 0):
			NHT.append(NH.sum())
			NHI.append(np.sum(E*NH))
                        #print 'for one absorber'
	for i in range(len(NHT)):
		print  int(id_in), D_env, NHT[i], NHI[i]
	        
"""
