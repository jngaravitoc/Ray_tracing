import numpy as np
import sys
import warnings
from scipy.special import gamma
from scipy.integrate import quad
#from scipy import integrate



warnings.simplefilter('ignore')
sys.stderr.write("Usage: " + sys.argv[0] + "test.txt 10000 1(Number of random halos) 1(NUmber of random directions) 4(Number of nearest neighboors: env)\n")



data = sys.argv[1]
P = float(sys.argv[2])
N = int(sys.argv[3])
K = int(sys.argv[4])
env = int(sys.argv[5])


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


#print "The number of halos of this catalogue is: ", len(X)

# Finding halos in a catalogue for a mass between 7^10 - 2^11 Msun

def host_halos(M, x, y, z, ids):
    host_halo = np.where( (M<20) & (M>7))
    M_hh = M[host_halo]
    x_hh = x[host_halo]
    y_hh = y[host_halo]
    z_hh = z[host_halo]
    id_hh = ids[host_halo]
    D3_mean = []
    for i in range(len(M_hh)):
	D = []
	for j in range(len(M_hh)):
		dist  = np.sqrt( (x_hh[i] - x_hh[j])**2 + (y_hh[i]-y_hh[j])**2 + (z_hh[i] - z_hh[j])**2  )
		D.append(dist)
	D3 = np.sort(D)
	D3_mean.append(D3[2])
    return id_hh, x_hh, y_hh, z_hh, np.mean(D3_mean)

# Finding a random emmiter hal, env):

def random_halo(id_hh, x, y, z, env, D3):
    N = len(id_hh)
    ran_halo = np.random.randint(0, N)
    id_h = id_hh[ran_halo]
    x_h = x[ran_halo]
    y_h = y[ran_halo]
    z_h = z[ran_halo]
    D = []
    for i in range(len(x)):
    	d = np.sqrt((x_h-x[i])**2 + (y_h-y[i])**2 + (z_h-z[i])**2)
	D.append(d)
    Dist = np.sort(D)
    r3 = Dist[2]
    delta3 = D3**3 * (1/(r3**3) - 1/(D3**3))   
    return x_h, y_h, z_h, id_h,  delta3

# Selecting a cube around the emmiter galaxy!

def selecting_halos(x_in, y_in, z_in, r, x, y, z, R, M, ids):
    x_max = x_in + r
    x_min = x_in - r
    y_max = y_in + r
    y_min = y_in - r
    z_max = z_in + r
    z_min = z_in - r
    cube = np.where( (x < x_max) & (x > x_min) & ( y < y_max )  & ( y > y_min) & (z < z_max ) & (z > z_min))
    
    x_cube = x[cube]
    y_cube = y[cube]
    z_cube = z[cube]
    R_cube = R[cube]
    M_cube = M[cube]
    id_cube = ids[cube]
      
    return x_cube, y_cube, z_cube, R_cube, M_cube, id_cube

# Setting a random direction

def random_direction(x, y, z, r):
    cos_ran = 2.0*np.random.random()-1.0
    sin_ran = np.sqrt(1.0-cos_ran**2)
    phi_ran = 2.0*np.pi*np.random.random()
    x_end = r*sin_ran*np.cos(phi_ran)
    y_end = r*sin_ran*np.sin(phi_ran)
    z_end = r*cos_ran
    x_out = x + x_end
    y_out = y + y_end
    z_out = z + z_end
    return x_out, y_out, z_out
# 

def impact_parameter(id_cube, x, y, z, L, x_in, y_in, z_in, R, M, x_out, y_out, z_out):
    d = np.sqrt((x - x_in)**2 + (y - y_in)**2 + (z - z_in)**2)
    dot_product = (x-x_in ) * (x_out - x_in) + (y -y_in)*(y_out - y_in) + (z - z_in)* (z_out - z_in)
    a_mag = L
    b_mag = d
    costheta = dot_product / (a_mag * b_mag)
    theta = np.arccos(costheta) 
    b = np.sin(theta)*d
    absorbers = np.where( (b<R) & (dot_product>0))
    b_abs = b[absorbers]
    x_abs = x[absorbers]
    y_abs = y[absorbers]
    z_abs = z[absorbers]
    R_abs = R[absorbers]
    M_abs = M[absorbers]
    id_abs = id_cube[absorbers]
    Vol = 4 / 3. * np.pi * R_abs**3
    rho =  M_abs / Vol
    return b_abs, x_abs, y_abs, z_abs, R_abs, M_abs, id_abs


def H(z, Omega0, H_0):
    Lambda0 = 1. - Omega0
    return H_0*(Omega0*(1+z)**3 - (Omega0+Lambda0-1)*(1+z)**2 + Lambda0)**0.5

def Omega_z(z,Omega0, H_0):
    return Omega0 * (1+z)**3 * (H_0/H(z,Omega0, H_0))**2

def model(x):
    b = 0.7
    return x**2 * ( (1+x)**(27*b/(2.*x)) )


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
    #zform = 6.
    d_c = 3000.0 * Omega_m * (1 + zform)**3
    #d_c = 200/3. * (c**3 / (log(1 + c) - c/(1+c)))
    return ( ( (f_gas * d_c * rho_c * Omega_b) / Omega_m ) * np.exp(27. * b / 2.) * (np.log(1+c) - c/(1+c)) )  /  integral(c)


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
    NH = (np.sqrt(np.pi)*(1/r_c**2)**(-3*beta/2.) * (B**2+ r_c**2)**(1/2. - 3*beta/2.) * gamma(-0.5 + 3*beta/2.) )/(2*gamma(3*beta/2.))
    return A*NH*rho/mp

# -----------------   NHI computation ---------------------------------------------------------------------

def LamdaT(T): # cm3/s Equation A6, Rahmati et al 2013
    return 1.17E-10 * T**(0.5) * np.exp(-157809/T) / ( 1 + np.sqrt(T/10**5))

def alphaA(T): # cm3 / s Equation A3
    l = 315614/T
    return 1.269E-13 * l**(1.503) / ( 1 + (l/0.522)**(0.47) )**1.923 

def Gammaphot(GammaUVB, nh, NHSSH):
    Gammap =  GammaUVB * ( 0.98 *  (  1 + ( nh  / NHSSH )**(1.64) )**(-2.28) + 0.22 * ( 1 + nh  / NSSH)**(-0.84) )
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
    eta = ( B - np.sqrt(B**2 - 4*A*C) ) / (2*A)
    #print "B", B
    #print "the left part", np.sqrt(B**2 - 4*C*A)
    return eta 

def tvir(M, z):
    T = 2554 * ( M / 1E6)**(2/3.) * ( (1 + z) / 31.0 )
    return T

idsh, x_hh, y_hh, z_hh, D3_mean = host_halos(M, x, y, z, ids)
print D3_mean

#print "#emmiter ID, Delta3 , Rvir(kpc), NH(1/cm2)"
print "#emmiter ID, Env, Total NH(1/cm2)"
for i in range(N):
	x_in, y_in, z_in, id_in, D_env = random_halo(idsh, x_hh, y_hh, z_hh, env, D3_mean)
	x_cube, y_cube, z_cube, R_cube, M_cube, ids_cube = selecting_halos(x_in, y_in, z_in, P, x, y, z, R, M, ids)
	NHT = []
	NHI = []
	for j in range(K):
		x_out, y_out, z_out = random_direction(x_in, y_in, z_in, P)
		babs, xabs, yabs, zabs, R_abs, M_abs, ids_abs = impact_parameter(ids_cube, x_cube, y_cube, z_cube, P, x_in, y_in, z_in, R_cube, M_cube, x_out , y_out, z_out)
		NH = nh(R_abs, babs)
		A = np.zeros(len(NH))
		B = np.zeros(len(NH))
		C = np.zeros(len(NH))
  		E = np.zeros(len(NH))
		T = tvir(10*M_abs*np.exp(10), Z)
		#print T
		for k in range(len(A)):
			A[k], B[k], C[k] = ABC(T[k], GUVB, NH[k], NSSH )
    			E[k] = Eta(A[k], B[k], C[k])
			#print E[k], A[k], B[k], C[k], NH[k]
			#rint E[k]
		if (NH.sum() > 0):
			NHT.append(NH.sum())
			NHI.append(E.sum())
		#print E
	for i in range(len(NHT)):
		print  int(id_in), D3_mean, NHT[i], NHI[i]
	        #print int(id_in), D_env, NH


