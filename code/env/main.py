"""

"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy import constants
from astropy import units as u
import warnings
#from scipy.special import gamma
#from scipy.integrate import quad
from sklearn.neighbors import KDTree
from neutral_h  import *
from ray_tracing import *
from density_profile import *
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


P = 10000
# N = number of emmiters
# M = number of traced rays in a radom direction

def reading_data(data):
    halos = np.loadtxt(data)
    h = 0.7
    ids = halos[:,0]
    x = halos[:,1]/h
    y = halos[:,2]/h
    z = halos[:,3]/h
    M = halos[:,4]/h
    R = halos[:,5]/h
    return ids, x, y, z, M, R

def concentration(M):
    c = 9.6*(M/1E12)**(-0.075)
    return c

# Read the data
ids, x, y, z, M, R = reading_data(data)

# Found possible emitter halos
ids_hh, x_hh, y_hh, z_hh, D3_mean, M_hh = host_halos(M, x, y, z, ids)

print M_hh

#Choose a random halo, which would be the emitter of the Lya photon
#x_in, y_in, z_in, id_in, D_env = random_halo(ids_hh, x_hh, y_hh, z_hh, D3_mean)
enviroment()
#Select the possible absorbers halos, by making a cube of length P
#centered in the emitter halo.
x_cube, y_cube, z_cube, R_cube, M_cube, ids_cube = selecting_halos(x_in, y_in, z_in, P, x, y, z, R, M, ids)

#Found a random direction in which the Lya would be emitted
print D_env, M_hh
nh_1 = np.zeros(100)

for i in range(100):
    x_out, y_out, z_out = random_direction(x_in, y_in, z_in, P)
    #Compute the impact parameter (b) to all the possible halos, if b is
    #less than the virial radiues the photon might be emitted.
    b_abs, xabs, yabs, zabs, R_abs, M_abs, ids_abs = impact_parameter(ids_cube, x_cube, y_cube, z_cube, P, x_in, y_in,
    z_in, R_cube, M_cube, x_out, y_out, z_out)
    c_abs = concentration(M_abs)
    nh_1[i] = nh(R_abs, b_abs, c_abs)
    print nh_1[i]
