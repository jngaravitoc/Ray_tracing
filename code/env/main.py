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

# Read the data
#ids, x, y, z, M, R = reading_data(data)

# Found possible emitter halos
#ids_hh, x_hh, y_hh, z_hh, D3_mean = host_halos(M, x, y, z, ids)

#Choose a random halo, which would be the emitter of the Lya photon
#x_in, y_in, z_in, id_in, D_env = random_halo(ids_hh, x_hh, y_hh, z_hh, D3_mean)

#Select the possible absorbers halos, by making a cube of length P
#centered in the emitter halo.
#x_cube, y_cube, z_cube, R_cube, M_cube, ids_cube = selecting_halos(x_in, y_in, z_in, P, x, y, z, R, M, ids)

#Found a random direction in which the Lya would be emitted
#x_out, y_out, z_out = random_direction(x_in, y_in, z_in, P)

#Compute the impact parameter (b) to all the possible halos, if b is
#less than the virial radiues the photon might be emitted.
#babs, xabs, yabs, zabs, R_abs, M_abs, ids_abs = impact_parameter(ids_cube, x_cube, y_cube, z_cube, P, x_in, y_in,
#z_in, R_cube, M_cube, x_out, y_out, z_out)

bb = np.linspace(0, 260, 100)
print len(bb)
nh_b = np.zeros(100)
for i in range(100):
    nh_b[i] = nh(260, bb[i])

plt.plot(bb, nh_b)
plt.ylabel('$N_H(g cm^{-3})$')
plt.xlabel('$Impact parameter b(kpc)$')
plt.show()

#print D_env, len(babs)


#plt.scatter(x, y, s=0.1, c='k')
#plt.scatter(x_hh, y_hh, s=5, c='k')
#plt.scatter(x_cube, y_cube, s=1, c='g')
#plt.scatter(x_in, y_in, s=40, c='b')
#plt.scatter(x_out, y_out, s=60, c='r')
#plt.show()

    #NHT = []
    #NHI = []

        #print 'for one emitter'
"""
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
