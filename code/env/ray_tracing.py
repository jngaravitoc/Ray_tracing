import numpy as np
from sklearn.neighbors import KDTree


def host_halos(M, x, y, z, ids):
    host_halo = np.where((M<20) & (M>7))[0]
    M_hh = M[host_halo]
    x_hh = x[host_halo]
    y_hh = y[host_halo]
    z_hh = z[host_halo]
    id_hh = ids[host_halo]
    r3 = []
    D = np.array([x, y, z])
    D = D.T
    tree = KDTree(D, leaf_size=20000)
    for i in range(len(x_hh)):
        dist, ind = tree.query(D[host_halo], k=4)
        r3.append(max(dist[0]))
    return id_hh, x_hh, y_hh, z_hh, np.mean(r3), M_hh

# Finding the environment \Delta_3 of a given halo:

def environment(x_h, y_h, z_h, x, y, z, D3):
    DD = np.array([x, y, z])
    DD = DD.T
    tree = KDTree(DD, leaf_size=20000)
    index = np.where(x_h == x)[0]
    dist, ind = tree.query(DD[index], k=4)
    r3 = max(dist[0])
    #print r3
    delta3 = D3**3.0 * (1.0/(r3**3.0) - 1.0/(D3**3.0))
    return  delta3

# Function that choose a random halo.

def random_halo(id_hh, x, y, z, D3):
    N = len(id_hh)
    ran_halo = np.random.randint(0, N)
    id_h = id_hh[ran_halo]
    x_h = x[ran_halo]
    y_h = y[ran_halo]
    z_h = z[ran_halo]
    delta3 = environment(x_h, y_h, z_h, x, y, z, D3)
    return x_h, y_h, z_h, id_h, delta3

# Selecting a cube around the emmiter galaxy!

def selecting_halos(x_in, y_in, z_in, r, x, y, z, R, M, ids):
    x_max = x_in + r
    x_min = x_in - r
    y_max = y_in + r
    y_min = y_in - r
    z_max = z_in + r
    z_min = z_in - r
    cube = np.where((x<x_max) & (x>x_min) & (y<y_max) & (y>y_min) & (z<z_max) & (z>z_min))
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
    sin_ran = np.sqrt(1.0-cos_ran**2.0)
    phi_ran = 2.0*np.pi*np.random.random()
    x_end = r*sin_ran*np.cos(phi_ran)
    y_end = r*sin_ran*np.sin(phi_ran)
    z_end = r*cos_ran
    x_out = x + x_end
    y_out = y + y_end
    z_out = z + z_end
    return x_out, y_out, z_out


#Computing the impact parameter of the halos to identify the absorbers
# xout, yout, zout is the random direction of the ray. x_in, y_in,
# z_in are the positions of the emitter halo. x, y, z, all the
# possible abosorbers positions.

def impact_parameter(id_cube, x, y, z, L, x_in, y_in, z_in, R, M, x_out, y_out, z_out):
    d = np.sqrt((x-x_in)**2.0 + (y-y_in)**2.0 + (z-z_in)**2.0)
    dot_product = (x-x_in)*(x_out-x_in) + (y-y_in)*(y_out-y_in) + (z-z_in)*(z_out-z_in)
    a_mag = L
    b_mag = d
    costheta = dot_product / (a_mag * b_mag)
    theta = np.arccos(costheta)
    b = np.sin(theta)*d
    absorbers = np.where((b<R) & (dot_product>0))
    b_abs = b[absorbers]
    x_abs = x[absorbers]
    y_abs = y[absorbers]
    z_abs = z[absorbers]
    R_abs = R[absorbers]
    M_abs = M[absorbers]
    id_abs = id_cube[absorbers]
    Vol = 4.0 / 3.0 * np.pi * R_abs**3.0
    rho =  M_abs / Vol
    return b_abs, x_abs, y_abs, z_abs, R_abs, M_abs, id_abs
