import numpy as np
import sys

sys.stderr.write("Usage: " + sys.argv[0] + "data.txt 10000 \n")

data = sys.argv[1]
P  = float(sys.argv[2])
print "#Period", P

halos = np.loadtxt("../data/" + data)
h = 0.7
Ids = halos[:,0]
X = halos[:,1] / h
Y = halos[:,2] / h
Z = halos[:,3] / h
m = halos[:,7]
r = halos[:,8] / h

print "#This catalogue has ", len(X), "halos"

def bpc(ids, x, y, z, M, R, P):
    L = 75000  / h
    region1 = np.where(x<P)
    ids1 = ids[region1]
    x1 = x[region1] + L
    y1 = y[region1]
    z1 = z[region1] 
    M1 = M[region1]
    R1 = R[region1]
    
    region2 = np.where(x>(L-P))
    ids2 = ids[region2]
    x2 = x[region2] - L 
    y2 = y[region2]
    z2 = z[region2] 
    M2 = M[region2]
    R2 = R[region2]
    
    region3 = np.where(y<P)
    
    ids3 = ids[region3] 
    x3 = x[region3] 
    y3 = y[region3] + L
    z3 = z[region3] 
    M3 = M[region3]
    R3 = R[region3]
    
    region4 = np.where(y>(L-P))
    
    ids4 = ids[region4] 
    x4 = x[region4] 
    y4 = y[region4] - L
    z4 = z[region4] 
    M4 = M[region4]
    R4 = R[region4]
    
    region5 = np.where(z<P)
    
    ids5 = ids[region5] 
    x5 = x[region5] 
    y5 = y[region5] 
    z5 = z[region5] + L
    M5 = M[region5]
    R5 = R[region5]
    
    region6 = np.where(z>(L-P))
    
    ids6 = ids[region6] 
    x6 = x[region6] 
    y6 = y[region6] 
    z6 = z[region6] - L
    M6 = M[region6]
    R6 = R[region6]
    
    region7 = np.where( (x>L-P) & (y>L-P))
    
    region8 = np.where( (x>L-P) & (y<P))
    
    region9 = np.where( (x<P) & (y<P))
    
    region10 = np.where( (x<P) & (y> L-P) )
    
    N = len(x) + len(region1[0]) + len(region2[0]) + len(region3[0]) + len(region4[0]) + len(region5[0]) + len(region6[0])
    id_new = np.zeros(N)
    x_new = np.zeros(N)
    y_new = np.zeros(N)
    z_new = np.zeros(N)
    M_new = np.zeros(N)
    R_new = np.zeros(N)
    
    N1 = len(x) + len(region1[0])
    N2 = len(x) + len(region1[0]) + len(region2[0])
    N3 = len(x) + len(region1[0]) + len(region2[0]) + len(region3[0])
    N4 = len(x) + len(region1[0]) + len(region2[0]) + len(region3[0]) + len(region4[0])
    N5 = len(x) + len(region1[0]) + len(region2[0]) + len(region3[0]) + len(region4[0]) + len(region5[0])
    N6 = len(x) + len(region1[0]) + len(region2[0]) + len(region3[0]) + len(region4[0]) + len(region5[0]) + len(region6[0])

    for i in range(len(x)):
        id_new[i] = ids[i]
        x_new[i] = x[i]
        y_new[i] = y[i]
        z_new[i] = z[i]
        M_new[i] = M[i]
        R_new[i] = R[i]
        
    for i in range(len(x), N1):
        id_new[i] = ids1[i-len(x)]
        x_new[i] = x1[i-len(x)]
        y_new[i] = y1[i-len(x)]
        z_new[i] = z1[i-len(x)]
        M_new[i] = M1[i-len(x)]
        R_new[i] = R1[i-len(x)]
    
    for i in range(N1, N2):
        id_new[i] = ids2[i - N1]
        x_new[i] = x2[i- N1]
        y_new[i] = y2[i- N1]
        z_new[i] = z2[i- N1]
        M_new[i] = M2[i- N1]
        R_new[i] = R2[i- N1]
        
    for i in range(N2, N3):
        id_new[i] = ids3[i-N2]
        x_new[i] = x3[i-N2]
        y_new[i] = y3[i-N2]
        z_new[i] = z3[i-N2]
        M_new[i] = M3[i-N2]
        R_new[i] = R3[i-N2]
    
    for i in range(N3, N4):
        id_new[i] = ids4[i-N3]
        x_new[i] = x4[i-N3]
        y_new[i] = y4[i-N3]
        z_new[i] = z4[i-N3]
        M_new[i] = M4[i-N3]
        R_new[i] = R4[i-N3]
        
    for i in range(N4, N5):
        id_new[i] = ids5[i-N4]
        x_new[i] = x5[i-N4]
        y_new[i] = y5[i-N4]
        z_new[i] = z5[i-N4]
        M_new[i] = M5[i-N4]
        R_new[i] = R5[i-N4]
        
    for i in range(N5, N6):
        id_new[i] = ids6[i-N5]
        x_new[i] = x6[i-N5]
        y_new[i] = y6[i-N5]
        z_new[i] = z6[i-N5]
        M_new[i] = M6[i-N5]
        R_new[i] = R6[i-N5]
        
    return id_new, x_new, y_new, z_new, M_new, R_new

ids, x, y, z, M, R = bpc(Ids, X, Y, Z, m, r, P)

for i in range(len(x)):
	print ids[i], x[i], y[i], z[i], M[i], R[i]
