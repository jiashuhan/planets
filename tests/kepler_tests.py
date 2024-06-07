import sys, os
current = os.path.dirname(os.path.realpath(__name__)) # current directory
parent = os.path.dirname(current)
sys.path.append(parent)

import numpy as np
from kepler import *

Msun = 1.9884e30 # [kg]

#-------------------
# Eccentric anomaly
#-------------------

es = np.linspace(0, 0.99498, 20)
Ms = np.linspace(0.001, 359, 100)

for i, e in enumerate(es):
    for j, M in enumerate(Ms):
        E = ecc_anomaly(e, M, tolerance=0.002)
        M_ = mean_anomaly(e, E)
        diff = np.abs(M_/(M + np.finfo(float).eps) - 1)
        assert diff < 0.001

#-----------------------------
# Orbital elements conversion
#-----------------------------

# kep2cart
test_cases1 = [[0.2,  1e10, 10.,   48.,  29.,  174., 1e23], # Mercury-like
               [0.02, 1e11, 0.01, -11.,  100., 300., 1e25], # Earth-like
               [0.05, 1e12, 1.,    100., 300., 20.,  1e27], # Jupiter-like
               [0.3,  1e13, 20.,   120., 120., 20.,  1e22], # Pluto-like
               [0.97, 1e12, 160.,  60.,  110., 1.,   1e14]] # Halley-like

for n, param in enumerate(test_cases1):
    e, a, i, Om, wb, L, m = param 
    X, V = kep2cart(e, a, i, Om, m + Msun, ϖ=wb, L=L)
    e_, a_, i_, Om_, wb_, L_ = cart2kep(*X, *V, m + Msun)

    assert np.abs(e_/e - 1) < 0.001
    assert np.abs(a_/a - 1) < 0.001
    assert np.abs(i_/i - 1) < 0.001
    assert np.abs(Om_/((Om + 360)%360) - 1) < 0.001
    assert np.abs(wb_/wb - 1) < 0.001
    assert np.abs(L_/L - 1) < 0.001

# cart2kep
test_cases2 = [[-1e10, -1e9,   1e9,  -1e3, -1e5, -1e4, 1e23], # Mercury-like
               [ 5e10, -1e11, -1e7,   3e4,  2e4,  4e0, 1e25], # Earth-like
               [ 1e12, -4e11, -2e10, -5e3,  1e4,  5e1, 1e27], # Jupiter-like
               [ 1e13, -2e12, -3e12, -4e2,  3e3, -5e2, 1e22], # Pluto-like
               [-2e12, -2e11, -5e11,  4e3,  2e3,  8e2, 1e14]] # Halley-like

for n, param in enumerate(test_cases2):
    x, y, z, vx, vy, vz, m = param 
    e, a, i, Om, wb, L = cart2kep(x, y, z, vx, vy, vz, m + Msun)
    X, V = kep2cart(e, a, i, Om, m + Msun, ϖ=wb, L=L)

    assert np.abs(X[0]/x - 1) < 0.001
    assert np.abs(X[1]/y - 1) < 0.001
    assert np.abs(X[2]/z - 1) < 0.001
    assert np.abs(V[0]/vx - 1) < 0.001
    assert np.abs(V[1]/vy - 1) < 0.001
    assert np.abs(V[2]/vz- 1) < 0.001
