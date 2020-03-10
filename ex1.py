import numpy as np
import Matrizes as m
import matplotlib as mp
import math

EPS = 10e-10

def funcaoF(t, x):
    betha2 = 10 ** 2
    xc = 0.7
    c2 = 10
    pi2 = (np.pi) ** 2

    if np.absolute(x - xc) < EPS:
        return 0

    aux1 = 1000 * c2
    aux2 = 1 - 2*betha2*pi2*(t ** 2)
    aux3 = np.exp(- betha2 * pi2 * (t ** 2))

    return aux1 * aux2 * aux3

def EDO(u, f):
