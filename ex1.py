import numpy as np
import Matrizes as m
import matplotlib as mp
import math

EPS = 10e-10

def funcaoF(t, x, betha2, xc, c2):

    pi2 = (np.pi) ** 2

    if np.absolute(x - xc) < EPS:
        return 0

    aux1 = 1000 * c2
    aux2 = 1 - 2*betha2*pi2*(t ** 2)
    aux3 = np.exp(- betha2 * pi2 * (t ** 2))

    return aux1 * aux2 * aux3

# Solucao de u(ti, xj) para os parâmetros dados
def EDO(x, t, T, nt, deltaX, deltaT, betha2, ut0, ut1, i, j):
    if np.absolute(x) < EPS or np.absolute(1 - x) < EPS:
        return 0

    valores = m.criarMatriz(num_linhas, num_colunas)
    valores

    nx = 1 / deltaX # Definicao de deltaX

    for it1 in range(nt):
        for it2 in range(nx):
            pass


betha2 = 10 ** 2
xc = 0.7
c2 = 10
T = 1
nt = 350 #inicial
deltaX = 0.01
