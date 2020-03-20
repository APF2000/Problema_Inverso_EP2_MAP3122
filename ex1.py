import numpy as np
import Matrizes as m
import matplotlib as mp
import math

EPS = 10e-10

def funcaoF(t, x, betha2, xc, c2):

    pi2 = (np.pi) ** 2

    if np.absolute(x - xc) > EPS:
        return 0

    aux1 = 1000 * c2
    aux2 = 1 - 2*betha2*pi2*(t ** 2)
    aux3 = np.exp(- betha2 * pi2 * (t ** 2))

    return aux1 * aux2 * aux3

# Solucao de u(ti, xj) para os parâmetros dados
def EDO(i, j, **params):
    # params = (x, t, T, nt, deltaX, deltaT, betha2, ut0, ut1)

    #import pdb; pdb.set_trace()

    T = params['T']
    nt = params['nt']
    deltaT = T / nt

    deltaX = params['deltaX']
    nx = 1 / deltaX # Definicao de deltaX

    betha2 = params['betha2']
    xc = params['xc']
    c2 = params['c2']
    c = math.sqrt(c2)

    ti = i * deltaT
    xj = j * deltaX

    #if np.absolute(xj) < EPS or np.absolute(1 - xj) < EPS:
    #    return 0

    if j == 0 or j == nx or i <= 1 :
        return 0

    alpha = c * deltaT/deltaX

    term1  = -EDO(i-2, j, **params)
    term2  = 2 * (1 - alpha**2) * EDO(i-1, j, **params)
    term3A = EDO(i-1, j+1, **params)
    term3B = EDO(i-1, j-1, **params)
    term3  = (alpha**2) * (term3A + term3B)
    term4  = (deltaT**2) * funcaoF(ti, xj, betha2, xc, c2)

    return term1 + term2 + term3 + term4

    # valores = m.criarMatriz(nt, nx)
#
    # Se j = 0 e i = 1 ou 0, entao uij = 0
# valores[0][0] = 0 # u(i=0, j=0) = u(t0, x0)
    # valores[1][0] = 0 # u(i=1, j=0) = u(t1, x0)
    # valores[0][1] = "naosei"# u(i=0, j=1) = u(t0, x1)
    # valores[1][1] = "tambemnao"# u(i=1, j=1) = u(t1, x1)
#
    # for it1 in range(nt):
        # if(it1 > 1):
            #
betha2 = 10 ** 2
xc = 0.7
c2 = 10
T = 1
nt = 350 #inicial
deltaX = 0.01


print(EDO(nt=nt, betha2=betha2, deltaX=deltaX, c2=c2, xc=xc, T=T, i=2, j=3))
