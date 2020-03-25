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

    firstTime = params['firstTime']
    T = params['T']
    nt = params['nt']
    deltaX = params['deltaX']
    betha2 = params['betha2']
    xc = params['xc']
    c2 = params['c2']

    if firstTime :
        params['firstTime'] = False
        deltaT = T / nt; params['deltaT'] = deltaT
        nx = 1 / deltaX; params['nx'] = nx # Definicao de deltaX
        c = math.sqrt(c2); params['c'] = c

    ti = i * params['deltaT']
    xj = j * params['deltaX']

    #if np.absolute(xj) < EPS or np.absolute(1 - xj) < EPS:
    #    return 0
    matrix = params['matrix']

    if j == 0 or j == params['nx'] or i <= 1 :
        if matrix[i][j] == 'a':
            matrix[i][j] = 0
            params['matrix'] = matrix
        return 0, matrix

    alpha = c * deltaT/deltaX

    term1  = -EDO(i-2, j, **params)[0]
    term2  = 2 * (1 - alpha**2) * EDO(i-1, j, **params)[0]
    term3A = EDO(i-1, j+1, **params)[0]
    term3B = EDO(i-1, j-1, **params)[0]
    term3  = (alpha**2) * (term3A + term3B)
    term4  = (deltaT**2) * funcaoF(ti, xj, betha2, xc, c2)

    print("term4 = {:.2f}".format(term4))
    sum = term1 + term2 + term3 + term4
    if matrix[i][j] == 'a':
        matrix[i][j] = sum
    return sum, matrix

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

#import pdb; pdb.set_trace()

M = m.criarMatriz(10, 10, fill='a')


print(EDO(nt=nt, betha2=betha2, deltaX=deltaX, c2=c2, xc=xc, T=T, i=2, j=3, matrix=M, firstTime=True))
