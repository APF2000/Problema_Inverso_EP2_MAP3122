import numpy as np
import Matrizes as m
import matplotlib.pyplot as plt
import math

EPS = 10e-10

def funcaoF(t, x, c2, xc):

    betha2 = 10 ** 2

    pi2 = (np.pi) ** 2

    if np.absolute(x - xc) > EPS:
        return 0

    aux1 = 1000 * c2
    aux2 = 1 - 2*betha2*pi2*(t ** 2)
    aux3 = np.exp(- betha2 * pi2 * (t ** 2))

    return aux1 * aux2 * aux3

def criarMatrizBordas(nt, nx):
    matrix = m.criarMatriz(nt + 1, nx+1)

    for i in range(nx + 1):
        matrix[0][i] = 0
        matrix[1][i] = 0

    for i in range(nt + 1):
        matrix[i][0] = 0
        matrix[i][nx] = 0

    return matrix

def fillEDOmatrix(nt, nx, alpha, T, matrix, c2, xc):
    deltaT = T / nt
    deltaX = 1 / nx
    for line in range(2, nt+1):
        for col in range(1, nx):
            ti = line * deltaT
            xj = col * deltaX

            term1 = -matrix[line-2][col]
            term2  = 2 * (1 - alpha**2) * matrix[line-1][col]
            term3A = matrix[line-1][col+1]
            term3B = matrix[line-1][col-1]
            term3  = (alpha**2) * (term3A + term3B)
            term4  = (deltaT **2) * funcaoF(ti, xj, c2, xc)

            matrix[line][col] = term1 + term2 + term3 + term4
    return matrix

# Solucao de u(ti, xj) para os parâmetros dados
def EDO(i, j, matrix, nt, nx, c2, T, firstTime, xc):

    #import pdb; pdb.set_trace()

    deltaT = T / nt
    deltaX = 1 / nx
    alpha = math.sqrt(c2) * deltaT / deltaX

    if firstTime:
        matrix = criarMatrizBordas(nt, nx)
        matrix = fillEDOmatrix(nt, nx, alpha, T, matrix, c2, xc)

    return matrix[i][j], matrix


c2 = 20
T = 1
nt = 1500 #inicial
deltaT = 1 / nt
deltaX = 0.005
nx = int(1 / deltaX)

def plotArt(t):
    valoresx = np.array([])
    valoresy = np.array([])
    i = int(t*nt)

    matrix = []

    firstTime = True
    for j in range (nx):
        valoresx = np.append(valoresx, j*deltaX)

        aux = EDO(i, j, matrix, nt, nx, c2, T, firstTime, xc=0.9)
        valoresy = np.append(valoresy, aux[0])
        matrix = aux[1]

        firstTime = False

    plt.plot(valoresx, valoresy)
    plt.ylim(-0.6, 0.6)
    plt.show()
#Parte principal
#plotArt(0.1)
#plotArt(0.2)
#plotArt(0.3)
#plotArt(0.4)
#plotArt(0.5)
# plotArt(0.6)
#print(EDO(nt=nt, deltaX=deltaX, c2=c2, T=T, i=70, j=50, firstTime=True)[0])
