import Matrizes as m
import math
import numpy as np
import ex1

def integral(u1, u2, s0, sN, deltaT):
    n = int( (sN - s0) / deltaT )
    sum = 0
    #import pdb; pdb.set_trace()
    for i in range(1, n):
        si = s0 + i * deltaT
        sum += (u1(si) * u2(si))

    sum *= 2
    sum += (u1(s0) * u2(s0)) + (u1(sN) * u2(sN))

    return sum * deltaT / 2

def linearSystem(K, funcs, dr, deltaT):
    B = m.criarMatriz(K, K)
    c = m.criarMatriz(K, 1)

    for i in range(K):
        for j in range(i):
            B[i][j] = integral(funcs[i], funcs[j], 0, 1, deltaT)
            B[j][i] = B[i][j]
        B[i][i] = integral(funcs[i], funcs[i], 0, 1, deltaT)

    for i in range(K):
        c[i][0] = integral(funcs[i], dr, 0, 1, deltaT)

    return B, c

def readNPY(name):
    return np.load(name)

nt = 350
deltaT = 1 / nt
deltaX = 0.01
nx = int(1 / deltaX)

c2 = 10
T = 1

dr = readNPY("dr3.npy")
def funDr(s):
    pos = int(s * nt)
    return dr[pos]

def ukx(xr, xc):
    j = int( xr / deltaX )

    x = np.array([])
    matrix = [] # default, sem significado

    firstTime = True
    for i in range(nt + 1):
        aux = ex1.EDO(i, j, matrix, nt, nx, c2, T, firstTime, xc=xc)
        x = np.append(x, aux[0])
        matrix = aux[1]
        firstTime = False
    return x

#a = ukx(0.7)
#print(type(a))

# def ukt(t):
    # xr = 0.7
    # pos = int(t * nt)
    # return ukx(xr, xc)[pos]


# Cria uma lista de funcoes a partir de um dos parametros
# do vetor xc fornecido, com K posições
def createFuncs(xcs, K):
    xr = 0.7
    funcs = []
    for i in range(K):
        aux = ukx(xr, xcs[i])
        funcs.append( lambda t: aux[int(t * nt)] )

    return aux


xcs = [0.7]
aux = createFuncs(xcs, len(xcs))
# aux = []
#for i in range(nt):
#    aux.append(funcs[0](i * deltaT))
#print(aux)

matrix = []
matrix = (ex1.EDO(100, 70, matrix, nt, nx, c2, T, True, xc=0.7)[1])
#print(matrix[100])

aux2 = []
for i in range(nt+1):
    aux2.append(matrix[i][70])
#print(aux2)

EPS = 0.00001
cont = 0
for i in range(nx):
    if (np.absolute(aux2[i] - aux[i]) >= EPS):
        cont += 1
        print(aux2[i], aux[i])

print(cont)
