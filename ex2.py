import Matrizes as m
import math
import numpy as np
import ex1

def integral(u1, u2, s0, sN, deltaT):
    n = int( (sN - s0) / deltaT )
    sum = 0
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


dr = readNPY("dr3.npy")
def funDr(s):
    nt = 1000
    pos = int(s * nt)

    return dr[pos]

def createFuncs():
    return 0


nt = 100
deltaT = 1 / nt
deltaX = 0.01
nx = int(1 / deltaX)

c2 = 20
T = 1

def ukx(xr):
    j = int( xr / deltaX )

    x = np.array([])
    matrix = [] # default, sem significado

    firstTime = True
    for i in range(nt + 1):
        aux = ex1.EDO(i, j, matrix, nt, nx, c2, T, firstTime)
        x = np.append(x, aux[0])
        matrix = aux[1]
        firstTime = False
    return x

a = ukx(0.7)
print(type(a))



def ukt(t):
    nt = 1000
    pos = int(t * nt)
    return uk(x)[pos]
