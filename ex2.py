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

def linearSystem(K, funcs, dr, deltaT, ti, tf):
    B = m.criarMatriz(K, K)
    c = m.criarMatriz(K, 1)

    for i in range(K):
        for j in range(i):
            B[i][j] = integral(funcs[i], funcs[j], ti, tf, deltaT)
            B[j][i] = B[i][j]
        B[i][i] = integral(funcs[i], funcs[i], ti, tf, deltaT)

    for i in range(K):
        c[i][0] = integral(funcs[i], dr, ti, tf, deltaT)

    return B, c

def readNPY(name):
    return np.load(name)

nt = 1000
deltaT = 1 / nt
deltaX = 0.01
nx = int(1 / deltaX)

c2 = 20
T = 1

dr = readNPY("dr10.npy")
def funDr(s):
    pos = int(s * nt)
    return dr[pos]
#dr = (lambda )


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


# Cria uma lista de funcoes a partir de um dos parametros
# do vetor xc fornecido, com K posições
def createFuncs(xcs, K):
    porFora = []
    xr = 0.7
    funcs = []
    for i in range(K):
        porFora.append(ukx(xr=xr, xc=xcs[i]))
        funcs.append( lambda t, k = i: porFora[k][int(t * nt)] )

    return funcs


xcs = [0.03, 0.15, 0.17, 0.25, 0.33, 0.34, 0.40, 0.44, 0.51, 0.73]
funcs = createFuncs(xcs, len(xcs))

#dr = readNPY("dr3.npy")
# A * x = b
A, b = linearSystem(len(xcs), funcs, funDr, deltaT, 0.9, 1)
print(A)
print(" * ")
x = m.cholesky(A, b)
print(x)
print(" == ")
print(b)

print(m.testarRespostaSistemaLinear(A, b, 10, x))
print("\n-------------------------\n")
x = m.metodoSOR(A, b, 10)
print(x)
print(m.testarRespostaSistemaLinear(A, b, 10, x))
