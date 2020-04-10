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
        funcs.append(lambda t, k = i: porFora[k][int(t * nt)])
    return funcs

nt = 1000
deltaT = 1 / nt
nx = 100
deltaX = 1 / nx
c2 = 20
T = 1

if __name__ == "__main__":

    xcs = [0.2, 0.3, 0.9]
    #xcs = [.03, .15, .17, .25, .33, .34, .40, 0.44, .51, .73]
    #xcs = ([ 0.1 + 0.025 * (k) for k in range(20) ])
    #xcs = [0.15]

    # Dando problema nas linhas 3, 8, 9, 14, 15, 19, 20
    #xcs = [0.15, 0.275, 0.3, 0.425, 0.45, 0.55, 0.575]

    K = len(xcs)
    funcs = createFuncs(xcs, len(xcs))
    #for i in range(K):
        #print(ukx(0.7, 0.275))

    def readNPY(name):
        return np.load(name)

    dr = readNPY("dr" + str(K) + ".npy")
    def funDr(t):
        pos = int(t * nt)
        return dr[pos]

    # for i in range(K):
        # aux = []
        # for j in range(K):
            # t = 0.9 + j * deltaT
            # aux.append("{:.1f}".format(funcs[i](t)*1e5))
        # print(aux)



    # A * x = b
    A, b = linearSystem(len(xcs), funcs, funDr, deltaT, ti=0.9, tf=1)
    print(A)
    print(" * ")
    x = m.cholesky(A, b)
    print(x)
    print(" == ")
    print(b)
    print(m.testarRespostaSistemaLinear(A, b, len(xcs), x))


    print("\n-------------------------\n")

    x = m.metodoSOR(A, b, len(xcs))
    print(x)
    print(m.testarRespostaSistemaLinear(A, b, len(xcs), x))
