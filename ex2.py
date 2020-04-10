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

# Apresenta o sistema linear em forma matricial
def presentSystem(B, c, K):
    for i in range(K):
        line = ""
        bLine = " |"
        equals, times = " == ", " * "
        space3, space4 = "   ", "    "
        for j in range(K):
            bLine += " {:10.3e}".format(B[i][j])

        bLine += " | "
        aLine = " | a" + str(i) + " |"

        line += bLine

        if( i == int(K / 2) or i == int(K / 2 - 1) and K % 2 == 0):
            line += times
        else:
            line += space3
        line += aLine

        if( i == int(K / 2) or i == int(K / 2 - 1) and K % 2 == 0):
            line += equals
        else:
            line += space4

        cLine = "|" + " {:10.3e}".format(c[i][0]) + " |"
        line += cLine

        print(line)

# Apresenta a solução do sistema na forma matricial
def presentResult(a, K):
    for i in range(K):
        line = " | {:10.5f} |".format(a[i][0])
        print(line)


# Cria uma lista de funcoes a partir de um dos parametros
# do vetor xc fornecido, com K posições
def createFuncs(xcs, K):
    porFora = []
    funcs = []
    for i in range(K):
        porFora.append(ukx(xr=xr, xc=xcs[i]))
        funcs.append(lambda t, k = i: porFora[k][int(t * nt)])
    return funcs

if __name__ == "__main__":
    K = int(input("Para qual valor de K? "))

    if K == 3:
        c2, nt, nx = 20, 1000, 100
        ti, tf = 0.5, 1.0
    elif K == 10:
        c2, nt, nx = 20, 1000, 100
        ti, tf = 0.9, 1.0
    elif K == 20:
        c2, nt, nx = 20, 2000, 200
        ti, tf = 0.9, 1.0
    else:
        print("Este valor de K é inválido")
        import sys
        sys.exit()

    deltaT = 1 / nt
    deltaX = 1 / nx
    T = 1
    xr = 0.7

    xcs = {
        3 : [0.2, 0.3, 0.9],
        10 : [.03, .15, .17, .25, .33, .34, .40, 0.44, .51, .73],
        20 : [ 0.1 + 0.025 * (k) for k in range(20) ]
    }
    xcs = xcs[K]

    funcs = createFuncs(xcs, len(xcs))

    def readNPY(name):
        return np.load(name)

    dr = readNPY("dr" + str(K) + ".npy")
    def funDr(t):
        pos = int(t * nt)
        return dr[pos]

    # B * a = c
    B, c = linearSystem(K, funcs, funDr, deltaT, ti, tf)
    presentSystem(B, c, K)

    print("\nUsando o método Cholesky, a resposta foi:\n")
    a = m.cholesky(B, c)
    presentResult(a, K)

    if K != 20:
        print("\n-------------------------\n")
        print("Usando o método SOR, a resposta foi:")
        a = m.metodoSOR(B, c, len(xcs))
        presentResult(a, K)
