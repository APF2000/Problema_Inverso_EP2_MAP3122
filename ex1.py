import numpy as np
import Matrizes as m
import matplotlib.pyplot as plt
import math
import os

EPS = 10e-10

# Função da fonte
# Recebe os indices i e j onde se deseja calcular f
def funcaoF(i, j, c2, xc, nx, nt):

    t = i / nt

    betha2 = 10 ** 2
    pi2 = (np.pi) ** 2

    if j != int(xc * nx):
        return 0

    aux1 = 1000 * c2
    aux2 = 1 - 2 * betha2 * pi2 * (t ** 2)
    aux3 = np.exp(- betha2 * pi2 * (t ** 2))

    return aux1 * aux2 * aux3

# Completa a matriz que contém os termos u(xi, tj)
def fillEDOmatrix(nt, nx, alpha, T, matrix, c2, xc):
    deltaT = T / nt
    deltaX = 1 / nx

    # Pelas condições iniciais, temos que :
    # u(0, xj) = u(1, xj) = u(ti, 0) = u(ti, nx) = 0
    for line in range(2, nt+1):
        for col in range(1, nx):
            term1 = -matrix[line-2][col]
            term2  = 2 * (1 - alpha**2) * matrix[line-1][col]
            term3A = matrix[line-1][col+1]
            term3B = matrix[line-1][col-1]
            term3  = (alpha**2) * (term3A + term3B)
            term4  = (deltaT **2) * funcaoF(line - 1, col, c2, xc, nx, nt)

            matrix[line][col] = term1 + term2 + term3 + term4
    return matrix

# Soluciona u(ti, xj) para os parâmetros dados
def EDO(i, j, matrix, nt, nx, c2, T, firstTime, xc):

    deltaT = T / nt
    deltaX = 1 / nx
    alpha = math.sqrt(c2) * (deltaT / deltaX)

    if firstTime:
        matrix = m.criarMatriz(nt + 1, nx + 1)
        matrix = fillEDOmatrix(nt, nx, alpha, T, matrix, c2, xc)

    return matrix[i][j], matrix

# Plota o gráfico de u(x, t) para um t dado
# e depende dos parâmetros nt, nx e c^2 (únicos que podem variar)
def plotArt(t):
    valoresx = np.array([])
    valoresy = np.array([])
    i = int(t*nt)

    matrix = []

    firstTime = True
    for j in range (nx):
        valoresx = np.append(valoresx, j*deltaX)

        aux = EDO(i, j, matrix, nt, nx, c2, T, firstTime, xc=0.7)
        valoresy = np.append(valoresy, aux[0])
        matrix = aux[1]

        firstTime = False

    plt.plot(valoresx, valoresy)
    nome = "Ex1 : c^2 = " + str(c2) + \
                    ", nt = " + str(nt) + \
                    ", nx = " + str(nx) + \
                    ", t = " + str(t)
    plt.title(nome,fontweight="bold")
    plt.grid(True)
    plt.xlabel("x")
    plt.ylabel("u(x, " + str(t) + ")")

    plt.savefig("Imagens/" + nome + ".jpeg")
    plt.clf()


# Execução principal deste arquivo
# Cria uma pasta para armazenar as imagens dos gráficos
if __name__ == "__main__":
    c2, nx, nt = 10, 100, 350
    T = 1
    deltaT = 1 / nt
    deltaX = 1 / nx

    try:
        os.mkdir("Imagens")
    except OSError:
        print("Não foi possível criar a pasta Imagens")
    else:
        print("As imagens do exercício 1 estão na pasta Imagens")

    while(nt > 300):
        plotArt(0.5)
        nt -= 10
    nt = 318
    while(nt > 314):
        plotArt(0.5)
        nt -= 1

    c2, nt = 20, 500
    while(nt > 430):
        plotArt(0.5)
        nt -= 10
    nt = 450
    while(nt > 440):
        plotArt(0.5)
        nt -= 1

    c2, nx, nt = 20, 200, 2500
    for i in range(6):
        plotArt((i + 1) * 0.1)
