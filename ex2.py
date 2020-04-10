import Matrizes as m
import math
import numpy as np
import ex1

# Calcula a integral do produto das funções u1 e u2
# Os limites da integral são s0 e sN, e o passo entre
# as áreas somadas é deltaT
# Aqui a aproximação é feita usando a Fórmula dos Trapézios Composta
def integral(u1, u2, s0, sN, deltaT):
    n = int( (sN - s0) / deltaT )
    sum = 0
    for i in range(1, n):
        si = s0 + i * deltaT
        sum += (u1(si) * u2(si))

    sum *= 2
    sum += (u1(s0) * u2(s0)) + (u1(sN) * u2(sN))

    return sum * deltaT / 2

# Cria um sistema linear a partir das integrais
# das funções Ui obtidas e o dr fornecido
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

# Retorna a coluna correspondente ao xr, obtida a partir
# da matriz gerada pela resolução da EDO para Ui(t, x)
# Esse vetor-coluna será utilizado para o cálculo
# das integrais (ver função createFuncs)
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

# Retorna o erro da solução encontrada
def erroReconstrucao(Ak_esperados, Ak_obtidos, K):
    sum = 0
    for i in range(K):
        sum += (Ak_obtidos[i] - Ak_esperados[i]) ** 2
    return math.sqrt(sum)

# Calcula a função custo para a solução encontrada
def residuoReconstrucao(ti, tf, Aks, dr, funcs, K):
    f = (lambda t : sum([Aks[j] * funcs[j](t) for j in range(K)]) - dr(t))
    result = integral(f, f, ti, tf, deltaT)
    return float( result / ( 2 * (tf - ti) ) )

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
def presentResult(vec, K):
    for i in range(K):
        line = " | {:10.5f} |".format(vec[i][0])
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

# Ordena um numpyndarray
# Fizemos esta função porque usar simplesmente
# o método sort() não estava funcionando
def numpyToList(vec, K):
    new =[]
    for i in range(K):
        new.append(vec[i])
    new.sort()
    return new

if __name__ == "__main__":
    K = int(input("Para qual valor de K? "))

    if K == 3:
        c2, nt, nx = 20, 1000, 100
        ti, tf = 0.5, 1.0
        esperado = [0.1, 40.0, 7.5]
    elif K == 10:
        c2, nt, nx = 20, 1000, 100
        ti, tf = 0.9, 1.0
        esperado = [7.3, 2.4, 5.7, 4.7, 0.1, 20.0, 5.1, 6.1, 2.8, 15.3]
    elif K == 20:
        c2, nt, nx = 20, 2000, 200
        ti, tf = 0.9, 1.0
    else:
        print("Este valor de K é inválido")
        import sys
        sys.exit()

    # Parâmetros para a EDO e o sistema linear
    deltaT = 1 / nt
    deltaX = 1 / nx
    T = 1
    xr = 0.7

    # Valores de xc (mudam a fonte) para cada uma das 3 situações requeridas
    xcs = {
        3 : [0.2, 0.3, 0.9],
        10 : [.03, .15, .17, .25, .33, .34, .40, 0.44, .51, .73],
        20 : [ 0.1 + 0.025 * (k) for k in range(20) ]
    }
    xcs = xcs[K]

    # Aqui temos um vetor com as funções uk(t, xr=0.7)
    funcs = createFuncs(xcs, len(xcs))

    # Lê o arquivo dr correspondente
    def readNPY(name):
        return np.load(name)

    # Definimos aqui a função dr
    dr = readNPY("dr" + str(K) + ".npy")
    def funDr(t):
        pos = int(t * nt)
        return dr[pos]

    # B * a = c; método Cholesky
    B, c = linearSystem(K, funcs, funDr, deltaT, ti, tf)
    presentSystem(B, c, K)

    print("\nUsando o método Cholesky, a resposta foi:")
    a = m.cholesky(B, c)
    presentResult(a, K)
    obtido = a

    print("\nO resíduo do Cholesky foi ")
    residuo = (residuoReconstrucao(ti, tf, obtido, funDr, funcs, K))
    print(residuo)

    if K != 20:
        print("\nO erro de reconstrução do Cholesky foi ")
        print(erroReconstrucao(esperado, obtido, K))

        print("\n-------------------------")
        # B * a = c; método SOR
        print("\nUsando o método SOR, a resposta foi:")
        a = m.metodoSOR(B, c, len(xcs))
        presentResult(a, K)

        print("\nO resíduo do SOR foi ")
        residuo = (residuoReconstrucao(ti, tf, a, funDr, funcs, K))
        print(residuo)

        print("\nO erro de reconstrução do SOR foi ")
        print(erroReconstrucao(esperado, a, K))
    else:
        print("\nAs raízes obtidas, em ordem crescente são")
        obtido = np.array(numpyToList(obtido, K))
        presentResult(obtido, K)
