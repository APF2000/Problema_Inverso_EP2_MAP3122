import numpy as np
import math

def criarMatriz(num_linhas,num_colunas, fill):
    dimensoes = [num_linhas,num_colunas]
    matriz = np.zeros(dimensoes)
    if fill :
        matriz = []
        for i in range(num_linhas):
            matriz.append([])
            for j in range(num_colunas):
                matriz[i].append(fill)
    return matriz


def somaMatrizes(matrizA,matrizB):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    linB = matrizB.shape[0]
    colB = matrizB.shape[1]
    if(linA != linB or colA != colB):
        return
    else:
        soma = criarMatriz(linA,colA, False)
        for i in range(linA):
            for j in range (colA):
                soma[i][j] = matrizA[i][j] + matrizB[i][j]
    return soma


def multiplicacaoMatrizes(matrizA,matrizB): # C = A*B
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    linB = matrizB.shape[0]
    colB = matrizB.shape[1]
    if(colA != linB ):
        return
    else:
        p = linB
        C = criarMatriz(linA,colB, False)
        for i in range(linA):
            for j in range(colB):
                soma = 0
                for k in range(p):
                    soma = soma + matrizA[i][k]*matrizB[k][j]
                C[i][j] = soma
    return C

def transpostaMatriz(matrizA):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    transposta = criarMatriz(colA,linA, False)
    for i in range(colA):
            for j in range (linA):
                transposta[i][j] = matrizA[j][i]
    return transposta

def multiplicacaoRealMatriz(matrizA,alfa):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    matriz = criarMatriz(colA,linA, False)
    for i in range(linA):
            for j in range (colA):
                matriz[i][j] = alfa*matrizA[j][i]
    return matriz

def matrixCholesky(A):
    n = len(A)
    # Cria matriz de zeros
    ch = criarMatriz(n, n, False)

    # Primeiro termo
    ch[0][0] = math.sqrt(A[0][0])

    # Primeira coluna
    auxSum = 0
    for i in range(1, n):
        ch[i][0] = A[i][0] / ch[0][0]

    #import pdb; pdb.set_trace()
    for i in range(1, n):
        # Resto dos elementos
        for j in range(1, i):
            auxSum = 0
            for k in range(j):
                auxSum += (ch[i][k] * ch[j][k])
            ch[i][j] = (A[i][j] - auxSum) / ch[j][j]
        auxSum = 0
        for k in range(i):
            auxSum += ( ch[i][k] * ch[i][k] )

        # Diagonal
        ch[i][i] = math.sqrt(A[i][i] - auxSum)
    return ch

# A * x = b
# isL determina se a matriz A é :
#   - L : triangular inferior
#   - U : triangular superior
def sistemaLouU(A, b, isL):
    n = len(A)
    x = criarMatriz(n, 1, False)

    #import pdb; pdb.set_trace()

    for i in range(n):
        if isL :
            line = i
        else :
            line = (n - 1) - i

        aux = b[line][0]
        for j in range(i):
            if isL :
                col = j
            else :
                col = (n - 1) - j

            if line != col:
                aux -= ( A[line][col] * x[col][0] )

        if A[line][line] == 0:
            print("Sistema nao pode ter 0 na diagonal")
            return False

        x[line][0] = ( aux / A[line][line] )

    return x


# A * x = b
# (ch * chT) * x = b
#
# ch * y  = b
# chT * x = y
def cholesky(A, b):
    ch = matrixCholesky(A)
    chT = transpostaMatriz(ch)

    # ch * y = b
    #print(ch)
    y = sistemaLouU(ch, b, isL=True)
    #print("\nsolucao = ", y)

    # chT * x = y
    print(chT)
    x = sistemaLouU(chT, y, isL=False)
    print("\nsolucao = ", x)

    return x

A = criarMatriz(3, 3, False)
A[0] = [   4,  12, -16] #    | 2, 0, 0 |   | 2, 6, -8 |
A[1] = [  12,  37, -43] # == | 6, 1, 0 | * | 0, 1,  5 |
A[2] = [ -16, -43,  98] #    |-8, 5, 3 |   | 0, 0,  3 |

b = criarMatriz(3, 1, False)
b[0][0] = 1
b[1][0] = 2
b[2][0] = 3

x = cholesky(A, b)
print(A)
print("\n         * ")
print(x)
print("\n         ==")
print(multiplicacaoMatrizes(A, x))

# ch = cholesky(A, b)
# chT = transpostaMatriz(ch)
#
# print(multiplicacaoMatrizes(ch, chT))
# print()
# print(multiplicacaoMatrizes(chT, ch))
