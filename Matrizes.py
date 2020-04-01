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
# (ch * chT) * x = b
#
# ch * y  = b
# chT * x = y
def cholesky(A, b):
    ch = matrixCholesky(A)
    chT = transpostaMatriz(ch)

    # ch * y = b
    return ch

A = criarMatriz(3, 3, False)
A[0] = [   4,  12, -16] #    | 2, 0, 0 |   | 2, 6, -8 |
A[1] = [  12,  37, -43] # == | 6, 1, 0 | * | 0, 1,  5 |
A[2] = [ -16, -43,  98] #    |-8, 5, 3 |   | 0, 0,  3 |

b = criarMatriz(3, 1, False)
b = [ 1, 2, 3 ]

ch = cholesky(A, b)
chT = transpostaMatriz(ch)

print(multiplicacaoMatrizes(ch, chT))
print()
print(multiplicacaoMatrizes(chT, ch))
