import numpy as np
import math

def criarMatriz(num_linhas,num_colunas):
    dimensoes = [num_linhas,num_colunas]
    matriz = np.zeros(dimensoes)
    return matriz


def somaMatrizes(matrizA,matrizB):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    linB = matrizB.shape[0]
    colB = matrizB.shape[1]
    if(linA != linB or colA != colB):
        return
    else:
        soma = criarMatriz(linA,colA)
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
        C = criarMatriz(linA, colB)
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
    transposta = criarMatriz(colA, linA)
    for i in range(colA):
            for j in range (linA):
                transposta[i][j] = matrizA[j][i]
    return transposta

def multiplicacaoRealMatriz(matrizA,alfa):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    matriz = criarMatriz(colA, linA)
    for i in range(linA):
            for j in range (colA):
                matriz[i][j] = alfa*matrizA[j][i]
    return matriz

def matrixCholesky(A):
    n = len(A)
    # Cria matriz de zeros
    ch = criarMatriz(n, n)

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
    x = criarMatriz(n, 1)

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
    y = sistemaLouU(ch, b, isL=True)

    # chT * x = y
    x = sistemaLouU(chT, y, isL=False)


    return x

def metodoSOR (matrizA,matrizb,n):
    numIteracoes = 1000
    omega = 1.6
    iteracoes = 0
    respostas = criarMatriz(n, 1)
    print(respostas)

    while(iteracoes < 10):
        for i in range(n):
            somatorias = 0
            for j in range(i):
                somatorias = matrizA[i][j]*respostas[j][0] + somatorias
            for j in range(i+1,n):
                somatorias = matrizA[i][j]*respostas[j][0] + somatorias
            #print(somatorias)
            respostas[i][0] = (1/matrizA[i][i])*(matrizb[i][0]-somatorias)
        iteracoes += 1

    iteracoes = 0

    while(iteracoes < 1000):
        for i in range(n):
            somatorias = 0
            for j in range(i):
                somatorias = matrizA[i][j]*respostas[j][0] + somatorias
            for j in range(i+1,n):
                somatorias = matrizA[i][j]*respostas[j][0] + somatorias
            #print(somatorias)
            respostas[i][0] = (1-omega)*respostas[i][0] + (omega/matrizA[i][i])*(matrizb[i][0]-somatorias)
        iteracoes += 1
    return respostas

def testarRespostaSistemaLinear (matrizA,matrizb,n,respostas):
    precisao = 0.01
    for i in range(n):
        soma = 0
        for j in range(n):
            soma = matrizA[i][j]*respostas[j][0] + soma
        if(np.absolute(soma - matrizb[i][0]) > precisao):
            return False
    return True
