
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


def criarMatriz(num_linhas,num_colunas):
    dimensoes = [num_linhas,num_colunas]
    matriz = np.ones(dimensoes)
    return matriz


# In[3]:


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


# In[4]:


def multiplicacaoMatrizes(matrizA,matrizB): # C = A*B
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    linB = matrizB.shape[0]
    colB = matrizB.shape[1]
    if(colA != linB ):
        return
    else:
        p = linB
        C = criarMatriz(linA,colB)
        for i in range(linA):
            for j in range(colB):
                soma = 0
                for k in range(p):
                    soma = soma + matrizA[i][k]*matrizB[k][j]
                C[i][j] = soma
    return C
        


# In[5]:


def transpostaMatriz(matrizA):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    transposta = criarMatriz(colA,linA)
    for i in range(colA):
            for j in range (linA):
                transposta[i][j] = matrizA[j][i]
    return transposta


# In[6]:


def multiplicacaoRealMatriz(matrizA,alfa):
    linA = matrizA.shape[0]
    colA = matrizA.shape[1]
    matriz = criarMatriz(colA,linA)
    for i in range(linA):
            for j in range (colA):
                matriz[i][j] = alfa*matrizA[j][i]
    return matriz


# In[7]:


A = criarMatriz(2,2)
B = criarMatriz(2,2)
B[0][1] = 0
B[1][0] = 0
print(B)
print(A)
C = multiplicacaoMatrizes(A,B)
print(C)
C = multiplicacaoRealMatriz(C,2)
print(C)
