import ex1
import ex2
import Matrizes as m
import numpy as np
import matplotlib.pyplot as plt

# Função que calcula o dr com ruído.
# Precisa que o dr carregue o arquivo certo
# no arquivo do exercicio 2
def drRuido(ni):
    dr10 = np.load("dr10.npy")
    ruido = np.array([])
    v = (lambda t : (2 * np.random.rand() - 1) )
    maximo = max(dr10)

    delta = (lambda t, ni=ni: ni * maximo * v(t))

    for i in range(len(dr10)):
        ruido = np.append(ruido, (1 + delta(i / nt)) * dr10[i])

    return ruido

def plotSismogramas(comRuido, semRuido, t):
    plt.subplot(1, 2, 2)
    plt.title("Sismograma com ruído")
    plt.ylabel("dr")
    plt.xlabel("t")
    plt.plot(t, comRuido)

    plt.subplot(1, 2, 1)
    plt.title("Sismograma sem ruído")
    plt.ylabel("dr")
    plt.xlabel("t")
    plt.plot(t, semRuido)
    plt.show()

def nivelRuido(ni, ti, tf):
    um = (lambda t : 1)
    ruido = drRuido(ni)
    dr10 = np.load("dr10.npy")

    numerador = (lambda t : abs(dr10[int(t * nt)] -  ruido[int(t * nt)]))
    denominador = (lambda t : dr10[int(t * nt)])

    num = ex2.integral(numerador, um, ti, tf, 1 / nt)
    den = ex2.integral(denominador, um, ti, tf, 1 / nt)

    return 100 * num / den

def drRuidoLambda(ni):
    ruido = drRuido(ni)
    return (lambda t : ruido[int(t * nt)])

nt = 1000
nx = 100
deltaT = 1 / nt
ni = 1e-4
K = 10
ti, tf = 0.9, 1.0

print("foi0")
ruido = (drRuido(ni))
print("foi1")
t = [ i / nt for i in range(nt + 1) ]
print("foi2")
semRuido = np.load("dr10.npy")
print("foi3")

# plotSismogramas(ruido, semRuido, t)
#print(nivelRuido(1e-3, 0.5, 1))


ruido = drRuidoLambda(ni)

xcs = [.03, .15, .17, .25, .33, .34, .40, .44, .51, .73]
funcs = ex2.createFuncs(xcs, K)

# A * x = b
A, b = ex2.linearSystem(K, funcs, ruido, deltaT, ti, tf)

print(A)
print(" * ")
x = m.cholesky(A, b)
print(x)
print(" == ")
print(b)
print(m.testarRespostaSistemaLinear(A, b, len(xcs), x))
