import ex1
import ex2
import Matrizes as m
import numpy as np
import matplotlib.pyplot as plt
import os

# Função que calcula o dr com ruído.
# Precisa que o dr carregue o arquivo dr10
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

    plt.savefig("ImagensEx3/ni={:.0e}.jpeg".format(ni))

    plt.clf()

def nivelRuido(ni, ti, tf):
    um = (lambda t : 1)
    ruido = drRuido(ni)
    dr10 = np.load("dr10.npy")

    numerador = (lambda t : abs(dr10[int(t * nt)] -  ruido[int(t * nt)]))
    denominador = (lambda t : dr10[int(t * nt)])

    num = ex2.integral(numerador, um, ti, tf, T / nt)
    den = ex2.integral(denominador, um, ti, tf, T / nt)

    return 100 * num / den

def drRuidoLambda(ni):
    ruido = drRuido(ni)
    return (lambda t : ruido[int(t * nt)])

if __name__ == "__main__":
    try:
        os.mkdir("ImagensEx3")
    except OSError:
        print("A pasta ImagensEx3 ou já existe ou não pôde ser criada")
    else:
        print("As imagens do exercício 3 estão na pasta ImagensEx3")

    nt = 1000
    nx = 100
    T = 1
    deltaT = T / nt
    K = 10
    ti, tf = 0.9, 1.0

    for j in range(3):
        ni = 10 ** (- j - 3)
        print("Para ni = ", ni)

        ruido = (drRuido(ni))
        t = [ i / nt for i in range(nt + 1) ]
        semRuido = np.load("dr10.npy")

        plotSismogramas(ruido, semRuido, t)
        print("\nO nível de ruído foi")
        print(nivelRuido(1e-3, 0.5, 1))

        ruido = drRuidoLambda(ni)

        xcs = [.03, .15, .17, .25, .33, .34, .40, .44, .51, .73]
        esperado = [7.3, 2.4, 5.7, 4.7, 0.1, 20.0, 5.1, 6.1, 2.8, 15.3]
        funcs = ex2.createFuncs(xcs, K)

        # B * a = c; método Cholesky
        B, c = ex2.linearSystem(K, funcs, ruido, deltaT, ti, tf)
        ex2.presentSystem(B, c, K)

        print("\nUsando o método Cholesky, a resposta foi:")
        a = m.cholesky(B, c)
        ex2.presentResult(a, K)
        obtido = a

        print("\nO erro de reconstrução do Cholesky foi ")
        print(ex2.erroReconstrucao(esperado, obtido, K))

        print("\nO resíduo do Cholesky foi ")
        residuo = (ex2.residuoReconstrucao(ti, tf, obtido, ruido, funcs, K))
        print(residuo)

        print("\n-------------------------")

        # B * a = c; método SOR
        print("\nUsando o método SOR, a resposta foi:")
        a = m.metodoSOR(B, c, len(xcs))
        ex2.presentResult(a, K)

        print("\nO resíduo do SOR foi ")
        residuo = (ex2.residuoReconstrucao(ti, tf, a, ruido, funcs, K))
        print(residuo)

        print("\nO erro de reconstrução do SOR foi ")
        print(ex2.erroReconstrucao(esperado, a, K))
