import Matrizes as m
import math

def integral(u1, u2, s0, sN, deltaT):
    n = int( (sN - s0) / deltaT )
    sum = 0
    for i in range(1, n):
        si = s0 + i * deltaT
        sum += (u1(si) * u2(si))

    sum *= 2
    sum += (u1(s0) * u2(s0)) + (u1(sN) * u2(sN))

    return sum * deltaT / 2

def y1(x):
    return 1

def y2(x):
    return math.exp(x)

print(integral(y1, y2, math.log(1.2345), math.log(2), 0.0000001))
