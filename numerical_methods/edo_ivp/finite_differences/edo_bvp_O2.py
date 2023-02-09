# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 08/02/23

import numpy as np
import matplotlib.pyplot as plt

"""
Este programa resuelve ED de orden 2 con valores en la frontera del tipo

        y''(x) + P(x)y'(x) + Q(x)y(x) = R(x) ; y(xo)=yo  &   y(xn)=yn
"""

# Matrices
def matrix(xo, xn, yo, yn, N, P, Q, R):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N)) ; b = np.zeros((N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    b[0] = yo   ;  b[N-1] = yn
    for i in range(1, N-1):
        A[i][i]   = -4 + 2*Q(X[i])*dx**2
        A[i][i+1] = 2 + dx*P(X[i])
        A[i][i-1] = 2 - dx*P(X[i])
        b[i]      = 2*R(X[i])*dx**2
    return A, b, X

def solve_ecns(A, b):
    Y = np.linalg.solve(A, b)
    return Y

"""
Parametros de la ecuación diferenecial
"""

# Parametro P(x) de la ecuación diferencial
def P(x):
    return 2

# Parametro Q(x) de la ecuación diferencial
def Q(x):
    return 1

# Parametro R(x) de la ecuación diferencial
def R(x):
    return 0

# Condiciones frontera
xo = 0 ; xn = 1
yo = 1 ; yn = 3

# Número de puntos que conforman la malla de [xo, xn]
N = 1001

"""
Llamamos el método para resolver la edo
"""

A, b, X = matrix(xo, xn, yo, yn, N, P, Q, R)
Y = solve_ecns(A, b)

plt.plot(X, Y)
plt.title(r"$y(x)$")
plt.show()