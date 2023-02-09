# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 08/02/23

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

"""
Este programa resuelve ED de orden 2 con valores en la frontera del tipo

        Y(x) = Psi(x);

        Q(x) = (2m/hbar^2)(E-V(x));

        Y''(x) + Q(x)y(x) = 0 ; y(xo)=yo  &   y(xn)=yn
"""

# Constantes

#Constante de planck h barra
hbar = 1.05457e-34   # J*s
# masa del electron
m = 9.10939e-31     # kg

# Matrices
def matrix_to_solve(xo, xn, yo, yn, N, P, Q, R, E):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N)) ; b = np.zeros((N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    b[0] = yo   ;  b[N-1] = yn
    for i in range(1, N-1):
        A[i][i]   = -4 + 2*Q(X[i], E)*dx**2
        A[i][i+1] = 2 + dx*P(X[i], E)
        A[i][i-1] = 2 - dx*P(X[i], E)
        b[i]      = 2*R(X[i])*dx**2
    return A, b, X

# Matrices
def matrix_eig(xo, xn, N, Q):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    for i in range(1, N-1):
        A[i][i]   = (-2 + Q(X[i], 0)*dx**2)
        A[i][i+1] = 1
        A[i][i-1] = 1
    return A

# Resolución de la matriz
def solve_ecns(A, b):
    Y = np.linalg.solve(A, b)
    return Y

"""
Parametros de la ecuación diferenecial
"""

# Parametro P(x) de la ecuación diferencial
def P(x):
    return 0

# Parametro Q(x) de la ecuación diferencial
def Q(x, E):
    hbar = 1.05457e-34   # J*s
    # masa del electron
    m = 9.10939e-31     # kg
    V = 0
    return (2*m/hbar**2)*(E-V)

# Parametro R(x) de la ecuación diferencial
def R(x):
    return 0

# Condiciones frontera
xo = 0 ; xn = 1
yo = 0 ; yn = 0

# Número de puntos que conforman la malla de [xo, xn]
N = 1001

dx = (xn-xo)/(N-1)

"""
Llamamos el método para resolver la edo
"""
A = matrix_eig(xo, xn, N, Q)
print(A)

energies = -((hbar**2)/(2*m*dx**2))*linalg.eigvals(A)
exact_energies = (((np.pi*hbar)**2)/(2*m*xn**2))*(np.arange(1,10,1))**2

print(energies)
print(exact_energies)

#A, b, X = matrix_to_solve(xo, xn, yo, yn, N, P, Q, R)
#Y = solve_ecns(A, b)

#plt.plot(X, Y)
#plt.title(r"$y(x)$")
#plt.show()