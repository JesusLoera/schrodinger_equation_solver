# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 08/02/23

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from scipy import integrate

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

def V(x):
    if ( x>=0 and x<=1):
        return 0
    else:
        return 10e10

# Matrices
def matrix_to_solve(xo, xn, yo, yn, N, P, Q, R, E):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N)) ; b = np.zeros((N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    b[0] = yo   ;  b[N-1] = yn
    for i in range(1, N-1):
        A[i][i]   = -4 + 2*Q(X[i], E, V)*dx**2
        A[i][i+1] = 2 - dx*P(X[i])
        A[i][i-1] = 2 + dx*P(X[i])
        b[i]      = 2*R(X[i])*dx**2
    return A, b, X

# Matrices
def matrix_to_solve2(xo, xn, yo, yn, N, P, Q, R, E):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N)) ; b = np.zeros((N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    b[0] = yo   ;  b[N-1] = yn
    for i in range(1, N-1):
        A[i][i]   = -2 + Q(X[i], 6.024642468562235170e-38, V)*dx**2
        A[i][i+1] = 1
        A[i][i-1] = 1
        b[i]      = 0
    return A, b, X

# Matrices
def matrix_eig(xo, xn, N, Q):
    X = np.linspace(xo, xn, N)
    dx = X[1] - X[0]
    A = np.zeros((N, N))
    A[0][0] = 1 ;  A[N-1][N-1] = 1
    for i in range(1, N-1):
        A[i][i]   = (-2 + Q(X[i], 0, V)*dx**2)
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
def Q(x, E, V):
    hbar = 1.05457e-34   # J*s
    # masa del electron
    m = 9.10939e-31     # kg
    V = V(x)
    return (2*m/hbar**2)*(E-V)

# Parametro R(x) de la ecuación diferencial
def R(x):
    return 0

# Condiciones frontera
xo = 0 ; xn = 1
yo = -0.00000001 ; yn = -0.000000001

# Número de puntos que conforman la malla de [xo, xn]
N = 1001

dx = (xn-xo)/(N-1)

"""
Llamamos el método para resolver la edo
"""
A = matrix_eig(xo, xn, N, Q)

eigenvalues = -((hbar**2)/(2*m*dx**2))*linalg.eigvals(A)
energies = []
for eigenvalue in eigenvalues:
    if (np.imag(eigenvalue) == 0):
        energies.append(np.real(eigenvalue))
energies = np.sort(energies)
exact_energies = (((np.pi*hbar)**2)/(2*m*xn**2))*(np.arange(1,len(energies)+1,1))**2

np.savetxt('energies.dat', np.column_stack([exact_energies, energies]))

for i in range(2, 101):
    A, b, X = matrix_to_solve(xo, xn, yo, yn, N, P, Q, R, energies[i])
    Y = solve_ecns(A, b)
    Y_squared = np.array(Y)**2

    A = integrate.simpson(Y_squared, X)
    Y = Y*(1/np.sqrt(A))

    Y_exact = np.sqrt(2)*np.sin( (i-1)*np.pi*X )

    plt.plot(X, Y_exact, label="Solución teórica", c="orange")
    plt.plot(X, Y, label="Solución numérica", c="blue")
    plt.legend()
    plt.grid()
    plt.title(r"$\psi(x)$, $E = %1.3e$" %energies[i])
    plt.savefig("graphs/energy_E"+str(i-1)+".png")
    #plt.save("psi_E"+str(i)+".png")
    plt.clf()