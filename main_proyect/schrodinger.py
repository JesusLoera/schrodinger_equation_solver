# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día dd/mm/aa
 
import numpy as np
from scipy import constants as const
from scipy.sparse.linalg import eigs

"""
Este programa resuelve la ecuación radial de Schrodinger para
potenciales arbitrarios V(r) [problema de dos cuerpos]
"""

# Matrices

def tridiagonal(N):
    matrix = np.zeros((N, N))
    for i in range(N):
        matrix[i][i] = -2
        if (i==0):
            matrix[i][i+1] = 1
        elif (i==(N-1)):
            matrix[i][i-1] = 1
        else:
            matrix[i][i+1] = 1
            matrix[i][i-1] = 1
    return matrix

def radial_interval(a, b, N):
    # R[0] = a
    # R[N]= b
    return np.linspace(a, b, N+1)

def potential_vector(R):
    V = e**2 / (4.0 * pi * epsilon_0) / R
    return V

def potential_matrix(R):
    V = e**2 / (4.0 * pi * epsilon_0) / R
    return np.diag(V)

def radialterm_matrix(R):
    radialterm = 1/R**2
    return np.diag(radialterm)

def radial_function(N):
    return np.zeros(N)

"""
Declaración de constantes físicas
"""
hbar = const.hbar               # Joule*Segundo
e = const.e                     # Coulomb
m_e = const.m_e                 # kilogramos  
mu = m_e                        # kilogramos                  
epsilon_0 = const.epsilon_0     # Coulomb^2 / Newton*m^2
pi = const.pi 

"""
Declaración del intervalo de integración
"""
a = 2e-10
b = 2e-8
N = 3000
R = radial_interval(a, b, N)
h = (b-a)/N

"""
Definición del número cuántico l
"""
l = 0

# Tridiagonal Matrix
tridiag = (1/h**2)*tridiagonal(N+1)

# Potential matrix
V_term = ((2*mu)/(hbar**2))*potential_matrix(R)

# Radial square term matrix
R_square = l*(l+1)*radialterm_matrix(R)

# Radial wave function array
U = radial_function(N+1)

# Matriz A
# A*U = -(2*mu*E/hbar^2)*U
A = tridiag - V_term - R_square

# Eigenvalues
number_of_eigenvalues = 30
eigenvalues, eigenvectors = eigs(A, k=number_of_eigenvalues, which='SM')
eigenvalues = -(hbar**2)*eigenvalues/(2*mu)
eigenvalues = eigenvalues*(6.242e+18)         # Joules -> eV                          

print(eigenvalues)
