# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 05/01/23
 
"""
Se resuelve las edo con condciones en la frontera del tipo

    y''(x) = k  ;  y(xo)=yo , y(xn)=yn

"""

# importamos las librerías requeridas
import numpy as np
import matplotlib.pyplot as plt

"""
Funciones que ejecuta el método numérico
"""

def finite_diff(k, xo, xn, yo, yn, N):
    X = np.linspace(xo, xn, N)
    dx = (xn-xo)/(N-1)
    b = np.zeros(N)
    b[0] = yo ; b[N-1] = yn
    for i in range(1, N-1):
        b[i] = k*dx**2
    A = A_matrix(N)
    # Resolvemos la ecuación AY = b con numpy
    Y = np.linalg.solve(A, b)
    return X, Y

def A_matrix(N):
    matrix = np.zeros((N, N))
    for i in range(N):
        matrix[i][i] = -2
        if (i==0):
            matrix[i][i] = 1
        elif (i==(N-1)):
            matrix[i][i] = 1
        else:
            matrix[i][i+1] = 1
            matrix[i][i-1] = 1
    return matrix
    

"""
Parametros de la EDO
"""
# Intervalo de integración
xo = 0 ; xn = 5
# Condiciones en la frontera
yo = 0 ; yn = 50
# Constante de la edo
k = -9.8
# Puntos en la malla de integración (N) {xo,x1,...,xN}
N = 11

"""
Llamamos nuestro método con una edo partícular
"""

X, Y = finite_diff(k, xo, xn, yo, yn, N)

# Graficamos nuestro ejemplo
plt.plot(X, Y, label = "y(x)")
plt.xlabel("Eje X")
plt.ylabel("Eje Y")
plt.show()