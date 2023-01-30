# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 30/01/23
 
import numpy as np
import matplotlib.pyplot as plt

"""
En este programa se implementa el método de Euler para resolver
ecuaciones diferenciales de la forma: 

        y'(x) = f(y(x), x)  con y(xo)=yo
"""   

"""
Función f(y(x), x)
"""
def funcion(y,x):
    return x*y

""""
Condiciones iniciales
"""
xo = 0
yo = 2

"""
Intervalo de integración
"""
xmin = xo
xmax = 1
N = 1000
X = np.linspace(xmin, xmax, N)

"""
Método de Euler
"""
def euler(X, xo, yo):
    Y = [yo]
    dx = X[1] - X[0]
    x = xo
    y = yo
    for i in range(1, len(X)):
        xi = x + dx
        yi = y + dx*funcion(y,x)
        Y.append(yi)
        x = xi
        y = yi
    return X, Y

Sol = euler(X, xo, yo)
Y = 2*np.exp((X**2)/2)


plt.scatter(Sol[0], Sol[1], label="Método de Euler", c= "blue")
plt.plot(X, Y, label ="Solución Analítica", c="orange")
plt.legend()
plt.grid()
plt.show()