# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 30/01/23
 
import numpy as np
import matplotlib.pyplot as plt

"""
En este programa se implementa el método de Leap-Frog para resolver
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

"""
Método de Leap Frog
"""
def leapfrog(X, xo, yo):
    Y = [yo]
    dx = X[1] - X[0]
    Y.append(yo + dx*funcion(yo, xo))
    for i in range(1, len(X) - 1):
        y = Y[i-1] + 2*dx*funcion(X[i],Y[i])
        Y.append(y)
    return X, Y

"""
Método de Runge-kutta orden 1 de 2 pasos
"""
def rk2o1(X, xo, yo):
    Y = np.zeros(len(X))
    Y[0] = yo
    dx = X[1] - X[0]
    for i in range(0, len(X)-1):
        k1 = funcion(Y[i], X[i])
        k2 = funcion( Y[i]+0.5*k1*dx, X[i]+0.5*dx)
        Y[i+1] = Y[i] + k2*dx
    return X, Y

Sol1 = euler(X, xo, yo)
Sol2 = leapfrog(X, xo, yo)
Sol3 = rk2o1(X, xo, yo)
Y = 2*np.exp((X**2)/2)


plt.plot(Sol1[0], Sol1[1], label="Método de Euler", c= "blue")
plt.plot(Sol2[0], Sol2[1], label="Método de Leap Frog", c= "green")
plt.plot(Sol3[0], Sol3[1], label="Método de Runge Kutta 2 pasos orden 1", c= "red")
plt.plot(X, Y, label ="Solución Analítica", c="orange")
plt.legend()
plt.grid()
plt.show()
