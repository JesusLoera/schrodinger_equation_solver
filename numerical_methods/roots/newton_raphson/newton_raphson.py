# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 29/02/23
 
import numpy as np

"""
En este programa se implementa el método de Newton-Raphson para 
encontrar las raíces de una función al dar una Xo de su dominio y 
después buscar la raíz más cercana a ella.
"""

"""
Función continua
"""
def f(x):
    return np.sqrt(x+1)*( (np.cos(x/2.0))**3 )

"""
Valor de Xo sobre alrededor del cuál buscar.
"""
xo = 1

"""
Primera derivada central
"""
def derivada(f, x):
    # El error escala con h^2
    h = 0.00001
    return (f(x+h)-f(x-h))/(2*h)

"""
Método de Newton-Raphson
"""
def newton_raphson(f, xo):
    roots = []
    tolerancia = 1e-8
    xi = xo
    step = 0
    limitsteps = 1000
    while(True):
        step += 1
        x = xi - f(xi)/derivada(f,xi)
        if(abs(x-xi) <= tolerancia ):
            roots.append(x)
            break
        elif(step > limitsteps):
            print("No se encontraron raíces en ", limitsteps)
            break
        xi = x
    return roots, step

roots = newton_raphson(f, xo)
print(roots)
            
