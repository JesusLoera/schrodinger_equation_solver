# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día dd/mm/aa
 
# importamos librerías 
import numpy as np 

"""
En este programa implementamos el método de bisección para determinar 
las raíces de una función continua f(x) en algún intervalo [a,b]
"""

"""
Función continua
"""
def f(x):
    return np.sin(x)

"""
Intervalo de inspección
"""
a = 0.1
b = 2*np.pi

"""
Método de bisección
"""
def bisection(f, a, b):
    # si la raíz está entre a y b, esperamos que b-a < tolerancia
    tolerancia = 0.0001            
    dx = b - a
    # contaremos los pasos en los que se ecnuetra la primera raíz
    step = 0
    while(0 < 1):
        roots = []
        step += 1
        c = 0.5*(a+b)
        dx = b - a
        if(dx <= tolerancia):
            roots.append(c)
            break
        if (f(a)*f(c) > 0):
            a = c
        elif(f(a)*f(c) < 0):
            b = c
        else:
            if(f(c)==0):
                roots.append(c)
            if(f(a)==0):
                roots.append(a)
            break
        if (step > 10000):
            break
    return roots, step

roots = bisection(f, a, b)
print(roots)



