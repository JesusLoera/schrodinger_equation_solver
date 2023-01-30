# Programa elaborado por Jesús Eduardo Loera Casas
# Elaborado el día 29/01/23
 
import numpy as np
from scipy import optimize

def function(x):
    return np.sqrt(x+1)*( (np.cos(x/2.0))**3 )

def derivada(x):
    # El error escala con h^2
    h = 0.00001
    return (function(x+h)-function(x-h))/(2*h)


sol = optimize.root_scalar(function, method="newton", x0 = 1, fprime=derivada)
print(sol.root, sol.iterations, sol.converged, sol.flag)

sol = optimize.root_scalar(function, method="bisect", bracket=[0, 6])
print(sol.root, sol.iterations, sol.converged, sol.flag)
