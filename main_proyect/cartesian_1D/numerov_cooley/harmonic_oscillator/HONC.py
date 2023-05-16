"""
Programado por Jesús Eduardo Loera Casas

Este programa resulve la ecn. 1D de Schrodinger
Con V = V(x), Phi(xo) = 0 y Phi(xn) = 0.

* En unidades atómicas *
me = 1, hbar = 1, qe = 1

"""

# importación de librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Configuracion de parametros de matplotlib
plt.style.use('ggplot')

# Potencial del oscilador armónico
def Potencial(x):
    me = 1.0 ; omega = 1.0
    return (0.5)*(me*omega**2)*(x**2)

# P(x)
def P(x, E):
    me = 1.0 ; hbar=1
    return 2*me*(E-Potencial(x))/(hbar**2)

# R(x)
def R(x, E):
    return 0

"""
Parametros iniciales
"""
# Definimos los parametros de la malla de integración
xo = -5.0 ; xn = 5.0 ; dx = 0.001
yo = 0.0  ; yn = 0.0
N_size = int((xn - xo)/dx + 1)
# Nivel de energía a buscar
n=0
# Parametro libre de disparo
s1 = 1.0e-5*(-1)**n   ;   s2 = 1.0e-5*(-1)**n
# Estimación inicial de la energía En
E = 1.2
# Unidades atómicas
me = 1 ; hbar = 1
# Parametro del potencial
omega = 1
# X matching
xm = np.sqrt(abs((2*E)/(me*omega**2)))

"""
Definición de la malla
"""
X = np.linspace(xo, xn, N_size)
Y = np.zeros(N_size)
Y[1] = s1
Y[N_size-2] = s2
V = Potencial(X)

"""
Integración de Numerov-Cooley
"""

def xm_index(xo, xm):
    x = xo
    index = 0
    condicion = True
    while(condicion):
        if(x>=xm):
            condicion = False
        x = x + dx
        index = index + 1
    return index-2

def y_Numerov(x, y, E):
    return y*(1 + (1/12)*(dx**2)*P(x, E))

def y_wavefunc(x, yknum, E):
    return yknum*(1/((1 + (1/12)*(dx**2)*P(x, E))))

def numerov_inward(X, Y, E):
    Ynum = y_Numerov(X, Y, E)
    for i in range(1, N_size-2):
        Ynum[i+1] = 2*Ynum[i] - Ynum[i-1] - (2*me/(hbar**2))*(dx**2)*(E-V[i])*Y[i]
        Y[i+1] = y_wavefunc(X[i+1], Ynum[i+1], E)
    # Normalizamos la función de onda
    Y = np.array(Y)
    A = np.sqrt(1.0/integrate.simpson(Y**2, X))
    Y = A*Y
    return Y
    
def numerov_backward(X, Y, E):
    Ynum = y_Numerov(X, Y, E)
    for i in range(N_size-2, 1,-1):
        Ynum[i-1] = 2*Ynum[i] - Ynum[i+1] - (2*me/(hbar**2))*(dx**2)*(E-V[i])*Y[i]
        Y[i-1] = y_wavefunc(X[i+1], Ynum[i+1], E)
    # Normalizamos la función de onda
    Y = np.array(Y)
    A = np.sqrt(1.0/integrate.simpson(Y**2, X))
    Y = A*Y
    return Y
    
def MetodoNumerov(X, Y, E):
    xmindex = xm_index(xo, xm)
    Yin = numerov_inward(X, Y, E)
    Ybw = numerov_backward(X, Y, E)
    match_constant = Yin[xmindex]/Ybw[xmindex]
    Ybw= match_constant*Ybw
    Y=[]
    for i in range(N_size):
        if(i<=xmindex):
            Y.append(Yin[i])
        else:
            Y.append(Ybw[i])
    # Normalizamos la función de onda
    Y = np.array(Y)
    A = np.sqrt(1.0/integrate.simpson(Y**2, X))
    Y = A*Y
    gE = Y[xmindex+1] + Y[xmindex-1] -2*Y[xmindex]
    return Y

def CooleyCorrection(X, Y, E):   
    m = xm_index(xo, xm) 
    Y = MetodoNumerov(X, Y, E)
    c1 = y_Numerov(X[m+1], Y[m+1], E) + y_Numerov(X[m-1], Y[m-1], E) - 2*y_Numerov(X[m], Y[m], E)
    c2 = -(hbar**2/(2*me))*(c1/(dx**2)) + (V[m] - E)*Y[m]
    dE = (1/integrate.simpson(Y**2, X))*(Y[m]*dx*(c2))
    return dE

def NumerovCooley(X, Y, E):
    tolerancia = 0.0000001
    search_energy = True
    Ek = E
    Energies = []
    iteration = 0
    while(search_energy):
        iteration = iteration + 1
        Energies.append(Ek)
        Yk = MetodoNumerov(X, Y, Ek)
        dE = CooleyCorrection(X, Yk, Ek)
        Ek = Ek + dE
        if (dE<=tolerancia):
            Energies.append(Ek)
            search_energy = False
    Yk = MetodoNumerov(X, Y, Ek)
    print("Se detiene la busqueda tras", iteration, "iteraciones")
    return Yk, Ek, Energies

"""
Prueba de los integradores de Numerov inward y backward

Energies = [0.5, 1.5, 2.5, 3.5, 4.5]
disp = 0

for energy in Energies:
    Yin = numerov_inward(X, Y, energy)
    plt.plot(X, Yin + disp, label=("E="+str(energy)))
    disp = disp + 0.05
plt.title((r"Primeros eigenestados para $V=\frac{1}{2}m\omega^{2}$"))
plt.ylabel(f"$\psi(x)$")
plt.xlabel(f"$X$")
plt.legend()
plt.show()

"""

Y1 = MetodoNumerov(X, Y, E)
plt.plot(X, Y1)
plt.show()

Y2, En, Energies = NumerovCooley(X, Y, E)
print(En)
print(Energies)
plt.plot(X, Y2)
plt.show()