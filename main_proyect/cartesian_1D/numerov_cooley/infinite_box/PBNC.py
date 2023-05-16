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
import scipy.special

# Configuracion de parametros de matplotlib
plt.style.use('ggplot')

# Potencial del oscilador armónico
def Potencial(x):
    v=x*0
    return v

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
xo = 0 ; xn = 10.0 ; dx = 0.001
yo = 0.0  ; yn = 0.0
N_size = int((xn - xo)/dx + 1)
# Nivel de energía a buscar
n = 1
# Parametro libre de disparo
s1 = 1.0e-5*(-1)**(n-1)   ;   s2 = 1.0e-5*(-1)**(n)
# Estimación inicial de la energía En
E = 0.0005
# Unidades atómicas
me = 1 ; hbar = 1
# Parametro del potencial
omega = 1
# X matching
xm = 7

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

def count_nodes(Y):
    nodes = 0
    for i in range (len(Y)-1):
        if (Y[i]*Y[i+1] < 0):
            nodes = nodes + 1
    return nodes

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
    nodes = count_nodes(Yk)
    print("Se detiene la busqueda tras", iteration, "iteraciones")
    print("Se contaron", nodes, "nodos")
    print("El eigenvalor de la energía es E=", Ek)
    return Yk, Ek, nodes, Energies


"""
Utilizamos el método
"""

def Y_exact(X, n):
    L = max(X) - min(X)
    Y = np.sqrt(2/L)*np.sin((n*np.pi*X)/L)
    return Y

def plot_solution(X, Y, E, n):
    Y1 = MetodoNumerov(X, Y, E)
    plt.title(r"Numerov-Cooley Method $E{}$".format(n))
    plt.plot(X, Y1, label=r"$E={:.2f}$ inicial".format(E))
    Y2, En, N, Energies = NumerovCooley(X, Y, E)
    plt.plot(X, Y2, label=r"$E{}={:.2f}$ final".format(N+1, En))
    Yexact = Y_exact(X, n)
    E_exact = ((n**2)*(np.pi**2)*(hbar**2))/(2*me*(max(X) - min(X))**2)
    plt.plot(X, Yexact, label=r"$E{}={:.2f}$ Sol. Exacta".format(n,E_exact), c="orange")
    plt.legend()
    plt.show()


plot_solution(X, Y, E, n)
