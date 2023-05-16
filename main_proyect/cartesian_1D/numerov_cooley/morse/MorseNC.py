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
    Co = 6.4 ; De = 0.3 ; a = 1 ; xe = 1.6
    return Co + De*(1-np.exp(-a*(x-xe)))**2

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
xo = 0.01 ; xn = 12 ; dx = 0.001
yo = 0.0  ; yn = 0.0
N_size = int((xn - xo)/dx + 1)
# Nivel de energía a buscar
n = 3
# Parametro libre de disparo
s1 = 1.0e-5*(-1)**n   ;   s2 = 1.0e-5*(-1)**n
# Estimación inicial de la energía En
E = 6.2
# Unidades atómicas
me = 1 ; hbar = 1
# X matching
xm = 3.3

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
Aplicación del método de Numerov-Cooley

"""

def Y_exact(X, n):
    c1 = 0.886463
    c2 = 0.774597
    Y = c1*np.exp(-c2*np.exp(1.6-X))*(np.exp(1.6-X))**(0.274597)
    return Y

def plot_solution(X, Y, En, i):
    Y, En, N, Energies = NumerovCooley(X, Y, E)
    plt.title(r"Numerov-Cooley Method $E{}$".format(N))
    plt.plot(X, Y, label=r"$E{}={:.2f}$ final".format(N, En))
    plt.legend()
    filename = "imagenes/"+str(i)+"_E"+ str(N) + ".png"
    plt.show()
    #plt.savefig(filename)

energies = np.linspace(0, 10, 5)

for i in range(len(energies)):
    E = energies[i]
    plot_solution(X, Y, E, i)
