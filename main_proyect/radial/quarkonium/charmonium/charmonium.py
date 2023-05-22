"""
Programado por Jesús Eduardo Loera Casas

Este programa resulve la ecn. 1D de Schrodinger
Con V = V(r), Phi(xo) = 0 y Phi(xn) = 0.

hbar = 1, c=1

"""

# importación de librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special

# Configuracion de parametros de matplotlib
plt.style.use('ggplot')

# Potencial del átomo de hidrógeno
def Potencial(x):
    alpha_s = 0.5461
    cf = 4/3
    sigma = 0.1425
    return -cf*alpha_s/x + sigma*x

# P(x)
def P(x, E, l):
    mu = 1.0 ; hbar=1
    return 2*mu*( (E-Potencial(x)) - (l*(l+1)*hbar**2)/(2*mu*x**2) )/(hbar**2)

# R(x)
def R(x, E, l):
    return 0

"""
Parametros iniciales
"""
# Definimos los parametros de la malla de integración
xo = 10e-15 ; xn = 20 ; dx = 0.001
yo = 0.0  ; yn = 0.0
N_size = int((xn - xo)/dx + 1)
# Nivel de energía a buscar
n = 0 ; l = 0
# Parametro libre de disparo
s1 = 1.0e-6   ;   s2 = 1.0e-6
# Estimación inicial de la energía En
E =  2.8
# Unidades atómicas
mc = 1.4830
mu = mc/2; hbar = 1 

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

def xm_index(Yin):
    for i in range(1, len(Yin)):
        # set it as the match point
        if Yin[i] < Yin[i-1]:
            m = i+1
            break
    return m

def y_Numerov(x, y, E, l):
    return y*(1 + (1/12)*(dx**2)*P(x, E, l))

def y_wavefunc(x, yknum, E, l ):
    return yknum*(1/((1 + (1/12)*(dx**2)*P(x, E, l))))

def numerov_inward(X, Y, E, l):
    Ynum = y_Numerov(X, Y, E, l)
    for i in range(1, N_size-2):
        y_Numerov(X, Y, E, l)
        Ynum[i+1] = 2*Ynum[i] - Ynum[i-1] - (dx**2)*P(X[i], E, l)*Y[i]
        Y[i+1] = y_wavefunc(X[i+1], Ynum[i+1], E, l)
    # Normalizamos la función de onda
    Y = np.array(Y)
    A = np.sqrt(1.0/integrate.simpson(Y**2, X))
    Y = A*Y
    return Y
    
def numerov_backward(X, Y, E, l):
    Ynum = y_Numerov(X, Y, E, l)
    for i in range(N_size-2, 1,-1):
        Ynum[i-1] = 2*Ynum[i] - Ynum[i+1] - (dx**2)*P(X[i], E, l)*Y[i]
        Y[i-1] = y_wavefunc(X[i+1], Ynum[i+1], E, l)
    # Normalizamos la función de onda
    Y = np.array(Y)
    A = np.sqrt(1.0/integrate.simpson(Y**2, X))
    Y = A*Y
    return Y
    
def MetodoNumerov(X, Y, E, l):
    Yin = numerov_inward(X, Y, E, l)
    Ybw = numerov_backward(X, Y, E, l)
    xmindex = xm_index(Yin)
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

def CooleyCorrection(X, Y, E, l):   
    m = xm_index(Y) 
    Y = MetodoNumerov(X, Y, E, l)
    c1 = y_Numerov(X[m+1], Y[m+1], E, l) + y_Numerov(X[m-1], Y[m-1], E, l) - 2*y_Numerov(X[m], Y[m], E, l)
    c2 = -(hbar**2/(2*mu))*(c1/(dx**2)) + (V[m] - E - (l*(l+1)*hbar**2)/(2*mu*(X[m]**2)) )*Y[m]
    dE = (1/integrate.simpson(Y**2, X))*(Y[m]*dx*(c2))
    return dE

def NumerovCooley(X, Y, E, l):
    tolerancia = 0.0001
    search_energy = True
    Ek = E
    Energies = []
    iteration = 0
    while(search_energy):
        iteration = iteration + 1
        Energies.append(Ek)
        Yk = MetodoNumerov(X, Y, Ek, l)
        dE = CooleyCorrection(X, Yk, Ek, l)
        Ek = Ek + dE
        if (dE<=tolerancia):
            Energies.append(Ek)
            search_energy = False
    Yk = MetodoNumerov(X, Y, Ek, l)
    nodes = count_nodes(Yk)
    Yk2 = Yk**2
    Ek = Ek*27.21
    print("Se detiene la busqueda tras", iteration, "iteraciones")
    print("Se contaron", nodes, "nodos")
    print("El eigenvalor de la energía es E=", Ek)
    return Yk, Yk2, Ek, nodes, Energies


"""
Aplicación del método de Numerov-Cooley

"""
Unl, Unl2, Enl, nodes, Energies = NumerovCooley(X, Y, E, l)

fig, ax = plt.subplots()

ax.plot(X, Unl, linestyle="--", c="#E24A33", label=f"$u(r)$")
ax.set_xlabel(r"r [$a_{o}$]")
ax.set_ylabel(f"$u_{{{n},{l}}}(r)$")
ax2 = ax.twinx()
ax2.plot(X, Unl2, c= "#348ABD", label=f"$|u(r)|^2$") 
ax2.set_ylabel(f'$|u_{{{n},{l}}}(r)|^{2}$')     
fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.title(f"Solución n={n}, l={l}, $E_{{{n},{l}}}={round(Enl, 2)}$ eV")
plt.savefig(f"imagenes/E{n}{l}.png")
np.savetxt(f'solutions/E{n}{l}.txt', np.array([X, Unl, Unl2]).T, delimiter='\t', fmt="%s")
