import numpy as np
from scipy import interpolate
import scipy.special as spe
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Configuracion de parametros de matplotlib
plt.style.use('ggplot')

# Quantum numbers
m = 0
l = 0
n = 2

# Cargamos la función radial de schrodinger calculada
#  para hacer una inteprolación 
Ur = np.loadtxt(f"solutions/E{n}{l}.txt", usecols=(1))
R = np.loadtxt(f"solutions/E{n}{l}.txt", usecols=(0))
model = interpolate.interp1d(R, Ur)

# Plotereamos en el plano z=0
X = np.linspace(-20, 20, 500)
Y = np.linspace(-20, 20, 500)
XX, YY = np.meshgrid(X, Y)
Z = np.zeros((len(X), len(Y)))
phi = np.pi/2

def wave_function(XX, YY, interpolation):
    r = np.sqrt(XX**2 + YY**2)
    theta = np.arctan2(YY, XX)
    sphHarm = spe.sph_harm(m,l,phi,theta)
    ZZ = np.real((interpolation(r)*sphHarm)/(r))
    return ZZ

WaveFunction = wave_function(XX, YY, model)
plt.title(f"$|\psi_{{{n}{l}{m}}}(r, \\theta, \phi)|^{2}$ ; $z=0$")
plt.xlabel("Eje X")
plt.ylabel("Eje Y")
plt.contourf(YY, XX, WaveFunction**2, 100, cmap="inferno")
plt.colorbar()
plt.show()
plt.clf()