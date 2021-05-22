# %matplotlib inline
import numpy as np
from numpy.fft import rfft,irfft
import matplotlib.pyplot as plt
from math import pi, sin, cos
# define constants and parameters
L = 1e-8 # length of the box
sigma = 1e-10 # width of the wave packet
k = 5e10 #wave number
x0 = L/2 # midpoint of the box
M = 9.109e-31 # mass of the electron in kg
hbar = 1.0546e-34
N = 1000

#Discrete sine transform
def dst(y):
    N = len(y)
    y2 = np.empty(2*N,float)
    y2[0] = y2[N] = 0.0
    y2[1:N] = y[1:]
    y2[:N:-1] = -y[1:]
    a = -np.imag(rfft(y2))[:N]
    a[0] = 0.0

    return a

#Inverse discrete sine transform
def idst(a):
    N = len(a)
    c = np.empty(N+1,complex)
    c[0] = c[N] = 0.0
    c[1:N] = -1j*a[1:]
    y = irfft(c)[:N]
    y[0] = 0.0

    return y

#Initial state of the wave packet (Gaussian)
import cmath
from math import exp
def psi0(x):
    return exp(-0.5*((x - x0)/sigma)**2)*cmath.exp(k*x*1j)
    # re = exp(-((x - x0) ** 2) / (2 * sigma ** 2)) * cos(k * x)
    # imag = exp(-((x - x0) ** 2) / (2 * sigma ** 2)) * sin(k * x)
    # return re, imag

# m_r = np.zeros(N + 1,float)
# m_i = np.zeros(N + 1,float)
# xgrid = range(N+1)
m = np.zeros(N,complex)
xgrid = range(N)
from statistics import mean

for K in xgrid:
    # m_r[K], m_i[K] = (psi0(K*L/N))
    m[K] = (psi0(K*L/N))

# print(max(m_r), min(m_r), mean(m_r))
# print(max(m_i), min(m_i), mean(m_i))

# r = dst(m_r) #alpha
# i = dst(m_i) #eta
r = dst(m.real) #alpha
i = dst(m.imag)
# print(max(r))
# print(max(i))

def Psi(T):
    Q = np.zeros(N, complex)
    for j in range(N):
        p = ((pi**2)*hbar*(j**2))/(2*M*(L**2))
        Q[j] = r[j] * cos(p * T) - sin(p * T) * i[j]
    # print(Q)
    return idst(Q)

T = 1.0e-16
plt.plot(xgrid, Psi(T), label=str(T))
plt.legend()
# plt.ylim(-1, 1)
plt.show()

from matplotlib import animation

fig, ax = plt.subplots()
plt.rc('animation', html='jshtml')

x = xgrid
line, = ax.plot(x,Psi(1e-17))


def animate(i):
    line.set_ydata(Psi(i))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, 100e-16,1e-17), interval=20, blit=True, save_count=50)
ani.save("pset 6 num 2 final.gif", writer="imagemagick")
