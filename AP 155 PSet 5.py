# PROBLEM 1
from numpy import array, arange
from pylab import plot, xlabel, show


def f(r, t):
    x = r[0]
    y = r[1]
    fx = alpha*x - B*x*y
    fy = Y*x*y - d*y
    return array([fx, fy], float)


alpha = 1
B = 0.5
Y = 0.5
d = 2

a = 0
b = 30
N = 1000
h = (b-a)/N

tpoints = arange(a, b, h)
xpoints = []
ypoints = []

r = array([2.0, 2.0], float)
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h*f(r, t)
    k2 = h*f(r+0.5*k1, t+0.5*h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h*f(r+k3, t+h)
    r += (k1+2*k2+2*k3+k4)/6
plot(tpoints, xpoints)
plot(tpoints, ypoints)
xlabel("t")
show()

# PROBLEM 2
from math import sin, pi
from numpy import array, arange
from pylab import plot, xlabel, show
g = 9.81
l = 0.1


def f(r, t):
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -(g/l)*sin(theta)
    return array([ftheta, fomega], float)


a = 0
b = 10
N = 1000
h = (b-a)/N

tpoints = arange(a, b, h)
thetapoints = []


r = array([179*pi/180, 0], float)
for t in tpoints:
    thetapoints.append(r[0])
    k1 = h*f(r, t)
    k2 = h*f(r+0.5*k1, t+0.5*h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h*f(r+k3, t+h)
    r += (k1+2*k2+2*k3+k4)/6
plot(tpoints, thetapoints)
xlabel("t")
show()

# PROBLEM 3
from numpy import array, arange

m = 9.1094e-31
hbar = 1.0546e-34
e = 1.6022e-19
L = 5.2918e-11
N = 1000
h = L/N


def V(x):
    return V0*x**2/a**2


V0 = 50*e
a = 1.0e-11


def f(r, x, E):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x)-E)*psi
    return array([fpsi, fphi], float)


def solve(E):
    psi = 0.0
    phi = 1.0
    r = array([psi, phi], float)

    for x in arange(0, L, h):
        k1 = h * f(r, x, E)
        k2 = h * f(r + 0.5 * k1, x + 0.5 * h, E)
        k3 = h * f(r + 0.5 * k2, x + 0.5 * h, E)
        k4 = h * f(r + k3, x + h, E)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return r[0]


E1 = 0.0
E2 = e
psi2 = solve(E1)

target = e/1000
while abs(E1-E2) > target:
    psi1, psi2 = psi2, solve(E2)
    E1, E2 = E2, E2-psi2*(E2-E1)/(psi2-psi1)

print("E = ", E2/e, "eV")
