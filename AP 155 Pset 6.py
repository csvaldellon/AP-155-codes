from numpy import copy

def banded(Aa,va,up,down):

    # Copy the inputs and determine the size of the system
    A = copy(Aa)
    v = copy(va)
    N = len(v)

    # Gaussian elimination
    for m in range(N):

        # Normalization factor
        div = A[up,m]

        # Update the vector first
        v[m] /= div
        for k in range(1,down+1):
            if m+k<N:
                v[m+k] -= A[up+k,m]*v[m]

        # Now normalize the pivot row of A and subtract from lower ones
        for i in range(up):
            j = m + up - i
            if j<N:
                A[i,j] /= div
                for k in range(1,down+1):
                    A[i+k,j] -= A[up+k,m]*A[i,j]

    # Backsubstitution
    for m in range(N-2,-1,-1):
        for i in range(up):
            j = m + up - i
            if j<N:
                v[m] -= A[i,j]*v[j]

    return v

import numpy as np
import matplotlib.pyplot as plt
import cmath
from math import exp


# define constants
m = 9.109e-31  # units of kg
L = 1e-8  # units of m
x0 = L/2  # units of m
sigma = 1e-10  # units of m
kappa = 5e10  # units of m^-1
hbar = 1.0546e-34  # in units of m^2 kg / s

# define parameters
N = 1000 # spatial slices
a = L/N # length of spatial slice
h = 10e-18 # length of time space in units of s

# define matrix entries
a1 = 1 + (h*hbar*1j)/(2*m*a**2)
a2 = -(h*hbar*1j)/(4*m*a**2)
b1 = 1 - (h*hbar*1j)/(2*m*a**2)
b2 = (h*hbar*1j)/(4*m*a**2)

# define initial wavefunction
def psi0(x):
  return exp(-0.5*((x - x0)/sigma)**2)*cmath.exp(kappa*x*1j)

# initialize arrays
xgrid = np.linspace(0, L, N+1)
vgrid = np.zeros(N+1, complex)

#initialize wavefunction
import numpy as np
psi0 = np.array(list(map(psi0, xgrid)), complex)
psi = psi0

# define A_banded
A_banded = np.zeros([3, N+1], complex)
A_banded[0,:] = np.full(N+1, a2, complex)
A_banded[1,:] = np.full(N+1, a1, complex)
A_banded[2,:] = np.full(N+1, a2, complex)
# print(A_banded)

vgrid[1: N] = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1]) # calculate components of v
psi = banded(A_banded, vgrid, 1, 1) # solve linear system

# plot wavefunction after a single step
# plt.plot(xgrid, psi, label="0")

# psi = psi0
# vgrid = np.zeros(N+1, complex)
# print(A_banded)

# for i in range(25):
  # vgrid[1: N] = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1])
  # psi = banded(A_banded, vgrid, 1, 1)

# plt.plot(xgrid, psi, label="250")

# psi = psi0
# vgrid = np.zeros(N+1, complex)
# print(A_banded)
q = 500
for i in range(q):
  vgrid[1: N] = b1*psi[1:N] + b2*(psi[2:N+1] + psi[0:N-1])
  psi = banded(A_banded, vgrid, 1, 1)

plt.plot(xgrid, psi, label=str(q))
plt.xlabel("x")
plt.ylabel(r"$\psi (x)$")
plt.ylim(-0.75, 1.0)
plt.xlim(0, 1e-08)
plt.legend()
plt.show()
plt.clf()
