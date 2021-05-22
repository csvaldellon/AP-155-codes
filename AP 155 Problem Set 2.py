# PROBLEM 1
def f(x):
    return 2*x**4 + 3*x**2 + 4*x + 5


# Trapezoidal Rule
for N in [10, 100, 1000]:
    a = 0
    b = 5
    h = (b - a) / N

    s = 0.5 * f(a) + 0.5 * f(b)
    for k in range(1, N):
        s += f(a + k * h)

    print(h * s)

# Simpson's Rule
for N in [10, 100, 1000]:
    a = 0
    b = 5
    h = (b - a) / N

    s = f(a) + f(b)
    for k in range(1, int(N / 2 + 1)):
        s += 4 * f(a + (2 * k - 1) * h)
    for k in range(1, int(N / 2)):
        s += 2 * f(a + 2 * k * h)

    print(h * s / 3)

print("")

# PROBLEM 2
import math
import numpy as np
import matplotlib.pyplot as plt


def g(x):
    return (math.sin(math.sqrt(100 * x))) ** 2


# Adaptive Trapezoidal Rule
Int = []
i = 0
T = 0
N = 1
N_list = []

a = 0
b = 1
h = (b - a) / N

# to compute the integral using trapezoidal rule
s = 0.5 * g(a) + 0.5 * g(b)
for k in range(1, N):
    s += g(a + k * h)
T = h * s
# print("T", T)
print(N, T)
Int += [T]

N_new = 2 * N
# T_a = 0.5 * T
# print("T_a", T_a)
# s = 0.5 * g(a) + 0.5 * g(b)
error_list = []
N_list = []
while N_new < 4100:
    a = 0
    b = 1
    h_new = (b - a) / N_new
    # if i > 0:
    s = 0
    for k in range(1, N_new, 2):
        s += g(a + k * h_new)
        # print("s", s)
    T = 0.5 * T + h_new * s

    # print(N, "T", T, "T_a_old", T_a_old, "s", s, "T_a", T_a)
    # to store integral into a list for error computation
    Int += [T]
    # print(Int)
    # print("Integral list so far: " + str(Int))

    # error computation
    error = abs((Int[i + 1] - Int[i]) / 3)
    # print(N, " T = ", T, "T_a = ", T_a, "error = " + str(error))
    print(N_new, T, error)
    error_list += [error]
    N_list += [N_new]
    # N_list += [N, N_new]
    N_new *= 2
    i += 1
plt.plot(N_list, error_list)
plt.show()

print("")

# Romberg Integration
import matplotlib.pyplot as plt

maxim = 66
R = np.random.random_sample((maxim - 1, maxim - 1))
# print("orig:", R)
i = 0
N = 1
N_list = []
error_list = []
while N < maxim:
    a = 0
    b = 1
    h = (b - a) / N

    # application of trapezoidal rule
    s = 0.5 * g(a) + 0.5 * g(b)
    for k in range(1, N):
        s += g(a + k * h)

    R[i][0] = h * s
    # print(N, i + 1, 1, R[i][0])
    # error_list = []

    # Equations 5.51 and 5.49
    if i > 0:
        for m in range(0, i):
            R[i][m + 1] = R[i][m] + (R[i][m] - R[i - 1][m]) / (4**(m+1) - 1)
            # error = abs((R[i][m] - R[i - 1][m]) / (4 ** (m + 1) - 1))
            # error_list += [error]
            # print(N, i + 1, m + 2, R[i][m + 1], error)
        error = abs((R[i][i - 1] - R[i - 1][i - 1]) / (4 ** i - 1))
        N_list += [N]
        error_list += [error]
        print(N, i + 1, i, R[i][i - 1], error)
        # print(N, R[i][i], max(error_list))
    i += 1
    N = 2 * N
# print("error list", len(error_list))
# print("N list", len(N_list))
plt.plot(N_list, error_list)
plt.show()
print("")

# print("a")

# PROBLEM 3
from numpy import ones, copy, cos, tan, pi, linspace


# print("b")
def f(x):
    return 2*x**4 + 3*x**2 + 4*x + 5


# print("c")
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3, 4*N-1, N) / (4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x, w


# print("d")
def gaussxwab(N, a, b):
    x, w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w


# print("e")
N = 3
a = 0
b = 5

# x, w = gaussxw(N)
x, w = gaussxwab(N, a, b)

# print("f")
s = 0.0
for k in range(N):
    s += w[k]*f(x[k])
    # print("s", s)
    # print("k", k)
print(s)

print("")

# PROBLEM 4
from numpy import ones, copy, cos, tan, pi, linspace
import math
import matplotlib.pyplot as plt


# print("b")
def V(x):
    return x**4


def f(x):
    return math.sqrt(8/(V(b) - V(x)))


# print("c")
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3, 4*N-1, N) / (4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x, w


# print("d")
def gaussxwab(N, a, b):
    x, w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w


# print("e")
N = 20
a = 0
b_array = linspace(0.1, 2, 10**4)
s_list = []

# x, w = gaussxw(N)
for b in b_array:
    x, w = gaussxwab(N, a, b)
    s = 0.0
    for k in range(N):
        s += w[k] * f(x[k])
        # print("s", s)
        # print("k", k)
    s_list += [s]
plt.plot(list(b_array), s_list)
# print('b array', len(b_array))
# print('s list', len(s_list))
plt.show()

print("")

# PROBLEM 5
import math
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return 1 + (math.tanh(2 * x)) / 2


h = 10 ** (-5)
x_array = np.linspace(-2, 2, 100)
derivative_list = []
analytic_list = []
for x in x_array:
    derivative = (f(x + h / 2) - f(x - h / 2)) / h
    derivative_list += [derivative]
    analytic_list += [(1 - (math.tanh(2 * x)) ** 2)]
# print(derivative_list)
# print(analytic_list)
plt.plot(x_array, analytic_list, linestyle='--', color='red', label='analytic')
plt.scatter(x_array, derivative_list, label='numerical')
plt.legend()
plt.show()
plt.clf()
