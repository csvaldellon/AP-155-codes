# PROBLEM 1
from numpy import zeros, sin
from cmath import exp, pi
import matplotlib.pyplot as plt

N = 1000


def dft(y):
    c = zeros(N // 2 + 1, complex)
    for k in range(N // 2 + 1):
        for n in range(N):
            c[k] += y[n] * exp(-2j * pi * k * n / N)
    return c


y = [sin(pi * n / N) * sin(20 * pi * n / N) for n in range(N)]
c = dft(y)
c = [abs(r) for r in c]

plt.plot(range(N // 2 + 1), c)
# plt.xlim(0, 11)
plt.scatter([10], [max(c)], color="red")
plt.show()

# PROBLEM 2
from numpy import loadtxt
from pylab import plot, xlim, show
from numpy.fft import rfft, irfft

dow = loadtxt("C:/Users/Val/Desktop/dow.txt", float)
plt.plot(dow)
plt.show()

# print(dow)

c = rfft(dow)
# c = [abs(r) for r in c]
# print(c)
# print(len(c))

# c_inv = irfft(c)
# print(c_inv)
# plt.plot(dow, c="blue")
# plt.plot(c_inv, c="red")
# plt.show()

first_10 = int(round(0.1 * len(c), 0))
# c_sorted = sorted(c)
i = len(c) - 1
while i > first_10:
    # c_sorted[i] = 0
    c[i] = 0
    i -= 1
# print(c_sorted)
# print(c)
# print(len(c))
c_inv = irfft(c)
plt.plot(dow, c="red")
plt.plot(c_inv, c="black")
plt.show()
plt.clf()

from numpy import loadtxt, exp, array
from pylab import imshow, show

# from matplotlib.pyplot import hist2d
blur = loadtxt("C:/Users/Val/Desktop/blur.txt")
print(blur.shape)
print(blur)
imshow(blur)
show()
sig = 25


def f(x, y):
    return exp(-((x ** 2) + (y ** 2)) / (2 * (sig ** 2)))


# def psf(x0, y0):
# psf_list = []
# for x in range(1024):
# row = []
# for y in range(1024):
# row += [f(x-x0, y-y0)]
# psf_list += [row]
# return psf_list


psf = []
for y in range(1024):
    row = []
    for x in range(1024):
        # row += [f(x, y) + f(x - 1023, y) + f(x, y - 1023) + f(x - 1023, y - 1023)]
        # row += [f(x, y)]
        row += [f(x, y) + f(x - 1024, y) + f(x, y - 1024) + f(x - 1024, y - 1024)]
    psf += [row]

# print(array(psf).shape)
imshow(psf)
print(array(psf).shape)
# show()
# imshow(psf(1023, 0))
# show()
# imshow(psf(0, 1023))
# show()
# imshow(psf(1023, 1023))
# show()
# imshow(psf(0, 0) + psf(0, 1023) + psf(1023, 0) + psf(1023, 1023))
show()

from numpy.fft import rfft2, irfft2
from numpy import zeros, complex_

b = rfft2(blur)
print(b.shape)
# print(b)
p = rfft2(psf)
print(p.shape)
# print(p)
i = 0
div = zeros((1024, 513), dtype=complex_)
while i < 1024:
    k = 0
    while k < 513:
        if p[i][k].real < 10**(-3) and p[i][k].imag < 10**(-3):
            div[i][k] = b[i][k]
        else:
            div[i][k] = b[i][k]/p[i][k]
        # print(div)
        k += 1
    i += 1
# while i < 1024:
    # k = 0
    # while k < 513:
       # if p[i][k].real < 10 ** (-3) and p[i][k].imag < 10 ** (-3):
           #  p[i][k] = 1
        # else:
          #   p[i][k] = p[i][k]
        # # print(div)
       #   k += 1
    # i += 1
# div = b / p
new = irfft2(div)
imshow(new)
show()
