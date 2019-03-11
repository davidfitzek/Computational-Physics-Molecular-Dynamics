import numpy as np
import matplotlib.pyplot as plt
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

data = np.genfromtxt(dir_path + '/plotg.dat', delimiter=' ')


n = 256.0 / 4700.0

Y = data[1:,1]
X = data[1:,0]
s_q = [1 + 4 * np.pi * n * np.sum(np.power(X,2)*(Y-1)*np.sin(q*X)/(q*X))*X[0] for q in np.arange(0.1,20,0.1)]
q = [q for q in np.arange(0.1,20,0.1)]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(q,s_q, label="Integral")


data = np.genfromtxt(dir_path + '/sq.dat', delimiter=' ')

Y = data[:,1:]
X = data[:,0]

ax1.plot(X,Y,label="Direct calculation")

plt.xlabel(r"q $[\AA^{-1}]$")
plt.ylabel("S(q)")
plt.legend()

plt.savefig('plot.png')

if "show" in sys.argv:
    plt.show()
