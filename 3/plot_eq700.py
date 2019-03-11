import numpy as np
import matplotlib.pyplot as plt
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

data1 = np.genfromtxt(dir_path + '/plottp_eq1.dat', delimiter=' ')
data2 = np.genfromtxt(dir_path + '/plottp_eq2.dat', delimiter=' ')

Y1 = data1[:,1:]
X1 = data1[:,0]
Y2 = data2[:,1:]
X2 = data2[:,0]

Y = np.concatenate((Y1,Y2))
X = np.concatenate((X1,X2+X1[-1]))

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.plot(X, Y[:,0],color='b')
ax1.set_xlabel(r"t $[pS]$")
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r"T $[K]$")
ax1.tick_params('y')

ax2 = ax1.twinx()
ax2.plot(X, Y[:,1], 'r')
ax2.set_ylabel(r"pressure $[bar]$", color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()

plt.savefig('plot_eq700.png')

if "show" in sys.argv:
    plt.show()
