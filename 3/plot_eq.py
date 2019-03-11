import numpy as np
import matplotlib.pyplot as plt
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
filename = sys.argv[1]

data = np.genfromtxt(dir_path + '/' + filename+'.dat', delimiter=' ')

Y = data[:,1:]
X = data[:,0]

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

plt.savefig(filename+'.png')

if "show" in sys.argv:
    plt.show()
