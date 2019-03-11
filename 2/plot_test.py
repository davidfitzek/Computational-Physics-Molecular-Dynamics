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
ax1.set_xlabel(r"t $[pS]$")
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel(r"$E [eV]$")
ax1.tick_params('y')
ax1.plot(X,Y[:,1],label=r"$E_{pot}$")
ax1.plot(X,Y[:,2],label=r"$E_{tot}$")

plt.legend()

ax2 = ax1.twinx()
ax2.set_ylabel(r"$E [eV]$", color='r')
ax2.tick_params('y', colors='r')
ax2.plot(X,Y[:,0]-800,label=r"$E_{kin}$", color='r')

fig.tight_layout()
plt.legend()

plt.savefig(filename+'.png')

if "show" in sys.argv:
    plt.show()
