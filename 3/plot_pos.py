import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from mpl_toolkits.mplot3d import Axes3D

dir_path = os.path.dirname(os.path.realpath(__file__))
filename = sys.argv[1]

data = np.genfromtxt(dir_path + '/' + filename+'.dat', delimiter=' ')

Y = data[:,1:]
X = data[:,0]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(X,Y[:,-1])
plt.xlabel(r"t $[pS]$")
plt.ylabel(r"deviation from p0 $[\AA]$")
plt.savefig(filename+'_all.png')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range((Y.shape[1]/3)):
    x = Y[:,3*i]
    y = Y[:,3*i+1]
    z = Y[:,3*i+2]
    ax.plot(x, y, z)
plt.savefig(filename+'_trac.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlabel(r"t $[pS]$")
plt.ylabel(r"position component $[\AA]$")
plt.plot(X,Y[:,:-1])
plt.savefig(filename+'_trac_2d.png')

