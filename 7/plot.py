import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.signal import argrelextrema

dir_path = os.path.dirname(os.path.realpath(__file__))
filename = sys.argv[1]

data = np.genfromtxt(dir_path + '/' + filename+'.dat', delimiter=' ')

Y = data[:,1:]
X = data[:,0]


extr = argrelextrema(Y[:,0], np.less)
min1 = extr[0][1]


fig = plt.figure()

ax1 = fig.add_subplot(111)
plt.xlabel(r"r $[\AA]$")
plt.ylabel(r"g(r)")

ax1.plot(X,Y)
ax1.annotate(r"$g("+str(X[min1])+r")="+str(Y[min1,0])+r"$", xy=(X[min1], Y[min1,0]), 
            xytext=(X[min1]+1, Y[min1,0]-0.2),
            arrowprops=dict(facecolor='black', shrink=0.05))

# Integral for coordination number
n = 256.0 / 4940.0
I = n * np.sum(Y[:min1,0]*4*np.pi*np.power(X[:min1],2)) * X[1]
print(I)
plt.savefig(filename+'.png')

if "show" in sys.argv:
    plt.show()
