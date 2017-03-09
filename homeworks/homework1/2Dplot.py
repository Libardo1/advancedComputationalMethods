import numpy as np
import matplotlib.pyplot as plt

N = 64
V = 200
delta_t = 5E-3
data = np.loadtxt('position.dat')
data = data.reshape(V,N)
plt.imshow(data, cmap='hot', interpolation='nearest')
plt.savefig('evolution.png')
