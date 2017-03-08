import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('serial.txt')
data = data.reshape(3, 1000)

plt.plot(data.T)
plt.savefig('serial.png')
