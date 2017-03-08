import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('serial.txt')
size = int(np.sqrt(data.shape[0]))
data = data.reshape(size, size)
plt.imshow(data, cmap='hot', interpolation='nearest')
plt.savefig('serial.png')

