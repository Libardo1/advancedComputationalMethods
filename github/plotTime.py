import numpy as np
import matplotlib.pyplot as plt

# read data
processes = [1,2,4]
for pro in processes:

    time.append(np.loadtxt('time'+pro+'.dat'))

plt.plot(time )
