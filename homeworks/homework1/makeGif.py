import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
import shutil

N = 64
V = 50

data = np.loadtxt('position.dat')
data = np.reshape(data,(V,N))

os.mkdir('temp')
with imageio.get_writer('./movimiento.gif', mode='I') as writer:
    for i in range(len(data[:,0])):
        fig=plt.figure()
        plt.grid()
        plt.plot(data[i,:], 'bo')
        plt.xlim((0,63))
        plt.ylim((-1,1))
        plt.savefig('./temp/'+str(i)+'.png', format='png')
#        image = imageio.imread('./temp/'+str(i)+'.png')
#        writer.append_data(image)

#shutil.rmtree('temp')
