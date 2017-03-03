import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

coor = np.loadtxt('coordenadas.dat')

node = coor[0:,0]
x    = coor[0:,1]
y    = coor[0:,2]

conectivity = np.loadtxt('conectividad.dat')

elements  = conectivity[0:,1]
nodes    = conectivity[:,1:4]

X = np.zeros(4)
Y = np.zeros(4)

plt.scatter(x,y)

for i in range(0,8):
    for j in range(0,3): 
        position = np.where(node == nodes[i,j])
        X[j] = x[position]
        Y[j] = y[position]
        
    position = np.where(node == nodes[i,0])
    X[3] = x[position]
    Y[3] = y[position]
    plt.plot(X,Y, linestyle='-')
     
plt.show()
