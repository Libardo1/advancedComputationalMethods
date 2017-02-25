import numpy as np 

delta_x = 0.01
delta_t = 0.01

#abs(delta_t/delta_x*u_max)=0.5

def initial():

L = 1.0
N = 100
rho = np.zeros(N)
x = np.linspace(0,1,N)
rho[:] = 0.125
rho[x<=0.5] = 1.0
p = np.zeros(N)
p[:] = 0.1
p[x<=0.5] = 1.0

