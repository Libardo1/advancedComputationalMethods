import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol, solve

gam=1.4
delx=0.01
L=1.0
N=int(L/delx)
mu = np.sqrt( (gam-1)/(gam+1) )

rl = 1.0
pl = 1.0
ul = 0.0
rr = 0.125
pr = 0.1
ur = 0.0
x0 = 0.5

pos=np.linspace(0,L,N)
rho=np.zeros(N)
pres=np.zeros(N)
vel=np.zeros(N)
U1=np.zeros(N)
U2=np.zeros(N)
U3=np.zeros(N)
U1_star=np.zeros(N)
U2_star=np.zeros(N)
U3_star=np.zeros(N)
F1=np.zeros(N)
F2=np.zeros(N)
F3=np.zeros(N)
F1_star=np.zeros(N)
F2_star=np.zeros(N)
F3_star=np.zeros(N)

# analytical solution
def sod_func(P):
    return (P - pr)*(np.sqrt( ((1 - mu**2)**2)*((rr*(P + mu*mu*pr))**(-1)) )**0.5) - 2*(np.sqrt(gam)/(gam - 1))*(1 - P**((gam - 1)/(2*gam)))

#PP = Symbol("PP")
#P_post = solve(sod_func(PP))

#print(P_post)
x = np.linspace(-10,10,100)
plt.plot(x,sod_func(x))
plt.show()
