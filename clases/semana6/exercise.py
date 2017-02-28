import numpy as np
import matplotlib.pyplot as plt

N = 100
x = np.linspace(0,1,N)

def function_fit(coef,x):
    return coef[0]*x +coef[1]*np.power(x,2) + coef[2]*np.power(x,3)
def function(x):
    return np.exp(-x)*np.power(x,1.0/3.0)

Y = np.array([0.21829,0.14478,0.103715])

# coefficients a1 x, a2 x2, a3 x3

X = np.zeros((3,3))

for i in range(0,3): 
    for j in range(0,3): 
        X[i,j] = 1.0/(i+1+ j+1 +1)
       # print i,j 

#print X
a = np.linalg.solve(X,Y)

plt.plot(x,function(x),color = 'blue', label='Original')
plt.plot(x,function_fit(a,x), linestyle = '--', color ='red', label='fit function')
#plt.label
plt.savefig('Figure.png')

 
