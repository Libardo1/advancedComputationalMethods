import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('data.dat')
x = data[:,0]
r = data[:,1]
v = data[:,2]
p = data[:,3]
e = data[:,4]

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, r)
plt.ylim((np.min(r)-0.1, np.max(r)+0.1))
plt.xlabel(r"Posicion ($x$)")
plt.ylabel(r"Densidad ($\rho$)")
plt.title(r"Densidad contra posicion")
plt.axvline(0.5, color='k', linestyle='solid')
plt.savefig('Densidad.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, p)
plt.ylim((np.min(p)-0.1, np.max(p)+0.1))
plt.xlabel(r"Posicion ($x$)")
plt.ylabel(r"Presion ($P$)")
plt.title(r"Presion contra posicion")
plt.axvline(0.5, color='k', linestyle='solid')
plt.savefig('Presion.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, v)
plt.ylim((np.min(v)-0.1, np.max(v)+0.1))
plt.xlabel(r"Posicion ($x$)")
plt.ylabel(r"Velocidad ($v$)")
plt.title(r"Velocidad contra posicion")
plt.axvline(0.5, color='k', linestyle='solid')
plt.savefig('Velocidad.pdf', format='pdf')
plt.close()

fig=plt.figure()
ax=plt.axes()
plt.grid()
plt.plot(x, e)
plt.ylim((np.min(e)-0.1, np.max(e)+0.1))
plt.xlabel(r"Posicion ($x$)")
plt.ylabel(r"Energia ($e$)")
plt.title(r"Energia contra posicion")
plt.axvline(0.5, color='k', linestyle='solid')
plt.savefig('Energia.pdf', format='pdf')
plt.close()
