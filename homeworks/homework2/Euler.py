import numpy as np
import matplotlib.pyplot as plt

gam=1.4
deltat=0.00001
delx=0.01
L=1.0
N=int(L/delx)

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

def max_velocity(v,p,r):
    return max(np.sqrt(gam*abs(p)/r) + abs(v))
def presion(e, r, v):
    return (gam-1)*(e - r*v**2/2)
def ener(p,r,v):
    return p/(gam-1) + r*v**2/2
def init(r, p, v, x, u1, u2, u3, f1, f2, f3):
    for i in range(N):
        v[i]=0.0
        if x[i]<=0.5:
            r[i]=1.0
            p[i]=1.0
        elif x[i]>0.5:
            r[i]=0.125
            p[i]=0.1
        e=ener(p[i], r[i], v[i])
        U1[i]=r[i]
        U2[i]=r[i]*v[i]
        U3[i]=e
        F1[i]=v[i]*r[i]
        F2[i]=r[i]*v[i]**2 + p[i]
        F3[i]=v[i]*(e+p[i])
def graph(r, p, v, e, x):
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
def fromU2F(u1, u2, u3):
    f1=np.zeros(N)
    f2=np.zeros(N)
    f3=np.zeros(N)
    for i in range(N):
        v=u2[i]/u1[i]
        p=presion(u3[i], u1[i], v)
        f1[i]=u2[i]
        f2[i]=u1[i]*v**2 + p
        f3[i]=v*(u3[i]+p)
    return f1, f2, f3
def Ustar_update(u1, u2, u3, f1, f2, f3, delt):
    u1_star=u1
    u2_star=u2
    u3_star=u3
    for i in range(1, N-1):
        u1_star[i]=u1[i]-delt/delx*(f1[i+1]-f1[i])
        u2_star[i]=u2[i]-delt/delx*(f2[i+1]-f2[i])
        u3_star[i]=u3[i]-delt/delx*(f3[i+1]-f3[i])
    return u1_star, u2_star, u3_star
def U_update(u1, u2, u3, u1_star, u2_star, u3_star, f1_star, f2_star, f3_star, delt):
    u1_new=u1
    u2_new=u2
    u3_new=u3
    for i in range(1,N-1):
        u1_new[i]=0.5*(u1[i]+u1_star[i] - delt/delx*(f1_star[i]-f1_star[i-1]))
        u2_new[i]=0.5*(u2[i]+u2_star[i] - delt/delx*(f2_star[i]-f2_star[i-1]))
        u3_new[i]=0.5*(u3[i]+u3_star[i] - delt/delx*(f3_star[i]-f3_star[i-1]))
    return u1_new, u2_new, u3_new
def fromU2var(u1, u2, u3):
    r=np.zeros(N)
    v=np.zeros(N)
    p=np.zeros(N)
    e=np.zeros(N)
    for i in range(N):
        r[i]=u1[i]
        v[i]=u2[i]/u1[i]
        p[i]=presion(u3[i], u1[i], v[i])
        e[i]=u3[i]
    return r, p, v, e
def Ustar_update_l(u1, u2, u3, f1, f2, f3, delt):
    u1_star=u1
    u2_star=u2
    u3_star=u3
    for i in range(1, N-1):
        u1_star[i]=0.5* ((u1[i+1] + u1[i]) - delt/delx*(f1[i+1]-f1[i]))
        u2_star[i]=0.5* ((u2[i+1] + u2[i]) - delt/delx*(f2[i+1]-f2[i]))
        u3_star[i]=0.5* ((u3[i+1] + u3[i]) - delt/delx*(f3[i+1]-f3[i]))
    return u1_star, u2_star, u3_star
def U_update_l(u1, u2, u3, u1_star, u2_star, u3_star, f1_star, f2_star, f3_star, delt):
    u1_new=u1
    u2_new=u2
    u3_new=u3
    for i in range(1,N-1):
        u1_new[i] = u1[i] - delt/delx*(f1_star[i]-f1_star[i-1])
        u2_new[i] = u2[i] - delt/delx*(f2_star[i]-f2_star[i-1])
        u3_new[i] = u3[i] - delt/delx*(f3_star[i]-f3_star[i-1])
    return u1_new, u2_new, u3_new

init(rho, pres, vel, pos, U1, U2, U3, F1, F2, F3)

t=0.0
T=0.1

while t<T:
    R,P,V,E = fromU2var(U1, U2, U3)
    deltat = (0.5*delx)/max_velocity(V,P,R)

    F1, F2, F3 = fromU2F(U1, U2, U3)
    print ('time',t)
    print ('delta_t',deltat)
    #for i in range(0,N):
    #    print(str(U1[i])+'     '+str(U2[i])+'     '+str(U3[i])+'   '+str(F1[i])+'     '+str(F2[i])+'     '+str(F3[i]))

    U1_star, U2_star, U3_star = Ustar_update_l(U1, U2, U3, F1, F2, F3, deltat)
    F1_star, F2_star, F3_star = fromU2F(U1_star, U2_star, U3_star)
    U1, U2, U3=U_update_l(U1, U2, U3, U1_star, U2_star, U3_star, F1_star, F2_star, F3_star, deltat)

    t+=deltat

R,P,V,E =fromU2var(U1, U2, U3)
graph(R, P, V, E, pos)
