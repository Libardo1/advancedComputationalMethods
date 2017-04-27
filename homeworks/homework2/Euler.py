import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol, solve
from scipy.optimize import fsolve

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
    return (P - pr)*(( ((1 - mu**2)**2)*((rr*(P + mu*mu*pr))**(-1)) )**0.5) - 2*(np.sqrt(gam)/(gam - 1))*(1 - P**((gam - 1)/(2*gam)))
def fun(x):
    return x**2-1
def analytic_sod(t):
    c_l = np.sqrt(gam*pl/rl)
    c_r = np.sqrt(gam*pr/rr)
    P_post = fsolve(sod_func,0.31)
    v_post = 2*(np.sqrt(gam)/(gam - 1))*(1 - P_post**((gam - 1)/(2*gam)))
    rho_post = rr*(( (P_post/pr) + mu**2 )/(1 + mu*mu*(P_post/pr)))
    v_shock = v_post*((rho_post/rr)/( (rho_post/rr) - 1))
    rho_middle = (rl)*(P_post/pl)**(1/gam)

    #Key Positions
    x1 = x0 - c_l*t
    x3 = x0 + v_post*t
    x4 = x0 + v_shock*t
    #determining x2
    c_2 = c_l - ((gam - 1)/2)*v_post
    x2 = x0 + (v_post - c_2)*t

    #start setting values
    #boundaries (can be set)
    x_min = 0.00
    x_max = L

    x = pos
    rho = np.zeros(N)
    P = np.zeros(N)
    u = np.zeros(N)
    e = np.zeros(N)

    for index in range(0,N):
        #print("index  :", str(index), str(x[index]))
        if x[index] < x1:
            #Solution b4 x1
            rho[index] = rl
            P[index] = pl
            u[index] = ul
        elif (x1 <= x[index] and x[index] <= x2):
            #Solution b/w x1 and x2
            c = mu*mu*((x0 - x[index])/t) + (1 - mu*mu)*c_l
            rho[index] = rl*(c/c_l)**(2/(gam - 1))
            P[index] = pl*(rho[index]/rl)**gam
            u[index] = (1 - mu*mu)*( (-(x0-x[index])/t) + c_l)
        elif (x2 <= x[index] and x[index] <= x3):
            # Solution b/w x2 and x3
            rho[index] = rho_middle
            P[index] = P_post
            u[index] = v_post
        elif (x3 <= x[index] and x[index] <= x4):
            #Solution b/w x3 and x4
            rho[index] = rho_post
            P[index] = P_post
            u[index] = v_post
        elif x4 < x[index]:
            #Solution after x4
            rho[index] = rr
            P[index] = pr
            u[index] = ur
        e[index] = P[index]/((gam - 1)*rho[index])
        #e[index] = ener(P[index],rho[index],u[index])

    return rho,P,u,e
def sound_velocity(p,r):
    return np.sqrt(p*gam/r)
def gas_constant(a,t):
    return a**2/(t*gam)
def max_velocity(v,p,r):
    return max(np.sqrt(gam*(p)/r) + (v))
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
def graph(r, p, v, e, x, r_a,p_a,v_a,e_a):

    fig=plt.figure()
    ax=plt.axes()
    plt.grid()
    plt.plot(x, r)
    plt.plot(x, r_a)
    #plt.ylim((np.min(r)-0.1, np.max(r)+0.1))
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
    plt.plot(x, p_a)
    #plt.ylim((np.min(p)-0.1, np.max(p)+0.1))
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
    plt.plot(x, v_a)
    #plt.ylim((np.min(v)-0.1, np.max(v)+0.1))
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
    plt.plot(x, e_a)
    #plt.ylim((np.min(e)-0.1, np.max(e)+0.1))
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
        #e[i] = p[i]/((gam - 1)*r[i])
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
T=0.2

while t<T:
    R,P,V,E = fromU2var(U1, U2, U3)
    deltat = (0.5*delx)/max_velocity(V,P,R)
    F1, F2, F3 = fromU2F(U1, U2, U3)
    U1_star, U2_star, U3_star = Ustar_update_l(U1, U2, U3, F1, F2, F3, deltat)
    F1_star, F2_star, F3_star = fromU2F(U1_star, U2_star, U3_star)
    U1, U2, U3=U_update_l(U1, U2, U3, U1_star, U2_star, U3_star, F1_star, F2_star, F3_star, deltat)

    t+=deltat

R,P,V,E =fromU2var(U1, U2, U3)
R_a,P_a,V_a,E_a = analytic_sod(T)
graph(R, P, V, E, pos, R_a,P_a,V_a,E_a)
