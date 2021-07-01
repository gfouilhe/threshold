import numpy as np
import matplotlib.pyplot as plt

a = 0.1
d = 0.2
J = 2**13
g = lambda u: u - u*(1-u)*(u-a)/d
u = np.linspace(0, 1, J)

plt.figure(1)
plt.plot(u,g(u), 'b', label = "g(u,a,d)")
plt.legend(loc = 'upper right')


if (a+1)**2-3*(a+d)>=0 :
    beta = (a+1 - np.sqrt((a+1)**2-3*(a+d)))/3
    gamma = (a+1 + np.sqrt((a+1)**2-3*(a+d)))/3    
    
    ### Recherche de g_^-1 ###
    
    def ginv_m(a, d, J):
    
        g = lambda u: u - u*(1-u)*(u-a)/d
        g_p = lambda u: 3/d * u**2 - 2*(a+1)/d * u + (1+ a/d)
        beta = (a+1 - np.sqrt((a+1)**2-3*(a+d)))/3
        m = g(beta)
        w = np.linspace(0,m, J)
        h = np.zeros(J)
        
        for i in range(J-1):
            V = h[i]
            for j in range(2):
                V =  V -(g(V)-w[i+1])/g_p(V)
            h[i+1] = V
        return h
    
    def ginv_m2(a, d, x):
        g = lambda u: u - u*(1-u)*(u-a)/d
        g_p = lambda u: 3/d * u**2 - 2*(a+1)/d * u + (1+ a/d)
        V = 0
        for j in range(10):
            V =  V -(g(V)-x)/g_p(V)
        return V
    
    m = g(beta)
    w1 = np.linspace(0, m, J)
    h1 = ginv_m2(a, d, w1)        
    # plt.figure(2)
    # plt.plot(w1,h1)
    
    ### Portrait de phase partie I ###
    def G(x,y):
        return -1/2 * x**2 +3/(4*d) * y**4 - 2/(3*d) *(1+a) * y**3 + 1/2 * (1+a/d)*y**2
    G_m = G(w1,h1)
        
    v1 = np.sqrt(-2*G_m)
    plt.figure(3)
    plt.plot(w1, v1, 'b')
    plt.plot(w1, -v1, 'b')
    plt.xlim(-0, m+0.1)
    plt.ylim(-.3,.3)
    
    ### Recherche g+^-1 ###
    
    def ginv_p(a, d, J):
        g = lambda u: u - u*(1-u)*(u-a)/d
        g_p = lambda u: 3/d * u**2 - 2*(a+1)/d * u + (1+ a/d)
        gamma = (a+1 + np.sqrt((a+1)**2-3*(a+d)))/3    
        M = g(gamma)
        w = np.linspace(M,1, J)
        h = np.ones(J)
        
        for i in reversed(range(1,J)):
            V = h[i]
            for j in range(2):
                V =  V -(g(V)-w[i-1])/g_p(V)
            h[i-1] = V
        return h
    
    def ginv_p2(a, d, x):
        g = lambda u: u - u*(1-u)*(u-a)/d
        g_p = lambda u: 3/d * u**2 - 2*(a+1)/d * u + (1+ a/d)
        V = 1
        for j in range(2):
            V =  V -(g(V)-x)/g_p(V)
        return V
 
    M = g(gamma)
    w2 = np.linspace(M, 1, J)
    h2 = ginv_p2(a, d, w2)        
    # plt.figure(5)
    # plt.plot(w2,h2)
    
    ### Portrait de phase partie 2 ###
    G_p = G(w2,h2)
    w, v = np.meshgrid(w2, np.linspace(-0.4, 0.4, J))
    z = 1/2*v**2+  G_p
    
    t = np.linspace(M,m, 20)
    vt = np.sqrt(-2*G(t,ginv_m2(a, d, t)))
    A = 1/2*vt**2 + G(t,ginv_p2(a, d, t))
    
    plt.figure(3)
    plt.contour(w, v, z, A, colors = 'r', vmax = 0.3)
    plt.vlines(M, -.3, .3, 'black', linestyle = 'dashed', linewidth = 0.8)
    plt.vlines(m, -.3, .3, 'black', linestyle = 'dashed', linewidth = 0.8)
    plt.hlines(0, 0, m+0.2, 'black', linewidth = 0.5)
    plt.xlabel('w')
    plt.ylabel("w'")
    plt.title("Portrait de phase Ground State", fontsize = 15)

else: 
    def ginv(a, d, x):
        g = lambda u: u - u*(1-u)*(u-a)/d
        g_p = lambda u: 3/d * u**2 - 2*(a+1)/d * u + (1+ a/d)
        for i in range(J-1):
            V = 0
            for j in range(10):
                V =  V -(g(V)-x)/g_p(V)
        return V
    def G(x,y):
        return -1/2 * x**2 +3/(4*d) * y**4 - 2/(3*d) *(1+a) * y**3 + 1/2 * (1+a/d)*y**2
    h = ginv(a, d, u)
    plt.figure(1)
    plt.plot(u,h)
    
    Ginv = G(u,h)
    
    plt.figure(6)
    plt.plot(u,np.sqrt(-2*Ginv),'b')
    plt.plot(u,-np.sqrt(-2*Ginv),'b')
    plt.hlines(0,0,1, "black", linewidth = .9)
    plt.xlim(0,.7)
    plt.title("Portrait de phase du Ground State")


### Regions ###

A = np.linspace(0,0.5,1000)
D = (A**2-A+1)/3
G = 1/3*(1-A+A**2 - np.sqrt(1-2*A))

plt.figure(4)
plt.plot(A,D, 'b', label = 'Sharp threshold')
plt.plot(A, G, 'r', label = 'Pinning')
plt.legend()
