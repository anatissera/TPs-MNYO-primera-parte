import numpy as np
import matplotlib.pyplot as plt

from punto2 import runge_kutta4_system

# def dN_dt (N, P, r, alpha):
#     return r*N - alpha*N*P

# def dP_dt (N, P, beta, q):
#     return beta*N*P - q*P


def dN_dt (N, P, r, alpha, K):
    return r*N*(1 - N/K) - alpha*N*P

def dP_dt (N, P, beta, q):
    return beta*N*P - q*P


h = 0.1
tf = 100
t0 = 0
N1_0 = 10
N2_0 = 10

# Caso a: 
r1_a = 0.1
r2_a = 0.1
K1_a = 4000
alpha = 0.3
beta = 2
N1_0_a = 10
N2_0_a = 10

# resolver con runge kutta 4
# definir uno nuevo
def runge_kutta4_system2 (f1, f2, t0, tf, h, N1_0, N2_0, r1, r2, alpha, beta, K):
    N = [N1_0]
    P = [N2_0]
    t = [t0]
    while t[-1] < tf:
        N.append(N[-1] + h*f1(N[-1], P[-1], r1, alpha, K))
        P.append(P[-1] + h*f2(N[-1], P[-1], beta, 1))
        t.append(t[-1] + h)
    return t, N, P

t, N1, N2 = runge_kutta4_system2(dN_dt, dP_dt, t0, tf, h, N1_0_a, N2_0_a, r1_a, r2_a, alpha, beta, K1_a)

plt.figure(figsize=(10, 10))
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.plot(t, N1, label='N1')
plt.plot(t, N2, label='N2')
plt.legend()
plt.show()
