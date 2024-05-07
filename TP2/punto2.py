import numpy as np
import matplotlib.pyplot as plt
from punto1 import runge_kutta_4



# runge-kutta 4 para ecuaciones diferenciales de primer orden (matriz)

# dN1/dt = r1N1 ((K1 − N1 − α12N2) / K1)
# dN2/dt = r2N2 ((K2 − N2 − α21N1)/K2)

# Donde para la especie i tenemos que Ni representa su población, Ki su capacidad de carga en el sistema,
# ri su tasa intrínseca de crecimiento y αij es el coeficiente de competencia interespecífica de 
# la especie j sobre la especie i. Y no necesariamente αij = αji.

# aproximen las soluciones N1(t) y N2(t) utilizando
# distintos valores de los parámetros (N1(t = 0), N2(t = 0), r1, r2, K1, K2, α12, α21) de forma 
# que cubran todos los casos relevantes.
# Para obtener más información de las soluciones es importante estudiar los puntos de equilibrio de cada
# especie, donde dN1/dt = 0 y dN2/dt = 0. Estos puntos en el grafico de N2 vs N1 definen curvas 
# denominadas isoclinas cero. 


def lotka_volterra(t, y, r1, r2, K1, K2, alpha12, alpha21):
    N1, N2 = y
    dN1dt = r1 * N1 * (K1 - N1 - alpha12 * N2) / K1
    dN2dt = r2 * N2 * (K2 - N2 - alpha21 * N1) / K2
    return np.array([dN1dt, dN2dt])

def runge_kutta4_system(f, y0, t0, tf, h, *args):
    t_values = np.arange(t0, tf + h, h)
    n = len(t_values)
    y_values = np.zeros((n, len(y0)))
    y_values[0] = y0
    
    for i in range(1, n):
        k1 = h * f(t_values[i-1], y_values[i-1], *args)
        k2 = h * f(t_values[i-1] + h/2, y_values[i-1] + k1/2, *args)
        k3 = h * f(t_values[i-1] + h/2, y_values[i-1] + k2/2, *args)
        k4 = h * f(t_values[i-1] + h, y_values[i-1] + k3, *args)
        y_values[i] = y_values[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6
    
    return t_values, y_values

# # Parámetros del sistema
# r1 = 0.1
# r2 = 0.1
# K1 = 4000
# K2 = 3800
# alpha12 = 0.3
# alpha21 = 2

# # Condiciones iniciales
# N1_0 = 10
# N2_0 = 10
# y0 = np.array([N1_0, N2_0])

# # Tiempo de integración
# t0 = 0
# tf = 100
# h = 0.1

# # Resolver el sistema de ecuaciones
# t_values, y_values = runge_kutta4_system(lotka_volterra, y0, t0, tf, h, r1, r2, K1, K2, alpha12, alpha21)

# # Graficar las soluciones
# plt.figure(figsize=(10, 6))
# plt.plot(t_values, y_values[:, 0], label='N1(t)')
# plt.plot(t_values, y_values[:, 1], label='N2(t)')
# plt.xlabel('Tiempo')
# plt.ylabel('Población')
# plt.title('Evolución temporal de las poblaciones N1 y N2')
# plt.legend()
# plt.grid(True)
# plt.show()


# Caso a y b: competidores intersp más fuertes eliminan a los más débiles

# Caso a: 
# K1 > K2/α21 & K2 < K1/α12
r1 = 0.1
r2 = 0.1
K1 = 4000
K2 = 3800
alpha12 = 0.3
alpha21 = 2
N1_0 = 10
N2_0 = 10
y0 = np.array([N1_0, N2_0])
t0 = 0
tf = 100
h = 0.1
t_values_a, y_values_a = runge_kutta4_system(lotka_volterra, y0, t0, tf, h, r1, r2, K1, K2, alpha12, alpha21)

# Caso b: 
# K1 /α12 < K2 & K1 < K2/α21
r1 = 0.1
r2 = 0.1
K1 = 4500
K2 = 5000  # Capacidad de carga de la especie 2 es menor
alpha12 = 2
alpha21 = 0.3
N1_0 = 10
N2_0 = 10
y0 = np.array([N1_0, N2_0])
t_values_b, y_values_b = runge_kutta4_system(lotka_volterra, y0, t0, tf, h, r1, r2, K1, K2, alpha12, alpha21)

# caso c: Puntos de equilibrio estables
# k1 > k2/α21 & k2 > k1/α12
r1 = 0.3
r2 = 0.6
K1 = 1500
K2 = 1400
alpha12 = 0.4  # Msima competencia interespecífica
alpha21 = 0.4
N1_0 = 10
N2_0 = 10
y0 = np.array([N1_0, N2_0])
t_values_c, y_values_c = runge_kutta4_system(lotka_volterra, y0, t0, tf, h, r1, r2, K1, K2, alpha12, alpha21)

# Caso d: Puntos de equilibrio inestables 
# K1 > K2 α12 & K2 > K1 α21
r1 = 0.2
r2 = 0.2
K1 = 1000
K2 = 1400
alpha12 = 1.7
alpha21 = 2  # Mayor competencia interespecífica
N1_0 = 10
N2_0 = 10
y0 = np.array([N1_0, N2_0])
t_values_d, y_values_d = runge_kutta4_system(lotka_volterra, y0, t0, tf, h, r1, r2, K1, K2, alpha12, alpha21)


plt.figure(figsize=(10, 10))

# Caso a
plt.subplot(2, 2, 1)
plt.plot(t_values_a, y_values_a[:, 0], label='N1(t)')
plt.plot(t_values_a, y_values_a[:, 1], label='N2(t)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Caso a: N1 sobrevive (es más fuerte), N2 se extingue')
plt.legend()
plt.grid(True)

# Caso b
plt.subplot(2, 2, 2)
plt.plot(t_values_b, y_values_b[:, 0], label='N1(t)')
plt.plot(t_values_b, y_values_b[:, 1], label='N2(t)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Caso b: N2 sobrevive, N1 se extingue')
plt.legend()
plt.grid(True)

# Caso c
plt.subplot(2, 2, 3)
plt.plot(t_values_c, y_values_c[:, 0], label='N1(t)')
plt.plot(t_values_c, y_values_c[:, 1], label='N2(t)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Caso c: Equilibrio estable entre N1 y N2')
plt.legend()
plt.grid(True)

# Caso d
plt.subplot(2, 2, 4)
plt.plot(t_values_d, y_values_d[:, 0], label='N1(t)')
plt.plot(t_values_d, y_values_d[:, 1], label='N2(t)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Caso d: Equilibrio inestable entre N1 y N2')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


