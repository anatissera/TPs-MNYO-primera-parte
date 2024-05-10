import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from punto2 import runge_kutta4_system

# def dN_dt (N, P, r, alpha):
#     return r*N - alpha*N*P

# def dP_dt (N, P, beta, q):
#     return beta*N*P - q*P


def dN_dt (N, P, r, alpha, K):
    return r * N * (1 - N / K) - alpha * N * P

def dP_dt (N, P, beta, q):
    return beta * N * P - q * P

# resolver con runge kutta 4
# definir uno nuevo
def lotka_volterra(t, y0, r1, r2, K1, alpha, beta, q):
    N, P = y0
    return np.array([dN_dt(N, P, r1, alpha, K1), dP_dt(N, P, beta, q)])

# sirve el mismo runge kutta 4

# Parámetros del sistema
h = 0.1
tf = 250
t0 = 0
N1_0_conejos = 2000
N2_0_zorros = 10

r1 = 0.1  # Tasa de crecimiento de los conejos (conejos/mes)
alpha = 0.005  # Eficiencia de captura (conejos*zorros/mes)
r2 = 0.04  # Tasa de crecimiento de los zorros (zorros/mes)
beta = 0.00004  # Eficiencia para convertir presas en nuevos depredadores (conejos*zorros/mes)
q = 0.05  # Tasa de mortalidad per cápita de los zorros (zorros/mes)
K1 = 10000  # Capacidad de carga del ambiente para los conejos


# Resolución del sistema de ecuaciones diferenciales
t_values, y_values = runge_kutta4_system(lotka_volterra, t0, [N1_0_conejos, N2_0_zorros], tf, h, r1, r2, K1, alpha, beta, q)

# Graficar las soluciones
plt.figure(figsize=(10, 6))
plt.plot(t_values, y_values[:, 0], label='Conejos')
plt.plot(t_values, y_values[:, 1], label='Zorros')
plt.xlabel('Tiempo (meses)')
plt.ylabel('Población')
plt.title('Dinámica de la población de conejos y zorros')
plt.legend()
plt.grid(True)
plt.show()

h = 0.1
tf = 250
t0 = 0
N1_0_conejos = 2000
N2_0_zorros = 10

r1 = 0.1  # Tasa de crecimiento de los conejos (conejos/mes)
alpha = 0.02  # Eficiencia de captura (conejos*zorros/mes)
r2 = 0.08  # Tasa de crecimiento de los zorros (zorros/mes)
beta = 0.0001  # Eficiencia para convertir presas en nuevos depredadores (conejos*zorros/mes) - Aumentada
q = 0.05  # Tasa de mortalidad per cápita de los zorros (zorros/mes)
K1 = 10000  # Capacidad de carga del ambiente para los conejos

# Resolución del sistema de ecuaciones diferenciales
t_values, y_values = runge_kutta4_system(lotka_volterra, t0, [N1_0_conejos, N2_0_zorros], tf, h, r1, r2, K1, alpha, beta, q)

# Graficar las soluciones
plt.figure(figsize=(10, 6))
plt.plot(t_values, y_values[:, 0], label='Conejos')
plt.plot(t_values, y_values[:, 1], label='Zorros')
plt.xlabel('Tiempo (meses)')
plt.ylabel('Población')
plt.title('Dinámica de la población de conejos y zorros')
plt.legend()
plt.grid(True)
plt.show()

# def calcular_isoclinas_y_graficar_contour_color(N_values, P_values, r, alpha, K, beta, q, legend_loc):
#     N, P = np.meshgrid(N_values, P_values)
#     dN = dN_dt(N, P, r, alpha, K)
#     dP = dP_dt(N, P, beta, q)
    
#     plt.figure(figsize=(10, 8))
#     plt.contour(N, P, dN, levels=[0], colors='blue', linewidths=2)
#     plt.contour(N, P, dP, levels=[0], colors='red', linewidths=2)
#     plt.quiver(N[::10, ::10], P[::10, ::10], dN[::10, ::10], dP[::10, ::10], scale=300, color='black', cmap='jet', label='Campo Vectorial')

#     # Encontrar los puntos de equilibrio
#     def f(x):
#         return [dN_dt(x[0], x[1], r, alpha, K), dP_dt(x[0], x[1], beta, q)]
#     punto_eq = fsolve(f, [0.5 * K, 0.5 * K])
#     plt.plot(punto_eq[0], punto_eq[1], 'o', color='green', markersize=8, label='Punto de Equilibrio')

#     plt.xlabel('Población de Presas (N)')
#     plt.ylabel('Población de Depredadores (P)')
#     plt.title('Isoclinas y Campo Vectorial')
#     plt.legend()
#     plt.show()

# # Definir parámetros
# r = 0.1  # Tasa de crecimiento de las presas
# q = 0.05  # Tasa de mortalidad per cápita de los depredadores
# alpha = 0.02  # Eficiencia de captura
# beta = 0.03  # Eficiencia para convertir presas en nuevos depredadores
# K = 1000  # Capacidad de carga del ambiente

# # Definir rangos de valores para las poblaciones de presas y depredadores
# N_values = np.linspace(0, K, 100)
# P_values = np.linspace(0, K, 100)

# # Graficar isoclinas y campo vectorial para los parámetros dados
# calcular_isoclinas_y_graficar_contour_color(N_values, P_values, r, alpha, K, beta, q, 'upper right')