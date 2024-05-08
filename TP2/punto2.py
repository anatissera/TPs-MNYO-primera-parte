import numpy as np
import matplotlib.pyplot as plt
from punto1 import runge_kutta_4

# diagrama de fases -> para comparar las Ns
# flechas

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




def dN1dt(N1, N2, r1, K1, alpha12):
    return r1 * N1 * ((K1 - N1 - alpha12 * N2) / K1)

def dN2dt(N1, N2, r2, K2, alpha21):
    return r2 * N2 * ((K2 - N2 - alpha21 * N1) / K2)

def lotka_volterra(t, y0, r1, r2, K1, K2, alpha12, alpha21):
    N1, N2 = y0
    return np.array([dN1dt(N1, N2, r1, K1, alpha12), dN2dt(N1, N2, r2, K2, alpha21)])

def runge_kutta4_system(f, t0, y0, tf, h, *args):
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

# Caso a y b: competidores intersp más fuertes eliminan a los más débiles

h = 0.1
tf = 100
t0 = 0

# Caso a: 
# K1 > K2/α21 & K2 < K1/α12
r1_a = 0.1
r2_a = 0.1
K1_a = 4000
K2_a = 3800
alpha12_a = 0.3
alpha21_a = 2
N1_0_a = 10
N2_0_a = 10

t_values_a, y_values_a = runge_kutta4_system(lotka_volterra, t0, [N1_0_a, N2_0_a], tf, h, r1_a, r2_a, K1_a, K2_a, alpha12_a, alpha21_a)

# Caso b: 
# K1 /α12 < K2 & K1 < K2/α21
r1_b = 0.1
r2_b = 0.1
K1_b = 4500
K2_b = 5000  # Capacidad de carga de la especie 2 es menor
alpha12_b = 2
alpha21_b = 0.3
N1_0_b = 10
N2_0_b = 10
y0_b = np.array([N1_0_b, N2_0_b])
t_values_b, y_values_b = runge_kutta4_system(lotka_volterra, t0, [N1_0_b, N2_0_b], tf, h, r1_b, r2_b, K1_b, K2_b, alpha12_b, alpha21_b)

# caso c: Puntos de equilibrio estables
# k1 > k2/α21 & k2 > k1/α12
r1_c = 0.3
r2_c = 0.6
K1_c = 1500
K2_c = 1400
alpha12_c = 0.4  # Msima competencia interespecífica
alpha21_c = 0.4
N1_0_c = 10
N2_0_c = 10
t_values_c, y_values_c = runge_kutta4_system(lotka_volterra, t0, [N1_0_c, N2_0_c], tf, h, r1_c, r2_c, K1_c, K2_c, alpha12_c, alpha21_c)

# Caso d: Puntos de equilibrio inestables 
# K1 > K2 α12 & K2 > K1 α21
r1_d = 0.2
r2_d = 0.2
K1_d = 1000
K2_d = 1400
alpha12_d = 1.7
alpha21_d = 2  # Mayor competencia interespecífica
N1_0_d = 10
N2_0_d = 10

t_values_d, y_values_d = runge_kutta4_system(lotka_volterra, t0, [N1_0_d, N2_0_d], tf, h, r1_d, r2_d, K1_d, K2_d, alpha12_d, alpha21_d)


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


# # Parámetros del sistema
# r1 = 0.1
# r2 = 0.1
# K1 = 4000
# K2 = 3800
# alpha12 = 0.3
# alpha21 = 2

# # Definir el rango de valores para N1 y N2
# N1_range = np.linspace(0, 5000, 100)
# N2_range = np.linspace(0, 5000, 100)

# # Calcular las derivadas en cada punto del plano
# dN1 = np.zeros((len(N1_range), len(N2_range)))
# dN2 = np.zeros((len(N1_range), len(N2_range)))

# for i in range(len(N1_range)):
#     for j in range(len(N2_range)):
#         dN1[i, j] = dN1dt(N1_range[i], N2_range[j], r1, K1, alpha12)
#         dN2[i, j] = dN2dt(N1_range[i], N2_range[j], r2, K2, alpha21)

# # Graficar las isoclinas y los campos vectoriales
# plt.figure(figsize=(10, 8))

# # Graficar las isoclinas
# plt.plot(N1_range, K1 - alpha12 * N2_range, 'r--', label='Isoclina N1')
# plt.plot(K2 - alpha21 * N1_range, N2_range, 'b--', label='Isoclina N2')

# # Graficar los campos vectoriales
# plt.quiver(N1_range, N2_range, dN1, dN2, scale=1000, color='g', alpha=0.5)

# # Etiquetas y leyenda
# plt.xlabel('N1')
# plt.ylabel('N2')
# plt.title('Isoclinas y Campos Vectoriales para el Sistema Lotka-Volterra')
# plt.legend()
# plt.grid(True)
# plt.show()


# import numpy as np
# import matplotlib.pyplot as plt

# # Define the range of values for N1 and N2
# N1_values = np.linspace(0, K1, 20)
# N2_values = np.linspace(0, K2, 20)

# # Create a grid of points
# N1, N2 = np.meshgrid(N1_values, N2_values)

# # Calculate the rate of change of each population at each point
# dN1, dN2 = lotka_volterra(0, [N1, N2], r1, r2, K1, K2, alpha12, alpha21)

# # Normalize the arrows so their size represents their speed
# M = np.hypot(dN1, dN2)
# dN1 /= M
# dN2 /= M

# plt.figure()

# # Draw the vector field
# plt.quiver(N1, N2, dN1, dN2, M, cmap='jet')

# # Draw the isoclines
# plt.contour(N1, N2, dN1, levels=[0], colors='r')
# plt.contour(N1, N2, dN2, levels=[0], colors='b')

# plt.xlabel('N1')
# plt.ylabel('N2')
# plt.title('Isoclines and vector field for case a')
# plt.show()

# puntos de equilibrio 
# Para la especie 1: N1 = k1 - α12 * N2

# Cuando N1= 0, N2= K1/ α12 y cuando N2= 0 N1= K1

# Para la especie 2:  N2 = k2 - α21 * N1 

# Cuando N2= 0 N1= K2 /α21 y cuando N1= 0 N2= K2

# dN2/dt = 0 tiene como ordenada al origen K2 y como raíz K2/α21
# dN1/dt = 0 tiene como ordenada al origen K1/α12 y como raíz K1.
# El eje de las absisas (X) es N1 y el eje de las ordenadas (Y) es N2.

# Parámetros del sistema


# Función para calcular las isoclinas
def isoclinas(N1_values, N2_values, K1, K2, alpha12, alpha21):
    isoclinas_N1 = np.where(alpha12 != 0, (K1 - alpha12 * N2_values) / alpha12, K1)
    isoclinas_N2 = np.where(alpha21 != 0, (K2 - alpha21 * N1_values) / alpha21, K2)
    return isoclinas_N1, isoclinas_N2

# def isoclinas(N1_values, N2_values, K1, K2, alpha12, alpha21):
#     isoclinas_N1 = (K1 - alpha12 * N2_values) / alpha12
#     isoclinas_N2 = (K2 - alpha21 * N1_values) / alpha21
#     return isoclinas_N1, isoclinas_N2

# # Definimos el rango de valores para N1 y N2
# N1_values_a = np.linspace(0, K1_a, 100)
# N2_values_a = np.linspace(0, K2_a, 100)

# # Calculamos las isóclinas
# isoclinas_N1, isoclinas_N2 = isoclinas(N1_values_a, N2_values_a, K1_a, K2_a, alpha12_a, alpha21_a)

# # Graficamos las isóclinas
# plt.plot(N1_values_a, isoclinas_N2, 'r--', label='Isóclina N2')
# plt.plot(isoclinas_N1, N2_values_a, 'b--', label='Isóclina N1')

# # Graficamos el punto de equilibrio
# plt.plot(K1_a/alpha12_a, K2_a/alpha21_a, 'ko', label='Punto de Equilibrio')

# # Calculamos el campo vectorial
# N1_a, N2_a = np.meshgrid(N1_values_a, N2_values_a)
# dN1dtvar = dN1dt(N1_a, N2_a, r1_a, K1_a, alpha12_a)
# dN2dtvar = dN2dt(N1_a, N2_a, r2_a, K2_a, alpha21_a)
# norm = np.sqrt(dN1dtvar**2 + dN2dtvar**2)
# dN1dt_normalized = dN1dtvar / norm
# dN2dt_normalized = dN2dtvar / norm

# # Graficamos el campo vectorial
# plt.quiver(N1_a, N2_a, dN1dt_normalized, dN2dt_normalized, color='g', label='Campo Vectorial')

# # Ajustamos los límites de los ejes
# plt.xlim(0, K1_a)
# plt.ylim(0, K2_a)

# # Añadimos leyenda y etiquetas
# plt.xlabel('N1')
# plt.ylabel('N2')
# plt.title('Isóclinas, Punto de Equilibrio y Campo Vectorial')
# plt.legend()
# plt.grid(True)
# plt.show()


def isoclinas(N1_values, N2_values, K1, K2, alpha12, alpha21):
    isoclinas_N1 = np.where(alpha12 != 0, (K1 - alpha12 * N2_values) / alpha12, K1)
    isoclinas_N2 = np.where(alpha21 != 0, (K2 - alpha21 * N1_values) / alpha21, K2)
    return isoclinas_N1, isoclinas_N2

def campo_vectorial(N1_values, N2_values, r1, r2, K1, K2, alpha12, alpha21):
    N1, N2 = np.meshgrid(N1_values, N2_values)
    dN1 = dN1dt(N1, N2, r1, K1, alpha12)
    dN2 = dN2dt(N1, N2, r2, K2, alpha21)
    norm = np.sqrt(dN1**2 + dN2**2)
    return N1, N2, dN1 / norm, dN2 / norm

# Definición de parámetros para cada caso
cases = [
    {'r1': 0.1, 'r2': 0.1, 'K1': 4000, 'K2': 3800, 'alpha12': 0.3, 'alpha21': 2, 'title': 'Caso a'},
    {'r1': 0.1, 'r2': 0.1, 'K1': 4500, 'K2': 5000, 'alpha12': 2, 'alpha21': 0.3, 'title': 'Caso b'},
    {'r1': 0.3, 'r2': 0.6, 'K1': 1500, 'K2': 1400, 'alpha12': 0.4, 'alpha21': 0.4, 'title': 'Caso c'},
    {'r1': 0.2, 'r2': 0.2, 'K1': 1000, 'K2': 1400, 'alpha12': 1.7, 'alpha21': 2, 'title': 'Caso d'}
]

# Configuración de la grilla
N1_values = np.linspace(0, 5000, 100)
N2_values = np.linspace(0, 5000, 100)
N1, N2 = np.meshgrid(N1_values, N2_values)

plt.figure(figsize=(12, 12))

for i, case in enumerate(cases, start=1):
    plt.subplot(2, 2, i)
    plt.title(case['title'])
    
    # Calcular isoclinas
    isoc_N1, isoc_N2 = isoclinas(N1_values, N2_values, case['K1'], case['K2'], case['alpha12'], case['alpha21'])
    plt.plot(N1_values, isoc_N2, 'r--', label='Isoclina N2')
    plt.plot(isoc_N1, N2_values, 'b--', label='Isoclina N1')
    
    # Calcular campo vectorial
    _, _, dN1_norm, dN2_norm = campo_vectorial(N1_values, N2_values, case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
    plt.quiver(N1, N2, dN1_norm, dN2_norm, scale=50, color='g', label='Campo Vectorial')
    
    plt.xlabel('N1')
    plt.ylabel('N2')
    plt.xlim(0, case['K1'])
    plt.ylim(0, case['K2'])
    plt.legend()

plt.tight_layout()
plt.show()



cases = {
    'a': {'r1': 0.1, 'r2': 0.1, 'K1': 4000, 'K2': 3800, 'alpha12': 0.3, 'alpha21': 2, 'title': 'Caso a'},
    'b': {'r1': 0.1, 'r2': 0.1, 'K1': 4500, 'K2': 5000, 'alpha12': 2, 'alpha21': 0.3, 'title': 'Caso b'},
    'c': {'r1': 0.3, 'r2': 0.6, 'K1': 1500, 'K2': 1400, 'alpha12': 0.4, 'alpha21': 0.4, 'title': 'Caso c'},
    'd': {'r1': 0.2, 'r2': 0.2, 'K1': 1000, 'K2': 1400, 'alpha12': 1.7, 'alpha21': 2, 'title': 'Caso d'}
}

for case, params in cases.items():
    N1_values = np.linspace(0, params['K1'], 100)
    N2_values = np.linspace(0, params['K2'], 100)

    isoc_N1, isoc_N2 = isoclinas(N1_values, N2_values, params['K1'], params['K2'], params['alpha12'], params['alpha21'])
    _, _, dN1_norm, dN2_norm = campo_vectorial(N1_values, N2_values, params['r1'], params['r2'], params['K1'], params['K2'], params['alpha12'], params['alpha21'])

    plt.figure()

    plt.quiver(N1_values, N2_values, dN1_norm, dN2_norm, scale=50, cmap='jet')

    plt.plot(N1_values, isoc_N2, 'r--', label='Isoclina N2')
    plt.plot(isoc_N1, N2_values, 'b--', label='Isoclina N1')

    plt.xlabel('N1')
    plt.ylabel('N2')
    plt.title(f'Isoclinas y campo vectorial para el caso {params["title"]}')
    plt.legend()
    plt.show()