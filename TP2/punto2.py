import numpy as np
import matplotlib.pyplot as plt
# from punto1 import runge_kutta_4
from scipy.optimize import fsolve

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

def punto_equilibrio(r1, r2, K1, K2, a12, a21):
    # Para encontrar los puntos de equilibrio, primero igualamos las derivadas a cero y resolvemos para N1 y N2.
    # Busco donde dN1/dt = 0 y dN2/dt = 0
    def f(x):
        return [r1 * x[0] * (K1 - x[0] - a12 * x[1]) / K1, r2 * x[1] * (K2 - x[1] - a21 * x[0]) / K2]
    
    return fsolve(f, [K1/2, K2/2])


# Caso a y b: competidores intersp más fuertes eliminan a los más débiles

h = 0.1
tf = 100
t0 = 0
N1_0 = 10
N2_0 = 10

# Caso a: 
# K1 > K2/α21 & K2 < K1/α12

# Caso b: 
# K1 /α12 < K2 & K1 < K2/α21

# caso c: Puntos de equilibrio estables
# k1 > k2/α21 & k2 > k1/α12

# Caso d: Puntos de equilibrio inestables 
# K1 > K2 α12 & K2 > K1 α21


# # Parámetros del sistema
# r1 = 0.1
# r2 = 0.1
# K1 = 4000
# K2 = 3800
# alpha12 = 0.3
# alpha21 = 2


# puntos de equilibrio 
# Para la especie 1: N1 = k1 - α12 * N2

# Cuando N1= 0, N2= K1/ α12 y cuando N2= 0 N1= K1

# Para la especie 2:  N2 = k2 - α21 * N1 

# Cuando N2= 0 N1= K2 /α21 y cuando N1= 0 N2= K2

# dN2/dt = 0 tiene como ordenada al origen K2 y como raíz K2/α21
# dN1/dt = 0 tiene como ordenada al origen K1/α12 y como raíz K1.
# El eje de las absisas (X) es N1 y el eje de las ordenadas (Y) es N2.






# # Define the range of values for N1 and N2
# N1_values = np.linspace(0, K1, 20)
# N2_values = np.linspace(0, K2, 20)

# # Create a grid of points
# N1, N2 = np.meshgrid(N1_values, N2_values)

# # Calculate the rate of change of each population at each point
# dN1, dN2 = lotka_volterra(0, [N1, N2], r1, r2, K1, K2, alpha12, alpha21)

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




cases = {
    'a': {'r1': 0.1, 'r2': 0.1, 'K1': 4000, 'K2': 3800, 'alpha12': 0.3, 'alpha21': 2, 'title': 'Caso a', 'case': 'N1 sobrevive (es más fuerte), N2 se extingue'},
    'b': {'r1': 0.1, 'r2': 0.1, 'K1': 4500, 'K2': 5000, 'alpha12': 2, 'alpha21': 0.3, 'title': 'Caso b', 'case': 'N2 sobrevive, N1 se extingue'},
    'c': {'r1': 0.3, 'r2': 0.6, 'K1': 1500, 'K2': 1400, 'alpha12': 0.4, 'alpha21': 0.4, 'title': 'Caso c', 'case': 'Equilibrio estable entre N1 y N2'},
    'd': {'r1': 0.2, 'r2': 0.2, 'K1': 1000, 'K2': 1400, 'alpha12': 1.7, 'alpha21': 2, 'title': 'Caso d', 'case': 'Equilibrio inestable entre N1 y N2'}
}
# def campo_vectorial(r1, r2, K1, K2, alpha12, alpha21):
#     N1, N2 = np.meshgrid(np.linspace(0, K1, 20), np.linspace(0, K2, 20))
#     dN1 = dN1dt(N1, N2, r1, K1, alpha12)
#     dN2 = dN2dt(N1, N2, r2, K2, alpha21)
#     # Normalize the arrows so their size represents their speed
#     norm = np.sqrt(np.max(dN1)**2 + np.max(dN2)**2)
#     return dN1 / norm, dN2 / norm

# def isoclinas(N1_values, N2_values, K1, K2, alpha12, alpha21):
#     isoclinas_N1 = (K1 - alpha12 * N2_values) / alpha12
#     isoclinas_N2 = (K2 - alpha21 * N1_values) / alpha21
#     return isoclinas_N1, isoclinas_N2

# plt.figure(figsize=(12, 12))

# for i, case in enumerate(cases.values(), start=1):
#     N1_values = np.linspace(0, case['K1'], 20)
#     N2_values = np.linspace(0, case['K2'], 20)
#     N1, N2 = np.meshgrid(N1_values, N2_values)
    
#     plt.subplot(2, 2, i)
#     plt.title(case['title'])
    
#     # Calcular isoclinas
#     isoc_N1, isoc_N2 = isoclinas(N1_values, N2_values, case['K1'], case['K2'], case['alpha12'], case['alpha21'])
    
#     # Graficar isóclinas
#     plt.plot(N1_values, isoc_N2, 'r--', label='Isoclina N2')
#     plt.plot(isoc_N1, N2_values, 'b--', label='Isoclina N1')
    
#     # Calcular campo vectorial
#     dN1_norm, dN2_norm = campo_vectorial(case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
#     plt.quiver(N1, N2, dN1_norm, dN2_norm, scale=50, color='g', label='Campo Vectorial', cmap = 'jet')
    
#     plt.xlabel('N1')
#     plt.ylabel('N2')
#     plt.xlim(0, case['K1'])
#     plt.ylim(0, case['K2'])
#     plt.legend()

# plt.tight_layout()
# plt.show()

def graficar_soluciones_rk(t0, N1_0, N2_0, tf, h, cases):
    plt.figure(figsize=(10, 10))

    for i, case in enumerate(cases.values(), start=1):
        t_values, y_values = runge_kutta4_system(lotka_volterra, t0, [N1_0, N2_0], tf, h, case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        plt.subplot(2, 2, i)        
        plt.plot(t_values, y_values[:, 0], label='N1(t)')
        plt.plot(t_values, y_values[:, 1], label='N2(t)')
        plt.xlabel('Tiempo')
        plt.ylabel('Población')
        plt.title(case['title'] + ': ' + case['case'])
        plt.legend()

    plt.tight_layout()
    plt.show()
    
def calcular_isoclinas_y_graficar_contour_color(title, r1, r2, K1, K2, alpha12, alpha21):
    N1 = np.linspace(0, K1, 1000)
    N2 = np.linspace(0, K2, 1000)
    N1, N2 = np.meshgrid(N1, N2)
    N1_isocline = dN1dt(N1, N2, r1, K1, alpha12)
    N2_isocline = dN2dt(N1, N2, r2, K2, alpha21)
    plt.figure(figsize=(10, 10))
    plt.contour(N1, N2, N1_isocline, levels=[0], colors='darkseagreen', linewidths=2) # label='dN1/dt = 0'
    plt.contour(N1, N2, N2_isocline, levels=[0], colors='darkred', linewidths=2) # label='dN2/dt = 0'
    punto_eq = punto_equilibrio(r1, r2, K1, K2, alpha12, alpha21)
    plt.plot(punto_eq[0], punto_eq[1], 'o', color = 'lightseagreen', markersize = 10, label='Punto de equilibrio')
    
    speed = np.sqrt(N1_isocline**2 + N2_isocline**2)    
    # Graficar el campo vectorial con un mapa de colores basado en la velocidad
    strm = plt.streamplot(N1, N2, N1_isocline, N2_isocline, color=speed, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    
    plt.xlabel('N1')
    plt.ylabel('N2')
    plt.xlim(0, K1)
    plt.ylim(0, K2)

    plt.legend()
    plt.title('Isoclinas: ' + title)
    
    # Agregar la barra de colores
    plt.colorbar(strm.lines, label='Velocidad')
    
    plt.tight_layout()
    plt.show()

    return N1, N2, N1_isocline, N2_isocline

def calcular_isoclinas_y_graficar_contour(title, r1, r2, K1, K2, alpha12, alpha21):
    N1 = np.linspace(0, K1, 1000)
    N2 = np.linspace(0, K2, 1000)
    N1, N2 = np.meshgrid(N1, N2)
    N1_isocline = dN1dt(N1, N2, r1, K1, alpha12)
    N2_isocline = dN2dt(N1, N2, r2, K2, alpha21)
    plt.figure(figsize=(10, 10))
    plt.contour(N1, N2, N1_isocline, levels=[0], colors='blue', linewidths=2) # label='dN1/dt = 0'
    plt.contour(N1, N2, N2_isocline, levels=[0], colors='red', linewidths=2) # label='dN2/dt = 0'
    plt.quiver(N1[::75, ::75], N2[::75, ::75], N1_isocline[::75, ::75], N2_isocline[::75, ::75], scale=85**2, label='Campo Vectorial', color='black', cmap = 'jet')
    punto_eq = punto_equilibrio(r1, r2, K1, K2, alpha12, alpha21)
    # agrandar el tamaño del punto de equilibrio
    plt.plot(punto_eq[0], punto_eq[1], 'o', color='violet', markersize=10, label='Punto de equilibrio')
    plt.xlabel('N1')
    plt.ylabel('N2')
    plt.xlim(0, K1)
    plt.ylim(0, K2)

    plt.legend()
    plt.title('Isoclinas: ' + title)
    plt.tight_layout()
    plt.show()

    return N1, N2, N1_isocline, N2_isocline


def graficar_soluciones_rk_separadas_informe(t0, N1_0, N2_0, tf, h, case):

        t_values, y_values = runge_kutta4_system(lotka_volterra, t0, [N1_0, N2_0], tf, h, case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        plt.figure(figsize=(10, 10))
        plt.plot(t_values, y_values[:, 0], label='N1(t)')
        plt.plot(t_values, y_values[:, 1], label='N2(t)')
        plt.xlabel('Tiempo')
        plt.ylabel('Población')
        plt.title(case['title'])
        plt.legend()
        plt.tight_layout()
        plt.show()
        
def calcular_isoclinas_y_graficar_contour_color_varios(cases):
    plt.figure(figsize=(15, 15))

    for i, case in enumerate(cases.values(), start=1):
        N1 = np.linspace(0, case['K1'], 1000)
        N2 = np.linspace(0, case['K2'], 1000)
        N1, N2 = np.meshgrid(N1, N2)
        N1_isocline = dN1dt(N1, N2, case['r1'], case['K1'], case['alpha12'])
        N2_isocline = dN2dt(N1, N2, case['r2'], case['K2'], case['alpha21'])
        
        plt.subplot(2, 2, i)
        plt.contour(N1, N2, N1_isocline, levels=[0], colors='blue', label='dN1/dt = 0')
        plt.contour(N1, N2, N2_isocline, levels=[0], colors='red', label='dN2/dt = 0')

        speed = np.sqrt(N1_isocline**2 + N2_isocline**2)    
        # Graficar el campo vectorial con un mapa de colores basado en la velocidad
        strm = plt.streamplot(N1, N2, N1_isocline, N2_isocline, color=speed, linewidth=2, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
       
        punto_eq = punto_equilibrio(case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        plt.plot(punto_eq[0], punto_eq[1], 'ro', label='Punto de equilibrio', markersize=10, color = 'violet')
        
        plt.xlabel('N1')
        plt.ylabel('N2')
        plt.xlim(0, case['K1'])
        plt.ylim(0, case['K2'])

        plt.title('Isoclinas: ' + case['title'])
        
        # Agregar la barra de colores
        plt.colorbar(strm.lines, label='Velocidad')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3) 
    plt.show()
        
def main():
    # graficar_soluciones_rk(t0, N1_0, N2_0, tf, h, cases)
    # calcular_isoclinas_y_graficar_contour_color_varios(cases)
    # # para el informe se puede ver separado en -> (después lo borramos)
    # for i, case in enumerate(cases.values(), start=1):
    #     graficar_soluciones_rk_separadas_informe(t0, N1_0, N2_0, tf, h, case)
    
    
    for i, case in enumerate(cases.values(), start=1):
        calcular_isoclinas_y_graficar_contour(case['title'], case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        calcular_isoclinas_y_graficar_contour_color(case['title'], case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])

if __name__ == '__main__':
    main()