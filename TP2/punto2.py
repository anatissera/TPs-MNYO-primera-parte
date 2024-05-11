import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.lines import Line2D

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

# Parámetros del sistema    
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

def graficar_soluciones_rk_varias(t0, N1_0, N2_0, tf, h, cases):
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

    # plt.tight_layout()
    plt.subplots_adjust(hspace=0.4)
    plt.show()

    
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
        
def isoclinas_cero(r1, r2, k1, k2, alpha12, alpha21, title, legend_loc):
    n1 = np.linspace(0, k1, 100)
    n2 = np.linspace(0, k2, 100)
    
    isocline1 = k1 - alpha12 * n2
    isocline2 = k2 - alpha21 * n1
    
    punto_eq = punto_equilibrio(r1, r2, k1, k2, alpha12, alpha21)
    
    vn1 = np.linspace(0, k1, 50)
    vn2 = np.linspace(0, k2, 50)
    VN1, VN2 = np.meshgrid(vn1, vn2)
    
    dN1 = dN1dt(VN1, VN2, r1, k1, alpha12)
    dN2 = dN2dt(VN1, VN2, r2, k2, alpha21)
    magnitude = np.sqrt(dN1**2 + dN2**2)
    
    plt.figure()
    plt.plot(n1, isocline2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
    plt.plot(isocline1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
    plt.plot(punto_eq[0], punto_eq[1], 'o', color='teal', markersize=10, label='Punto de equilibrio')
    
    strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('N1', fontsize = 17)
    plt.ylabel('N2', fontsize = 17)
    plt.xlim(0, k1)
    plt.ylim(0, k2)

    plt.title('Isoclinas: ' + title, fontsize = 20)
    
    plt.legend(loc=legend_loc)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=14)
    
    plt.show()
        
def isoclinas__cero_y_graficar_varios(cases):

    plt.figure()
    for i, case in enumerate(cases.values(), start=1):
        
        plt.subplot(2, 2, i)
        
        n1 = np.linspace(0, case['K1'], 100)
        n2 = np.linspace(0, case['K2'], 100)
        
        isocline1 = case['K1'] - case['alpha12'] * n2
        isocline2 = case['K2'] - case['alpha21'] * n1
        
        punto_eq = punto_equilibrio(case['K1'], case['r2'], case['K1'], case['r2'], case['alpha12'], case['alpha21'])
        
        vn1 = np.linspace(0, case['K1'], 50)
        vn2 = np.linspace(0, case['K2'], 50)
        VN1, VN2 = np.meshgrid(vn1, vn2)
        
        dN1 = dN1dt(VN1, VN2, case['r1'], case['K1'], case['alpha12'])
        dN2 = dN2dt(VN1, VN2, case['r2'], case['K2'], case['alpha21'])
        magnitude = np.sqrt(dN1**2 + dN2**2)
        
        plt.plot(n1, isocline2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
        plt.plot(isocline1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
        plt.plot(punto_eq[0], punto_eq[1], 'o', color='teal', markersize=10, label='Punto de equilibrio')
        
        strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
        cbar = plt.colorbar(strm.lines)
        cbar.set_label(label='Magnitud del campo vectorial', fontsize=10)
        
 
        plt.xlabel('N1', fontsize = 14)
        plt.ylabel('N2', fontsize = 14)
        plt.xlim(0, case['K1'])
        plt.ylim(0, case['K2'])

        plt.title('Isoclinas: ' + case['title'], fontsize = 20)
        
        plt.legend(loc=case['legend_loc'])
      
    plt.subplots_adjust(hspace=0.4) 
    plt.show()

cases = {
    'a': {'r1': 0.1, 'r2': 0.1, 'K1': 4000, 'K2': 3800, 'alpha12': 0.3, 'alpha21': 2, 'title': 'Caso a', 'case': 'N1 sobrevive (es más fuerte), N2 se extingue', 'legend_loc': 'upper center'},
    'b': {'r1': 0.1, 'r2': 0.1, 'K1': 4500, 'K2': 5000, 'alpha12': 2, 'alpha21': 0.3, 'title': 'Caso b', 'case': 'N2 sobrevive, N1 se extingue', 'legend_loc': 'center right'},
    'c': {'r1': 0.3, 'r2': 0.6, 'K1': 1500, 'K2': 1400, 'alpha12': 0.4, 'alpha21': 0.4, 'title': 'Caso c', 'case': 'Equilibrio estable entre N1 y N2', 'legend_loc': 'lower left'},
    'd': {'r1': 0.2, 'r2': 0.2, 'K1': 1000, 'K2': 1400, 'alpha12': 1.7, 'alpha21': 2, 'title': 'Caso d', 'case': 'Equilibrio inestable entre N1 y N2', 'legend_loc': 'upper right'}
}

def main():
    graficar_soluciones_rk_varias(t0, N1_0, N2_0, tf, h, cases)
    isoclinas__cero_y_graficar_varios(cases)
    
    # para el informe se puede ver separado en -> (después lo borramos)
    for i, case in enumerate(cases.values(), start=1):
        graficar_soluciones_rk_separadas_informe(t0, N1_0, N2_0, tf, h, case)
    
    
    for i, case in enumerate(cases.values(), start=1):
        # calcular_isoclinas_y_graficar_contour(case['title'], case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        isoclinas_cero(case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'], case['title'], case['legend_loc'])
        
        
if __name__ == '__main__':
    main()