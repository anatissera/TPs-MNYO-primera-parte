import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from punto2 import runge_kutta4_system
from scipy.integrate import odeint
from punto1 import runge_kutta_4

# def dN_dt (N, P, r, alpha):
#     return r*N - alpha*N*P

# def dP_dt (N, P, beta, q):
#     return beta*N*P - q*P

def dN_dt (N, P, r, alpha, K):
    return r * N * (1 - N / K) - alpha * N * P

def dP_dt (N, P, beta, q):
    return beta * N * P - q * P

def lotka_volterra(t, y0, r1, r2, K1, alpha, beta, q):
    N, P = y0
    return np.array([dN_dt(N, P, r1, alpha, K1), dP_dt(N, P, beta, q)])


# # sirve el mismo runge kutta 4

# h = 0.1
# tf = 200
# t0 = 0
# N0_conejos = 2000
# P0_zorros = 10

# r1 = 0.1  # Tasa de crecimiento de los conejos (conejos/mes)
# alpha = 0.005  # Eficiencia de captura (conejos*zorros/mes)
# r2 = 0.04  # Tasa de crecimiento de los zorros (zorros/mes)
# beta = 0.00004  # Eficiencia para convertir presas en nuevos depredadores (conejos*zorros/mes)
# q = 0.05  # Tasa de mortalidad per cápita de los zorros (zorros/mes)
# K1 = 10000  # Capacidad de carga del ambiente para los conejos


def punto_equilibrio(r, K, alpha, beta, q):
    # Para encontrar los puntos de equilibrio, primero igualamos las derivadas a cero y resolvemos para N y P.
    # Busco donde dN/dt = 0 y dP/dt = 0
    
    # x[0] = N, x[1] = P
    def f(x):
        return [r * x[0] * (1 - x[0] / K ) - alpha * x[0] * x[1], beta * x[0] * x[1] - q * x[1]]
    
    return fsolve(f, [K/2, K/2]) # punto de convergencia

def isoclinas_y_campo_vectorial(r, K, alpha, beta, q, title, legend_loc):
    N = np.linspace(0, K, 100)
    P = np.linspace(0, K, 100)
    
    isocline_N = q / beta * np.ones_like(N)
    isocline_P = (r * (1 - N / K)) / alpha * np.ones_like(P)
    
    punto_eq = punto_equilibrio(r, K, alpha, beta, q)
    
    VN, VP = np.meshgrid(N, P)
    
    dN = r * VN * (1 - VN / K) - alpha * VN * VP
    dP = beta * VN * VP - q * VP
    magnitude = np.sqrt(dN**2 + dP**2)
    
    plt.figure()
    plt.plot(isocline_N, P, label='dN/dt = 0', color='limegreen', linewidth=2)
    plt.plot(N, isocline_P, label='dP/dt = 0', color='firebrick', linewidth=2)
    plt.plot(punto_eq[0], punto_eq[1], 'o', color='teal', markersize=10, label='Punto de equilibrio')
    
    strm = plt.streamplot(VN, VP, dN, dP, color=magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('Población de Presas (N)', fontsize=17)
    plt.ylabel('Población de Depredadores (P)', fontsize=17)
    plt.title('Isoclinas y Campo Vectorial: ' + title, fontsize=20)
    
    plt.legend(loc=legend_loc, fontsize=17, handlelength=0.75)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=18)
    
    plt.show()

def isoclinas(r, K, alpha, beta, q, title, legend_loc):
    N = np.linspace(0, K, 100)
    P = np.linspace(0, K, 100)
    
    isocline_N = q / beta * np.ones_like(N)
    isocline_P = (r * (1 - N / K)) / alpha * np.ones_like(P)
    
    punto_eq = punto_equilibrio(r, K, alpha, beta, q)
    
    VN, VP = np.meshgrid(N, P)
    
    dN = dN_dt(VN, VP, r, alpha, K)
    dP = dP_dt(VN, VP, beta, q)
    magnitude = np.sqrt(dN**2 + dP**2)
    
    plt.figure()
    plt.plot(isocline_N, P, label='dN/dt = 0', color='limegreen', linewidth=2)
    plt.plot(N, isocline_P, label='dP/dt = 0', color='firebrick', linewidth=2)
    plt.plot(punto_eq[0], punto_eq[1], 'o', color='teal', markersize=10, label='Punto de equilibrio')
    
    strm = plt.streamplot(VN, VP, dN, dP, color=magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('Población de Presas (N)', fontsize=17)
    plt.ylabel('Población de Depredadores (P)', fontsize=17)
    plt.title('Isoclinas y Campo Vectorial: ' + title, fontsize=20)
    
    plt.legend(loc=legend_loc, fontsize=17, handlelength=0.75)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=18)
    
    plt.show()

def graficar_soluciones_rk_varias(t0, N1_0, N2_0, tf, h, cases):
    plt.figure(figsize=(10, 10))

    for i, case in enumerate(cases.values(), start=1):
        t_values, y_values = runge_kutta4_system(lotka_volterra, t0, [N1_0, N2_0], tf, h, case['r1'], case['r2'], case['K1'], case['alpha12'], case['beta'], case['q'])
        plt.subplot(2, 2, i)        
        plt.plot(t_values, y_values[:, 0], label='N1(t)')
        plt.plot(t_values, y_values[:, 1], label='N2(t)')
        plt.xlabel('Tiempo', fontsize=15)
        plt.ylabel('Población', fontsize=15)
        plt.title(case['title'] + ': ' + case['case'], fontsize=16)
        plt.legend(fontsize = 15)
        
    plt.subplots_adjust(hspace=0.55, wspace=0.3)
    plt.show()

def graficar_sol_rk (t0, y0, tf, h, r, K, alpha, beta, q, title):
    t_values, y_values = runge_kutta4_system(lotka_volterra, t0, y0, tf, h, r, r, K, alpha, beta, q)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, y_values[:, 0], label='presas')
    plt.plot(t_values, y_values[:, 1], label='depredadores')
    plt.xlabel('Tiempo (meses)')
    plt.ylabel('Población')
    plt.title('Dinámica de la población de presas y depredadores (' + title +')')
    plt.legend()
    plt.grid(True)
    plt.show()

 
def lotka_q_r_K_constant(y, t, alpha, beta): # r, q = 1, K = infinito
    yp = (1 - alpha * y[1]) * y[0] # r * (1 - N / K) está siendo = 1
    yp = np.append(yp, (-1 + beta * y[0]) * y[1]) # q = 1
    return yp

def plano_de_fases(alpha, beta, y0_values, r, K, q):
 
    t = np.linspace(0, 15, 1000)

    plt.figure(figsize=(8, 6))
    for y0 in y0_values:
        y = odeint(lotka_q_r_K_constant, y0, t, args=(alpha, beta))
        plt.plot(y[:, 0], y[:, 1])

    plt.title('Plano de Fases: Poblaciones de Presas vs. Depredadores')
    plt.xlabel('Población de Presas')
    plt.ylabel('Población de Depredadores')
    plt.grid(True)
    plt.show()
 
h = 0.1
tf = 300
t0 = 0
   
r = 1 # Tasa de crecimiento de las presas
alpha = 0.05 # éxito en la caza, afecta al crecimiento de las presas
beta = 0.01 # éxito en la caza, afecta al crecimiento de los depredadores
q = 0.1 # tasa de mortalidad de los depredadores
K = 1000 # Capacidad de carga del ambiente para las presas
N0 = 100
P0 = 10
y0 = [N0, P0]

# caso a: ciclo
# caso b: depredadores extinguen a las presas
# # r > q and alpha > beta
# # r < q and alpha < beta
# # r > q and alpha < beta
# # r < q and alpha > beta

# 

cases = {
    'a': {'r': 1, 'alpha': 0.05, 'K': 1000, 'beta': 0.01, 'q': 0.1, 'title': 'Caso a', 'legend_loc': 'upper center'},
    'b': {'r': 0.1, 'alpha': 0.09, 'K': 1000, 'beta': 0.004, 'q': 0.005, 'title': 'Caso b', 'legend_loc': 'center right'},
    'c': {'r': 2, 'alpha': 0.005, 'K': 5000, 'beta': 0.5, 'q': 0.05, 'title': 'Caso c', 'legend_loc': 'lower left'},
    'd': {'r': 0.1, 'alpha': 0.02, 'K': 1000, 'beta': 0.005, 'q': 0.05, 'title': 'Caso d', 'legend_loc': 'upper right'},
    'e': {'r': 0.01, 'alpha': 0.5, 'K': 1000, 'beta': 0.0005, 'q': 0.0005, 'title': 'Caso d', 'legend_loc': 'upper right'}
}

alpha_1 = 0.01
beta_1 = 0.02
y0_values = ([1, 1], [10, 10], [20, 20], [30, 30], [40, 40])

def main():
    graficar_sol_rk(t0, y0, tf, h, r, K, alpha, beta, q, 'Modelo Lotka-Volterra Presa-Depredador')
    isoclinas_y_campo_vectorial(r, K, alpha, beta, q, 'Modelo Lotka-Volterra Presa-Depredador', 'upper right')
    
    plano_de_fases(alpha_1, beta_1, y0_values ) # el centro de esto es el punto de equilibrio
    isoclinas_y_campo_vectorial(r, K, alpha_1, beta_1, q, 'Modelo Lotka-Volterra Presa-Depredador', 'upper right')
    
    for i, case in enumerate(cases.values(), start=1):
        isoclinas_y_campo_vectorial(case['r'], case['K'], case['alpha'], case['beta'], case['q'], 'Modelo Lotka-Volterra Presa-Depredador ' + case['title'], case['legend_loc'])
        graficar_sol_rk(t0, y0, tf, h, case['r'], case['K'], case['alpha'], case['beta'], case['q'], case['title'])
        isoclinas_y_campo_vectorial(case['r'], case['K'], case['alpha'], case['beta'], case['q'], 'Modelo Lotka-Volterra Presa-Depredador ' + case['title'], 'upper right')
        
if __name__ == '__main__':
    main()