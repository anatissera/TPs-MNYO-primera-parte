import numpy as np
import matplotlib.pyplot as plt
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

def punto_equilibrio(r1, r2, K1, K2, a12, a21, punto=[0, 0]):
    # Para encontrar los puntos de equilibrio, primero igualamos las derivadas a cero y resolvemos para N1 y N2.
    # Busco donde dN1/dt = 0 y dN2/dt = 0
    def f(x):
        return [r1 * x[0] * (K1 - x[0] - a12 * x[1]) / K1, r2 * x[1] * (K2 - x[1] - a21 * x[0]) / K2]
    
    return fsolve(f, punto)

# Parámetros del sistema    
h = 0.1
tf = 100
t0 = 0
N1_0 = 10
N2_0 = 10


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
        plt.xlabel('Tiempo', fontsize=15)
        plt.ylabel('Población', fontsize=15)
        plt.title(case['title'] + ': ' + case['case'], fontsize=16)
        plt.legend(fontsize = 15)
        
    plt.subplots_adjust(hspace=0.55, wspace=0.3)
    plt.show()
    
        
def isoclinas_cero(r1, r2, k1, k2, alpha12, alpha21, title, legend_loc):
    n1 = np.linspace(0, k1, 100)
    n2 = np.linspace(0, k2, 100)
    
    isocline1 = k1 - alpha12 * n2
    isocline2 = k2 - alpha21 * n1
    
    puntos_eq = calcular_todos_los_equililbrios(r1, r2, k1, k2, alpha12, alpha21)
    puntos_eq_x = [p[0] for p in puntos_eq]
    puntos_eq_y = [p[1] for p in puntos_eq]


    vn1 = np.linspace(0, k1, 50)
    vn2 = np.linspace(0, k2, 50)
    VN1, VN2 = np.meshgrid(vn1, vn2)
    
    dN1 = dN1dt(VN1, VN2, r1, k1, alpha12)
    dN2 = dN2dt(VN1, VN2, r2, k2, alpha21)
    magnitude = np.sqrt(dN1**2 + dN2**2)
    
    
    plt.figure()
    plt.plot(n1, isocline2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
    plt.plot(isocline1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
   
    plt.scatter(puntos_eq_x, puntos_eq_y, color='teal', s=100, label='Puntos de equilibrio', zorder=3)
    
    strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('N1', fontsize = 17)
    plt.ylabel('N2', fontsize = 17)
    plt.xlim(0, k1)
    plt.ylim(0, k2)

    plt.title('Isoclinas: ' + title, fontsize = 20)
    
    plt.legend(loc=legend_loc, fontsize=17, handlelength=0.75)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=18)
    
    plt.show()
    
def calcular_todos_los_equililbrios(r1, r2, k1, k2, alpha12, alpha21):
    n1 = np.linspace(0, k1, 100)
    n2 = np.linspace(0, k2, 100)
    
    puntos = []
    
    for i in range(0, 100):
        punto_eq = punto_equilibrio(r1, r2, k1, k2, alpha12, alpha21, [n1[i], n2[i]])
        if not es_punto_repetido(puntos, punto_eq, 1e-6):
            puntos.append(punto_eq)
    
    return puntos

def es_punto_repetido(lista_puntos, nuevo_punto, tolerancia):
    for punto in lista_puntos:
        if np.all(np.abs(np.array(punto) - np.array(nuevo_punto)) < tolerancia):
            return True
    return False

def casoab(r1, r2, k1, k2, alpha12, alpha21, title, legend_loc, punto_eq_inestable, punto_eq_estable):
    n1 = np.linspace(0, k1, 100)
    n2 = np.linspace(0, k2, 100)
    
    isocline1 = k1 - alpha12 * n2
    isocline2 = k2 - alpha21 * n1
    
    plt.scatter(punto_eq_inestable[0], punto_eq_inestable[1], edgecolor='teal', facecolor='lightcyan', s=100, zorder=3)
    plt.scatter(punto_eq_estable[0], punto_eq_estable[1], color='teal', s=100, label='Punto de equilibrio', zorder=3)


    vn1 = np.linspace(0, k1, 50)
    vn2 = np.linspace(0, k2, 50)
    VN1, VN2 = np.meshgrid(vn1, vn2)
    
    dN1 = dN1dt(VN1, VN2, r1, k1, alpha12)
    dN2 = dN2dt(VN1, VN2, r2, k2, alpha21)
    magnitude = np.sqrt(dN1**2 + dN2**2)
    
    plt.plot(n1, isocline2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
    plt.plot(isocline1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
   
    strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('N1', fontsize = 17)
    plt.ylabel('N2', fontsize = 17)
    plt.xlim(0, k1)
    plt.ylim(0, k2)

    plt.title('Isoclinas: ' + title, fontsize = 20)
    
    plt.legend(loc=legend_loc, fontsize=9, handlelength=0.75)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=12)
    
def casocd(r1, r2, k1, k2, alpha12, alpha21, title, legend_loc, punto_eq_inestable, punto_eq_estable, bool):
    n1 = np.linspace(0, k1, 100)
    n2 = np.linspace(0, k2, 100)
    
    isocline1 = k1 - alpha12 * n2
    isocline2 = k2 - alpha21 * n1
    
   
    puntos_eq_x = [p[0] for p in punto_eq_inestable]
    puntos_eq_y = [p[1] for p in punto_eq_inestable]
    
    if bool:
        plt.scatter(puntos_eq_x, puntos_eq_y,  color='teal', s=100, zorder=3)
        plt.scatter(punto_eq_estable[0], punto_eq_estable[1], edgecolor='teal', facecolor='lightcyan', label='Punto de equilibrio', s=100, zorder=3)
  
    else:
        plt.scatter(punto_eq_estable[0], punto_eq_estable[1], color='teal', s=100, label='Punto de equilibrio', zorder=3)
        


    vn1 = np.linspace(0, k1, 50)
    vn2 = np.linspace(0, k2, 50)
    VN1, VN2 = np.meshgrid(vn1, vn2)
    
    dN1 = dN1dt(VN1, VN2, r1, k1, alpha12)
    dN2 = dN2dt(VN1, VN2, r2, k2, alpha21)
    magnitude = np.sqrt(dN1**2 + dN2**2)
    
    
    plt.plot(n1, isocline2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
    plt.plot(isocline1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
   
    strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
    plt.grid()
    
    plt.xlabel('N1', fontsize = 17)
    plt.ylabel('N2', fontsize = 17)
    plt.xlim(0, k1)
    plt.ylim(0, k2)

    plt.title('Isoclinas: ' + title, fontsize = 20)
    
    plt.legend(loc=legend_loc, fontsize=9, handlelength=0.75)
    cbar = plt.colorbar(strm.lines)
    cbar.set_label(label='Magnitud del campo vectorial', fontsize=12)
    
def isoclinas_cero_y_graficar_varios_con_estabilidad(cases):
    
    plt.figure()
    
    plt.subplot(2, 2, 1)
    casoab(cases['a']['r1'], cases['a']['r2'], cases['a']['K1'], cases['a']['K2'], cases['a']['alpha12'], cases['a']['alpha21'], cases['a']['title'], cases['a']['legend_loc'], [7.79903684e-11, 3.60000000e+03], [4.200000e+03, 1.789256e-10])
  
    
    plt.subplot(2, 2, 2)
    casoab(cases['b']['r1'], cases['b']['r2'], cases['b']['K1'], cases['b']['K2'], cases['b']['alpha12'], cases['b']['alpha21'], cases['b']['title'], cases['b']['legend_loc'], [4.40000000e+03, 4.06958441e-10], [4.57001831e-11, 5.00000000e+03])
    
  
    plt.subplot(2, 2, 3)
    casocd(cases['c']['r1'], cases['c']['r2'], cases['c']['K1'], cases['c']['K2'], cases['c']['alpha12'], cases['c']['alpha21'], cases['c']['title'], cases['c']['legend_loc'], [[ 1.00000000e+03, -1.07384375e-11], [1.43044112e-11, 1.35000000e+03]], [527.27272727, 295.45454545], True)
    

    plt.subplot(2, 2, 4)
    casocd(cases['d']['r1'], cases['d']['r2'], cases['d']['K1'], cases['d']['K2'], cases['d']['alpha12'], cases['d']['alpha21'], cases['d']['title'], cases['d']['legend_loc'], [[1.60000000e+03, 2.00774581e-09], [4.00278193e-12, 1.50000000e+03]], [1133.33333333,  933.33333333], False)
   
    plt.subplots_adjust(hspace=0.55, wspace=0.26) 
    plt.show()
    
    

def isoclinas__cero_y_graficar_varios(cases):

    plt.figure()
    for i, case in enumerate(cases.values(), start=1):
        
        plt.subplot(2, 2, i)
        
        n1 = np.linspace(0, case['K1'], 100)
        n2 = np.linspace(0, case['K2'], 100)
        
        isocline_N1 = case['K1'] - case['alpha12'] * n2
        isocline_N2 = case['K2'] - case['alpha21'] * n1
        
        puntos_eq = calcular_todos_los_equililbrios(case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'])
        print(puntos_eq)
        puntos_eq_x = [p[0] for p in puntos_eq]
        puntos_eq_y = [p[1] for p in puntos_eq]

        vn1 = np.linspace(0, case['K1'], 50)
        vn2 = np.linspace(0, case['K2'], 50)
        VN1, VN2 = np.meshgrid(vn1, vn2)
        
        dN1 = dN1dt(VN1, VN2, case['r1'], case['K1'], case['alpha12'])
        dN2 = dN2dt(VN1, VN2, case['r2'], case['K2'], case['alpha21'])
        magnitude = np.sqrt(dN1**2 + dN2**2)
        
        plt.plot(n1, isocline_N2, label='dN2/dt = 0', color ='limegreen', linewidth=2)
        plt.plot(isocline_N1, n2, label='dN1/dt = 0', color = 'firebrick', linewidth=2)
    
        plt.scatter(puntos_eq_x, puntos_eq_y, color='teal', s=100, label='Puntos de equilibrio', zorder=3)
        
        strm = plt.streamplot(VN1, VN2, dN1, dN2, color= magnitude, linewidth=1, cmap='CMRmap', arrowstyle='->', arrowsize=1.5)
        cbar = plt.colorbar(strm.lines)
        cbar.set_label(label='Magnitud del campo vectorial', fontsize=10)
        
        plt.xlabel('N1', fontsize = 14)
        plt.ylabel('N2', fontsize = 14)
        plt.xlim(0, case['K1'])
        plt.ylim(0, case['K2'])

        plt.title('Isoclinas: ' + case['title'], fontsize = 20)
        
        plt.legend(loc=case['legend_loc'], handlelength=0.75)
      
    plt.subplots_adjust(hspace=0.4) 
    plt.show()

# exclusión competitiva
# Caso a: 
# K1 > K2/α21 & K2 < K1/α12

# Caso b: 
# K1 < K2/α21 & K2 > K1 /α12

# Caso c: Dominancia indeterminada
# K1 < K2 * α12 & K2 < K1 * α21

# Caso d: Coexistencia 
# K1 > K2 * α12 & K2 > K1 * α21

cases = {
    'a': {'r1': 0.1, 'r2': 0.1, 'K1': 4200, 'K2': 3600, 'alpha12': 0.35, 'alpha21': 2.2, 'title': 'Caso a', 'case': 'Exclusión competitiva por parte de N1', 'legend_loc': 'upper center'},
    'b': {'r1': 0.1, 'r2': 0.1, 'K1': 4400, 'K2': 5000, 'alpha12': 2.1, 'alpha21': 0.7, 'title': 'Caso b', 'case': 'Exclusión competitiva por parte de N2', 'legend_loc': 'center right'},
    'c': {'r1': 0.2, 'r2': 0.2, 'K1': 1000, 'K2': 1350, 'alpha12': 1.6, 'alpha21': 2, 'title': 'Caso c', 'case': 'Dominancia indeterminada', 'legend_loc': 'upper right'},
    'd': {'r1': 0.4, 'r2': 0.8, 'K1': 1600, 'K2': 1500, 'alpha12': 0.5, 'alpha21': 0.5, 'title': 'Caso d', 'case': 'Coexistencia', 'legend_loc': 'lower left'}
}

def main():
   
    graficar_soluciones_rk_varias(t0, N1_0, N2_0, tf, h, cases)
    isoclinas_cero_y_graficar_varios_con_estabilidad(cases)
    isoclinas__cero_y_graficar_varios(cases)
    
    # para el informe se puede ver separado en -> (después lo borramos)
    for i, case in enumerate(cases.values(), start=1):
        graficar_soluciones_rk_separadas_informe(t0, N1_0, N2_0, tf, h, case)
    
    
    for i, case in enumerate(cases.values(), start=1):
        isoclinas_cero(case['r1'], case['r2'], case['K1'], case['K2'], case['alpha12'], case['alpha21'], case['title'], case['legend_loc'])
        
        
if __name__ == '__main__':
    main()