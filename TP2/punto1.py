import numpy as np
import matplotlib.pyplot as plt

# introducción, el código de rk4 y euler. 
# método elegido, el h que eleigmos y por qué


# definiciones de consigna
# N(t) =

def exponential_solution(t, N0, h):
    """
    Calcula la solución exponencial de la ecuación diferencial de crecimiento poblacional.

    Parámetros:
    t: arreglo unidimensional de tiempo
    N0: población inicial en t=0
    r: tasa de crecimiento intrínseca de la población
   

    Retorna:
    N: arreglo unidimensional de tamaño poblacional en cada instante de tiempo t
    """
    return N0 * np.exp(h * t)

def exponential_dNdt(N0, h):
    dNdt = h * N0
    return dNdt

def logistic_solution(t, N0, h, K):
    """
    Mismo que exponencial pero calcula la solución logística
    Nuevo parámetro:
    K: capacidad de carga de la población 
    """
    exponent = -h * t
    return K / (1 + (K / N0 - 1) * np.exp(exponent)) #sol chequeada en internet

def logistic_dNdt(N0, h, K):
    dNdt = h * N0 * (1-N0/K)
    return dNdt


def euler_method(f, t_span, y0, N):
    a, b = t_span
    h = (b - a) / N
    t_values = [a]
    y_values = [y0]

    for i in range(1, N + 1):
        t = a + i * h
        t_values.append(t)
        y = y_values[-1] + h * f(t_values[-1], y_values[-1])
        y_values.append(y)

    return t_values[:-1], y_values[:-1]

def runge_kutta_4(f, t_span, y0, N):
    a, b = t_span
    h = (b - a) / N
    t_values = [a]
    y_values = [y0]

    for i in range(1, N + 1):
        t = a + i * h
        t_values.append(t)

        K1 = h * f(t_values[i - 1], y_values[i - 1])
        K2 = h * f(t_values[i - 1] + h / 2, y_values[i - 1] + K1 / 2)
        K3 = h * f(t_values[i - 1] + h / 2, y_values[i - 1] + K2 / 2)
        K4 = h * f(t_values[i - 1] + h, y_values[i - 1] + K3)

        y = y_values[i - 1] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        y_values.append(y)

    return t_values[:-1], y_values[:-1]   # el [:-1] es para que ambos tengan la misma dimensión

# consigna
# 1. Obtener la solución exacta de la ecuación logística y la ecuación exponencial. 
# 2. Con métodos de Euler y Runge-Kutta de orden 4 para aproximar las soluciones de ambas ecuaciones. 
# 3. Graficar el tamaño poblacional en función del tiempo y la variación poblacional en función del 
# tamaño poblacional para ambos modelos, utilizando diferentes valores de los parámetros 𝑁0, 𝑟 y 𝐾 
# 4. Comparar las soluciones numéricas con las soluciones exactas (error).


# ¿Qué es K?  La forma más simple de f(N) es la linea recta. Cuando la densidad poblacional es baja, 
# la tasa de crecimiento per cápita es similar a r (tasa intrínseca de crecimiento). 
# Pero cercanos a un determinado tamaño poblacional, f(N) se hace 0 y no hay más crecimiento poblacional. 
# Esta tamaño poblacional se llama capacidad de carga de la población y se nota como K.

# Propiedades

# 1) La tasa de crecimiento per cápita no es aquí constante, sino que disminuye con la densidad 
# dN/dt 1/N= r(K-N)/N.

# 2) La curva logística difiere de la curva geométrica en dos puntos: 
# tiene una asíntota superior, y se acerca a esta asíntota suavemente, no bruscamente.

# 3) La curva predice un equilibrio dinámico estable de la población cuando N=K.

# 4) Hay dos atributos de la curva logística que la hacen muy atractiva: 
# su simplicidad matemática, y su aparente realidad. Sólo contiene dos constantes, K y r.

# 5) La curva es simétrica respecto a su punto central= K/2.


# hacer gráfico que demuestre qué son los parámetros. cambio k, la población max va a ser distinto. 
# en un eje un parámetro k y en el otro cómo cambia la población
# diagrama de Fases


def calculate_exact_solutions(t, N0, h, K):
    N_exact_logistic = logistic_solution(t, N0, h, K)
    N_exact_exponential = exponential_solution(t, N0, h)
    return N_exact_logistic, N_exact_exponential
    
def plot_solutions_exact(t, N0, r, K, title):
    N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, N0, r, K)
   
    plt.figure(figsize=(10, 5))

    plt.plot(t, N_exact_exponential, label='Solución Exponencial', color='mediumseagreen', linewidth=2.5)
    plt.plot(t, N_exact_logistic, label='Solución Logística', color='lightcoral', linewidth=2.5)

    plt.xlabel('Tiempo', fontsize = 17.5)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 17.5)
    plt.title(title, fontsize = 20)
    plt.legend(fontsize = 17.5, handlelength=1.1)
    plt.grid(True)

    plt.show()

def plot_solutions_exact_varios(variables):
    plt.figure(figsize=(12, 18))
    
    for i, (title, casos) in enumerate(variables.items(), start=1):
        plt.subplot(3, 2, i)
        plt.title(title + ' - Solución Logística')
        t = np.linspace(0, 100, 1000)
         
        for caso in casos:
           
            N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, caso['N0'], caso['r'], caso['K'])
            plt.plot(t, N_exact_logistic, label=f"N0={caso['N0']}, r={caso['r']}, K={caso['K']}")
        
        plt.xlabel('Tiempo')
        plt.ylabel('Tamaño Poblacional (N)')
        plt.legend()
        plt.grid(True)
    
    for i, (title, casos) in enumerate(variables.items(), start=4):
        plt.subplot(3, 2, i)
        plt.title(title + ' - Solución Exponencial')
        
        for caso in casos:
        
            N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, caso['N0'], caso['r'], caso['K'])
            plt.plot(t, N_exact_exponential, label=f"N0={caso['N0']}, r={caso['r']}, K={caso['K']}")
        
        plt.xlabel('Tiempo')
        plt.ylabel('Tamaño Poblacional (N)')
        plt.legend()
        plt.grid(True)
    
    plt.subplots_adjust(hspace=0.8)
    plt.show()
    

def plot_solutions_exact_pares(casos, title):
    fig = plt.figure(figsize=(12, 15))  
        
    for i, caso in enumerate(casos, start=1):
        if title == 'N0':
            label_t = f"N0={caso['N0']}"
            label_generic = f"r={caso['r']}, K={caso['K']}"
        elif title == 'r':
            label_t = f"r={caso['r']}"
            label_generic = f"N0={caso['N0']}, K={caso['K']}"
        else:
            label_t = f"K={caso['K']}"
            label_generic = f"N0={caso['N0']}, r={caso['r']}"


        plt.subplot(2, 2, 1)
        plt.title('Solución Logística', fontsize=15)
    
        t = np.linspace(0, caso['time'], 1000)
        N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, caso['N0'], caso['r'], caso['K'])
        plt.plot(t, N_exact_logistic, label=label_t)
    
        plt.xlabel('Tiempo', fontsize=13)
        plt.ylabel('Tamaño Poblacional (N)', fontsize=13)
        plt.legend(fontsize = 13)
        plt.grid(True)
        
        plt.subplot(2, 2, 2)
        plt.title('Solución Exponencial', fontsize=15)
        
        plt.plot(t, N_exact_exponential, label=label_t)
    
        plt.xlabel('Tiempo', fontsize=13)
        plt.ylabel('Tamaño Poblacional (N)', fontsize=13)
        plt.legend(fontsize = 13)
        plt.grid(True)

        plt.subplot(2, 2, 3)
        plt.title('Variación Logística', fontsize=15)
        dN_exact_logistic = logistic_dNdt(N_exact_logistic, h, caso['K'])
        plt.plot(N_exact_logistic, dN_exact_logistic, label=label_t)
         
        plt.xlabel('Tamaño Poblacional (N)', fontsize=12)
        plt.ylabel('Variación Poblacional (dN/dt)', fontsize=12)
        plt.legend(fontsize = 13)
        plt.grid(True)
        
        plt.subplot(2, 2, 4)
        plt.title('Variación Exponencial', fontsize=15)
        dN_exact_exponential = exponential_dNdt(N_exact_exponential, h)
        plt.plot(N_exact_exponential, dN_exact_exponential, label=label_t)
        
        plt.xlabel('Tamaño Poblacional (N)', fontsize=12)
        plt.ylabel('Variación Poblacional (dN/dt)', fontsize=12)
        plt.legend(fontsize = 13)
        plt.grid(True)
    
    
    plt.suptitle(f"varío {title}  -  {label_generic}", fontsize=18)
    plt.subplots_adjust(hspace=0.5,  wspace=0.3)

    plt.show()

    
def plot_population_variation_exact(variables):
    plt.figure(figsize=(12, 10))  
    
    for i, (title, casos) in enumerate(variables.items(), start=1):
        plt.subplot(3, 1, i)
        plt.title(title)
        
        for caso in casos:
            t = np.linspace(0, caso['time'], 1000)
            N_exact_logistic, _ = calculate_exact_solutions(t, caso['N0'], caso['r'], caso['K'])
            dN_exact_logistic = logistic_dNdt(N_exact_logistic, h, caso['K'])
            
            plt.plot(N_exact_logistic, dN_exact_logistic, label=f"N0={caso['N0']}, r={caso['r']}, K={caso['K']}")
        
        plt.xlabel('Tamaño Poblacional (N)')
        plt.ylabel('Variación Poblacional (dN/dt)')
        plt.legend()
        plt.grid(True)

    plt.subplots_adjust(hspace=0.8)
    plt.show()


def plot_solutions_numerical(t, N0, h, K, time, space, title):
    plt.figure(figsize=(10, 5))
    N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, N0, h, K)

    _, N_logistic_euler = euler_method(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)  
    _, N_logistic_rk4 = runge_kutta_4(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)
    _, N_exponential_euler = euler_method(lambda t, N: exponential_dNdt(N, h), (0, time), N0, space) 
    _, N_exponential_rk4 = runge_kutta_4(lambda t, N: exponential_dNdt(N, h) , (0, time), N0, space)

    plt.subplot(1, 2, 1)
    
    plt.plot(t, N_exact_logistic, linestyle ='dashed', label='Exacta', c= 'darkred')
    plt.plot(t, N_logistic_euler, label='Euler', color='mediumseagreen')
    plt.plot(t, N_logistic_rk4, label='RK4', color='rebeccapurple')
    
    plt.xlabel('Tiempo', fontsize = 18)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 18)
    plt.title('Crecimiento Logístico', fontsize = 20)
    plt.legend(fontsize = 17, handlelength=1)
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(t, N_exact_exponential, linestyle ='dashed', label='Exacta', c='darkgreen')
    plt.plot(t, N_exponential_euler, label='Euler', color='lightcoral')
    plt.plot(t, N_exponential_rk4, label='RK4', color='skyblue')

    plt.xlabel('Tiempo', fontsize = 18)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 18)
    plt.title('Crecimiento Exponencial', fontsize = 20)
    plt.legend(fontsize = 17, handlelength=1)
    plt.grid(True)

    plt.subplots_adjust(wspace=0.3)
    plt.show()
    
def plot_error_relativo (t, N0, h, K, time, space):
    N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, N0, h, K)

    _, N_logistic_euler = euler_method(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)  
    _, N_logistic_rk4 = runge_kutta_4(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)
    _, N_exponential_euler = euler_method(lambda t, N: exponential_dNdt(N, h), (0, time), N0, space) 
    _, N_exponential_rk4 = runge_kutta_4(lambda t, N: exponential_dNdt(N, h) , (0, time), N0, space)

    error_r_euler = np.abs(N_exact_exponential - N_exponential_euler) / N_exact_exponential
    error_r_rk4 = np.abs(N_exact_logistic - N_logistic_rk4) / N_exact_logistic
    error_r_euler_log = np.abs(N_exact_logistic - N_logistic_euler) / N_exact_logistic
    error_r_rk4_exp = np.abs(N_exact_exponential - N_exponential_rk4) / N_exact_exponential
    
    plt.figure(figsize=(10, 5))
    plt.plot(t, error_r_euler, label='ER Euler Exponencial', color='skyblue', linewidth=2.5)
    plt.plot(t, error_r_rk4, label='ER RK4 Logística', color='rebeccapurple', linewidth=2.5)
    plt.plot(t, error_r_euler_log, label='ER Euler Logística', color='mediumseagreen', linewidth=2.5)
    plt.plot(t, error_r_rk4_exp, label='ER RK4 Exponencial', color='lightcoral', linewidth=2.5)
    
    
    plt.xlabel('Tiempo', fontsize = 18)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 18)
    plt.title('Error Relativo de los Métodos Numéricos', fontsize = 20)
    plt.legend(fontsize = 18, handlelength=1)
    plt.grid(True)

    plt.show()

def plot_error_abs (t, N0, h, K, time, space):
    
    N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, N0, h, K)

    _, N_logistic_euler = euler_method(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)  
    _, N_logistic_rk4 = runge_kutta_4(lambda t, N: logistic_dNdt(N, h, K), (0, time), N0, space)
    _, N_exponential_euler = euler_method(lambda t, N: exponential_dNdt(N, h), (0, time), N0, space) 
    _, N_exponential_rk4 = runge_kutta_4(lambda t, N: exponential_dNdt(N, h) , (0, time), N0, space)

    error_a_euler = np.abs(N_exact_exponential - N_exponential_euler)
    error_a_rk4 = np.abs(N_exact_logistic - N_logistic_rk4) 
    error_a_euler_log = np.abs(N_exact_logistic - N_logistic_euler)
    error_a_rk4_exp = np.abs(N_exact_exponential - N_exponential_rk4)
    
    plt.figure(figsize=(10, 7))
    plt.subplot(1, 2, 1)
    
    plt.plot(t, error_a_euler, label='EA Euler Exponencial', color='skyblue', linewidth=2.5)
    plt.plot(t, error_a_rk4_exp, label='EA RK4 Exponencial', color='lightcoral', linewidth=2.5)
    plt.xlabel('Tiempo', fontsize = 16)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 16)
    plt.title('Error Absoluto de los Métodos Numéricos', fontsize = 18)
    plt.legend(fontsize = 16, handlelength=0.75)
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(t, error_a_euler_log, label='EA Euler Logística', color='mediumseagreen', linewidth=2.5)
    plt.plot(t, error_a_rk4, label='EA RK4 Logística', color='rebeccapurple', linewidth=2.5)
    plt.xlabel('Tiempo', fontsize = 16)
    plt.ylabel('Tamaño Poblacional (N)', fontsize = 16)
    plt.title('Error Absoluto de los Métodos Numéricos', fontsize = 18)
    plt.legend(fontsize = 15, handlelength=0.75)
    plt.grid(True)

    plt.subplots_adjust(wspace=0.3)
    plt.show()
  
def plot_population_variation(t, N0, h, K, title):
    
    N_exact_logistic, N_exact_exponential = calculate_exact_solutions(t, N0, h, K)
    dN_exact_logistic = logistic_dNdt(N_exact_logistic, h, K)
    dN_exact_exponential = exponential_dNdt(N_exact_exponential, h)
    
    plt.figure(figsize=(10, 5))

    if title == 'Variación Exponencial':
        plt.plot(N_exact_exponential, dN_exact_exponential, label = title, color = 'mediumseagreen')
    else:
        plt.plot(N_exact_logistic, dN_exact_logistic,  label = title, color = 'lightcoral')

    plt.xlabel('Tamaño Poblacional (N)')
    plt.ylabel('Variación Poblacional (dN/dt)')
    plt.title(title)
    plt.grid(True)

    plt.show()


N0 = 10
h = 0.1 
K_values = [100, 200, 100000] 


t1 = np.linspace(0, 200, 1000)
t2 = np.linspace(0, 250, 1000)

variables = { 'K': [{'N0': 10, 'r': 0.1, 'K': 100, 'time': 90, 'space': 1000}, {'N0': 10, 'r': 0.1, 'K': 150, 'time': 90, 'space': 1000}, {'N0': 10, 'r': 0.1, 'K': 300, 'time': 90, 'space': 1000} ], 
             'N0': [{'N0': 10, 'r': 0.1, 'K': 150, 'time': 100, 'space': 1000}, {'N0': 50, 'r': 0.1, 'K': 150, 'time': 100, 'space': 1000}, {'N0': 100, 'r': 0.1, 'K': 150, 'time': 100, 'space': 1000}], 
             'r': [{'N0': 50, 'r': -0.2, 'K': 150, 'time': 25, 'space': 1000}, {'N0': 50, 'r': 0.5, 'K': 150, 'time': 25, 'space': 1000}, {'N0': 50, 'r': 1, 'K': 150, 'time': 25, 'space': 1000}]
             }

def main():
    plot_solutions_exact(np.linspace(0, 110, 1000), N0, h, K_values[2], 'Soluciones Exactas')
    
    plot_solutions_exact_pares(variables['N0'], 'N0')
    plot_solutions_exact_pares(variables['r'], 'r')
    plot_solutions_exact_pares(variables['K'], 'K')
    
    plot_solutions_exact_varios(variables)
    
    plot_solutions_numerical(t1, N0, h, K_values[2], 200, 1000, 'Soluciones Numéricas')
    plot_error_relativo(t1, N0, h, K_values[2], 200, 1000)
    plot_error_abs(t1, N0, h, K_values[2], 200, 1000)

    # plot_population_variation(t2, N0, h, K_values[2], 'Variación Exponencial')
    # plot_population_variation(t2, N0, h, K_values[2], 'Variación Logística')
    
    
if __name__ == '__main__':
    main()