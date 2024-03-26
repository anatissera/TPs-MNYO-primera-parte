from punto1a import equispaced_points, chebyshev_points
import matplotlib.pyplot as plt # para el a
from mpl_toolkits.mplot3d import Axes3D # para el b
from scipy.interpolate import RectBivariateSpline
import numpy as np
from scipy.interpolate import griddata


# sp.interpolate.RegularGridInterpolator con método lineal, 5x5, 25x25 puntos.
# Observamos que a medida que aumentamos la cantidad de
# puntos equidistantes el error disminuye. Los resultados tambien muestran que, al menos para esta funci ´ on, los puntos de ´
# Chebyshev no mejoran el error.

# para la interpolación por splines
# cúbicos utilizar puntos equidistantes resulta en una estabilidad
# mayor en la estimacion de la función

# ejemplos
# poner cuál fue la más rápida y cuál fue la más precisa
# Por otro lado, los metodos de interpolación nearest neighbor
# e interpolacion lineal destacan por su rapidez en la ejecución.
# Estos metodos son más rápidos y parecieran ser más adecuados
# para aplicaciones en tiempo real o para conjuntos de datos muy
# grandes. Sin embargo, proporcionan una interpolacion menos
# precisa en comparacion con los métodos más complejos. En
# conclusion, la elección del método de interpolación dependerá
# del caso de uso específico.
# Ademas, se encontró que aumentar el número de puntos de interpolacion 
# no siempre resulta en una disminución del error (ej lagrange).


# Definir la función fb(x1, x2) (de consigna)
def fb(x1, x2):
    return 0.75 * np.exp(-((10 * x1 - 2)**2 / 4) - ((9 * x2 - 2)**2 / 4)) + \
           0.65 * np.exp(-((9 * x1 + 1)**2 / 9) - ((10 * x2 + 1)**2 / 2)) + \
           0.55 * np.exp(-((9 * x1 - 6)**2 / 4) - ((9 * x2 - 3)**2 / 4)) - \
           0.01 * np.exp(-((9 * x1 - 7)**2 / 4) - ((9 * x2 - 3)**2 / 4))
xb_min = -1
xb_max = 1

# variables globales
# malla de puntos para evaluar la interpolación
x1_grid = np.linspace(xb_min, xb_max, 100)
x2_grid = np.linspace(xb_min, xb_max, 100)
X1_grid, X2_grid = np.meshgrid(x1_grid, x2_grid)

def definir_puntos(num_points):
    x1_equispaced_fb = equispaced_points(xb_min, xb_max, num_points)
    x2_equispaced_fb = equispaced_points(xb_min, xb_max, num_points)

    x1_nonequispaced_fb = np.sort(chebyshev_points(xb_min, xb_max, num_points))
    x2_nonequispaced_fb = np.sort(chebyshev_points(xb_min, xb_max, num_points))
   
    # Evaluar la función fb en la malla de puntos
    X1_equigrid, X2_equigrid = np.meshgrid(x1_equispaced_fb, x2_equispaced_fb)
    z_equispaced_fb = fb(X1_equigrid, X2_equigrid)

    X1_nonequigrid, X2_nonequigrid = np.meshgrid(x1_nonequispaced_fb, x2_nonequispaced_fb)
    z_nonequispaced_fb = fb(X1_nonequigrid, X2_nonequigrid)
    
    return x1_equispaced_fb, x2_equispaced_fb, x1_nonequispaced_fb, x2_nonequispaced_fb, X1_equigrid, X2_equigrid, z_equispaced_fb, X1_nonequigrid, X2_nonequigrid, z_nonequispaced_fb

def splines_cubicos(num_points):
    x1_equispaced_fb, x2_equispaced_fb, x1_nonequispaced_fb, x2_nonequispaced_fb, X1_equigrid, X2_equigrid, z_equispaced_fb, X1_nonequigrid, X2_nonequigrid, z_nonequispaced_fb = definir_puntos(num_points)

    spline_equispaced_fb = RectBivariateSpline(x1_equispaced_fb, x2_equispaced_fb, z_equispaced_fb)
    spline_nonequispaced_fb = RectBivariateSpline(x1_nonequispaced_fb, x2_nonequispaced_fb, z_nonequispaced_fb)

    Y_interp_equispaced_fb = spline_equispaced_fb(x1_grid, x2_grid)
    Y_interp_nonequispaced_fb = spline_nonequispaced_fb(x1_grid, x2_grid)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), subplot_kw={'projection': '3d'})



    axes[0].plot_surface(X1_grid, X2_grid, Y_interp_equispaced_fb, cmap='viridis', alpha=0.8)
    axes[0].set_xlabel('$x_1$')
    axes[0].set_ylabel('$x_2$')
    axes[0].set_title(f"Interpolación de $f_b(x_1, x_2)$ con Splines Cúbicos (P. Equiespaciados: {num_points})")

    axes[1].plot_surface(X1_grid, X2_grid, Y_interp_nonequispaced_fb, cmap='viridis', alpha=0.8)
    axes[1].set_xlabel('$x_1$')
    axes[1].set_ylabel('$x_2$')
    axes[1].set_title(f"Interpolación de $f_b(x_1, x_2)$ con Splines Cúbicos (P. No Equiespaciados: {num_points})")

    axes[0].plot_wireframe(X1_grid, X2_grid, fb(X1_grid, X2_grid), alpha=0.3, color= "navy")
    axes[1].plot_wireframe(X1_grid, X2_grid, fb(X1_grid, X2_grid), alpha=0.3, color= "navy")
    
    plt.tight_layout()
    plt.show()

def graficar_error_por_nodos(nodes_max, func_text):
    x = np.linspace(xb_min, xb_max, 100)
    y = np.linspace(xb_min, xb_max, 100)
    X, Y = np.meshgrid(x, y)
    Z = fb(X, Y) 
    
    error_equispaced_median = np.zeros(nodes_max)
    error_nonequispaced_median = np.zeros(nodes_max)
    error_equispaced_max = np.zeros(nodes_max)
    error_nonequispaced_max = np.zeros(nodes_max)
   
    error_equispaced_median[0] = np.nan
    error_nonequispaced_median[0] = np.nan
    error_equispaced_max[0] = np.nan
    error_nonequispaced_max[0] = np.nan
    
    error_equispaced_median[1] = np.nan
    error_nonequispaced_median[1] = np.nan
    error_equispaced_max[1] = np.nan
    error_nonequispaced_max[1] = np.nan
   
    for node in range(2, nodes_max + 1):
        x_eq = np.linspace(xb_min, xb_max, node) 
        y_eq = np.linspace(xb_min, xb_max, node)
        X_eq, Y_eq = np.meshgrid(x_eq, y_eq)
        z_eq = fb(X_eq, Y_eq)
        
        x_noneq = np.sort(chebyshev_points(xb_min, xb_max, node))
        y_noneq = np.sort(chebyshev_points(xb_min, xb_max, node))
        X_noneq, Y_noneq = np.meshgrid(x_noneq, y_noneq)
        z_noneq = fb(X_noneq, Y_noneq)
        
        z_eq_interp = griddata((X_eq.ravel(), Y_eq.ravel()), z_eq.ravel(), (X, Y), method='cubic')
        z_noneq_interp = griddata((X_noneq.ravel(), Y_noneq.ravel()), z_noneq.ravel(), (X, Y), method='cubic')
        error_equispaced = np.abs(z_eq_interp - fb(X, Y))
        error_nonequispaced = np.abs(z_noneq_interp - fb(X, Y))
        error_equispaced_max[node - 1] = np.max(error_equispaced)
        error_equispaced_median[node - 1] = np.median(error_equispaced)
        error_nonequispaced_max[node - 1] = np.max(error_nonequispaced)
        error_nonequispaced_median[node - 1] = np.median(error_nonequispaced)

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot()
    ax.plot(np.arange(1, nodes_max + 1), error_equispaced_max, label='Máximo del error equiespaciado', color="red")
    ax.plot(np.arange(1, nodes_max + 1), error_equispaced_median, label='Mediana del error equiespaciado', color="salmon")
    ax.plot(np.arange(1, nodes_max + 1), error_nonequispaced_max, label='Máximo del error no equiespaciado', color="green")
    ax.plot(np.arange(1, nodes_max + 1), error_nonequispaced_median, label='Mediana del error no equiespaciado', color="seagreen")
    ax.set_xlabel('Cantidad de nodos')
    ax.set_xticks(np.arange(0, nodes_max + 1, 5))
    ax.set_ylabel('Error')
    ax.legend()
    plt.suptitle(f"Error de la interpolación usando {func_text}", fontsize=16)
    plt.tight_layout()
    plt.show()

# ERROR  --> investigar para hacerlo en 2d
def main():
    graficar_error_por_nodos(50, "Splines Cúbicos")
    splines_cubicos(15)
    splines_cubicos(25)
    

if __name__ == "__main__":
    main()