import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from scipy.interpolate import CubicSpline

def fa(x):
    return 0.3 ** abs(abs((x) ))* np.sin(4 * x) - np.tanh(2 * x) + 2

xa_min = -4
xa_max = 4

def equispaced_points(a, b, num_points):
    return np.linspace(a, b, num_points)

def chebyshev_points(a, b, num_points):
    return (a + b) / 2 + (b - a) / 2 * np.cos((2 * np.arange(1, num_points + 1) - 1) * np.pi / (2 * num_points))

def graficar_interpol_ambos_puntos(f_interpol_equispaced, f_interpol_nonequispaced, x_equispaced_fa, y_equispaced_fa, x_nonequispaced_fa, y_nonequispaced_fa, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points):
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    
    points_to_study_function = equispaced_points(-3.98, 3.98, 150)
    points_to_study_equifunction = equispaced_points(-3.46, 3.46, 150)
    
    axs[0].plot(points_to_study_function, fa(points_to_study_function), label='$f_a(x)$', linestyle='--', color='black')  # Graficar la función fa(x)
    axs[0].plot(x_equispaced_fa, y_equispaced_fa, 'o', label='$Puntos de Colocación$ (Equispaciado)', color = 'blue')
    axs[0].plot(points_to_study_equifunction, f_interpol_equispaced(points_to_study_equifunction), label='$Interpolación$ (Equispaciado)', color = 'green')
    axs[0].plot(x_nonequispaced_fa, y_nonequispaced_fa, 'o', label='$Puntos de Colocación$ (No equiespaciado)', color = 'orange')
    axs[0].plot(points_to_study_function, f_interpol_nonequispaced(points_to_study_function), label='$Interpolación$ (No equiespaciado)', color = 'red')
    axs[0].set_xlabel('$x$')
    axs[0].set_ylabel('$f_a(x)$')
    axs[0].set_title(f"Interpolación de $f_a(x)$ con {method} con {q_points} Puntos de Colocación")
    axs[0].legend()
    axs[0].grid(True)
    
    error_equispaced_w_midpoints = np.abs(fa(x_compare_equipoints_fa) - f_interpol_equispaced(x_compare_equipoints_fa))
    error_nonequispaced_w_midpoints = np.abs(fa(x_compare_nonequipoints_fa) - f_interpol_nonequispaced(x_compare_nonequipoints_fa))

    axs[1].plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, 'o')
    axs[1].plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, 'o')
    axs[1].plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, label='$Error$ (Equispaciado)')
    axs[1].plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, label='$Error$ (No equiespaciado)')
    axs[1].set_xlabel('$x$')
    axs[1].set_ylabel('$Error$')
    axs[1].set_title(f"Error de Interpolación con {method} de $f_a(x)$ con {q_points} puntos")
    axs[1].legend()
    axs[1].grid(True)
    
    plt.tight_layout()
    plt.show()
    
    
# def graficar_error(f_interpol_equispaced, f_interpol_nonequispaced, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points):
#     error_equispaced_w_midpoints = np.abs(fa(x_compare_equipoints_fa) - f_interpol_equispaced(x_compare_equipoints_fa))
#     error_nonequispaced_w_midpoints = np.abs(fa(x_compare_nonequipoints_fa) - f_interpol_nonequispaced(x_compare_nonequipoints_fa))

#     plt.figure(figsize=(10, 6))
#     plt.plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, 'o')
#     plt.plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, 'o')
#     plt.plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, label='$Error$ (Equispaciado)')
#     plt.plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, label='$Error$ (No equiespaciado)')
#     plt.xlabel('$x$')
#     plt.ylabel('$Error$')
#     plt.title(f"Error de Interpolación con {method} de $f_a(x)$ con {q_points} puntos")
#     plt.legend()
#     plt.grid(True)
#     plt.show()

def generate_midpoints(lst): # -> solo los midpoints
        midpoints = []
        for i in range(len(lst) - 1):
            midpoint = (lst[i] + lst[i+1]) / 2.0
            midpoints.append(midpoint)
        return np.array(midpoints)

    # def generate_midpoints(lst): # -> los midpoints y los puntos de interpolación
    #     midpoints = []
    #     for i in range(len(lst) - 1):
    #         midpoint = (lst[i] + lst[i+1]) / 2.0
    #         midpoints.extend([lst[i], midpoint])
    #     midpoints.append(lst[-1])  # Agregar el último elemento de la lista original
    #     return np.array(midpoints)

def lagrange_interpol(num_points, method, q_points):
            
    x_equispaced_fa = equispaced_points(xa_min, xa_max, num_points)
    x_nonequispaced_fa = chebyshev_points(xa_min, xa_max, num_points)

    y_equispaced_fa = fa(x_equispaced_fa)
    y_nonequispaced_fa = fa(x_nonequispaced_fa)

    fa_lagrange_equispaced = lagrange(x_equispaced_fa, y_equispaced_fa)
    fa_lagrange_nonequispaced = lagrange(x_nonequispaced_fa, y_nonequispaced_fa)

    # buscar puntos para comparar -> para el error
    x_compare_equipoints_fa = generate_midpoints(x_equispaced_fa)
    x_compare_nonequipoints_fa = generate_midpoints(x_nonequispaced_fa)
    
    graficar_interpol_ambos_puntos(fa_lagrange_equispaced, fa_lagrange_nonequispaced, x_equispaced_fa, y_equispaced_fa, x_nonequispaced_fa, y_nonequispaced_fa, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points)

    # graficar_error(fa_lagrange_equispaced, fa_lagrange_nonequispaced, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points)


def splines_interpol(num_points, method, q_points):
    x_equispaced_fa = equispaced_points(xa_min, xa_max, num_points)
    x_nonequispaced_fa = chebyshev_points(xa_min, xa_max, num_points)

    y_equispaced_fa = fa(x_equispaced_fa)
    y_nonequispaced_fa = fa(x_nonequispaced_fa)

    spline_equispaced_fa = CubicSpline(x_equispaced_fa, y_equispaced_fa)

    sorted_few_indices = np.argsort(x_nonequispaced_fa)
    x_nonequispaced_fa_sorted = x_nonequispaced_fa[sorted_few_indices]
    y_nonequispaced_fa_sorted = y_nonequispaced_fa[sorted_few_indices]
    
    x_compare_equipoints_fa = generate_midpoints(x_equispaced_fa)
    x_compare_nonequipoints_fa = generate_midpoints(x_nonequispaced_fa_sorted)

    spline_nonequispaced_fa = CubicSpline(x_nonequispaced_fa_sorted, y_nonequispaced_fa_sorted)
    
    graficar_interpol_ambos_puntos(spline_equispaced_fa, spline_nonequispaced_fa, x_equispaced_fa, y_equispaced_fa, x_nonequispaced_fa, y_nonequispaced_fa, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points)

    # graficar_error(spline_equispaced_fa, spline_nonequispaced_fa, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points)

# plt.plot(x_compare_equipoints_fa, y_compare_equipoints_fa, 'o', label='$Puntos de Colocación$ (Equispaciado)')
# plt.plot(x_compare_nonequipoints_fa, y_nonequispaced_fa, 'o', label='$Puntos de Colocación$ (No equiespaciado)')

# graficar fa en los midpoints
# plt.plot(x_compare_equipoints_fa, fa(x_compare_equipoints_fa), label='$f_a(x)$', linestyle='--', color='black')  # Graficar la función fa(x)
# plt.plot(x_compare_equipoints_fa, fa_equispaced(x_compare_equipoints_fa), label='$Interpolación$ (Equispaciado)')
# plt.plot(x_compare_nonequipoints_fa, fa_nonequispaced(x_compare_nonequipoints_fa), label='$Interpolación$ (No equiespaciado)')

# # Calcular los errores absolutos de las interpolaciones equiespaciadas y no equiespaciadas
# error_equispaced = np.abs(y_equispaced_fa - fa_interpolation_equispaced(x_equispaced_fa))
# error_nonequispaced = np.abs(y_nonequispaced_fa - fa_interpolation_nonequispaced(x_nonequispaced_fa))

# # Histograma del error
# plt.figure(figsize=(10, 6))
# plt.hist(error_equispaced, bins=20, alpha=0.5, label='$Equispaciado$')
# plt.hist(error_nonequispaced, bins=20, alpha=0.5, label='$No\ equiespaciado$')
# plt.xlabel('$Error$')
# plt.ylabel('$Frecuencia$')
# plt.title('Histograma del Error de Interpolación de $f_a(x)$')
# plt.legend()
# plt.grid(True)
# plt.show()

# # Boxplot del error
# plt.figure(figsize=(10, 6))
# plt.boxplot([error_equispaced, error_nonequispaced], labels=['Equispaciado', 'No equiespaciado'])
# plt.xlabel('Tipo de Puntos')
# plt.ylabel('Error Absoluto')
# plt.title('Boxplot del Error de Interpolación de $f_a(x)$')
# plt.grid(True)
# plt.show()


# # Graficar las bases de Lagrange
# for i in range(len(x_equispaced_fa)):
#     plt.plot(x_equispaced_fa, fa_interpolation_equispaced(x_equispaced_fa[i]) * lagrange(x_equispaced_fa, np.eye(len(x_equispaced_fa))[i])(x_equispaced_fa), linestyle='--', color='gray', alpha=0.5)
#     plt.plot(x_nonequispaced_fa, fa_interpolation_nonequispaced(x_nonequispaced_fa[i]) * lagrange(x_nonequispaced_fa, np.eye(len(x_nonequispaced_fa))[i])(x_nonequispaced_fa), linestyle='--', color='gray', alpha=0.5)

# plt.xlabel('$x$')
# plt.ylabel('$f_a(x)$')
# plt.title('Interpolación de $f_a(x)$ y Bases de Lagrange')
# plt.legend()
# plt.grid(True)
# plt.show()


# # Calcular estadísticas del error
# median_error_equispaced = np.median(error_equispaced)
# max_error_equispaced = np.max(error_equispaced)
# min_error_equispaced = np.min(error_equispaced)

# median_error_nonequispaced = np.median(error_nonequispaced)
# max_error_nonequispaced = np.max(error_nonequispaced)
# min_error_nonequispaced = np.min(error_nonequispaced)

def graficar():
    lagrange_interpol(13, "Lagrange", "pocos")
    lagrange_interpol(20, "Lagrange", "muchos")
    splines_interpol(13, "Splines", "pocos")
    splines_interpol(20, "Splines", "muchos")

def main():
    graficar()
    
if __name__ == "__main__":
    main()