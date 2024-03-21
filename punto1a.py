
# latex formula generator -> te facilita el código

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import lagrange
from scipy.interpolate import CubicSpline

# Consigna 1
# Estudiar el desempeño de distintos esquemas de interpolación en las funciones
# a)
# fa(x) = 0.3 **∥x∥ * sin (4x) − tanh (2x) + 2
# con x ∈ [−4, 4].

#b )
# fb(x) = 0.75 * exp ( − (10 * x1 − 2)**2 / 4 − (9 * x2 − 2)**2 / 4 ) + 0.65 * exp ( − (9 * x1 + 1)**2 / 9 − (10 * x2 + 1)** / 2 ) + 0.55 * exp ( − (9 * x1 − 6)**2 / 4 − (9 * x2 − 3)**2 / 4 ) − 0.01 * exp (− (9 *x1 − 7)**2 / 4 − (9 * x2 − 3)**2 / 4 )
# con x1, x2 ∈ [−1, 1].

# Hágalo primero tomando puntos de colocación equiespaciados. Luego proponga (al menos) una regla
# para elegir puntos no equiespaciados. Compare los resultados.

# qué hace exp()?: np.exp() calculates e^x for each value of x in your input array.

# escribir en código la consigna
def fa(x):
    return 0.3 ** abs(abs((x) ))* np.sin(4 * x) - np.tanh(2 * x) + 2

xa_min = -4
xa_max = 4

# cómo calcular puntos equiespaciados? con np.linspace() ->
def equispaced_points(a, b, num_points):
    return np.linspace(a, b, num_points)

def chebyshev_points(a, b, num_points):
    return (a + b) / 2 + (b - a) / 2 * np.cos((2 * np.arange(1, num_points + 1) - 1) * np.pi / (2 * num_points))

def graficar_interpol_ambos_puntos(f_interpol_equispaced, f_interpool_nonequispaced, x_equispaced_fa, y_equispaced_fa, x_nonequispaced_fa, y_nonequispaced_fa, method, q_points):
    plt.figure(figsize=(12, 8))
    plt.plot(points_to_study_function, fa(points_to_study_function), label='$f_a(x)$', linestyle='--', color='black')  # Graficar la función fa(x)
    plt.plot(x_equispaced_fa, y_equispaced_fa, 'o', label='$Puntos de Colocación$ (Equispaciado)', color = 'blue')
    plt.plot(points_to_study_equifunction, f_interpol_equispaced(points_to_study_equifunction), label='$Interpolación$ (Equispaciado)', color = 'green')
    plt.plot(x_nonequispaced_fa, y_nonequispaced_fa, 'o', label='$Puntos de Colocación$ (No equiespaciado)', color = 'orange')
    plt.plot(points_to_study_function, f_interpool_nonequispaced(points_to_study_function), label='$Interpolación$ (No equiespaciado)', color = 'red')

    plt.xlabel('$x$')
    plt.ylabel('$f_a(x)$')
    plt.title(f"Interpolación de $f_a(x)$ con {method} con {q_points} Puntos de Colocación")
    plt.legend()
    plt.grid(True)
    plt.show()
    
def graficar_error(f_interpol_equispaced, f_interpol_nonequispaced, x_compare_equipoints_fa, x_compare_nonequipoints_fa, method, q_points):
    error_equispaced_w_midpoints = np.abs(fa(x_compare_equipoints_fa) - f_interpol_equispaced(x_compare_equipoints_fa))
    error_nonequispaced_w_midpoints = np.abs(fa(x_compare_nonequipoints_fa) - f_interpol_nonequispaced(x_compare_nonequipoints_fa))

    plt.figure(figsize=(10, 6))
    plt.plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, 'o')
    plt.plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, 'o')
    plt.plot(x_compare_equipoints_fa, error_equispaced_w_midpoints, label='$Error$ (Equispaciado)')
    plt.plot(x_compare_nonequipoints_fa, error_nonequispaced_w_midpoints, label='$Error$ (No equiespaciado)')
    plt.xlabel('$x$')
    plt.ylabel('$Error$')
    plt.title(f"Error de Interpolación con {method} de $f_a(x)$ con {q_points} puntos")
    plt.legend()
    plt.grid(True)
    plt.show()
    
x_few_equispaced_fa = equispaced_points(xa_min, xa_max, 13) # con Lagrange es mejor menos puntos, si pongo más se arranca a disparar
x_few_nonequispaced_fa = chebyshev_points(xa_min, xa_max, 14)

x_more_equispaced_fa = equispaced_points(xa_min, xa_max, 20) 
x_more_nonequispaced_fa = chebyshev_points(xa_min, xa_max, 20)

y_few_equispaced_fa = fa(x_few_equispaced_fa)
y_few_nonequispaced_fa = fa(x_few_nonequispaced_fa)

y_more_equispaced_fa = fa(x_more_equispaced_fa)
y_more_nonequispaced_fa = fa(x_more_nonequispaced_fa)

fa_lagrange_few_equispaced = lagrange(x_few_equispaced_fa, y_few_equispaced_fa)
fa_lagrange_few_nonequispaced = lagrange(x_few_nonequispaced_fa, y_few_nonequispaced_fa)

fa_lagrange_more_equispaced = lagrange(x_more_equispaced_fa, y_more_equispaced_fa)
fa_lagrange_more_nonequispaced = lagrange(x_more_nonequispaced_fa, y_more_nonequispaced_fa)

# Graficar la interpolación de fa(x) con ambos puntos y fa(x)
points_to_study_function = equispaced_points(-3.98, 3.98, 150)
points_to_study_equifunction = equispaced_points(-3.46, 3.46, 150)

# Graficar la interpolación de fa(x) con ambos puntos y fa(x)

# buscar puntos para comparar -> para el error

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

x_compare_few_equipoints_fa = generate_midpoints(x_few_equispaced_fa)
x_compare_few_nonequipoints_fa = generate_midpoints(x_few_nonequispaced_fa)

x_compare_more_equipoints_fa = generate_midpoints(x_more_equispaced_fa)
x_compare_more_nonequipoints_fa = generate_midpoints(x_more_nonequispaced_fa)

# hacer lo mismo con splines 
spline_few_equispaced_fa = CubicSpline(x_few_equispaced_fa, y_few_equispaced_fa)
spline_more_equispaced_fa = CubicSpline(x_more_equispaced_fa, y_more_equispaced_fa)

sorted_few_indices = np.argsort(x_few_nonequispaced_fa)
sorted_more_indices = np.argsort(x_more_nonequispaced_fa)
x_few_nonequispaced_fa_sorted = x_few_nonequispaced_fa[sorted_few_indices]
y_few_nonequispaced_fa_sorted = y_few_nonequispaced_fa[sorted_few_indices]
x_more_nonequispaced_fa_sorted = x_more_nonequispaced_fa[sorted_more_indices]
y_more_nonequispaced_fa_sorted = y_more_nonequispaced_fa[sorted_more_indices]

# Interpolación con splines cúbicos
spline_few_nonequispaced_fa = CubicSpline(x_few_nonequispaced_fa_sorted, y_few_nonequispaced_fa_sorted)
spline_more_nonequispaced_fa = CubicSpline(x_more_nonequispaced_fa_sorted, y_more_nonequispaced_fa_sorted)

points_to_study_function = equispaced_points(-3.98, 3.98, 150)
points_to_study_equifunction = equispaced_points(-3.46, 3.46, 150)

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
    
    graficar_interpol_ambos_puntos(fa_lagrange_few_equispaced, fa_lagrange_few_nonequispaced, x_few_equispaced_fa, y_few_equispaced_fa, x_few_nonequispaced_fa, y_few_nonequispaced_fa, "Lagrange", "pocos")

    graficar_error(fa_lagrange_few_equispaced, fa_lagrange_few_nonequispaced, x_compare_few_equipoints_fa, x_compare_few_nonequipoints_fa, "Lagrange", "pocos")

    graficar_interpol_ambos_puntos(fa_lagrange_more_equispaced, fa_lagrange_more_nonequispaced, x_more_equispaced_fa, y_more_equispaced_fa, x_more_nonequispaced_fa, y_more_nonequispaced_fa, "Lagrange", "más")

    graficar_error(fa_lagrange_more_equispaced, fa_lagrange_more_nonequispaced, x_compare_more_equipoints_fa, x_compare_more_nonequipoints_fa, "Lagrange", "más")

    # Graficar la interpolación de fa(x) con ambos puntos y fa(x)
    graficar_interpol_ambos_puntos(spline_few_equispaced_fa, spline_few_nonequispaced_fa, x_few_equispaced_fa, y_few_equispaced_fa, x_few_nonequispaced_fa, y_few_nonequispaced_fa, "Splines Cúbicos", "pocos")

    graficar_error(spline_few_equispaced_fa, spline_few_nonequispaced_fa, x_compare_few_equipoints_fa, x_compare_few_nonequipoints_fa, "Splines Cúbicos", "pocos")

    graficar_interpol_ambos_puntos(spline_more_equispaced_fa, spline_more_nonequispaced_fa, x_more_equispaced_fa, y_more_equispaced_fa, x_more_nonequispaced_fa, y_more_nonequispaced_fa, "Splines Cúbicos", "más")

    graficar_error(spline_more_equispaced_fa, spline_more_nonequispaced_fa, x_compare_more_equipoints_fa, x_compare_more_nonequipoints_fa, "Splines Cúbicos", "más")

def main():
    graficar()
    
if __name__ == "__main__":
    main()