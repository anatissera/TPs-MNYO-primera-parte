# latex formula generator -> te facilita el código


import numpy as np
import matplotlib.pyplot as plt # para el a
from mpl_toolkits.mplot3d import Axes3D # para el b
# comparar la interpolación
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import lagrange


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
    return 0.3 ** np.linalg.norm(x) * np.sin(4 * x) - np.tanh(2 * x) + 2

xa_min = -4
xa_max = 4

# cómo calcular puntos equiespaciados? con np.linspace() ->
def equispaced_points(a, b, num_points):
    return np.linspace(a, b, num_points)

def chebyshev_points(a, b, num_points):
    return (a + b) / 2 + (b - a) / 2 * np.cos((2 * np.arange(1, num_points + 1) - 1) * np.pi / (2 * num_points))

# pensar cómo puedo calcular puntos no equiespaciados, puede ser con la fórmula de Chebyshev (está en internet)??


# definir los puntos
x_equispaced_fa = equispaced_points(xa_min, xa_max, 100) # operación vectorizada en NumPy
x_nonequispaced_fa = chebyshev_points(xa_min, xa_max, 100)

y_equispaced_fa = fa(x_equispaced_fa)
y_nonequispaced_fa = fa(x_nonequispaced_fa)


# Interpolar fa(x) con lagrange, haciendo un for iterando por puntos:
# (es lo mismo que usar interpld())
# from scipy.interpolate import Lagrangre también
# teórica dice: fa = lagrange (x, y)
# pero hay que justificar por qué y cómo
# scypy.interpolate.lagrange también me deja obtener los intervalos de error, las bases, revisar documentación, se puede usar para el informe

fa_equispaced = lagrange(x_equispaced_fa, y_equispaced_fa)
fa_nonequispaced = lagrange(x_nonequispaced_fa, y_nonequispaced_fa)

# Graficar la interpolación de fa(x) con matplotlib contra fa(x)
# se pueden también graficar las bases de lagrange
# además graficar el error en función del dominio (función del error sería, en los puntos evaluados debe ser 0)
# graficar el error en función de n también es idea
# histograma para el error -> plt.hist(array)

# Mediana -> np.median(array)
# boxplot -> plt.boxplot(array), te da una idea de la distribución. te da el que tiene un cuarto arriba y un cuarto abajo

# el error se puede calcular restando mi ground truth a la interpolación o también acotándola empíricamente en a y b porque conozco la función
# si no tenés el ground truth, al error lo podés calcular tomando un subconjunto de los puntos que tenés, interpolo, y mido el error contra los restantes
# puedo agarrar 2 puntos al azar con numpy.uniform y repito muchas veces el proceso y calculo el error contra lo otro.
# puedo agarrar primero 2, después 3, después 4, y calculo el error todas las veces.
# puedo agarrar 1 punto sí y uno no. puedo explicitar los errorees en cadaa punto
# se puede sacar el mínimo y el máximo, y la mediana (los ordenas de mayor a menor y elegis el del medio) MEDIANA MEJOR QUE PROMEDIO, NO SE USA PROMEDIO -> para tener idea de la densidad

# también se puede hacer una tabla de error absoluto para splines y lagrange ponele

# buscar puntos para comparar -> para el error
def generate_midpoints(lst):
    midpoints = []
    for i in range(len(lst) - 1):
        midpoint = (lst[i] + lst[i+1]) / 2.0
        midpoints.append(midpoint)
    return midpoints

x_compare_equipoints_fa = generate_midpoints(x_equispaced_fa)
x_compare_nonequipoints_fa = generate_midpoints(x_nonequispaced_fa)

y_compare_equipoints_fa = fa(x_compare_equipoints_fa) #ver que onda, f recibe una lista y funciona igual
y_compare_nonequipoints_fa = fa(x_compare_nonequipoints_fa)


# Graficar la interpolación de fa(x) con ambos puntos y fa(x)
plt.figure(figsize=(12, 8))
plt.plot(x_equispaced_fa, y_equispaced_fa, 'o', label='$Puntos de Colocación$ (Equispaciado)')
plt.plot(x_nonequispaced_fa, y_nonequispaced_fa, 'o', label='$Puntos de Colocación$ (No equiespaciado)')

plt.plot(x_compare_equipoints_fa, y_compare_equipoints_fa, 'o', label='$Puntos de Colocación$ (Equispaciado)')
plt.plot(x_compare_nonequipoints_fa, y_nonequispaced_fa, 'o', label='$Puntos de Colocación$ (No equiespaciado)')

plt.plot(x_equispaced_fa, fa(x_equispaced_fa), label='$f_a(x)$', linestyle='--', color='black')  # Graficar la función fa(x)
plt.plot(x_equispaced_fa, fa_equispaced(x_equispaced_fa), label='$Interpolación$ (Equispaciado)')
plt.plot(x_nonequispaced_fa, fa_nonequispaced(x_nonequispaced_fa), label='$Interpolación$ (No equiespaciado)')

# graficar fa en los midpoints
plt.plot(x_compare_equipoints_fa, fa(x_compare_equipoints_fa), label='$f_a(x)$', linestyle='--', color='black')  # Graficar la función fa(x)
plt.plot(x_compare_equipoints_fa, fa_equispaced(x_compare_equipoints_fa), label='$Interpolación$ (Equispaciado)')
plt.plot(x_compare_nonequipoints_fa, fa_nonequispaced(x_compare_nonequipoints_fa), label='$Interpolación$ (No equiespaciado)')


plt.xlabel('$x$')
plt.ylabel('$f_a(x)$')
plt.title('Interpolación de $f_a(x)$ con Ambos Puntos de Colocación')
plt.legend()
plt.grid(True)
plt.show()

# Graficar las bases de Lagrange
for i in range(len(x_equispaced_fa)):
    plt.plot(x_equispaced_fa, fa_equispaced(x_equispaced_fa[i]) * lagrange(x_equispaced_fa, np.eye(len(x_equispaced_fa))[i])(x_equispaced_fa), linestyle='--', color='gray', alpha=0.5)
    plt.plot(x_nonequispaced_fa, fa_nonequispaced(x_nonequispaced_fa[i]) * lagrange(x_nonequispaced_fa, np.eye(len(x_nonequispaced_fa))[i])(x_nonequispaced_fa), linestyle='--', color='gray', alpha=0.5)

plt.xlabel('$x$')
plt.ylabel('$f_a(x)$')
plt.title('Interpolación de $f_a(x)$ y Bases de Lagrange')
plt.legend()
plt.grid(True)
plt.show()

# Calcular el error en función del dominio
error_equispaced = np.abs(y_equispaced_fa - fa_equispaced(x_equispaced_fa))
error_nonequispaced = np.abs(y_nonequispaced_fa - fa_nonequispaced(x_nonequispaced_fa))


# Graficar el error en función del dominio
plt.figure(figsize=(10, 6))
plt.plot(x_equispaced_fa, error_equispaced, label='$Error$ (Equispaciado)')
plt.plot(x_nonequispaced_fa, error_nonequispaced, label='$Error$ (No equiespaciado)')
plt.xlabel('$x$')
plt.ylabel('$Error$')
plt.title('Error de Interpolación de $f_a(x)$')
plt.legend()
plt.grid(True)
plt.show()

# Calcular estadísticas del error
median_error_equispaced = np.median(error_equispaced)
max_error_equispaced = np.max(error_equispaced)
min_error_equispaced = np.min(error_equispaced)

median_error_nonequispaced = np.median(error_nonequispaced)
max_error_nonequispaced = np.max(error_nonequispaced)
min_error_nonequispaced = np.min(error_nonequispaced)

print("Estadísticas del error para puntos equiespaciados:")
print("Mediana:", median_error_equispaced)
print("Máximo:", max_error_equispaced)
print("Mínimo:", min_error_equispaced)

print("Estadísticas del error para puntos no equiespaciados:")
print("Mediana:", median_error_nonequispaced)
print("Máximo:", max_error_nonequispaced)
print("Mínimo:", min_error_nonequispaced)

# Graficar el histograma del error
plt.figure(figsize=(10, 6))
plt.hist(error_equispaced, bins=20, alpha=0.5, label='$Equispaciado$')
plt.hist(error_nonequispaced, bins=20, alpha=0.5, label='$No equiespaciado$')
plt.xlabel('$Error$')
plt.ylabel('$Frecuencia$')
plt.title('Histograma del Error de Interpolación de $f_a(x)$')
plt.legend()
plt.grid(True)
plt.show()


# Definir la función fb(x1, x2) (de consigna)
def fb(x1, x2):
    return 0.75 * np.exp(-((10 * x1 - 2)**2 / 4) - ((9 * x2 - 2)**2 / 4)) + \
           0.65 * np.exp(-((9 * x1 + 1)**2 / 9) - ((10 * x2 + 1)**2 / 2)) + \
           0.55 * np.exp(-((9 * x1 - 6)**2 / 4) - ((9 * x2 - 3)**2 / 4)) - \
           0.01 * np.exp(-((9 * x1 - 7)**2 / 4) - ((9 * x2 - 3)**2 / 4))
xb_min = -1
xb_max = 1

# hacer lo mismo que en fa(x)
# Generar los puntos según el intervalo dado

x1_equispaced_fb = equispaced_points(xb_min, xb_max, 100)
x2_equispaced_fb = equispaced_points(xb_min, xb_max, 100)

# la derivada se debe explicitar, pero también hay librerías, las podés calcular con librerías pero hay que ponerlas. la derivada la calculé (la puedo dejar en el apéndice si es muy larga)

x1_nonequispaced_fb = chebyshev_points(-1, 1, 100)
x2_nonequispaced_fb = chebyshev_points(-1, 1, 100)

# Evalúa fb(x1, x2) en lo de R2 
# z_equispaced_fb = fb(...)
# z_nonequispaced_fb = fb(...)

# Interpolar fb(x1, x2) 

# graficar