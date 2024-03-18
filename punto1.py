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

# pensar cómo puedo calcular puntos no equiespaciados, puede ser con la fórmula de Chebyshev (está en internet)??


# definir los puntos
x_equispaced_fa = equispaced_points(xa_min, xa_max, 100)
x_nonequispaced_fa = 

y_equispaced_fa = fa(x_equispaced_fa)
y_nonequispaced_fa = fa(x_nonequispaced_fa)

# Interpolar fa(x) con lagrange, haciendo un for iterando por puntos:
# (es lo mismo que usar interpld())
# from scipy.interpolate import Lagrangre también
# teórica dice: fa = lagrange (x, y)
# pero hay que justificar por qué y cómo
# scypy.interpolate.lagrange también me deja obtener los intervalos de error, las bases, revisar documentación, se puede usar para el informe

fa_equispaced = 
fa_nonequispaced = 

# Graficar la interpolación de fa(x) con matplotlib
# se pueden también graficar las bases de lagrange
# además graficar el error contra mi ground truth
# el error se puede calcular restando mi ground truth a la interpolación o también acotándola empíricamente en a y b porque conozco la función
# si no tenés el ground truth, al error lo podés calcular tomando un subconjunto de los puntos que tenés, interpolo, y mido el error contra los restantes
# puedo agarrar 2 puntos al azar con numpy.uniform y repito muchas veces el proceso y calculo el error contra lo otro.
# puedo agarrar primero 2, después 3, después 4, y calculo el error todas las veces.
# puedo agarrar 1 punto sí y uno no. puedo explicitar los errorees en cadaa punto y luego hacer el promedio

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

x1_nonequispaced_fb =
x2_nonequispaced_fb = 

# Evalúa fb(x1, x2) en lo de R2 
z_equispaced_fb = fb(...)
z_nonequispaced_fb = fb(...)

# Interpolar fb(x1, x2) 

# graficar