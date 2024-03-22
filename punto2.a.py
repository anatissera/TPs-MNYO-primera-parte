# se puede usar "splines" -> interpolar a trozos
# es una función partida (de grado 1, modelo a trozos)
# splines cúbico me da c2

# se puede hacer una tabla de error absoluto comparando splines y lagrange ponele

# splines cúbicos, splines de 6to orden, Newton, Gauss-Newton

import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

mediciones_1 = pd.read_csv('myno_mediciones.csv')
ground_truth = pd.read_csv('myno_ground_truth.csv')

x1_registered = mediciones_1.iloc[0].values
x2_registered = mediciones_1.iloc[1].values

x1 = ground_truth.iloc[0].values
x2 = ground_truth.iloc[1].values

# X1_equigrid, X2_equigrid = np.meshgrid(x1, x1)

splines_func = CubicSpline(x1_registered)

y_interp = splines_func (x1, x2)

plt.figure(figsize=(8, 6))
plt.plot(x1, y_interp, label='Interpolación', color='red')
plt.scatter(x1_registered, x2_registered, label='Mediciones Registradas', color='blue')
plt.plot(x1, x2, label='Ground Truth', linestyle='--', color='green')
plt.xlabel('x1')
plt.ylabel('x2')
plt.title('Interpolación con Splines Cúbicos')
plt.legend()
plt.grid(True)
plt.show()

