# se puede usar "splines" -> interpolar a trozos
# es una función partida (de grado 1, modelo a trozos)
# splines cúbico me da c2

# se puede hacer una tabla de error absoluto comparando splines y lagrange ponele

# splines cúbicos, splines de 6to orden, Newton, Gauss-Newton

import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline

mediciones_1 = pd.read_csv('myno_mediciones.csv')

x1 = mediciones_1.iloc[0].values
x2 = mediciones_1.iloc[1].values

# X1_equigrid, X2_equigrid = np.meshgrid(x1, x1)

splines = CubicSpline(x1, x2)


