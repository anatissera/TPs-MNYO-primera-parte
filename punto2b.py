import pandas as pd
from punto2a import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

mediciones_2_df = pd.read_csv("mnyo_mediciones2.csv", sep=" ", header=None, names=["x1", "x2"])

points_to_evaluate_v2 = np.linspace(mediciones_2_df.index.min(), mediciones_2_df.index.max(), 1000)

interpol_x1_vehic2 = CubicSpline(mediciones_2_df.index, mediciones_2_df["x1"]) # uso los índices de las filas del DataFrame como los puntos de referencia para la interpolación
interpol_x2_vehic2 = CubicSpline(mediciones_2_df.index, mediciones_2_df["x2"])

x1_vehiculo1, x2_vehiculo1 = mediciones_1_df["x1"], mediciones_1_df["x2"]
x1_vehiculo2, x2_vehiculo2 = mediciones_2_df["x1"], mediciones_2_df["x2"]


# # Interpolación de la trayectoria del segundo vehículo con Lagrange
# polinomio_interpolado = CubicSpline(x1_vehiculo2, x2_vehiculo2)
# x1_interpolados = np.linspace(min(x1_vehiculo2), max(x1_vehiculo2), len(x1_vehiculo1))
# x2_interpolados = polinomio_interpolado(x1_interpolados)

# def diff(x):
#     return np.interp(x, x1_vehiculo1, x2_vehiculo1) - np.interp(x, x1_interpolados, x2_interpolados)

# # Usar el método de Newton para encontrar la intersección
# x_interseccion = newton(diff, x0=0)

# # Calcular las coordenadas de la intersección
# y_interseccion = np.interp(x_interseccion, x1_vehiculo1, x2_vehiculo1)

# # Mostrar las coordenadas de la intersección
# print(f"Las coordenadas de la intersección son: ({x_interseccion}, {y_interseccion})")



# TERMINA SECCIÓN INTERSECCIÓN


# Graficar trayectorias de ambos vehículos
plt.figure(figsize=(10, 6))
plt.plot(ground_truth_df["x1"], ground_truth_df["x2"], label='Ground Truth', linestyle='-.', color='black')

plt.scatter(mediciones_1_df["x1"], mediciones_1_df["x2"], label="Mediciones", color='orange')
plt.plot(interpol_x1_vehic1(points_to_evaluate_v1), interpol_x2_vehic1(points_to_evaluate_v1), label='Interpolación con Splines Cúbicos v1', color='red')
plt.plot(interpol_x1_vehic2(points_to_evaluate_v2), interpol_x2_vehic2(points_to_evaluate_v2), label='Interpolation con Splines Cúbicos v2', color='blue')
# plt.plot(x1_vehiculo1, x2_vehiculo1, label="Trayectoria Vehículo 1", marker='o')
# plt.plot(x1_vehiculo2, x2_vehiculo2, label="Trayectoria Vehículo 2", marker='o')
# plt.plot(x1_interpolados, x2_interpolados, label="Trayectoria Interpolada Vehículo 2")
# plt.plot(x_interseccion, y_interseccion, 'ro', label="Intersección")
plt.plot(x_interseccion, y_interseccion, 'ro', label="Intersección")

# Configuración de la gráfica
plt.xlabel("Coordenada x1")
plt.ylabel("Coordenada x2")
plt.title("Trayectorias de Vehículos")
plt.legend()
plt.grid()
plt.show()
