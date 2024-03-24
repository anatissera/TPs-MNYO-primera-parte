# se puede hacer una tabla de error absoluto comparando splines y lagrange ponele

# splines cúbicos, splines de 6to orden, Newton, Gauss-Newton

import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

ground_truth_df = pd.read_csv("mnyo_ground_truth.csv", sep=" ", header=None, names=["x1", "x2"])
mediciones_1_df = pd.read_csv("mnyo_mediciones.csv", sep=" ", header=None, names=["x1", "x2"])

splines3_x1_vehic1 = CubicSpline(mediciones_1_df.index, mediciones_1_df["x1"]) # uso los índices de las filas del DataFrame como los puntos de referencia para la interpolación
splines3_x2_vehic1 = CubicSpline(mediciones_1_df.index, mediciones_1_df["x2"])

points_to_evaluate_v1 = np.linspace(mediciones_1_df.index.min(), mediciones_1_df.index.max(), 1000)

def graficar_trayectoria_splines_3():
    plt.figure(figsize=(12, 5))

    plt.plot(ground_truth_df["x1"], ground_truth_df["x2"], label='Ground Truth', linestyle='-.', color='black')

    plt.scatter(mediciones_1_df["x1"], mediciones_1_df["x2"], label="Mediciones", color='orange')
    plt.plot(splines3_x1_vehic1(points_to_evaluate_v1), splines3_x2_vehic1(points_to_evaluate_v1), label='Interpolación con Splines Cúbicos', color='red')

    plt.title('Trayedoria Interpolada del primer vehículo')
    plt.xlabel('eje X')
    plt.ylabel('eje Y')

    plt.legend()
    plt.show()

points_to_compare = np.linspace(mediciones_1_df.index.min(), mediciones_1_df.index.max(), 100)

def graficar_error_abs():
    error_x1 = np.abs(ground_truth_df["x1"] - splines3_x1_vehic1(points_to_compare))
    error_x2 = np.abs(ground_truth_df["x2"] - splines3_x2_vehic1(points_to_compare))

    plt.figure(figsize=(10, 6))
    plt.plot(points_to_compare, error_x1, 'o')
    plt.plot(points_to_compare, error_x2, 'o')
    plt.plot(points_to_compare, error_x1, label='$Error x1$')
    plt.plot(points_to_compare, error_x2, label='$Error x2$')
    plt.xlabel('$x$')
    plt.ylabel('$Error$')
    plt.title(f"$Error absoluto$ de trayectoria interpolada con Splines Cúbicos contra el ground truth")
    plt.legend()
    plt.grid(True)
    plt.show()

def graficar_error_en_tiempo():
    error_x1 = np.abs(ground_truth_df["x1"] - splines3_x1_vehic1(ground_truth_df.index))
    error_x2 = np.abs(ground_truth_df["x2"] - splines3_x2_vehic1(ground_truth_df.index))

    plt.figure(figsize=(10, 6))
    plt.plot(ground_truth_df.index, error_x1, 'o')
    plt.plot(ground_truth_df.index, error_x2, 'o')
    plt.plot(ground_truth_df.index, error_x1, label='$Error x1$')
    plt.plot(ground_truth_df.index, error_x2, label='$Error x2$')
    plt.xlabel('$Tiempo$')
    plt.ylabel('$Error$')
    plt.title(f"Error de Interpolación de trayectoria en función del tiempo con Splines Cúbicos")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    graficar_trayectoria_splines_3()
    graficar_error_abs()
    graficar_error_en_tiempo()
    
if __name__ == "__main__":
    main()