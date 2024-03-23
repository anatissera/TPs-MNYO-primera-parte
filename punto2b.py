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

def interseccion(t, t2):
    x1_trayec_vehic1, x2_trayec_vehic1 = interpol_x1_vehic1(t), interpol_x2_vehic1(t)
    x1_trayec_vehic2, x2_trayec_vehic2 = interpol_x1_vehic2(t2), interpol_x2_vehic2(t2)
    return (x1_trayec_vehic1 - x1_trayec_vehic2, x2_trayec_vehic1 - x2_trayec_vehic2)

def newton_2d(f, P0, P1, tol=1e-6, max_iter=100):
    P = np.array([P0, P1])  # Aproximación inicial
    for _ in range(max_iter):
        f_val = f(*P)
        if np.linalg.norm(f_val) < tol:
            break
        J = np.array([[f(P[0]+tol, P[1])[0]-f_val[0], f(P[0], P[1]+tol)[0]-f_val[0]],
                      [f(P[0]+tol, P[1])[1]-f_val[1], f(P[0], P[1]+tol)[1]-f_val[1]]]) / tol
        delta_P = np.linalg.solve(J, np.array([-f_val[0], -f_val[1]]))
        P = P + delta_P
    return tuple(P)

t0 = mediciones_1_df.index.min()
t1 = mediciones_2_df.index.min()

t_interseccion, t2_interseccion = newton_2d(interseccion, t0, t1)

x_interseccion_1, y_interseccion_1 = interpol_x1_vehic1(t_interseccion), interpol_x2_vehic1(t_interseccion)
x_interseccion_2, y_interseccion_2 = interpol_x1_vehic2(t2_interseccion), interpol_x2_vehic2(t2_interseccion)

error_x1 = np.abs(x_interseccion_1 - x_interseccion_2)
error_x2 = np.abs(y_interseccion_1 - y_interseccion_2)

print(f"Coordenadas de la intersección-> vehículo 1:({x_interseccion_1}, {y_interseccion_1}) y vehículo 2:({x_interseccion_2}, {y_interseccion_2})")

def graficar_trayectorias_intersec():
    plt.figure(figsize=(10, 6))
    plt.plot(ground_truth_df["x1"], ground_truth_df["x2"], label='Ground Truth', linestyle='-.', color='dimgrey')

    plt.scatter(mediciones_1_df["x1"], mediciones_1_df["x2"], label="Mediciones", color='yellowgreen', s=20)
    plt.plot(interpol_x1_vehic1(points_to_evaluate_v1), interpol_x2_vehic1(points_to_evaluate_v1), label='Interpolación con Splines Cúbicos v1', color='teal')

    plt.plot(interpol_x1_vehic2(points_to_evaluate_v2), interpol_x2_vehic2(points_to_evaluate_v2), label='Interpolation con Splines Cúbicos v2', color='mediumpurple')

    plt.plot(x_interseccion_1, y_interseccion_1, 'x', label="Intersección", markersize=8.5, color='darkred',  markeredgewidth=4)

    plt.xlabel("Coordenada $x1$")
    plt.ylabel("Coordenada $x2$")
    plt.title("Trayectorias de los 2 Vehículos")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

def graficar_error_abs():
    plt.figure(figsize=(7, 6))
    index = np.arange(2)

    plt.bar(index, [error_x1, error_x2], 0.4, color=['khaki', 'darkseagreen'])
    plt.xlabel('Dimensiones')
    plt.ylabel('Error')
    plt.title('Error entre las Coordenadas de la Intersección $(Error bound: 1e-10)$')
    plt.xticks(index, ('$x1$', '$x2$'))
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    graficar_trayectorias_intersec()
    graficar_error_abs()