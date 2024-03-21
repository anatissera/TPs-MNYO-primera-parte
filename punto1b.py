from punto1a import equispaced_points, chebyshev_points
import matplotlib.pyplot as plt # para el a
from mpl_toolkits.mplot3d import Axes3D # para el b
from scipy.interpolate import RectBivariateSpline
import numpy as np
from scipy.interpolate import lagrange


# Definir la función fb(x1, x2) (de consigna)
def fb(x1, x2):
    return 0.75 * np.exp(-((10 * x1 - 2)**2 / 4) - ((9 * x2 - 2)**2 / 4)) + \
           0.65 * np.exp(-((9 * x1 + 1)**2 / 9) - ((10 * x2 + 1)**2 / 2)) + \
           0.55 * np.exp(-((9 * x1 - 6)**2 / 4) - ((9 * x2 - 3)**2 / 4)) - \
           0.01 * np.exp(-((9 * x1 - 7)**2 / 4) - ((9 * x2 - 3)**2 / 4))
xb_min = -1
xb_max = 1

# hacer lo mismo que en fa(x)

x1_equispaced_fb = equispaced_points(xb_min, xb_max, 15)
x2_equispaced_fb = equispaced_points(xb_min, xb_max, 15)

x1_nonequispaced_fb = np.sort(chebyshev_points(xb_min, xb_max, 15))
x2_nonequispaced_fb = np.sort(chebyshev_points(xb_min, xb_max, 15))
# x2_nonequispaced_fb = np.sort(np.cos(np.linspace(0, np.pi, 15))) esta forma tiene más error porque no tiene en cuenta la derivada segunda

# Evaluar la función fb en la malla de puntos
X1_equigrid, X2_equigrid = np.meshgrid(x1_equispaced_fb, x2_equispaced_fb)
z_equispaced_fb = fb(X1_equigrid, X2_equigrid)

X1_nonequigrid, X2_nonequigrid = np.meshgrid(x1_nonequispaced_fb, x2_nonequispaced_fb)
z_nonequispaced_fb = fb(X1_nonequigrid, X2_nonequigrid)

spline_equispaced_fb = RectBivariateSpline(x1_equispaced_fb, x2_equispaced_fb, z_equispaced_fb)
spline_nonequispaced_fb = RectBivariateSpline(x1_nonequispaced_fb, x2_nonequispaced_fb, z_nonequispaced_fb)

# malla de puntos para evaluar la interpolación
x1_grid = np.linspace(xb_min, xb_max, 100)
x2_grid = np.linspace(xb_min, xb_max, 100)
X1_grid, X2_grid = np.meshgrid(x1_grid, x2_grid)

Y_interp_equispaced_fb = spline_equispaced_fb(x1_grid, x2_grid)
Y_interp_nonequispaced_fb = spline_nonequispaced_fb(x1_grid, x2_grid)

fig, axes = plt.subplots(1, 2, figsize=(16, 6), subplot_kw={'projection': '3d'})

axes[0].plot_surface(X1_grid, X2_grid, fb(X1_grid, X2_grid), cmap='viridis', alpha=0.3)
axes[1].plot_surface(X1_grid, X2_grid, fb(X1_grid, X2_grid), cmap='viridis', alpha=0.3)

axes[0].plot_surface(X1_grid, X2_grid, Y_interp_equispaced_fb, cmap='viridis', alpha=0.8)
axes[0].set_xlabel('$x_1$')
axes[0].set_ylabel('$x_2$')
axes[0].set_title('Interpolación de $f_b(x_1, x_2)$ con Splines Cúbicos (P. Equiespaciados)')

axes[1].plot_surface(X1_grid, X2_grid, Y_interp_nonequispaced_fb, cmap='viridis', alpha=0.8)
axes[1].set_xlabel('$x_1$')
axes[1].set_ylabel('$x_2$')
axes[1].set_title('Interpolación de $f_b(x_1, x_2)$ con Splines Cúbicos (P. No Equiespaciados)')

plt.tight_layout()
plt.show()


# ERROR 
error_equispaced_fb = np.abs(fb(X1_grid, X2_grid) - Y_interp_equispaced_fb)
error_nonequispaced_fb = np.abs(fb(X1_grid, X2_grid) - Y_interp_nonequispaced_fb)

fig_error, axes_error = plt.subplots(1, 2, figsize=(16, 6), subplot_kw={'projection': '3d'})

axes_error[0].plot_surface(X1_grid, X2_grid, error_equispaced_fb, cmap='viridis', alpha=0.8)
axes_error[0].set_xlabel('$x_1$')
axes_error[0].set_ylabel('$x_2$')
axes_error[0].set_title('Error de Interpolación (P. Equiespaciados)')

axes_error[1].plot_surface(X1_grid, X2_grid, error_nonequispaced_fb, cmap='viridis', alpha=0.8)
axes_error[1].set_xlabel('$x_1$')
axes_error[1].set_ylabel('$x_2$')
axes_error[1].set_title('Error de Interpolación (P. No Equiespaciados)')

plt.tight_layout()
plt.show()




# def lagrange_interp_3d(x, y, z, data, xi, yi, zi):
#     interp_values = np.zeros_like(xi)
#     for i in range(len(xi)):
#         for j in range(len(yi)):
#             for k in range(len(zi)):
#                 x_interp = xi[i]
#                 y_interp = yi[j]
#                 z_interp = zi[k]
#                 interpolated_value = 0.0
#                 for m in range(len(x)):
#                     for n in range(len(y)):
#                         for o in range(len(z)):
#                             L_x = lagrange(x, [x_interp])
#                             L_y = lagrange(y, [y_interp])
#                             L_z = lagrange(z, [z_interp])
#                             interpolated_value += data[m, n, o] * L_x * L_y * L_z
#                 interp_values[i, j, k] = interpolated_value
#     return interp_values

# # Suponiendo que tienes los siguientes datos
# data = z_equispaced_fb

# # Definir puntos en los que se desea interpolar
# zi = np.linspace(xb_min, xb_max, 100)

# # Interpolación con Lagrange en 3D
# interp_values = lagrange_interp_3d(x1_equispaced_fb, x2_equispaced_fb, z_equispaced_fb, data, x1_grid, x2_grid, zi)
