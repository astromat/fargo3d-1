import numpy as np
import matplotlib.pyplot as plt

# Constantes (en unidades normalizadas)
G = 1  # Constante gravitacional
M = 1  # Masa del objeto central

def runge_kutta_4(func, y0, t):
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    for i in range(1, n):
        h = t[i] - t[i-1]
        k1 = func(y[i-1], t[i-1])
        k2 = func(y[i-1] + k1 * h/2., t[i-1] + h/2.)
        k3 = func(y[i-1] + k2 * h/2., t[i-1] + h/2.)
        k4 = func(y[i-1] + k3 * h, t[i-1] + h)
        y[i] = y[i-1] + (h/6.) * (k1 + 2*k2 + 2*k3 + k4)
    return y

def derivative(state, t):
    x, y, vx, vy = state
    r = np.sqrt(x**2 + y**2)
    ax = -G * M * x / r**3
    ay = -G * M * y / r**3
    return np.array([vx, vy, ax, ay])

def orbital_parameters_to_initial_state(e, v_r, v_phi):
    energy = 0.5 * (v_r**2 + v_phi**2) - G * M / 1
    p = v_phi**2 / (G * M)
    a = -G * M / (2 * energy)
    x0 = p / (1 + e)
    y0 = 0
    vx0 = v_r
    vy0 = v_phi
    return np.array([x0, y0, vx0, vy0]), a

# Parámetros orbitales
e = 1.0  # Excentricidad
v_r = 0  # Velocidad radial inicial
v_phi = np.sqrt((1 + e) * G * M / 1)  # Velocidad azimutal para la órbita deseada

initial_state, a = orbital_parameters_to_initial_state(e, v_r, v_phi)

# Cálculo del período orbital para órbitas elípticas
if e < 1:
    T = 2 * np.pi * np.sqrt(a**3 / (G * M))
    num_orbits = 2  # Número de órbitas completas
    t = np.linspace(0, num_orbits * T, 1000 * num_orbits)
else:
    # Para órbitas hiperbólicas, usamos un tiempo más largo
    t = np.linspace(0, 20, 2000)

solution = runge_kutta_4(derivative, initial_state, t)

x = solution[:, 0]
y = solution[:, 1]

plt.figure(figsize=(10, 10))
plt.plot(x, y)
plt.plot(0, 0, 'yo', markersize=10)
plt.title(f'Órbita del objeto (e = {e})')
plt.xlabel('Posición x (UA)')
plt.ylabel('Posición y (UA)')
plt.axis('equal')
plt.grid(True)

max_dist = max(np.max(np.abs(x)), np.max(np.abs(y)))
plt.xlim(-max_dist*1.1, max_dist*1.1)
plt.ylim(-max_dist*1.1, max_dist*1.1)

plt.show()