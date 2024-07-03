import numpy as np
import matplotlib.pyplot as plt

# Constantes (en unidades normalizadas)
G = 1  # Constante gravitacional
M = 1  # Masa del objeto central
tolerance = 1e-10  # Tolerancia para precisión numérica

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
    if r < 1e-8:  # Evitar singularidad
        r = 1e-8
    ax = -G * M * x / r**3
    ay = -G * M * y / r**3
    return np.array([vx, vy, ax, ay])

def calculate_orbit(x0, y0, v0, theta):
    vx0 = v0 * np.cos(theta)
    vy0 = v0 * np.sin(theta)
    
    r0 = np.sqrt(x0**2 + y0**2)
    
    epsilon = 0.5 * (vx0**2 + vy0**2) - G * M / r0
    h = x0 * vy0 - y0 * vx0  # Momento angular específico
    
    print(f"x0: {x0}, y0: {y0}, v0: {v0}, theta: {theta}")
    print(f"vx0: {vx0}, vy0: {vy0}, r0: {r0}")
    print(f"epsilon: {epsilon}, h: {h}")
    
    if abs(epsilon) < tolerance:
        # Órbita parabólica
        a = np.inf
        T_parabolica = 2 * r0 / v0
        t = np.linspace(0, 5 * T_parabolica, 2000)  # Aumentamos el tiempo para capturar mejor la trayectoria
        e = 1.0
    elif epsilon < 0:
        # Órbita elíptica
        a = -G * M / (2 * epsilon)
        T = 2 * np.pi * np.sqrt(a**3 / (G * M))
        t = np.linspace(0, 2 * T, 1000)  # Dos períodos orbitales
        e = np.sqrt(1 + (2 * epsilon * h**2) / (G**2 * M**2))
    else:
        # Órbita hiperbólica
        a = -G * M / (2 * epsilon)
        v_infinito = np.sqrt(2 * epsilon)
        T_hiperbolica = r0 / v_infinito
        t = np.linspace(0, 5 * T_hiperbolica, 2000)  # Aumentamos el tiempo para capturar mejor la trayectoria
        e = np.sqrt(1 + (2 * epsilon * h**2) / (G**2 * M**2))
    
    print(f"a: {a}, e: {e}")
    
    initial_state = np.array([x0, y0, vx0, vy0])
    
    solution = runge_kutta_4(derivative, initial_state, t)
    
    return solution, t, a, epsilon, e

# Parámetros iniciales
x0 = 30.  # Coordenada x inicial (en unidades astronómicas)
y0 = 0.  # Coordenada y inicial (en unidades astronómicas)
v0 = np.sqrt(2 * G * M / np.sqrt(x0**2 + y0**2))  # Velocidad inicial (velocidad de escape)
theta = np.pi/2  # Ángulo de lanzamiento (en radianes)

# Calcula la órbita
solution, t, a, epsilon, e = calculate_orbit(x0, y0, v0, theta)

# Extrae posiciones
x = solution[:, 0]
y = solution[:, 1]

# Gráfica de la órbita
factor_de_escala = x0
vx0_esc = factor_de_escala * v0 * np.cos(theta)
vy0_esc = factor_de_escala * v0 * np.sin(theta)


plt.figure(figsize=(10, 10))
plt.plot(x, y)
plt.plot(0, 0, 'yo', markersize=10)
plt.title(f'Órbita del objeto (e = {e:.2f})')
plt.xlabel('Posición x (UA)')
plt.ylabel('Posición y (UA)')
plt.axis('equal')
plt.grid(True)

max_dist = max(np.max(np.abs(x)), np.max(np.abs(y)))
plt.xlim(-max_dist*1.1, max_dist*1.1)
plt.ylim(-max_dist*1.1, max_dist*1.1)
plt.plot(x0, y0, 'ro') # posicion incial del planeta
plt.arrow(x0, y0, vx0_esc, vy0_esc, head_width=factor_de_escala/10, head_length=factor_de_escala/10, fc='k', ec='k')  # Ajusta los parámetros para la flecha

plt.show()
