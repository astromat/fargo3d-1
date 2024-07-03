import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation
from PIL import Image
import glob
import re

# Cargar los datos
path = "/Users/matiasmontesinos/Simulations/fargo3d/outputs/flyby2d/"

variables_par = np.genfromtxt(path+"/variables.par",dtype={'names': ("parametros","valores"),'formats': ("|S30","|S300")}).tolist()#Lee archivo var    iable.pary convierte a string       esita un int o float
parametros_par, valores_par = [],[]                                                                                                                                                                                                                         #Reparte entre parametros y valores
for posicion in variables_par:                                                                                                                                                                                                                              #
        parametros_par.append(posicion[0].decode("utf-8"))                                                                                                                                                                                  #
        valores_par.append(posicion[1].decode("utf-8"))                                                                                                                                                                                     #

def P(parametro):
        return valores_par[parametros_par.index(parametro)]                                                                                                                                                                                 #Retorna siempre str, recordad tranformar si necesita un int o float


# Identificar el número de outputs disponibles
# Patrón para archivos gasdens*.dat con un número entero
pattern = re.compile(r'gasdens(\d+)\.dat')
# Lista de archivos que coinciden con el patrón gasdens*.dat
files = glob.glob(path + "gasdens*.dat")
# Filtrar archivos que se ajustan al patrón correcto
valid_files = [f for f in files if pattern.match(os.path.basename(f))]
# Contar el número de archivos válidos
Ntot = len(valid_files)
print(Ntot)


# Datos de la grilla
Rmin = float(P("YMIN"))
Rmax = float(P("YMAX"))

delta = 10.0

# Cargar la grilla de la simulación
domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3]

NX = len(domain_x) - 1
NY = len(domain_y) - 1

dr = (Rmax - Rmin) / NY

# Función para cargar la densidad de gas
def cargar_densidad(output_number, path):
    dens_out = np.fromfile(path + f"gasdens{output_number}.dat").reshape(NY, NX)
    return dens_out

# Leer las coordenadas del planeta desde planet0.dat
def leer_coordenadas_planeta(path, snapshot):
    planet_file0 = path + "planet0.dat"
    planet0 = np.genfromtxt(planet_file0)
    xp0 = planet0[snapshot][1]  # Coordenada x del planeta
    yp0 = planet0[snapshot][2]  # Coordenada y del planeta

    planet_file1 = path + "planet1.dat"
    planet1 = np.genfromtxt(planet_file1)
    xp1 = planet1[snapshot][1]  # Coordenada x del planeta
    yp1 = planet1[snapshot][2]  # Coordenada y del planeta

    return xp0, yp0, xp1, yp1

def Grilla_XY():
    R = 0.5 * (domain_y[1:] + domain_y[:-1])
    Phi = 0.5 * (domain_x[1:] + domain_x[:-1])
    P, R = np.meshgrid(Phi, R)
    X, Y = R * np.cos(P), R * np.sin(P)
    return X, Y

X, Y = Grilla_XY()

# Establecer límites fijos para los ejes
x_limits = (-(Rmax+delta), (Rmax+delta))
y_limits = (-(Rmax+delta), (Rmax+delta))

# Crear la carpeta para guardar los PNGs si no existe
output_dir = os.path.join(path, "gas_png")
os.makedirs(output_dir, exist_ok=True)




# Calcular vmin y vmax globales
vmin, vmax = None, None
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)
    log_dens_out = np.log10(dens_out)
    if vmin is None or vmax is None:
        vmin = log_dens_out.min()
        vmax = log_dens_out.max()
    else:
        vmin = min(vmin, log_dens_out.min())
        vmax = max(vmax, log_dens_out.max())


### graficos de tests:
# Solicitar al usuario el snapshot a usar
snapshot = int(input("Ingrese el número de snapshot: "))
dens_out = cargar_densidad(snapshot, path)
xp0, yp0, xp1, yp1 = leer_coordenadas_planeta(path, snapshot)

plt.figure(figsize=(8, 8))
mesh = plt.pcolormesh(X, Y, np.log10(dens_out), cmap='jet', shading='nearest', vmin = -4, vmax=-3)
plt.scatter(xp0, yp0, c='blue', s=50, edgecolor='black', marker='x', label='Planet Position 1')
plt.scatter(xp1, yp1, c='black', s=50, edgecolor='black', marker='x', label='Planet Position 2')
cbar = plt.colorbar(mesh, label='Log Gas Density [gr / cm$^2$]')
#cbar.set_clim(-4.72, -3.24)  # Ajusta los límites de la barra de colores
plt.xlabel('X [AU]', fontsize=14)
plt.ylabel('Y [AU]', fontsize=14)
plt.title('Gas Density Distribution', fontsize=16)
plt.xlim(x_limits)
plt.ylim(y_limits)
plt.gca().set_aspect('equal', adjustable='datalim')
#plt.axis('equal')
plt.legend()
plt.tight_layout()
plt.show()



# Cargar el último archivo de salida para calcular vmin y vmax
ultimo_output = Ntot - 1
dens_out = cargar_densidad(ultimo_output, path)
log_dens_out = np.log10(dens_out)

#vmin = log_dens_out.min()
#vmax = log_dens_out.max()

# Generar los plots para cada output
for output_number in range(Ntot):
    dens_out = cargar_densidad(output_number, path)
    xp0, yp0, xp1, yp1 = leer_coordenadas_planeta(path, output_number)
    plt.figure(figsize=(8, 8))
    mesh = plt.pcolormesh(X, Y, np.log10(dens_out), cmap='jet', shading='nearest', vmin=vmin, vmax=vmax)
    plt.scatter(xp0, yp0, c='blue', s=50, edgecolor='black', marker='x', label='Planet Position 1')
    plt.scatter(xp1, yp1, c='blue', s=50, edgecolor='black', marker='x', label='Planet Position 2')
    plt.colorbar(mesh, label='Log Gas Density [gr / cm$^2$]')
    plt.xlabel('X [AU]', fontsize=14)
    plt.ylabel('Y [AU]', fontsize=14)
    plt.title(f'Gas Density Distribution - Output {output_number}', fontsize=16)
    plt.xlim(x_limits)
    plt.ylim(y_limits)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.legend()
    plt.tight_layout()
    
    # Guardar el PNG en la carpeta especificada
    file_path = os.path.join(output_dir, f"output_{output_number:04d}.png")
    plt.savefig(file_path)
    plt.close()

### Animar los PNGs
# Ruta donde están almacenados los archivos PNG
png_dir = os.path.join(path, "gas_png")

# Obtener una lista de todos los archivos PNG
png_files = [os.path.join(png_dir, f"output_{i:04d}.png") for i in range(Ntot)]

# Crear la figura para la animación
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

# Ocultar los ejes y eliminar los márgenes
ax.axis('off')
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

# Función para cargar las imágenes
def load_image(file):
    img = Image.open(file)
    return [plt.imshow(img, animated=True)]

# Crear la animación
frames = [load_image(file) for file in png_files]
ani = animation.ArtistAnimation(fig, frames, interval=100, repeat_delay=1000, blit=True)

# Guardar la animación en un archivo MP4
ani.save(png_dir + 'gas_density_animation_from_png.mp4', fps=5, extra_args=['-vcodec', 'libx264'])

plt.show()
