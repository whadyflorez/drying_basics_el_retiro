
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog

# Datos
proveedores = ['P1', 'P2', 'P3', 'P4', 'P5']
hornos = ['H1', 'H2']
costos = np.array([
    [4, 6],   # P1
    [5, 4],   # P2
    [3, 5],   # P3
    [6, 3],   # P4
    [7, 4]    # P5
])
oferta = [30, 20, 25, 15, 10]
demanda = [60, 40]

# Formulación
c = costos.flatten()
A_ub = np.zeros((len(proveedores), len(proveedores) * len(hornos)))
for i in range(len(proveedores)):
    for j in range(len(hornos)):
        A_ub[i, i * len(hornos) + j] = 1
b_ub = oferta
A_eq = np.zeros((len(hornos), len(proveedores) * len(hornos)))
for j in range(len(hornos)):
    for i in range(len(proveedores)):
        A_eq[j, i * len(hornos) + j] = 1
b_eq = demanda
bounds = [(0, None)] * (len(proveedores) * len(hornos))

# Optimización
res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
envios = res.x.reshape(len(proveedores), len(hornos))
print("Envíos óptimos (toneladas):\n", envios)
print("Costo total mínimo:", res.fun)

# Gráfica: Esquema del supply chain
plt.figure(figsize=(10, 6))
for i, proveedor in enumerate(proveedores):
    plt.plot([0, 1], [i, i], 'o-', color='brown')
    plt.text(-0.1, i, proveedor, ha='right', va='center', fontsize=10, fontweight='bold')
for j, horno in enumerate(hornos):
    plt.plot([2, 3], [j + 1, j + 1], 'o-', color='orange')
    plt.text(3.1, j + 1, horno, ha='left', va='center', fontsize=10, fontweight='bold')
for i in range(len(proveedores)):
    for j in range(len(hornos)):
        cantidad = envios[i, j]
        if cantidad > 0:
            plt.arrow(1, i, 1, (j + 1 - i) * 0.9, head_width=0.2, head_length=0.2,
                      length_includes_head=True, color='gray')
            plt.text(1.5, (i + j + 1) / 2, f"{cantidad:.0f}", ha='center', va='center', fontsize=9)
plt.xlim(-0.5, 3.5)
plt.ylim(-1, max(len(proveedores), len(hornos)) + 1)
plt.axis('off')
plt.title("Esquema de la cadena de suministro de madera")
plt.tight_layout()
plt.show()

# Gráfica: Diagrama de barras
bar_width = 0.35
indices = np.arange(len(proveedores))
plt.figure(figsize=(10, 6))
plt.bar(indices - bar_width / 2, envios[:, 0], bar_width, label='H1', color='saddlebrown')
plt.bar(indices + bar_width / 2, envios[:, 1], bar_width, label='H2', color='peru')
plt.xlabel('Proveedores')
plt.ylabel('Toneladas enviadas')
plt.title('Distribución de madera desde proveedores a hornos')
plt.xticks(indices, proveedores)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
