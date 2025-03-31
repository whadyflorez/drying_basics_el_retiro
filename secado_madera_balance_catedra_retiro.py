
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constantes
A = 0.2  # m²
m_s = 5  # kg
c_madera = 2500  # J/kg·K
L_v = 2260000  # J/kg
T_inf = 80  # °C
T_i = 25  # °C
h_m = 0.002  # kg/m²·s
t = 3600  # s (1 hora)
R = 461  # J/kg·K
X_i = 0.30  # humedad base seca inicial

# Aproximación de P_sat(T) en kPa: lineal entre 40°C y 80°C
def P_sat(T):
    return 0.996 * T - 32.45  # en kPa

# Concentración de vapor a la temperatura de la tabla
def C_s(T):
    return P_sat(T) * 1000 / (R * (T + 273))  # convertir kPa a Pa

# Función para resolver el balance de energía para una h dada
def calcular_T_Y_Xf(h):
    def balance(T):
        cs = C_s(T)
        m_evap = h_m * A * cs * t
        Q_evap = m_evap * L_v
        Q_calor = m_s * c_madera * (T - T_i)
        Q_aire = h * A * (T_inf - T) * t
        return Q_aire - (Q_calor + Q_evap)

    T_sol = fsolve(balance, 50)[0]
    cs = C_s(T_sol)
    m_evap = h_m * A * cs * t
    X_f = (X_i * m_s - m_evap) / m_s
    return T_sol, max(X_f, 0)  # evitar humedad negativa

# Rango de coeficientes de convección h (W/m²·K)
h_values = np.linspace(5, 50, 100)
T_results = []
Xf_results = []

for h in h_values:
    T, Xf = calcular_T_Y_Xf(h)
    T_results.append(T)
    Xf_results.append(Xf * 100)  # convertir a porcentaje

# Primera gráfica: temperatura y humedad vs h
fig, ax1 = plt.subplots()
ax1.set_xlabel("Coeficiente de convección térmica h [W/m²·K]")
ax1.set_ylabel("Temperatura final de la tabla [°C]", color='tab:red')
ax1.plot(h_values, T_results, color='tab:red', label="Temperatura final")
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid(True)

ax2 = ax1.twinx()
ax2.set_ylabel("Humedad final [% base seca]", color='tab:blue')
ax2.plot(h_values, Xf_results, color='tab:blue', linestyle='--', label="Humedad final")
ax2.tick_params(axis='y', labelcolor='tab:blue')

plt.title("Efecto del coeficiente de convección sobre el secado de la madera")
fig.tight_layout()
plt.show()

# Segunda gráfica: solo humedad vs h
plt.figure()
plt.plot(h_values, Xf_results, color='tab:blue')
plt.xlabel("Coeficiente de convección térmica h [W/m²·K]")
plt.ylabel("Humedad final [% base seca]")
plt.title("Humedad final de la madera vs coeficiente de convección")
plt.grid(True)
plt.tight_layout()
plt.show()
