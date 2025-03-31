
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI

# Constantes
A = 0.2  # m²
m_s = 5  # kg
c_madera = 2500  # J/kg·K
T_inf = 80  # °C
T_i = 25  # °C
h_m = 0.002  # kg/m²·s
t = 3600  # s (1 hora)
R = 461  # J/kg·K
X_i = 0.30  # humedad base seca inicial

# Función para obtener presión de saturación y entalpía de vaporización desde CoolProp
def obtener_propiedades_agua(T_C):
    T_K = T_C + 273.15
    P_sat = PropsSI('P', 'T', T_K, 'Q', 0, 'Water')  # Pa
    h_liq = PropsSI('H', 'T', T_K, 'Q', 0, 'Water')  # J/kg
    h_vap = PropsSI('H', 'T', T_K, 'Q', 1, 'Water')  # J/kg
    L_v = h_vap - h_liq
    return P_sat, L_v

# Concentración de vapor a la temperatura de la tabla
def C_s(T):
    P_sat, _ = obtener_propiedades_agua(T)
    return P_sat / (R * (T + 273.15))  # kg/m³

# Función para resolver el balance de energía para una h dada
def calcular_T_Y_Xf(h):
    def balance(T):
        cs = C_s(T)
        P_sat, L_v = obtener_propiedades_agua(T)
        m_evap = h_m * A * cs * t
        Q_evap = m_evap * L_v
        Q_calor = m_s * c_madera * (T - T_i)
        Q_aire = h * A * (T_inf - T) * t
        return Q_aire - (Q_calor + Q_evap)

    T_sol = fsolve(balance, 50)[0]
    cs = C_s(T_sol)
    P_sat, L_v = obtener_propiedades_agua(T_sol)
    m_evap = h_m * A * cs * t
    X_f = (X_i * m_s - m_evap) / m_s
    return T_sol, max(X_f, 0)

# Rango de coeficientes h
h_values = np.linspace(5, 50, 100)
T_results = []
Xf_results = []

for h in h_values:
    T, Xf = calcular_T_Y_Xf(h)
    T_results.append(T)
    Xf_results.append(Xf * 100)

# Gráficas
fig, ax1 = plt.subplots()
ax1.set_xlabel("h [W/m²·K]")
ax1.set_ylabel("Temperatura final [°C]", color='tab:red')
ax1.plot(h_values, T_results, color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax2 = ax1.twinx()
ax2.set_ylabel("Humedad final [% bs]", color='tab:blue')
ax2.plot(h_values, Xf_results, color='tab:blue', linestyle='--')
ax2.tick_params(axis='y', labelcolor='tab:blue')
plt.title("Temperatura y humedad final vs coeficiente h")
fig.tight_layout()
plt.show()
