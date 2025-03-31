
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve
import numpy as np

# Datos
A = 0.2               # m² (área superficial de la tabla)
m_s = 5.0             # kg (masa seca de la tabla)
X_i = 0.30            # humedad inicial base seca
T_inf = 80.0          # °C (temperatura del aire)
T_i = 25.0            # °C (temperatura inicial de la tabla)
t = 3600              # s (tiempo de secado)

h = 12.0              # W/m²·K (coeficiente de transferencia de calor)
h_m = 0.002           # kg/m²·s (coef. de transferencia de masa)
c_madera = 2500       # J/kg·K (calor específico promedio de la madera)
R = 461               # J/kg·K (constante del gas del vapor de agua)

# Función para obtener presión de saturación y entalpía de vaporización desde CoolProp
def obtener_propiedades_agua(T_C):
    T_K = T_C + 273.15
    P_sat = PropsSI('P', 'T', T_K, 'Q', 0, 'Water')          # Pa
    h_liq = PropsSI('H', 'T', T_K, 'Q', 0, 'Water')          # J/kg
    h_vap = PropsSI('H', 'T', T_K, 'Q', 1, 'Water')          # J/kg
    L_v = h_vap - h_liq                                      # J/kg
    return P_sat, L_v

# Función a resolver: balance energético
def balance(T_C):
    T_K = T_C + 273.15
    P_sat, L_v = obtener_propiedades_agua(T_C)
    Cs = P_sat / (R * T_K)
    m_evap = h_m * A * Cs * t

    Q_entregado = h * A * (T_inf - T_C) * t
    Q_necesario = m_s * c_madera * (T_C - T_i) + m_evap * L_v

    return Q_entregado - Q_necesario

# Resolver temperatura final de la tabla
T_tabla_sol = fsolve(balance, 60)[0]

# Calcular propiedades finales
T_K = T_tabla_sol + 273.15
P_sat, L_v = obtener_propiedades_agua(T_tabla_sol)
Cs = P_sat / (R * T_K)
m_evap = h_m * A * Cs * t
X_f = (X_i * m_s - m_evap) / m_s

# Resultados
print(f"Temperatura final de la tabla: {T_tabla_sol:.2f} °C")
print(f"Presión de saturación a T: {P_sat/1000:.2f} kPa")
print(f"Entalpía de vaporización a T: {L_v/1000:.2f} kJ/kg")
print(f"Masa de agua evaporada: {m_evap:.3f} kg")
print(f"Humedad final: {X_f*100:.2f} % base seca")
