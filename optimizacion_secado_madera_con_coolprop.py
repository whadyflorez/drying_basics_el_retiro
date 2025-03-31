
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve

# Constantes
A = 0.2  # m²
m_s = 5  # kg
c_madera = 2500  # J/kg·K
# L_v será calculado dinámicamente con CoolProp
T_i = 25  # °C
h_m = 0.002  # kg/m²·s
t = 3600  # s
R = 461  # J/kg·K
X_i = 0.30  # humedad base seca inicial

from CoolProp.CoolProp import PropsSI

def P_sat(T):
    T_K = T + 273.15
    return PropsSI('P', 'T', T_K, 'Q', 0, 'Water') / 1000  # kPa

def C_s(T):
    return P_sat(T) * 1000 / (R * (T + 273))  # kg/m³

def calcular_T_Y_Xf(h, T_inf):
    def balance(T):
        cs = C_s(T)
        m_evap = h_m * A * cs * t
        T_K = T + 273.15
        h_liq = PropsSI('H', 'T', T_K, 'Q', 0, 'Water')
        h_vap = PropsSI('H', 'T', T_K, 'Q', 1, 'Water')
        L_v = h_vap - h_liq
        Q_evap = m_evap * L_v
        Q_calor = m_s * c_madera * (T - T_i)
        Q_aire = h * A * (T_inf - T) * t
        return Q_aire - (Q_calor + Q_evap)

    T_tabla = fsolve(balance, 50)[0]
    cs = C_s(T_tabla)
    m_evap = h_m * A * cs * t
    X_f = (X_i * m_s - m_evap) / m_s
    return T_tabla, max(X_f, 0), m_evap

def objective(x):
    h, T_inf = x
    T_tabla, _, _ = calcular_T_Y_Xf(h, T_inf)
    Q_air = h * A * (T_inf - T_tabla) * t
    return Q_air / 1000  # kJ

def humidity_constraint(x):
    h, T_inf = x
    _, Xf, _ = calcular_T_Y_Xf(h, T_inf)
    return 0.25 - Xf

bounds = [(10, 50), (30, 90)]
x0 = [30, 70]
constraints = {'type': 'ineq', 'fun': humidity_constraint}

result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints)
opt_h, opt_Tinf = result.x
opt_Ttabla, opt_Xf, opt_mevap = calcular_T_Y_Xf(opt_h, opt_Tinf)
opt_Q = objective(result.x)

# Malla de alta resolución
h_vals = np.linspace(10, 50, 100)
Tinf_vals = np.linspace(40, 95, 100)
H, TINF = np.meshgrid(h_vals, Tinf_vals)
Q_vals = np.zeros_like(H)
Xf_vals = np.zeros_like(H)

for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        Ttabla, Xf, _ = calcular_T_Y_Xf(H[i, j], TINF[i, j])
        Xf_vals[i, j] = Xf
        if Xf <= 0.25:
            Q_vals[i, j] = objective([H[i, j], TINF[i, j]])
        else:
            Q_vals[i, j] = np.nan

# Gráfico de contorno mejorado
plt.figure(figsize=(12, 7))
cp = plt.contourf(H, TINF, Q_vals, levels=40, cmap='viridis')
plt.colorbar(cp, label='Energía total entregada (kJ)')
plt.scatter(opt_h, opt_Tinf, color='red', label='Óptimo', edgecolors='black')
plt.xlabel('Coef. de convección h [W/m²·K]')
plt.ylabel('Temperatura del aire T_inf [°C]')
plt.title('Mapa ampliado de energía entregada (restricción: humedad ≤ 25%)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Resultados
print(f"Óptimo encontrado:")
print(f"  h = {opt_h:.2f} W/m²·K")
print(f"  T_inf = {opt_Tinf:.2f} °C")
print(f"  T_tabla = {opt_Ttabla:.2f} °C")
print(f"  Humedad final = {opt_Xf*100:.2f} %")
print(f"  Agua evaporada = {opt_mevap:.3f} kg")
print(f"  Energía entregada = {opt_Q:.2f} kJ")


# Gráfico mejorado que cubre toda la región factible con colores continuos
import matplotlib.cm as cm
import matplotlib.colors as mcolors

Q_vals_filled = np.copy(Q_vals)
max_val = np.nanmax(Q_vals_filled)
Q_vals_filled[np.isnan(Q_vals_filled)] = max_val + 100

#mask = np.where(Xf_vals <= 0.25, 1, np.nan)

plt.figure(figsize=(12, 7))
cmap = cm.viridis
norm = mcolors.Normalize(vmin=np.nanmin(Q_vals), vmax=max_val)

plt.pcolormesh(H, TINF, Q_vals_filled, cmap=cmap, norm=norm, shading='auto')
plt.colorbar(label='Energía total entregada (kJ)')
plt.contour(H, TINF, Xf_vals, levels=[0.25], colors='white', linewidths=1.5, linestyles='--')
plt.scatter(opt_h, opt_Tinf, color='red', label='Óptimo', edgecolors='black')
plt.xlabel('Coef. de convección h [W/m²·K]')
plt.ylabel('Temperatura del aire T_inf [°C]')
plt.title('Mapa completo de energía entregada (humedad ≤ 25%)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
