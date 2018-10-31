# =============================================================================
# Modelling of RC-circuit
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Data
# =============================================================================

data = np.genfromtxt("forsøgsdata.csv", delimiter=",")

tid = data[1:, 0]
sqwave_data = data[1:, 1]
cap = data[1:, 2]

time_charge = (4*10**-6 - tid[0])

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================

R = 4732
C = 98.98 * 10**(-9)
tau = R*C
V0 = 2

K = -4

t_vc = np.linspace(0, time_charge, 5000)

t_vt = np.linspace(0*tau, time_charge, 5000)

def V_C(t):
    return np.exp(-t/tau)*V0 - V0/2

def V_C2(t):
    return V0 + np.exp(-t/tau)*(V0 + K) - V0/2

# =============================================================================
# Square wave
# =============================================================================

t_sqwave = np.linspace(0, time_charge*2, 1000)

sqwave = np.piecewise(t_sqwave,
                      [t_sqwave > time_charge,
                       t_sqwave < time_charge,
                       t_sqwave == 0],
                      [V0/2, -V0/2, V0/2]
                      )

# =============================================================================
# Plot of data
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(tid - tid[0], cap, 'b,')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Data of RC-circuit')
# plt.savefig('Data.pdf')

# =============================================================================
# Plot of mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t_vc, V_C(t_vc), 'tab:orange')
plt.plot(time_charge + t_vt, V_C2(t_vt), 'tab:orange')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Mathematical model')
# plt.savefig('Mathematical_model.pdf')

# =============================================================================
# Plot of data and mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t_vc, V_C(t_vc), 'tab:orange', label='Modelled voltage')
plt.plot(time_charge + t_vt, V_C2(t_vt), 'tab:orange')
plt.plot(tid - tid[0], cap, 'b,', label='Data points')
plt.plot(tau, np.exp(-1)*V0 - V0/2, 'kx', label='$V_C(tau)$')
plt.plot(t_sqwave, sqwave, 'k--', label='Square wave')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.legend()
plt.title('Voltage across capacitor')
# plt.savefig('Afladning_kapacitor.pdf')

# =============================================================================
# Plot of deviation between data and mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t_vc, cap[:5000] - V_C(t_vc), 'k,',
         time_charge + t_vt, cap[5000:] - V_C2(t_vt), 'k,')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_{data} - V_C$ [V]')
plt.title('Difference between data and model')
# plt.savefig('Forskel_data_model.pdf')
