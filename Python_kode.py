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

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================

R = 4732
C = 98.98 * 10**(-9)
tau = R*C
V0 = 2

V = 4

time_charge = (4*10**-6 - tid[0])

t_vc = np.linspace(0, time_charge, 5000)

V_C = np.exp(-t_vc/tau)*V0 - V0/2

t_vt = np.linspace(0*tau, time_charge, 5000)

V_t = V0 + (V0 - V)*np.exp(-t_vt/tau) - V0/2

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
# Plotting
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t_sqwave, sqwave, 'y-',
         t_vc, V_C, 'c-',
         time_charge + t_vt, V_t, 'c-',
         tid - tid[0], cap, 'k,',
         tau, np.exp(-1)*V0 - V0/2, 'mx')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Voltage across capacitor')
# plt.savefig('Afladning_kapacitor.pdf')

plt.figure(figsize=(12, 8))
plt.plot(t_vc, cap[:5000] - V_C, 'k,',
         time_charge + t_vt, cap[5000:] - V_t, 'k,')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_{data} - V_C$ [V]')
plt.title('Difference between data and model')
# plt.savefig('Forskel_data_model.pdf')
