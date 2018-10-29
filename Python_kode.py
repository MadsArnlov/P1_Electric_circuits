# =============================================================================
# Modelling of RC-circuit
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================

R = 4732
C = 98.98 * 10**(-9)
tau = R*C
V0 = 2

V = 4

t_vc  = np.linspace(0, 10*tau, 1000)

V_C = np.exp(-t/tau)*V0 - V0/2

V_t = V0 + (V0 - V)*np.exp(-t/tau) - V0/2

t_vt = np.linspace(10*tau, 20*tau, 1000)

# =============================================================================
# Square wave
# =============================================================================

t_sqwave = np.linspace(0, 20*tau, 1000)

sqwave = np.piecewise(t, [t > 10*tau, t < 10*tau, t == 0], [V0/2, -V0/2, V0/2])

# =============================================================================
# Plotting
# =============================================================================

plt.figure(figsize = (12, 8))
plt.plot(t_vc, V_C, 'r-', t_sqwave, sqwave, 'b-', t_vt, V_t, 'g-')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Voltage across capacitor')
# plt.savefig('Afladning_kapacitor.pdf')