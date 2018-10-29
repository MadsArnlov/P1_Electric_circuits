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

t = np.linspace(0, 10*tau, 1000)

V_C = np.exp(-t/tau)*V0 - V0/2

# =============================================================================
# Square wave
# =============================================================================

sqwave = np.piecewise(t, [t == 0, t > 0], [V0/2, -V0/2])

# =============================================================================
# Plotting
# =============================================================================

plt.plot(t, V_C, 'r-', t, sqwave, 'b-')
