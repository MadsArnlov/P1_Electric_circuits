# =============================================================================
# Modelling of RC-circuit
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np

R = 4732
C = 98.98 * 10**(-9)
tau = R*C
V0 = 2

t = np.linspace(0, 10*tau, 1000)

V_C = np.exp(-t/tau)*V0 - V0/2

plt.plot(t, V_C)