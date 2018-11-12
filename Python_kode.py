# =============================================================================
# Modelling of RC-circuit
# =============================================================================
"The program uses the modules matplotlib and numpy"
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Data
# =============================================================================
"The data is exported from 'Waveform'."
data = np.genfromtxt("forsøgsdata.csv", delimiter=",")

time = data[1:, 0]
sqwave_data = data[1:, 1]
cap = data[1:, 2]

"'time_charge' is the measured time before the capacitor starts charging."
time_charge = (4*10**-6 - time[0])

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================
"The values of the components used in the experiment are defined:"
R = 4732
C = 98.98 * 10**(-9)
tau = R*C
V0 = 2

K = -4

"""
The measurements for discharge/charge has the duration of 'time_charge',
and these measurements each have 5000 points. To simulate this accordingly,
an arary for time is made, t, which has the length of 5000, uniform spaced
values between 0 and 'time_charge'.
"""
t = np.linspace(0, time_charge, 5000)
s = np.linspace(1, 100/tau, 5000)


def V_C(t):
    """
    Returns the voltage of the capacitor as it discharges through the resistor.

    The formula is derived in chapter "RC-time analysis", section
    "Zero state":
        Vc(t) = e^(-t/tau)V0
    The data is offset by V0/2, for which the model is modified.

    Parameters
    ----------
    t: float
        The time after zero state
    """
    return np.exp(-t/tau)*V0 - V0/2


def V_C2(t):
    """
    Returns the voltage of the capacitor as it charges.

    The formula is derived in chapter "RC time analysis", section
    "Steady state":
        Vc(t) = V + e^(-t/tau)(V + K)
    The data is offset by V0/2, for which the model is modified.

    The value V, is the same as V0 from the earlier function.

    Parameters
    ----------
    t: float
        The time as the circuit enters steady state
    """
    return V0 + np.exp(-t/tau)*(V0 + K) - V0/2


def H_lp(s):
    """
    Returns the transfer function of a first order Low Pass filter.

    The formula is derived in chapter "Passive Analogue Filters", section
    "Bodeplot":
        H(s) = 1/sqrt(s^2tau^2 + 1)

    Parameters
    ----------
    s: float
        The complex frequency
    """
    return 1/np.sqrt(s**2*tau**2 + 1)


# =============================================================================
# Square wave
# =============================================================================
"A piecewise function is a simple way to create a square wave."
t_sqwave = np.linspace(0, time_charge*2, 1000)

sqwave = np.piecewise(t_sqwave,                 # The parameter
                      [t_sqwave > time_charge,  # The first condition
                       t_sqwave < time_charge,  # The second condition
                       t_sqwave == 0],          # The third condition
                      [V0/2, -V0/2, V0/2]       # The corresponding values
                      )

# =============================================================================
# Plot of data
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(time - time[0], cap, 'b,')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Data of RC-circuit')
plt.savefig('data.png')

# =============================================================================
# Plot of mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t, V_C(t), 'tab:orange')
plt.plot(time_charge + t, V_C2(t), 'tab:orange')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.title('Mathematical model')
plt.savefig('mathematical_model.png')

# =============================================================================
# Plot of data and mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t, V_C(t), 'tab:orange', label='Modelled voltage')
plt.plot(time_charge + t, V_C2(t), 'tab:orange')
plt.plot(time - time[0], cap, 'b,', label='Data points')
plt.plot(tau, V_C(tau), 'kx', label='$V_C(tau)$')
plt.plot(t_sqwave, sqwave, 'k--', label='Square wave')
plt.xlabel('$t$ [s]')
plt.ylabel('$V_C$ [V]')
plt.legend()
plt.title('Voltage across capacitor')
plt.savefig('data_vs_model.png')

# =============================================================================
# Plot of deviation between data and mathematical model
# =============================================================================

plt.figure(figsize=(12, 8))
plt.plot(t, cap[:5000] - V_C(t), 'k,',
         time_charge + t, cap[5000:] - V_C2(t), 'k,')
plt.xlabel('t [s]')
plt.ylabel('$V_{data} - V_C$ [V]')
plt.title('Difference between data and model')
plt.savefig('deviation.png')

# =============================================================================
# Bodeplot of RC LP filter
# =============================================================================

plt.figure(figsize=(12, 8))
plt.subplot(2, 1, 1)
plt.semilogx(s, 20*np.log10(H_lp(s)), 'b-', label='LP Transfer function')
plt.semilogx(1/tau, 20*np.log10(H_lp(1/tau)), 'kx', label='Gain at time=1/tau')
plt.legend()
plt.ylabel('Magnitude [dB]')
plt.title('Bodeplot of RC LP filter')

plt.subplot(2, 1, 2)
plt.semilogx(s, np.arctan(1/(s*tau))*180/np.pi, 'b-', label='LP Phase shift')
plt.yticks(np.arange(0, 105, step=15))
plt.legend()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [degrees]')

plt.show()
