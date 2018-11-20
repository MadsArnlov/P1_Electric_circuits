# =============================================================================
# Modelling of RC-circuit
# =============================================================================
"The program uses the modules matplotlib and numpy"
import matplotlib.pyplot as plt
import numpy as np
import sys

# =============================================================================
# Data
# =============================================================================
"The data is exported from 'Waveform'."
"Voltage across capacitor in time domain"
data_time = np.genfromtxt("forsøgsdata.csv", delimiter=",")

time = data_time[1:, 0]
sqwave_data = data_time[1:, 1]
cap = data_time[1:, 2]

"'time_charge' is the measured time before the capacitor starts charging."
time_charge = (time[5001] - time[0])

"Magnitude of capacitor/resistor in frequency domain"
data_magnitude_C = np.genfromtxt("bode_low_1v4753ohm96_73nf.csv", delimiter=",")

frequency_C = data_magnitude_C[1:, 0]
angular_frequency_C = frequency_C*2*np.pi
magnitude_C = data_magnitude_C[1:, 2]
phase_C = data_magnitude_C[1:, 3]

data_magnitude_R = np.genfromtxt("bode_high_1v4753ohm96_73nf.csv", delimiter=",")

frequency_R = data_magnitude_R[1:, 0]
angular_frequency_R = frequency_R*2*np.pi
magnitude_R = data_magnitude_R[1:, 2]
phase_R = data_magnitude_R[1:, 3]

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================
"The values of the components used in the experiment are defined:"
R = 4753                                            # Resistance
C = 96.73 * 10**(-9)                                # Capacitance
tau = R*C                                           # Time constant
omega_c = 1/tau                                     # Cutoff frequency
V0 = 2                                              # Initial voltage

R2 = R*10
tau2 = R2*C

"""
The measurements for discharge/charge has the duration of 'time_charge',
and these measurements each have 5000 points. To simulate this accordingly,
an arary for time is made, t, which has the length of 5000, uniform spaced
values between 0 and 'time_charge'.
"""
t = np.linspace(0, time_charge, 5000)
omega = np.linspace(angular_frequency_C[0], angular_frequency_C[-1], 5000)


def V_C(t):
    """
    Returns the voltage of the capacitor as it discharges through the resistor.

    The formula is derived in chapter "RC-time analysis", section
    "Zero-input response":
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
    "Complete response":
        Vc(t) = V + (V0 - V)e^(-t/tau)
    The data is offset by V0/2, for which the model is modified.

    The value V, is the same as V0 from V_C(t).

    Parameters
    ----------
    t: float
        The time starting from 0
    """
    return V0 + (0 - V0)*np.exp(-t/tau) - V0/2


def V_R(t):
    """
    Returns the voltage of the resistor in a complete response circuit.

    The formula is derived in chapter "RC time analysis", section
    "Complete response":
        VR(t) = -(V0 - V)e^(-t/tau)
    The data is offset by V0/2, for which the model is modified.

    The value V, is the same as V0 from V_C(t).

    Parameters
    ----------
    t: float
        The time starting from 0
    """
    return -(0 - V0)*np.exp(-t/tau) - V0/2


def H_lp(omega):
    """
    Returns the transfer function of a first order Low Pass filter.

    The formula is derived in chapter "Passive Analogue Filters", section
    "Low pass filter":
        H(s) = 1/sqrt(s^2tau^2 + 1)

    Parameters
    ----------
    s: float
        The complex frequency
    """
    return 1/np.sqrt(1 + omega**2 * tau**2)


def H_hp(omega):
    """
    Returns the transfer function of a first order High Pass Filter.

    The formula is derived in chapter "Passive Analogue Filters", section
    "High pass filter":
        H(s) = 1/sqrt(1 + (1/s*tau)^2)

    Parameters
    ----------
    s: float
        The complex frequency
    """
    return np.sqrt(omega**2*tau**2)/np.sqrt(1 + omega**2*tau**2)


def H_bp(omega):
    """
    Returns the transfer function of a first order Band Pass filter.

    The formula is derived in chapter "Passive Analogue Filters", section
    "Band pass filter":
        H(s) = 1/sqrt(s^2tau^2 + 1/(s^2tau^2) + 4)

    Parameters
    ----------
    s: float
        The complex frequency
    """
    return 1/np.sqrt(4 + omega**2 * tau**2 + 1/(omega**2 * tau2**2))


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
if len(sys.argv) == 2:
    if sys.argv[1] == 'data':
        plt.figure(figsize=(12, 8))
        plt.plot(time - time[0], cap, 'b,')
        plt.xlabel('$t$ [s]')
        plt.ylabel('$V_C$ [V]')
        plt.title('Data of RC-circuit')
        plt.savefig('data.png')

# =============================================================================
# Plot of mathematical model
# =============================================================================

    if sys.argv[1] == 'model':
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

    if sys.argv[1] == 'data_vs_model':
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

    if sys.argv[1] == 'deviation':
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

    if sys.argv[1] == 'RCLP':
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20*np.log10(H_lp(omega)),
                     'b-', label='LP Transfer function')
        plt.semilogx(omega_c, 20*np.log10(H_lp(omega_c)),
                     'kx', label='Gain at $\omega=\omega_c$')
        plt.plot(angular_frequency_C, magnitude_C,
                 'r-', label='Data')
        plt.legend()
        plt.grid(True)
        plt.ylabel('Magnitude $G(j\omega)$ [dB]')


        plt.subplot(2, 1, 2)
        plt.semilogx(omega, np.arctan(-omega*tau)*180/np.pi,
                     'b-', label='LP Phase shift')
        plt.semilogx(omega_c, np.arctan(-omega_c*tau)*180/np.pi,
                     'kx', label='Phase at $\omega=\omega_c$')
        plt.plot(angular_frequency_C, phase_C,
                     'r-', label='Data')
        plt.yticks(np.arange(0, -105, step=-15))
        plt.legend()
        plt.grid(True)
        plt.xlabel('Angular frequency $\omega$ [Hz]')
        plt.ylabel('Phase $\u03B8$ [degrees]')
        
        plt.savefig('data_bodeplots_rc_lp.png')
        plt.show()

# =============================================================================
# Bodeplot of RC HP filter
# =============================================================================

    if sys.argv[1] == 'RCHP':
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20*np.log10(H_hp(omega)),
                     'b-', label='HP Transfer function')
        plt.semilogx(1/tau, 20*np.log10(H_hp(1/tau)),
                     'kx', label='Gain at $\omega_c$')
        plt.plot(angular_frequency_R, magnitude_R,
                 'r-', label='Data')
        plt.legend()
        plt.grid(True)
        plt.ylabel('Magnitude $G(j\omega)$ [dB]')

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, np.arctan(1/(omega*tau))*180/np.pi,
                     'b-', label='HP Phase shift')
        plt.semilogx(omega_c, np.arctan(1/(omega_c*tau))*180/np.pi,
                     'kx', label='Phase at $\omega_c$')
        plt.plot(angular_frequency_R, phase_R,
                 'r-', label='Data')
        plt.yticks(np.arange(0, 105, step=15))
        plt.legend()
        plt.grid(True)
        plt.xlabel('Angular frequency $\omega$ [Hz]')
        plt.ylabel('Phase $\u03B8$ [degree]')

        plt.savefig('data_bodeplots_rc_hp.png')
        plt.show()
