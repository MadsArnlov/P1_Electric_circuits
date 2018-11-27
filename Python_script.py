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
data_magnitude_C = np.genfromtxt("bode_low_1v4753ohm96_73nf.csv",
                                 delimiter=",")

frequency_C = data_magnitude_C[1:, 0]
angular_frequency_C = frequency_C*2*np.pi
magnitude_C = data_magnitude_C[1:, 2]
phase_C = data_magnitude_C[1:, 3]

data_magnitude_R = np.genfromtxt("bode_high_1v4753ohm96_73nf.csv",
                                 delimiter=",")

frequency_R = data_magnitude_R[1:, 0]
angular_frequency_R = frequency_R*2*np.pi
magnitude_R = data_magnitude_R[1:, 2]
phase_R = data_magnitude_R[1:, 3]

"Output sine wave across capacitor for simpel input wave"
data_sim = np.genfromtxt("outputWave_1v4755ohm96_54nF.csv", delimiter=",")

f_sim = 1000
A_sim = 1
phi_sim = 0
k_sim = 0
time_sim = data_sim[1:, 0]
sine_sim = data_sim[1:, 1]

"Output sine wave across capacitor with phase shift and amplitude"
data_out = np.genfromtxt("outputWave_3_5v30phase4755ohm96_54nF500f.csv",
                         delimiter=",")

f_out = 500
A_out = 3.5
phi_out = 30*(np.pi/180)
k_out = 0
time_out = data_out[1:, 0]
sine_out = data_out[1:, 1]

# =============================================================================
# Voltage drop across capacitor as a function of time
# =============================================================================
"The values of the components used in the experiment are defined:"
R = 4753                                            # Resistance
C = 96.73 * 10**(-9)                                # Capacitance
tau = R*C                                           # Time constant
omega_c = 1/tau                                     # Cutoff frequency
V0 = 2                                              # Initial voltage

"Values for the input voltage V(t), which is a sine wave:"
A = 1                                               # Amplitude of sine wave
phi = 0                                             # Phase of sine wave
w = 500                                             # Angular frequency sine
k = 0.000                                           # Oscilliation constant


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
    return 0 + (V0 - 0)*np.exp(-t/tau) - V0/2


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


def V(t):
    """
    Returns a sine wave as the input voltage.

    A = Amplitude
    w = angular frequency
    phi = phase

    Parameters
    ----------
    t: float
        The time for the sine wave
    """
    return A*np.sin(w*t + phi)


def V_out(t):
    """
    Returns a sine wave as the output voltage.

    A = amplitude of input voltage

    Parameters
    ----------
    t: float
        The time for the sine wave
    """
    return A*H_lp(w)[0]*np.sin(w*t + phi + H_lp(w)[1]) + k*A/tau*np.exp(-t/tau)


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

    Returns
    -------
    modulus: float
        The length of the transfer function
    argument: float
        The phase of the transfer function in radians
    """
    modulus = 1/np.sqrt(1 + omega**2 * tau**2)
    argument = np.arctan(-omega*tau)
    return modulus, argument


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

    Returns
    -------
    modulus: float
        The length of the transfer function
    argument: float
        The phase of the transfer function in radians
    """
    modulus = np.sqrt(omega**2*tau**2)/np.sqrt(1 + omega**2*tau**2)
    argument = np.arctan(1/(omega*tau))
    return modulus, argument


# =============================================================================
# Plot of data
# =============================================================================
if len(sys.argv) >= 2:
    if sys.argv[1] == 'data':
        plt.figure(figsize=(12, 8))
        plt.plot(time - time[0], cap, 'b,', label='Data')
        plt.plot(time - time[0], sqwave_data, 'k--', label='Square wave')
        plt.xlabel('$t$ [s]')
        plt.ylabel('$V_C$ [V]')
        plt.legend()
        plt.savefig('data.pdf')
        plt.show()

# =============================================================================
# Plot of mathematical model
# =============================================================================

    if sys.argv[1] == 'model':
        plt.figure(figsize=(12, 8))
        plt.plot(t, V_C(t), 'tab:orange')
        plt.plot(time_charge + t, V_C2(t), 'tab:orange',
                 label='Mathematical model')
        plt.plot(time - time[0], sqwave_data, 'k--', label='Square wave')
        plt.xlabel('$t$ [s]')
        plt.ylabel('$V_C$ [V]')
        plt.legend()
        plt.savefig('mathematical_model.pdf')
        plt.show()

# =============================================================================
# Plot of data and mathematical model
# =============================================================================

    if sys.argv[1] == 'data_vs_model' or sys.argv[1] == 'model_vs_data':
        plt.figure(figsize=(12, 8))
        plt.plot(t, V_C(t), 'tab:orange', label='Modelled voltage')
        plt.plot(time_charge + t, V_C2(t), 'tab:orange')
        plt.plot(time - time[0], cap, 'b,', label='Data points')
        plt.plot(tau, V_C(tau), 'kx', label='$V_C(tau)$')
        plt.plot(time - time[0], sqwave_data, 'k--', label='Square wave')
        plt.xlabel('$t$ [s]')
        plt.ylabel('$V_C$ [V]')
        plt.legend()

        plt.savefig('data_vs_model.pdf')
        plt.show()

# =============================================================================
# Plot of deviation between data and mathematical model
# =============================================================================

    if sys.argv[1] == 'deviation':
        plt.figure(figsize=(12, 8))
        plt.plot(t, cap[:5000] - V_C(t), 'k,',
                 time_charge + t, cap[5000:] - V_C2(t), 'k,')
        plt.xlabel('t [s]')
        plt.ylabel('$V_{data} - V_C$ [V]')

        plt.savefig('deviation.pdf')
        plt.show()

# =============================================================================
# Fucked up error grafer - HJÆLP
# =============================================================================

    if sys.argv[1] == 'HELP':
        plt.figure(figsize=(12, 8))
        plt.plot(t, ((cap[:5000] - V_C(t))*100)/V_C(t), 'k.',
                 time_charge + t, ((cap[5000:] - V_C(t))*100)/V_C(t), 'k.')
        plt.plot(t[176], V_C(t[176]), 'bD')
        plt.xlabel('t [s]')
        plt.ylabel('% Difference')
        plt.title('% Deviation')
        plt.show()
        a = ((cap[:5000] - V_C(t))*100)/V_C(t)
        b = ((V_C(time_charge + t) - cap[5000:])*100)/V_C(time_charge + t)
        c = []
        c.append(np.argwhere(a >= 10))
        c.append(np.argwhere(a <= -10))
        sum_a = sum(abs(a)/len(a))
        sum_b = sum(abs(b)/len(b))
        print(sum_a, "\n", sum_b)

# =============================================================================
# Bodeplot of RC LP filter
# =============================================================================

    if sys.argv[1] == 'RCLP':
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20*np.log10(H_lp(omega)[0]),
                     'b-', label='LP Transfer function')
        plt.plot(angular_frequency_C, magnitude_C,
                 'r-', label='Data')
        plt.semilogx(omega_c, 20*np.log10(H_lp(omega_c)[0]),
                     'kx', label='Gain at $\omega=\omega_c$')
        plt.legend()
        plt.grid(True)
        plt.ylabel('Magnitude $G(j\omega)$ [dB]')

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, H_lp(omega)[1]*180/np.pi,
                     'b-', label='LP Phase shift')
        plt.plot(angular_frequency_C, phase_C,
                 'r-', label='Data')
        plt.semilogx(omega_c, H_lp(omega_c)[1]*180/np.pi,
                     'kx', label='Phase at $\omega=\omega_c$')
        plt.yticks(np.arange(0, -105, step=-15))
        plt.legend()
        plt.grid(True)
        plt.xlabel('Angular frequency $\omega$ [Hz]')
        plt.ylabel('Phase $\u03B8$ [degrees]')

        plt.savefig('data_bodeplots_rc_lp.pdf')
        plt.show()

# =============================================================================
# Bodeplot of RC HP filter
# =============================================================================

    if sys.argv[1] == 'RCHP':
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20*np.log10(H_hp(omega)[0]),
                     'b-', label='HP Transfer function')
        plt.plot(angular_frequency_R, magnitude_R,
                 'r-', label='Data')
        plt.semilogx(omega_c, 20*np.log10(H_hp(omega_c)[0]),
                     'kx', label='Gain at $\omega_c$')
        plt.legend()
        plt.grid(True)
        plt.ylabel('Magnitude $G(j\omega)$ [dB]')

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, H_hp(omega)[1]*180/np.pi,
                     'b-', label='HP Phase shift')
        plt.plot(angular_frequency_R, phase_R,
                 'r-', label='Data')
        plt.semilogx(omega_c, H_hp(omega_c)[1]*180/np.pi,
                     'kx', label='Phase at $\omega_c$')
        plt.yticks(np.arange(0, 105, step=15))
        plt.legend()
        plt.grid(True)
        plt.xlabel('Angular frequency $\omega$ [Hz]')
        plt.ylabel('Phase $\u03B8$ [degree]')

        plt.savefig('data_bodeplots_rc_hp.pdf')
        plt.show()

# =============================================================================
# Simulation of output sine wave
# =============================================================================

    if sys.argv[1] == 'sine':
        if len(sys.argv) == 3:
            w = eval(sys.argv[2])
        else:
            w = 100*2*np.pi
        t = np.linspace(0, 5/(w/(2*np.pi)), 5000)
        tmax = 3/(w/(2*np.pi))
        A = 1
        phi = 0
        k = 0
        plt.figure(figsize=(12, 8))
        plt.plot(t, V(t), 'b-', label='Input')
        plt.plot(t, V_out(t), 'r-', label='Output')
        plt.axhline(A*(1/np.sqrt(2)), label='Cutoff')

    if sys.argv[1] == 'sine_sim':
        w = f_sim*2*np.pi
        t = np.linspace(0, 10/(w/(2*np.pi)), 5000)
        tmax = 3/(w/(2*np.pi))
        A = A_sim
        phi = phi_sim
        k = k_sim
        plt.figure(figsize=(12, 8))
        plt.plot(t, V(t), 'b-', label='Input')
        plt.plot(t, V_out(t), 'r-', label='Output')
        plt.plot(time_sim, sine_sim, 'k-', label='Data')
        plt.axhline(A*(1/np.sqrt(2)), label='Cutoff')

    if sys.argv[1] == 'sine_hard':
        w = f_out*2*np.pi
        t = np.linspace(0, 5/(w/(2*np.pi)), 5000)
        tmax = 3/(w/(2*np.pi))
        A = A_out
        phi = phi_out
        k = k_out
        plt.figure(figsize=(12, 8))
        plt.plot(t, V(t), 'b-', label='Input')
        plt.plot(t, V_out(t), 'r-', label='Output')
        plt.plot(time_out, sine_out, 'k-', label='Data')
        plt.axhline(A*(1/np.sqrt(2)), label='Cutoff')

        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Voltage [V]')
        plt.title('Angular frequency $\omega = {:.2f}$'.format(w))

        plt.savefig('sine.pdf')
        plt.show()
