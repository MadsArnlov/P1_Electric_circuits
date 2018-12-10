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

"Output sine wave across capacitor with phase shift and amplitude"
data_out = np.genfromtxt("outputWave_3_5v30phase4755ohm96_54nF500f.csv",
                         delimiter=",")

f_out = 500
A_out = 3.5
phi_out = 30
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
t = np.linspace(0, abs(time[0]), 5000)

"""
The frequencies are distributed evenly across the logarithmic values,
using the start and end of the collected data. The data has 200 points, for
which reason the model also have 200 points.
"""
omega = np.logspace(np.log10(angular_frequency_C[0]),
                    np.log10(angular_frequency_C[-1]), num=200)


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


def V(t, w, phi, A):
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
    return A*np.sin(w*t + phi*(np.pi/180))


def V_out(t, w, phi, A, k):
    """
    Returns a sine wave as the output voltage.

    A = amplitude of input voltage

    Parameters
    ----------
    t: float
        The time for the sine wave
    """
    return A*H_lp(w)[0]*np.sin(w*t + (phi + H_lp(w)[2])*(np.pi/180)) + k*A/tau*np.exp(-t/tau)


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
    magnitude = 20*np.log10(modulus)
    argument = np.arctan(-omega*tau)*180/np.pi
    return modulus, magnitude, argument


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
    modulus = np.sqrt(omega**2 * tau**2)/np.sqrt(1 + omega**2 * tau**2)
    magnitude = 20*np.log10(modulus)
    argument = np.arctan(1/(omega*tau))*180/np.pi
    return modulus, magnitude, argument


# =============================================================================
# Plot of data
# =============================================================================


def data():
    plt.figure(figsize=(12, 8))
    plt.plot(time + abs(time[0]), cap, 'b,', label='Data')
    plt.plot(time + abs(time[0]), sqwave_data, 'k--', label='Square wave')
    plt.xlabel('$t$ [s]')
    plt.ylabel('$V_C$ [V]')
    plt.legend()
    plt.savefig('data.pdf')
    plt.show()

# =============================================================================
# Plot of mathematical model
# =============================================================================


def model():
    plt.figure(figsize=(12, 8))
    plt.plot(t, V_C(t), 'tab:orange')
    plt.plot(abs(time[0]) + t, V_C2(t), 'tab:orange',
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


def data_vs_model():
    plt.figure(figsize=(12, 8))
    plt.plot(t, V_C(t), 'tab:orange', label='Modelled voltage')
    plt.plot(abs(time[0]) + t, V_C2(t), 'tab:orange')
    plt.plot(time + abs(time[0]), cap, 'b,', label='Data points')
    plt.plot(tau, V_C(tau), 'kx', label='$V_C(tau)$')
    plt.plot(time + abs(time[0]), sqwave_data, 'k--', label='Square wave')
    plt.xlabel('$t$ [s]')
    plt.ylabel('$V_C$ [V]')
    plt.legend()

    plt.savefig('data_vs_model.pdf')
    plt.show()

# =============================================================================
# Plot of deviation between data and mathematical model
# =============================================================================


def deviation():
    plt.figure(figsize=(12, 8))
    plt.plot(t, cap[:5000] - V_C(t), 'k,',
             abs(time[0]) + t, cap[5000:] - V_C2(t), 'k,')
    plt.xlabel('t [s]')
    plt.ylabel('$V_{data} - V_C$ [V]')

    plt.savefig('deviation.pdf')
    plt.show()

# =============================================================================
# Relative percentage difference
# =============================================================================


def relative():
    discharging = ((cap[:5000] - V_C(t))*100)/abs(V_C(t))
    charging = ((cap[5000:] - V_C2(t))*100)/abs(V_C2(t))
    discharging1 = ((cap[:5000] - V_C(t))*100)/abs(V_C(t))
    charging1 = ((cap[5000:] - V_C2(t))*100)/abs(V_C2(t))
    for value in range(len(discharging)):
        if discharging[value] >= 100:
            discharging[value] = None
        elif discharging[value] <= -100:
            discharging[value] = None
    for value in range(len(charging)):
        if charging[value] <= -100:
            charging[value] = None
        elif charging[value] >= 100:
            charging[value] = None
    average_discharging = sum(abs(discharging)/len(discharging))
    average_charging = sum(abs(charging)/len(charging))
    print("The percentage difference of", "\n",
          "Discharging: {:.5f}%".format(average_discharging), "\n",
          "Charging:    {:.5f}%".format(average_charging))
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.plot(t, discharging1, 'k,',
             abs(time[0]) + t, charging1, 'k,')

    plt.subplot(2, 1, 2)
    plt.plot(t, discharging, 'k,',
             abs(time[0]) + t, charging, 'k,')
    plt.ylabel('% Difference')
    plt.xlabel('t [s]')

    plt.savefig("relative_percentage_difference.pdf")
    plt.show()

# =============================================================================
# Bodeplot of RC LP filter
# =============================================================================


def RCLP():
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, H_lp(omega)[1], 'tab:orange', label='LP Transfer function')
    plt.plot(angular_frequency_C, magnitude_C, 'b-', label='Data')
    plt.semilogx(omega_c, H_lp(omega_c)[1], 'kx', label='Gain at $\omega=\omega_c$')
    plt.legend()
    plt.grid(True)
    plt.ylabel('Magnitude $G(j\omega)$ [dB]')

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, H_lp(omega)[2], 'tab:orange', label='LP Phase shift')
    plt.plot(angular_frequency_C, phase_C, 'b-', label='Data')
    plt.semilogx(omega_c, H_lp(omega_c)[2], 'kx', label='Phase at $\omega=\omega_c$')
    plt.yticks(np.arange(0, -105, step=-15))
    plt.legend()
    plt.grid(True)
    plt.xlabel('Angular frequency $\omega$ [Hz]')
    plt.ylabel('Phase $\u03B8$ [degrees]')

    plt.savefig('data_bodeplots_rc_lp.pdf')
    plt.show()

    percent_difference_mag = (magnitude_C - H_lp(omega)[1])*100/H_lp(omega)[1]
    percent_difference_phase = (phase_C - H_lp(omega)[2])*100/H_lp(omega)[2]
    sum_percent_mag = sum(abs(percent_difference_mag))/len(magnitude_C)
    sum_percent_phase = sum(abs(percent_difference_phase))/len(magnitude_C)
    print("The percentage difference of", "\n",
          "Magnitude: {:.5f}%".format(sum_percent_mag), "\n",
          "Phase:     {:.5f}%".format(sum_percent_phase))

# =============================================================================
# Bodeplot of RC HP filter
# =============================================================================


def RCHP():
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, H_hp(omega)[1], 'tab:orange', label='HP Transfer function')
    plt.plot(angular_frequency_R, magnitude_R, 'b-', label='Data')
    plt.semilogx(omega_c, H_hp(omega_c)[1], 'kx', label='Gain at $\omega_c$')
    plt.legend()
    plt.grid(True)
    plt.ylabel('Magnitude $G(j\omega)$ [dB]')

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, H_hp(omega)[2], 'tab:orange', label='HP Phase shift')
    plt.plot(angular_frequency_R, phase_R, 'b-', label='Data')
    plt.semilogx(omega_c, H_hp(omega_c)[2], 'kx', label='Phase at $\omega_c$')
    plt.yticks(np.arange(0, 105, step=15))
    plt.legend()
    plt.grid(True)
    plt.xlabel('Angular frequency $\omega$ [Hz]')
    plt.ylabel('Phase $\u03B8$ [degree]')

    plt.savefig('data_bodeplots_rc_hp.pdf')
    plt.show()

    percent_difference_mag = (magnitude_R - H_hp(omega)[1])*100/H_hp(omega)[1]
    percent_difference_phase = (phase_R - H_hp(omega)[2])*100/H_hp(omega)[2]
    sum_percent_mag = sum(abs(percent_difference_mag))/len(magnitude_R)
    sum_percent_phase = sum(abs(percent_difference_phase))/len(magnitude_R)
    print("The percentage difference of", "\n",
          "Magnitude: {:.5f}%".format(sum_percent_mag), "\n",
          "Phase:     {:.5f}%".format(sum_percent_phase))

# =============================================================================
# Simulation of output sine wave
# =============================================================================


def sine():
    w = f_out*2*np.pi
    t = np.linspace(time_out[0], time_out[-1], 8192)
    A = A_out
    phi = phi_out
    k = k_out
    plt.figure(figsize=(12, 8))
    plt.plot(t, V(t, w, phi, A), 'k-', label='$V(t) = {}\cdot\sin(\omega t + {:.0f}\N{DEGREE SIGN})$'.format(A, phi))
    plt.plot(t, V_out(t, w, phi, A, k), 'tab:orange', label='$V_C(t) = {:.1f}\cdot\sin(\omega t {:.2f}\N{DEGREE SIGN})$'.format(A*H_lp(w)[0], (phi + H_lp(w)[2])))
    plt.axhline(A*(1/np.sqrt(2)), label='A of $V_C(t)$ at $\omega_c$')
    plt.plot(time_out, sine_out, 'b--', label='Data')
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Voltage [V]')
    plt.title('Angular frequency $\omega = {:.2f}$ and phase $\phi = {:.0f}\N{DEGREE SIGN}$'.format(w, phi))

    plt.savefig('sine.pdf')
    plt.show()

    percent_difference_sine = (sine_out - V_out(t, w, phi, A, k))*100/V_out(t, w, phi, A, k)
    sum_percent_sine = sum(abs(percent_difference_sine))/len(sine_out)
    print("The percentage difference is:", "\n",
          "{:.5f}%".format(sum_percent_sine))


if len(sys.argv) >= 2:
    if sys.argv[1].lower() == 'sine':
        sine()
    elif sys.argv[1].upper() == 'RCHP':
        RCHP()
    elif sys.argv[1].lower() == 'relative':
        relative()
    elif sys.argv[1].lower() == 'deviation':
        deviation()
    elif sys.argv[1].upper() == 'RCLP':
        RCLP()
    elif sys.argv[1].lower() == 'data_vs_model' or sys.argv[1] == 'model_vs_data':
        data_vs_model()
    elif sys.argv[1].lower() == 'model':
        model()
    elif sys.argv[1].lower() == 'data':
        data()
