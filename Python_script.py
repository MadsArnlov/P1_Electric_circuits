# =============================================================================
# Modelling of RC-circuit
# =============================================================================
"The program uses the modules matplotlib and numpy"
import matplotlib.pyplot as plt
import numpy as np
import sys

# =============================================================================
# How to use the script
# =============================================================================
"""
The following commands can be used to generate plots:
    data() - plots the data from the experiment of the time domain.
    model() - plots the modelled voltage of the capacitor in the time domain.
    data_vs_model() - plots both the data and model of the capacitor voltage.
    deviation() - plots the deviation between the data and model.
    relative() - plots the relative percentage difference in the time domain.
    RCLP() - plots the Bode plot of the low pass filter.
    RCHP() - plots the Bode plot of the high pass filter.
    sine() - plots the output sine wave across the capacitor.
"""
# =============================================================================
# Data
# =============================================================================
"The data is exported from 'Waveform'."
"Voltage across capacitor in time domain"
data_time = np.genfromtxt("charging_discharging.csv", delimiter=",")

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
    discharging = 0 + (V0 - 0)*np.exp(-t/tau) - V0/2
    charging = V0 + (0 - V0)*np.exp(-t/tau) - V0/2
    return discharging, charging


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
    Returns the voltage of the input wave to time t.

    Parameters
    ----------
    t: float
        The time for the sine wave
    w: float
        The angular frequency
    phi: float
        The phase
    """
    return A*np.sin(w*t + phi*(np.pi/180))


def V_out(t, w, phi, A, k):
    """
    Returns the voltage of the output wave to time t.

    Parameters
    ----------
    t: float
        The time for the sine wave
    w: float
        The angular frequency
    phi: float
        The phase
    k: float
        The oscilliation constant
    """
    return A*H_lp(w)[0]*np.sin(w*t + (phi + H_lp(w)[2])*(np.pi/180)) + k*A/tau*np.exp(-t/tau)


def H_lp(omega):
    """
    Returns modulus, magnitude and the argument of the transfer function
    for a first order low pass filter.

    The formula is derived in chapter "RC Frequency Analysis", section
    "Response and Transfer Function of the Capacitor":
         H(s)  = 1/sqrt(omega*tau + 1)
        |H(s)| = 1/sqrt(omega^2*tau^2 + 1)
         G(s)  = 20*log10(|H(s)|)
     arg(H(s)) = arctan(-omega*tau)

    Parameters
    ----------
    omega: float
        The angular frequency

    Returns
    -------
    modulus: float
        The absolute value of the transfer function
    magnitude: float
        Logarithmic scale of the modulus
    argument: float
        The phase of the transfer function in degrees
    """
    modulus = 1/np.sqrt(omega**2*tau**2 + 1)
    magnitude = 20*np.log10(modulus)
    argument = np.arctan(-omega*tau)*180/np.pi
    return modulus, magnitude, argument


def H_hp(omega):
    """
    Returns modulus, magnitude and the argument of the transfer function
    for a first order high pass filter.

    The formula is derived in chapter "RC Frequency Analysis", section
    "Response and transfer function of the Resistor":
         H(s)  = omega*tau/sqrt(1 + omega*tau)
        |H(s)| = sqrt(omega^2*tau^2)/sqrt(omega^2*tau^2 + 1)
         G(s)  = 20*log10(|H(s)|)
     arg(H(s)) = arctan(1/omega*tau)

    Parameters
    ----------
    omega: float
        The angular frequency

    Returns
    -------
    modulus: float
        The absolute value of the transfer function
    magnitude: float
        Logarithmic scale of the modulus
    argument: float
        The phase of the transfer function in degrees
    """
    modulus = np.sqrt(omega**2*tau**2)/np.sqrt(omega**2*tau**2 + 1)
    magnitude = 20*np.log10(modulus)
    argument = np.arctan(1/(omega*tau))*180/np.pi
    return modulus, magnitude, argument


# =============================================================================
# Plot of data
# =============================================================================
"The time is shifted, so that it starts from 0"


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
"To simulate discharging and charging, a piecewise function is used - V_C(t)."


def model():
    plt.figure(figsize=(12, 8))
    plt.plot(t, V_C(t)[0], 'tab:orange')
    plt.plot(abs(time[0]) + t, V_C(t)[1], 'tab:orange',
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
"The data and simulated model are plotted together with the square wave."


def data_vs_model():
    plt.figure(figsize=(12, 8))
    plt.plot(t, V_C(t)[0], 'tab:orange', label='Modelled voltage')
    plt.plot(abs(time[0]) + t, V_C(t)[1], 'tab:orange')
    plt.plot(time + abs(time[0]), cap, 'b,', label='Data points')
    plt.plot(tau, V_C(tau)[0], 'kx', label='$V_C(\u03C4)$')
    plt.plot(time + abs(time[0]), sqwave_data, 'k--', label='Square wave')
    plt.xlabel('$t$ [s]')
    plt.ylabel('$V_C$ [V]')
    plt.legend()

    plt.savefig('data_vs_model.pdf')
    plt.show()


# =============================================================================
# Plot of deviation between data and mathematical model
# =============================================================================
"The deviation between data and model are calculated and plottet."


def deviation():
    plt.figure(figsize=(12, 8))
    plt.plot(t, cap[:5000] - V_C(t)[0], 'k,',
             abs(time[0]) + t, cap[5000:] - V_C(t)[1], 'k,')
    plt.xlabel('t [s]')
    plt.ylabel('$V_{data} - V_C$ [V]')

    plt.savefig('deviation.pdf')
    plt.show()


# =============================================================================
# Relative percentage difference
# =============================================================================
"""
First the percentage difference for discharging and charging are calculated.
Afterwards all indices where the percentage difference is greater than 100%
are set to "None".

The average percentage difference is calculated and printed where all
indices are included.

The percentage difference with and without percentages greater than
100% are plotted.
"""


def relative():
    discharging = ((cap[:5000] - V_C(t)[0])*100)/abs(V_C(t)[0])
    charging = ((cap[5000:] - V_C(t)[1])*100)/abs(V_C(t)[1])
    discharging1 = ((cap[:5000] - V_C(t)[0])*100)/abs(V_C(t)[0])
    charging1 = ((cap[5000:] - V_C(t)[1])*100)/abs(V_C(t)[1])

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

    average_discharging = sum(abs(discharging1)/len(discharging1))
    average_charging = sum(abs(charging1)/len(charging1))

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
"""
Bodeplot of low pass filter consisting of the magnitude and phase.

After the plot the average percentage difference is also calculated.
"""

def RCLP():
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, H_lp(omega)[1], 'tab:orange',
                 label='LP transfer function')
    plt.plot(angular_frequency_C, magnitude_C, 'b-', label='Data')
    plt.semilogx(omega_c, H_lp(omega_c)[1], 'kx',
                 label='$G(j\omega_c)$')
    plt.legend()
    plt.grid(True)
    plt.ylabel('Magnitude $G(j\omega)$ [dB]')

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, H_lp(omega)[2], 'tab:orange', label='LP phase shift')
    plt.plot(angular_frequency_C, phase_C, 'b-', label='Data')
    plt.semilogx(omega_c, H_lp(omega_c)[2], 'kx',
                 label='Phase at $\omega=\omega_c$')
    plt.yticks(np.arange(0, -105, step=-15))
    plt.legend()
    plt.grid(True)
    plt.xlabel('Angular frequency $\omega$ [rad/s]')
    plt.ylabel('Phase $\u03B8$ [degrees]')

    plt.savefig('data_bodeplots_rc_lp.pdf')
    plt.show()

    amplitude_percent_mag = ((10**(magnitude_C/20) -
                              H_lp(angular_frequency_C)[0]))*100/(
                              H_lp(angular_frequency_C)[0])
    sum_amp = sum(abs(amplitude_percent_mag))/len(magnitude_C)

    percent_difference_phase = (phase_C - H_lp(omega)[2])*100/H_lp(omega)[2]
    sum_percent_phase = sum(abs(percent_difference_phase))/len(magnitude_C)

    print("The percentage difference of", "\n",
          "Amplitude: {:.5f}%".format(sum_amp), "\n",
          "Phase:     {:.5f}%".format(sum_percent_phase))


# =============================================================================
# Bodeplot of RC HP filter
# =============================================================================
"""
Bodeplot of high pass filter consisting of the magnitude and phase.

After the plot the average percentage difference is also calculated.
"""


def RCHP():
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, H_hp(omega)[1], 'tab:orange',
                 label='HP transfer function')
    plt.plot(angular_frequency_R, magnitude_R, 'b-', label='Data')
    plt.semilogx(omega_c, H_hp(omega_c)[1], 'kx', label='$G(j\omega_c)$')
    plt.legend()
    plt.grid(True)
    plt.ylabel('Magnitude $G(j\omega)$ [dB]')

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, H_hp(omega)[2], 'tab:orange', label='HP phase shift')
    plt.plot(angular_frequency_R, phase_R, 'b-', label='Data')
    plt.semilogx(omega_c, H_hp(omega_c)[2], 'kx', label='Phase at $\omega = \omega_c$')
    plt.yticks(np.arange(0, 105, step=15))
    plt.legend()
    plt.grid(True)
    plt.xlabel('Angular frequency $\omega$ [rad/s]')
    plt.ylabel('Phase $\u03B8$ [degree]')

    plt.savefig('data_bodeplots_rc_hp.pdf')
    plt.show()

    amplitude_percent_mag = ((10**(magnitude_R/20) -
                              H_hp(angular_frequency_R)[0]))*100/(
                              H_hp(angular_frequency_R)[0])
    sum_amp = sum(abs(amplitude_percent_mag))/len(magnitude_R)

    percent_difference_phase = (phase_R - H_hp(omega)[2])*100/H_hp(omega)[2]
    sum_percent_phase = sum(abs(percent_difference_phase))/len(magnitude_R)

    print("The percentage difference of", "\n",
          "Amplitude: {:.5f}%".format(sum_amp), "\n",
          "Phase:     {:.5f}%".format(sum_percent_phase))


# =============================================================================
# Simulation of output sine wave
# =============================================================================
"""
Plots the returning voltage across the capacitor for an input sine wave
with an amplitude of 3.5 and phase of 30 degrees.

The percentage difference is calculated.
"""


def sine():
    w = f_out*2*np.pi
    t = np.linspace(time_out[0], time_out[-1], 8192)
    A = A_out
    phi = phi_out
    k = k_out
    plt.figure(figsize=(12, 8))
    plt.plot(t, V(t, w, phi, A), 'k-',
             label='$V(t) = {}\cdot\sin(\omega t + {:.0f}\N{DEGREE SIGN})$'
             .format(A, phi))
    plt.plot(t, V_out(t, w, phi, A, k), 'tab:orange',
             label='$V_C(t) = {:.1f}\cdot\sin(\omega t {:.2f}\N{DEGREE SIGN})$'
             .format(A*H_lp(w)[0], (phi + H_lp(w)[2])))
    plt.axhline(A*(1/np.sqrt(2)), label='Amplitude of $V_C(t)$ at $\omega_c$')
    plt.plot(time_out, sine_out, 'b--', label='Data')
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Voltage [V]')
    plt.title("Angular frequency $\omega = {:.2f}$ and phase $\phi = {:.0f}\N{DEGREE SIGN}$"
              .format(w, phi))

    plt.savefig('sine.pdf')
    plt.show()

    percent_difference_sine = (sine_out - V_out(t, w, phi, A, k))*100/V_out(
            t, w, phi, A, k)
    sum_percent_sine = sum(abs(percent_difference_sine))/len(sine_out)
    print("The percentage difference is:", "\n",
          "{:.5f}%".format(sum_percent_sine))


"The program can be run from the command line."

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
