import matplotlib.pyplot as mtplot
import numpy as np
from scipy.integrate import odeint

# Model parameters
R = 8.314  # Joules/mole-Kelvin
T = 298  # Kelvin
V = 1.0  # liters
kappa = 1.0  # liters/mole-seconds
muA_standard = 1.0
muB_standard = 1.0   # moles/liter
muC_standard = 0.00001   # moles/liter
muD_standard = 0.00001   # moles/liter
# I changed these from 1 to 1719.0 and 1 to 0.00001 while examining k1 and k2 in the output.
#I basically did a manual bisection method to find the value of muA_standard that would give me a k1 of 2.0
# These numbers are more reasonable based on the output for k1 and k2 I get.
# I aimed to match k1 and k2 to the mass action model in massaction.py.

# Initial conditions
qA = 1.0  # moles
qB = 1.0  # moles
qC = 0.0  # moles
qD = 0.0  # moles
time_step = 0.01  # seconds
t_initial = 0.0 # seconds
t_final = 10.0  # seconds


def convert_qtoconc(q, V):
    return q / V
def mu(mu_standard, R, T, concq):
    return mu_standard + R * T * np.log(concq)


# define the model
def model(z, t):
    # unpack the state variables
    qA = z[0]
    qB = z[1]
    qC = z[2]
    qD = z[3]
    A = convert_qtoconc(qA, V)
    B = convert_qtoconc(qB, V)
    C = convert_qtoconc(qC, V)
    D = convert_qtoconc(qD, V)
    muA = mu(muA_standard, R, T, A)
    muB = mu(muB_standard, R, T, B)
    muC = mu(muC_standard, R, T, C)
    muD = mu(muD_standard, R, T, D)

    mu3 = muA + muB
    mu4 = muC + muD
    v3 = kappa*(np.exp((mu3/(R*T)))-np.exp((mu4/(R*T))))
    v1 = v3
    v2 = v3
    v5 = v3
    v6 = v3


    # define the ODEs
    dqA_dt = -v1
    dqB_dt = -v2
    dqC_dt = v5
    dqD_dt = v6

    # return the ODEs
    return [dqA_dt, dqB_dt, dqC_dt,dqD_dt]
#evalulate mass action model reaction rate equivalents
k1 = kappa*(np.exp((muA_standard+muB_standard)/(R*T)))/V #+np.exp(muB_standard/(R*T)))/V
k2 = kappa*(np.exp((muC_standard+muD_standard)/(R*T)))/V #+np.exp(muD_standard/(R*T)))/V
print('The equivalent k1 and k2 are: ', k1, k2)

# define the initial condition vector
z0 = [qA, qB, qC, qD]
# define the time vector
t = np.arange(t_initial, t_final, time_step)
# solve the ODEs
z = odeint(model, z0, t)
# plot the results
A_t = convert_qtoconc(z[:, 0], V)
B_t = convert_qtoconc(z[:, 1], V)
C_t = convert_qtoconc(z[:, 2], V)
D_t = convert_qtoconc(z[:, 3], V)
mtplot.plot(t, A_t, 'b-', label='[A]')
mtplot.plot(t, B_t, 'r--', label='[B]')
mtplot.plot(t, C_t, 'g:', label='[C]')
mtplot.plot(t, D_t, 'k-.', label='[D]')
mtplot.xlim([0, 10])
mtplot.ylabel('Concentration (M)')
mtplot.xlabel('Time (sec)')
mtplot.title('Bond graph derived')
mtplot.legend(loc='best')
mtplot.show()
