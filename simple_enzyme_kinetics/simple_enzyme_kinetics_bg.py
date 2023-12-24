import matplotlib.pyplot as mtplot
import numpy as np
from scipy.integrate import odeint

# Model parameters
R = 8.314  # Joules/mole-Kelvin
T = 298  # Kelvin
V = 1.0  # liters
kappa1 = 1.0  # liters/mole-seconds
kappa2 = 1.0  # liters/mole-seconds
muA_standard = 10.0
muB_standard = 10.0   # moles/liter
muC_standard = -50.0   # moles/liter
muD_standard = -50.0   # moles/liter
muE_standard = -100.0   # moles/liter


# Initial conditions
qA = 1.0  # moles
qB = 1.0  # moles
qC = 1.0  # moles
qD = 1.0  # moles
qE = 1.0  # moles
time_step = 0.0001  # seconds
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
    qE = z[4]
    A = convert_qtoconc(qA, V)
    B = convert_qtoconc(qB, V)
    C = convert_qtoconc(qC, V)
    D = convert_qtoconc(qD, V)
    E = convert_qtoconc(qE, V)

    muA = mu(muA_standard, R, T, A)
    muB = mu(muB_standard, R, T, B)
    muC = mu(muC_standard, R, T, C)
    muD = mu(muD_standard, R, T, D)
    muE = mu(muE_standard, R, T, E)

    mu4 = muA + muB + muE
    mu5 = muC
    mu7 = muC
    mu8 = muD+muE
    v4 = kappa1*(np.exp((mu4/(R*T)))-np.exp((mu5/(R*T))))
    v7 = kappa2*(np.exp((mu7/(R*T)))-np.exp((mu8/(R*T))))
    v1 = v4
    v2 = v4
    v3 = v4
    v5 = v4
    v6 = v4

    v10 = v7
    v8 = v7
    v9 = v7

    # define the ODEs
    dqA_dt = -v2
    dqB_dt = -v1
    dqC_dt = v6-v7
    dqD_dt = -v9
    dqE_dt = v10-v3
    # return the ODEs
    return [dqA_dt, dqB_dt, dqC_dt,dqD_dt,dqE_dt]
#evalulate mass action model reaction rate equivalents
#k1 = kappa*(np.exp((muA_standard+muB_standard)/(R*T)))/V #+np.exp(muB_standard/(R*T)))/V
#k2 = kappa*(np.exp((muC_standard+muD_standard)/(R*T)))/V #+np.exp(muD_standard/(R*T)))/V
#print('The equivalent k1 and k2 are: ', k1, k2)

# define the initial condition vector
z0 = [qA, qB, qC, qD,qE]
# define the time vector
t = np.arange(t_initial, t_final, time_step)
# solve the ODEs
z = odeint(model, z0, t)
# plot the results
A_t = convert_qtoconc(z[:, 0], V)
B_t = convert_qtoconc(z[:, 1], V)
C_t = convert_qtoconc(z[:, 2], V)
D_t = convert_qtoconc(z[:, 3], V)
E_t = convert_qtoconc(z[:, 4], V)

mtplot.plot(t, A_t, 'b-', label='[A]')
mtplot.plot(t, B_t, 'r--', label='[B]')
mtplot.plot(t, C_t, 'g:', label='[C]')
mtplot.plot(t, D_t, 'k-.', label='[D]')
mtplot.plot(t, E_t, 'm-.', label='[E]')

mtplot.xlim([0, 10])
mtplot.ylabel('Concentration (M)')
mtplot.xlabel('Time (sec)')
mtplot.title('Bond graph derived')
mtplot.legend(loc='best')
mtplot.show()
