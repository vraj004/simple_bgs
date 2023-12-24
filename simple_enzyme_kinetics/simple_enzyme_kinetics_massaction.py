# This is a sample Python script to show how mass action law can be used to solve chemical reacion dynamics.
#The reaction is: A+B+E->C->D+E
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
# import libraries
import matplotlib.pyplot as mtplot
import numpy as np
from scipy.integrate import odeint

# Model parameters

R = 8.314  # Joules/mole-Kelvin
T = 300  # Kelvin
V = 1.0  # liters

k1 = 1.0  # 1/seconds #forward reaction A+B+E->C
k2 = 1.0  # 1/seconds #forward reaction C->D+E
k3 = 1.0  # 1/seconds #reverse reaction D+E->C
k4 = 1.0  # 1/seconds #reverse reaction C->A+B+E
# Initial conditions
t = 0.0  # seconds
concA = 1.0  # moles/liter
concB = 1.0  # moles/liter
concC = 0.0  # moles/liter
concD = 0.0  # moles/liter
concE = 1.0  # moles/liter
time_step = 0.00001  # seconds
t_initial = 0.0 # seconds
t_final = 10.0  # seconds

#define the model
def model(z, t):
    # unpack the state variables
    concA = z[0]
    concB = z[1]
    concC = z[2]
    concD = z[3]
    concE = z[4]

    # define the ODEs
    dconcA_dt = -k1 * (concA * concB*concE) + k4 * (concC)
    dconcB_dt = -k1 * (concA * concB*concE) + k4 * (concC)
    dconcC_dt = k1 * (concA * concB*concE) - k2 * (concD*concE) - k4 * (concC)
    dconcD_dt = k2 * (concC) - k3 * (concE*concD)
    dconcE_dt = -k1 * (concA * concB*concE) + k2 * (concD*concE) - k3 * (concE*concD)

    # return the ODEs
    return [dconcA_dt, dconcB_dt, dconcC_dt,dconcD_dt,dconcE_dt]

# define the initial condition vector
z0 = [concA, concB, concC, concD, concE]
# define the time vector
t = np.arange(t_initial, t_final, time_step)
# solve the ODEs
z = odeint(model, z0, t)
# plot the results
mtplot.plot(t, z[:, 0], 'b-', label='concA')
mtplot.plot(t, z[:, 1], 'r--', label='concB')
mtplot.plot(t, z[:, 2], 'g:', label='concC')
mtplot.plot(t, z[:, 3], 'k-.', label='concD')
mtplot.plot(t, z[:, 4], 'm-.', label='concE')
mtplot.ylabel('Concentration (M)')
mtplot.xlim([0, 10])
mtplot.xlabel('Time (sec)')
mtplot.title('Mass Action Law')
mtplot.legend(loc='best')
mtplot.show()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
