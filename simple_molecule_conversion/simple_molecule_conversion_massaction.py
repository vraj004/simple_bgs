# This is a sample Python script to show how mass action law can be used to solve chemical reacion dynamics.

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
kappa = 1.0  # liters/mole-seconds
k1 = 2.0  # 1/seconds
k2 = 1.0  # 1/seconds

# Initial conditions
t = 0.0  # seconds
concA = 1.0  # moles/liter
concB = 1.0  # moles/liter
time_step = 0.01  # seconds
t_initial = 0.0 # seconds
t_final = 10.0  # seconds

#define the model
def model(z, t):
    # unpack the state variables
    concA = z[0]
    concB = z[1]

    # define the ODEs
    dconcA_dt = -k1 * concA+ k2 * concB
    dconcB_dt = -k2 * concB + k1 * concA

    # return the ODEs
    return [dconcA_dt, dconcB_dt]

# define the initial condition vector
z0 = [concA, concB]
# define the time vector
t = np.arange(t_initial, t_final, time_step)
# solve the ODEs
z = odeint(model, z0, t)
# plot the results
mtplot.plot(t, z[:, 0], 'b-', label='concA')
mtplot.xlim([0, 10])
mtplot.plot(t, z[:, 1], 'r--', label='concB')
mtplot.ylabel('Concentration (M)')
mtplot.xlabel('Time (sec)')
mtplot.title('Mass Action Law')
mtplot.legend(loc='best')
mtplot.show()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
