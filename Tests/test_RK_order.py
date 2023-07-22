import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def vector1_exact(t):
    sol0 = (-8.0/5.0)*np.exp(-t)*(-1) - (4.0/5.0)*np.exp(4*t)*2
    sol1 = (-8.0/5.0)*np.exp(-t)*1 - (4.0/5.0)*np.exp(4*t)*3
    arr = np.array([sol0, sol1])
    return arr

def Lambert_exact(t):
    sol0 = 2*np.exp(-t) + np.sin(t)
    sol1 = 2*np.exp(-t) + np.cos(t)
    arr = np.array([sol0, sol1])
    return arr


def order_plot(directory, steps, exact_sol, interval, order):
    tau = np.array([])
    errors = np.array([])
    for i in range(1,15):
        file_name = directory+"vector_output_" + str(int(steps)) + ".csv"
        df = pd.read_csv(file_name, header=None)

        # Scale the error norm by dividing with the step-size to the
        # power of the order of the method. Should see correct order if
        # graph levels off.
        step_size = interval/steps
        l2_error = np.linalg.norm(exact_sol-np.array(df.iloc[-1,:].values[1:-1]), np.inf) * step_size**(-order)

        tau = np.append(tau, step_size)
        errors = np.append(errors, l2_error)

        steps = steps*2
    # TODO: Add asteriks to the datapoints
    plt.loglog(tau, errors, '-*')

# Explicit Euler method
#order_plot(10, vector1_exact(1), 1-0, 1)
#order_plot("./runs/Lambert1/", 10, Lambert_exact(1), 1-0, 1)
order_plot("./runs/Lambert2/", 10, Lambert_exact(1), 1-0, 1)
plt.show()
