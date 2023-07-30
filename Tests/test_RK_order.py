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

def P62sin_exact(t, Psize):
    sol = np.array([])
    dx = (1.0-0.0)/(Psize+1);
    for i in range(0,Psize):
        x = dx * (i+1);
        sol = np.append(sol, 10*(1 - x)*x*(1+np.sin(t)) + 2)
    return sol


def order_plot(directory, num_runs, steps, exact_sol, interval, order):
    tau = np.array([])
    errors = np.array([])
    for i in range(1,num_runs):
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
    plt.loglog(tau, errors, '-*')

# Explicit Euler method
#order_plot(10, vector1_exact(1), 1-0, 1)
#order_plot("./runs/Lambert1/", 15, 10, Lambert_exact(1), 1-0, 1)
#order_plot("./runs/Lambert2/", 15, 10, Lambert_exact(1), 1-0, 4)
#order_plot("./runs/Lambert2Exp/", 15, 10, Lambert_exact(1), 1-0, 1)
order_plot("./runs/Problem_6_2_sin/", 9, 10, P62sin_exact(1,199), 1-0, 3)
#print(P62sin_exact(0,19))
plt.show()
