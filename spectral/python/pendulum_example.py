# The second order differential equation for the angle theta of a pendulum 
# acted on by gravity with friction can be written:

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



# solver
def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt


# constants
b = 0.25
c = 5.0

# initial conditions
y0 = [np.pi - 0.1, 0.0]

# We generate a solution 101 evenly spaced samples in the interval 0 <= t <= 10. 
# So our array of times is:
t = np.linspace(0, 10, 101)

sol = odeint(pend, y0, t, args=(b, c))

plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
