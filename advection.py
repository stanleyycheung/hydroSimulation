"""
Upstream advection code - from chapter 3 in http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/index.shtml
Section 3.10
"""

import numpy as np
import matplotlib.pyplot as plt

nx = 100
x, q, qnew = np.empty(nx), np.empty(nx), np.empty(nx)
dx = 1
for i in range(0, nx):
    x[i] = i
    if i < 30:
        q[i] = 1
    else:
        q[i] = 0
u = 1
dt = 2e-1
tend = 30

qleft = q[0]
time = 0

while time < tend:
    if time + dt < tend:
        dtused = dt
    else:
        dtused = tend-time
    for i in range(1, nx):
        qnew[i] = q[i] - u * (q[i] - q[i-1]) * dtused/dx
        q[i] = qnew[i]
    q[0] = qleft
    time = time + dtused

plt.plot(x, q, '.-')
plt.xlabel('x')
plt.ylabel('q')
plt.grid()
plt.show()
