"""
Upstream advection code - from chapter 4 in http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/index.shtml
Section 4.6
"""

import numpy as np
import matplotlib.pyplot as plt

nx = 100
# 100 grid points + 2 ghost cells : x
# 103 cell walls: xi
x, xi = np.empty(nx+2), np.empty(nx+3)
q, qnew, f = np.empty(nx+2), np.empty(nx+2), np.empty(nx+3)

dx = 1

for i in range(0, nx+2):
    x[i] = (i-1)*dx
for i in range(1, nx+2):
    xi[i] = 0.5 * (x[i] + x[i-1])
xi[0] = x[0] - 0.5 * (x[1]-x[0])
xi[nx+2] = x[nx+1] + 0.5 * (x[nx+1] - x[nx])
# wall i is left of cell i

# initial conditions
for i in range(0, nx):
    if x[i] < 30:
        q[i] = 1
    else:
        q[i] = 0
u = np.full(nx+3, 1)

dt = 2e-01
tend = 70
time = 0

while time < tend:
    if time + dt < tend:
        dtused = dt
    else:
        dtused = tend - time
    q[0] = q[nx]
    q[nx+1] = q[1]

    for i in range(1, nx+2):
        if u[i] > 0:
            f[i] = q[i-1] * u[i]
        else:
            f[i] = q[i] * u[i]
    for i in range(1, nx+1):
        qnew[i] = q[i] - dtused * (f[i+1] - f[i]) / (xi[i+1] - xi[i])
        q[i] = qnew[i]
    time = time + dtused

plt.plot(x[1:nx], q[1:nx])
plt.show()
