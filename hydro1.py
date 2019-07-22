"""
Basic 1-D hydro solver
Section 5.2
"""

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D


plt.rcParams["figure.figsize"] = [12, 8]


def hydroiso_cen(x, xi, rho, rhou, e, gamma, dt):
    nx = len(x)

    # compute velocities at the cell interfaces; eqn 5.10
    ui = np.zeros(nx+1)
    for ix in range(1, nx):
        ui[ix] = 0.5 * (rhou[ix]/rho[ix] + rhou[ix-1]/rho[ix-1])

    # compute flux for rho; eqn 5.11
    fluxrho = np.zeros(nx+1)
    for ix in range(1, nx):
        if ui[ix] > 0:
            fluxrho[ix] = rho[ix-1] * ui[ix]
        else:
            fluxrho[ix] = rho[ix] * ui[ix]

    # update density; eqn 5.12
    for ix in range(0, nx):
        rho[ix] = rho[ix] - (dt/(xi[ix+1]-xi[ix]))*(fluxrho[ix+1]-fluxrho[ix])

    # compute flux for rho u; repeat all for rho*u
    fluxrhou = np.zeros(nx+1)
    for ix in range(1, nx):
        if ui[ix] > 0:
            fluxrhou[ix] = rhou[ix-1]**2/rho[ix-1]
        else:
            fluxrhou[ix] = rhou[ix]**2/rho[ix]

    # update the momentum
    for ix in range(0, nx):
        rhou[ix] = rhou[ix] - (dt/(xi[ix+1]-xi[ix]))*(fluxrhou[ix+1]-fluxrhou[ix])

    # compute pressure
    p = (gamma-1.)*rho*e

    # now add the pressure force, for all cells except the ones near the boundary
    for ix in range(1, nx-1):
        rhou[ix] = rhou[ix] - dt*(p[ix+1]-p[ix-1])/(x[ix+1]-x[ix-1])

    # now do the boundary cells, assuming mirror symmetry in the boundaries
    rhou[0] = rhou[0] - 0.5*dt*(p[1]-p[0])/(x[1]-x[0])
    rhou[nx-1] = rhou[nx-1] - 0.5*dt*(p[nx-1]-p[nx-2])/(x[nx-1]-x[nx-2])


# define a test problem
nx = 100
nt = 1000
x0 = 0
x1 = 100
xmid = 0.5*(x0+x1)
dt = 0.25
cfl = 0.5
x = x0 + (x1-x0)*(np.arange(nx)/(nx-1.))
gamma = 7./5.
rho = np.zeros((nx, nt+1))
rhou = np.zeros((nx, nt+1))
e = np.ones(nx)
time = np.zeros(nt+1)
dg = 0.1*(x1-x0)
rho[:, 0] = 1 + 3*np.exp(-(x-xmid)**2/dg**2)

# additional arrays needed
xi = np.zeros(nx+1)
xi[1:nx] = 0.5 * (x[1:nx] + x[0:nx-1])
xi[0] = 2*xi[1] - xi[2]
xi[nx] = 2*xi[nx-1] - xi[nx-2]
dx = (xi[1:nx+1] - xi[0:nx])

# main loop
counter = 1
for it in range(1, nt+1):
    qrho = rho[:, it-1]
    qrhou = rhou[:, it-1]
    cs = np.sqrt(gamma*(gamma-1)*e)
    dum = dx/(cs+abs(qrhou/qrho))
    dt = cfl*min(dum)
    time[it] = time[it-1]+dt
    # print("Time step: {}, Time = {}, Dt = {}".format(it, time[it], dt))
    hydroiso_cen(x, xi, qrho, qrhou, e, gamma, dt)
    rho[:, it] = qrho
    rhou[:, it] = qrhou
    np.savetxt('output/output{:05d}.dat'.format(counter),
               np.c_[x, qrho, qrho], header='{}'.format(time[it-1]))
    counter += 1


'''
# 3D-plot:
# print(rhou.shape)
# print(time)

time_plot, x_plot = np.meshgrid(time, x)

fig = plt.figure()
ax = plt.axes(projection='3d')

print(x_plot.shape)
print(time_plot.shape)
print(rho.shape)

print(rho)

ax.plot_surface(x_plot, time_plot, rho, cmap='viridis')
plt.show()
'''
