"""
1-D Hydro simulation
Flux limiter options: donor cell/ superbee/ van Leer
BC options: mirror/ periodic/ mirror-free/ disc

Added spherical graviational field: -GM/R
Added energy equation
Added rotational velocity and force

Written by SC
"""

import numpy as np
import matplotlib.pyplot as plt


def advect(x, xi, q, ui, dt, fluxlim, nghost):
    nx = len(xi) - 1
    if len(x) != nx or len(xi) != nx + 1 or len(q) != nx or len(ui) != nx + 1 or nghost < 1:
        raise ValueError('Wrong input shape')

    # Determine the r_{i-1/2} for the flux limiter
    r = np.zeros(nx + 1)
    for i in range(2, nx - 1):
        dq = (q[i] - q[i - 1])
        if abs(dq) > 0:
            if ui[i] > 0:
                r[i] = (q[i - 1] - q[i - 2]) / dq
            else:
                r[i] = (q[i + 1] - q[i]) / dq

    #  Determine the flux limiter (many other flux limiters can be implemented here!)
    if fluxlim == 'donor-cell':
        phi = np.zeros(nx + 1)
    elif fluxlim == 'superbee':
        phi = np.zeros(nx + 1)
        for i in range(1, nx):
            a = min([1., 2. * r[i]])
            b = min([2., r[i]])
            phi[i] = max([0., a, b])
    elif fluxlim == 'van Leer':
        phi = (r + abs(r)) / (1 + abs(r))
    else:
        raise Exception("Choose a valid flux limiter")

    # Now construct the flux
    flux = np.zeros(nx + 1)
    for i in range(1, nx):
        if ui[i] > 0:
            flux[i] = ui[i] * q[i - 1]
        else:
            flux[i] = ui[i] * q[i]
        flux[i] = flux[i] + 0.5 * \
            abs(ui[i]) * (1 - abs(ui[i] * dt / (x[i] - x[i - 1]))) * \
            phi[i] * (q[i] - q[i - 1])

    # Update the cells, except the ghost cells
    for i in range(nghost, nx - nghost):
        q[i] = q[i] - dt * (flux[i + 1] - flux[i]) / (xi[i + 1] - xi[i])


def boundary(rho, rhou, rhov, bc):
    # Get the number of grid points including the ghost cells
    nx = len(rho)

    # If periodic, then install periodic BC, using two ghost cells on each side
    # (two are required for the non-lin flux limiter)
    if bc == 'periodic':
        rho[0] = rho[nx - 4]
        rho[1] = rho[nx - 3]
        rho[nx - 2] = rho[2]
        rho[nx - 1] = rho[3]
        rhou[0] = rhou[nx - 4]
        rhou[1] = rhou[nx - 3]
        rhou[nx - 2] = rhou[2]
        rhou[nx - 1] = rhou[3]

    # If mirror symmetry, then install mirror BC, using two ghost cells
    # on each side (two are required for the non-lin flux limiter)
    elif bc == 'mirror':
        rho[0] = rho[3]
        rho[1] = rho[2]
        rho[nx - 2] = rho[nx - 3]
        rho[nx - 1] = rho[nx - 4]
        rhou[0] = -rhou[3]
        rhou[1] = -rhou[2]
        rhou[nx - 2] = -rhou[nx - 3]
        rhou[nx - 1] = -rhou[nx - 4]
    elif bc == 'mirror/free':
        # Left BC
        rho[0] = rho[2]
        rho[1] = rho[3]
        rhou[0] = -rhou[2]
        rhou[1] = -rhou[3]
        rhov[0] = -rhov[2]
        rhov[1] = -rhov[3]
        # Right BC
        rho[nx - 2] = rho[nx - 4]
        rho[nx - 1] = rho[nx - 3]
        rhou[nx - 2] = rhou[nx - 4]
        rhou[nx - 1] = rhou[nx - 3]
        rhov[nx - 2] = rhov[nx - 4]
        rhov[nx - 1] = rhov[nx - 3]
    elif bc == 'disc':
        # Left BC
        rho[0] = N * mu  # equation 1.4
        rho[1] = N * mu
        rhou[0] = rhou[2]/rho[2] * N * mu
        rhou[1] = rhou[3]/rho[3] * N * mu
        rhov[0] = np.sqrt(G*M_star/(x[0]+starDist)) * rho[0]
        rhov[1] = np.sqrt(G*M_star/(x[1]+starDist)) * rho[1]
        # Right BC
        rho[nx - 2] = rho[nx - 4]
        rho[nx - 1] = rho[nx - 3]
        rhou[nx - 2] = rhou[nx - 4]
        rhou[nx - 1] = rhou[nx - 3]
        rhov[nx - 2] = rhov[nx - 4]
        rhov[nx - 1] = rhov[nx - 3]
    else:
        raise Exception("Choose valid BC")


def hydrostep(x, xi, rho, rhou, rhoe, rhov, gamma, dt, bc, fluxlim='donor-cell'):
    # Use 2 ghost cells on each side
    nghost = 2

    # Get the number of grid points including the ghost cells
    nx = len(x)
    # Impose boundary conditions
    boundary(rho, rhou, rhov, bc)

    # Compute the velocities at the cell interfaces
    ui = np.zeros(nx + 1)
    vi = np.zeros(nx + 1)
    for ix in range(1, nx):
        ui[ix] = 0.5 * (rhou[ix] / rho[ix] + rhou[ix - 1] / rho[ix - 1])
        vi[ix] = 0.5 * (rhov[ix] / rho[ix] + rhov[ix - 1] / rho[ix - 1])

    # Advect rho
    advect(x, xi, rho, ui, dt, fluxlim, nghost)

    # Advect rhou
    advect(x, xi, rhou, ui, dt, fluxlim, nghost)

    # Advect rhoe
    advect(x, xi, rhoe, ui, dt, fluxlim, nghost)

    # Advect rhov
    advect(x, xi, rhov, vi, dt, fluxlim, nghost)

    # Re-impose boundary conditions
    boundary(rho, rhou, rhov, bc)

    # Compute the pressure for each quantity
    u = rhou / rho
    v = rhov / rho
    e_tot = rhoe / rho
    e_kin = u**2 / 2
    e_th = e_tot - e_kin
    p = (gamma - 1.) * rho * e_th

    # Calculate graviational potiental
    V = G * M_star / (x + starDist)
    # Now add all external forces, for all cells except the ghost cells
    for ix in range(2, nx - 2):
        rhou[ix] = rhou[ix] - dt * \
            (p[ix + 1] - p[ix - 1] + rho[ix + 1] * V[ix + 1] -
             rho[ix - 1] * V[ix - 1]) / (2 * (x[ix + 1] - x[ix - 1]))
        rhoe[ix] = rhoe[ix] - dt * (p[ix + 1] * u[ix + 1] -
                                    p[ix - 1] * u[ix - 1] + rho[ix + 1] * V[ix + 1] * u[ix + 1] -
                                    rho[ix - 1] * V[ix - 1] * u[ix - 1]) / (2 * (x[ix + 1] - x[ix - 1]))
        # impliment rotational force
        rhou[ix] = rhou[ix] + dt * rho[ix]*v[ix]**2/x[ix]
    # Re-impose boundary conditions a last time (not strictly necessary)
    boundary(rho, rhou, rhov, bc)


def rho_initial(option):
    # Gaussian
    if option == 'gaussian':
        xmid = 0.5 * (x0 + x1)
        dg = 0.1 * (x1 - x0)
        result = 1 + 3 * np.exp(-(x - xmid)**2 / dg**2)
    # Sod shock tube
    elif option == 'sod':
        result = np.zeros(len(x))
        for i in range(0, len(x)):
            if i < int(len(x) / 2):
                result[i] = 1
            else:
                result[i] = 0.125
    # Atmosphere
    elif option == 'atmosphere':
        # Density of air is ~ 1kg/m^3
        result = np.full(len(x), 1)
    elif option == 'zero':
        result = np.zeros(nx)
    else:
        raise Exception("Choose a valid density")
    return result


def e_initial(option):
    # Flat
    if option == 'flat':
        result = np.ones(nx)
    # Sod shock tube
    elif option == 'sod':
        result = np.zeros(nx)
        for i in range(0, len(result)):
            if i < int(len(result) / 2):
                result[i] = 2.5
            else:
                result[i] = 2.
    # Atmosphere
    elif option == 'atmosphere':
        # typical pressure is 100kPa
        result = np.full(len(x), 2.5e5)
    elif option == 'zero':
        result = np.zeros(nx)
    else:
        raise Exception("Choose a valid initial e")
    return result


def v_initial(option):
    # Keplarian
    if option == 'kep':
        result = np.zeros(nx)
        result[0] = np.sqrt(G*M_star/(x[0]+starDist))
    elif option == 'step':
        result = np.zeros(nx)
        for i in range(0, len(result)):
            if i < int(len(result) / 2):
                result[i] = 10
            else:
                result[i] = 5
    elif option == 'zero':
        result = np.zeros(nx)
    else:
        raise Exception("Choose a valid initial v")
    return result


# physical constants - mass in solar masses, length in m
G = 6.67e-11
AU = 1.496e11
M_solar = 1.9891e30
# M_solar = M_star
M_star = 5.972e24
# starDist = 100 * AU
starDist = 6.371e6
m_h = 1.673e-27
mu = 1.3 * m_h  # following Facchini et al. 2016
N = 10**9.2
T = 20

# atmosphere test problem
nx = 2500
nt = 10000
x0 = 0
x1 = 1200*AU
dt = 0.5
cfl = 0.3
x = x0 + (x1 - x0) * (np.arange(nx) / (nx - 1.))
gamma = 7. / 5.
rho = np.zeros((nx, nt + 1))
rhou = np.zeros((nx, nt + 1))
rhoe = np.zeros((nx, nt + 1))
rhov = np.zeros((nx, nt + 1))
e = e_initial('atmosphere')
v = v_initial('step')
time = np.zeros(nt + 1)
rho[:, 0] = rho_initial('atmosphere')
rhoe[:, 0] = rho[:, 0] * e
rhov[:, 0] = rho[:, 0] * v

# plot initial conditions
# plt.plot(x, rho[:, 0])
# plt.grid()
# plt.show()

# additional arrays needed
xi = np.zeros(nx + 1)
xi[1:nx] = 0.5 * (x[1:nx] + x[0:nx - 1])
xi[0] = 2 * xi[1] - xi[2]
xi[nx] = 2 * xi[nx - 1] - xi[nx - 2]
dx = (xi[1:nx + 1] - xi[0:nx])

counter = 1
print_counter = 0
debug_counter = 1
time_interval = 10

debug = False

for it in range(1, nt + 1):
    qrho = rho[:, it - 1]
    qrhou = rhou[:, it - 1]
    qrhoe = rhoe[:, it - 1]
    qrhov = rhov[:, it - 1]

    # remove ghost cells in dt calculation
    cs = np.sqrt(gamma * (gamma - 1) * abs(qrhoe[2:-2]/qrho[2:-2]))
    dum = dx[2:-2] / (cs + abs(qrhou[2:-2] / qrho[2:-2]))
    dt = cfl * min(dum)

    time[it] = time[it - 1] + dt

    # debugging phase
    if debug:
        np.savetxt('debug/log{:03d}.dat'.format(debug_counter),
                   np.c_[cs, dum, qrhou[2:-2] / qrho[2:-2], qrhou[2:-2], qrho[2:-2]], header='{}'.format(dt))
        print("producing log{:03d}".format(debug_counter))
        debug_counter += 1
    # print("Time step: {}, Time = {}, dt = {}".format(it, time[it]/3.154e7, dt/3.154e7))
    # print(dt)
    print_counter += 1
    if print_counter == time_interval:
        print("Time step: {}, Time = {}, dt = {}".format(it, time[it]/3.154e7, dt/3.154e7))
        np.savetxt('outputdisc/output{:05d}.dat'.format(counter),
                   np.c_[x, qrho, qrhou, qrhoe, qrhov], header='{}'.format(time[it - 1]))
        print_counter = 0
        counter += 1
    hydrostep(x, xi, qrho, qrhou, qrhoe, qrhov, gamma, dt, 'disc', 'van Leer')
    rho[:, it] = qrho
    rhou[:, it] = qrhou
    rhoe[:, it] = qrhoe
    rhov[:, it] = qrhov
