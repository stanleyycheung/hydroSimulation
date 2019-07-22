"""
More advanced 1-D hydro solver
Section 5.2
"""

import numpy as np
import matplotlib.pyplot as plt


def advect(x, xi, q, ui, dt, fluxlim, nghost):
    nx = len(xi) - 1
    if len(x) != nx or len(xi) != nx+1 or len(q) != nx or len(ui) != nx+1 or nghost < 1:
        raise ValueError('Wrong input shape')

    # Determine the r_{i-1/2} for the flux limiter
    r = np.zeros(nx+1)
    for i in range(2, nx-1):
        dq = (q[i]-q[i-1])
        if abs(dq) > 0:
            if ui[i] > 0:
                r[i] = (q[i-1]-q[i-2])/dq
            else:
                r[i] = (q[i+1]-q[i])/dq

    #  Determine the flux limiter (many other flux limiters can be implemented here!)
    if fluxlim == 'donor-cell':
        phi = np.zeros(nx+1)
    elif fluxlim == 'superbee':
        phi = np.zeros(nx+1)
        for i in range(1, nx):
            a = min([1., 2.*r[i]])
            b = min([2., r[i]])
            phi[i] = max([0., a, b])
    else:
        raise Exception("Choose a valid flux limiter")

    # Now construct the flux
    flux = np.zeros(nx+1)
    for i in range(1, nx):
        if ui[i] > 0:
            flux[i] = ui[i] * q[i-1]
        else:
            flux[i] = ui[i] * q[i]
        flux[i] = flux[i] + 0.5 * \
            abs(ui[i]) * (1-abs(ui[i]*dt/(x[i]-x[i-1]))) * phi[i] * (q[i]-q[i-1])

    # Update the cells, except the ghost cells
    for i in range(nghost, nx-nghost):
        q[i] = q[i] - dt * (flux[i+1]-flux[i]) / (xi[i+1] - xi[i])


def boundary(rho, rhou, bc):
    # Get the number of grid points including the ghost cells
    nx = len(rho)

    # If periodic, then install periodic BC, using two ghost cells on each side
    # (two are required for the non-lin flux limiter)
    if bc == 'periodic':
        rho[0] = rho[nx-4]
        rho[1] = rho[nx-3]
        rho[nx-2] = rho[2]
        rho[nx-1] = rho[3]
        rhou[0] = rhou[nx-4]
        rhou[1] = rhou[nx-3]
        rhou[nx-2] = rhou[2]
        rhou[nx-1] = rhou[3]

    # If mirror symmetry, then install mirror BC, using two ghost cells
    # on each side (two are required for the non-lin flux limiter)
    elif bc == 'mirror':
        rho[0] = rho[3]
        rho[1] = rho[2]
        rho[nx-2] = rho[nx-3]
        rho[nx-1] = rho[nx-4]
        rhou[0] = -rhou[3]
        rhou[1] = -rhou[2]
        rhou[nx-2] = -rhou[nx-3]
        rhou[nx-1] = -rhou[nx-4]


def hydrostep(x, xi, rho, rhou, e, gamma, dt, bc, fluxlim='donor-cell'):
    # Use 2 ghost cells on each side
    nghost = 2

    # Get the number of grid points including the ghost cells
    nx = len(x)
    # Impose boundary conditions
    boundary(rho, rhou, bc)

    # Compute the velocity at the cell interfaces
    ui = np.zeros(nx+1)
    for ix in range(1, nx):
        ui[ix] = 0.5 * (rhou[ix]/rho[ix] + rhou[ix-1]/rho[ix-1])

    # Advect rho
    advect(x, xi, rho, ui, dt, fluxlim, nghost)

    # Re-impose boundary conditions
    boundary(rho, rhou, bc)

    # Compute the pressure
    p = (gamma-1.)*rho*e

    # Now add the pressure force, for all cells except the ghost cells
    for ix in range(2, nx-2):
        rhou[ix] = rhou[ix] - dt*(p[ix+1]-p[ix-1])/(x[ix+1]-x[ix-1])

    # Re-impose boundary conditions a last time (not strictly necessary)
    boundary(rho, rhou, bc)
