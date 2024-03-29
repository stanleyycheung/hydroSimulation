'''
1D Hydrodyanmics evolution code

Written by SC
'''

import numpy as np


class discSimulation:
    def __init__(self, nx, nt, x0, x1, cfl, init_rho, init_e, init_v, bc, fluxlim, debug=False):
        self._nx = nx
        self._nt = nt
        self._x0 = x0
        self._x1 = x1
        self._cfl = cfl
        self._init_e = init_e
        self._init_v = init_v
        self._debug = debug
        self._bc = bc
        self._fluxlim = fluxlim

        self._counter = 1
        self._print_counter = 0
        self._debug = debug
        if self._debug:
            self._debug_counter = 1

        self._x = self._x0 + (self._x1 - self._x0) * (np.arange(self._nx) / (self._nx - 1.))
        self._rho = np.zeros((self._nx, self._nt + 1))
        self._rhou = np.zeros((self._nx, self._nt + 1))
        self._rhoe = np.zeros((self._nx, self._nt + 1))
        self._rhov = np.zeros((self._nx, self._nt + 1))
        self._time = np.zeros(self._nt + 1)

        self._xi = np.zeros(self._nx + 1)
        self._xi[1:self._nx] = 0.5 * (self._x[1:self._nx] + self._x[0:self._nx - 1])
        self._xi[0] = 2 * self._xi[1] - self._xi[2]
        self._xi[self._nx] = 2 * self._xi[self._nx - 1] - self._xi[self._nx - 2]
        self._dx = (self._xi[1:self._nx + 1] - self._xi[0:self._nx])

        self._e = self.e_initial(init_e)
        self._v = self.v_initial(init_v)
        self._rho[:, 0] = self.rho_initial(init_rho)
        self._rhoe[:, 0] = self._rho[:, 0] * self._e
        self._rhov[:, 0] = self._rho[:, 0] * self._v

    def set_boundary(self):
        # Get the number of grid points including the ghost cells
        nx = len(self._qrho)  # maybe redudant
        # If periodic, then install periodic BC, using two ghost cells on each side
        # (two are required for the non-lin flux limiter)
        if self._bc == 'periodic':
            self._qrho[0] = self._qrho[nx - 4]
            self._qrho[1] = self._qrho[nx - 3]
            self._qrho[nx - 2] = self._qrho[2]
            self._qrho[nx - 1] = self._qrho[3]
            self._qrhou[0] = self._qrhou[nx - 4]
            self._qrhou[1] = self._qrhou[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[2]
            self._qrhou[nx - 1] = self._qrhou[3]
        # If mirror symmetry, then install mirror BC, using two ghost cells
        # on each side (two are required for the non-lin flux limiter)
        elif self._bc == 'mirror':
            self._qrho[0] = self._qrho[3]
            self._qrho[1] = self._qrho[2]
            self._qrho[nx - 2] = self._qrho[nx - 3]
            self._qrho[nx - 1] = self._qrho[nx - 4]
            self._qrhou[0] = -self._qrhou[3]
            self._qrhou[1] = -self._qrhou[2]
            self._qrhou[nx - 2] = -self._qrhou[nx - 3]
            self._qrhou[nx - 1] = -self._qrhou[nx - 4]
        elif self._bc == 'mirror/free':
            # Left BC
            self._qrho[0] = self._qrho[2]
            self._qrho[1] = self._qrho[3]
            self._qrhou[0] = -self._qrhou[2]
            self._qrhou[1] = -self._qrhou[3]
            self._qrhov[0] = -self._qrhov[2]
            self._qrhov[1] = -self._qrhov[3]
            # Right BC
            self._qrho[nx - 2] = self._qrho[nx - 4]
            self._qrho[nx - 1] = self._qrho[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[nx - 4]
            self._qrhou[nx - 1] = self._qrhou[nx - 3]
            self._qrhov[nx - 2] = self._qrhov[nx - 4]
            self._qrhov[nx - 1] = self._qrhov[nx - 3]
        elif self._bc == 'disc':
            # Left BC
            self._qrho[0] = N * mu  # equation 1.4
            self._qrho[1] = N * mu
            self._qrhou[0] = self._qrhou[2]/self._qrho[2] * N * mu
            self._qrhou[1] = self._qrhou[3]/self._qrho[3] * N * mu
            self._qrhov[0] = np.sqrt(G*M_star/(self._x[0]+starDist)) * self._qrho[0]
            self._qrhov[1] = np.sqrt(G*M_star/(self._x[1]+starDist)) * self._qrho[1]
            # Right BC
            self._qrho[nx - 2] = self._qrho[nx - 4]
            self._qrho[nx - 1] = self._qrho[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[nx - 4]
            self._qrhou[nx - 1] = self._qrhou[nx - 3]
            self._qrhov[nx - 2] = self._qrhov[nx - 4]
            self._qrhov[nx - 1] = self._qrhov[nx - 3]
        else:
            raise Exception("Choose valid BC")

    def e_initial(self, init_e):
        if init_e == 'flat':
            result = np.ones(self._nx)
        # Sod shock tube
        elif init_e == 'sod':
            result = np.zeros(self._nx)
            for i in range(0, len(result)):
                if i < int(len(result) / 2):
                    result[i] = 2.5
                else:
                    result[i] = 2.
        # Atmosphere
        elif init_e == 'atmosphere':
            # typical pressure is 100kPa
            result = np.full(len(self._x), 2.5e5)
        elif init_e == 'zero':
            result = np.zeros(self._nx)
        else:
            raise Exception("Choose a valid initial e")
        return result

    def v_initial(self, init_v):
        # Keplarian
        if init_v == 'kep':
            result = np.zeros(self._nx)
            result[0] = np.sqrt(G*M_star/(self._x[0]+starDist))
        elif init_v == 'step':
            result = np.zeros(self._nx)
            for i in range(0, len(result)):
                if i < int(len(result) / 2):
                    result[i] = 10
                else:
                    result[i] = 5
        elif init_v == 'zero':
            result = np.zeros(self._nx)
        else:
            raise Exception("Choose a valid initial v")
        return result

    def rho_initial(self, init_rho):
        # Gaussian
        if init_rho == 'gaussian':
            xmid = 0.5 * (self._x0 + self._x1)
            dg = 0.1 * (self._x1 - self._x0)
            result = 1 + 3 * np.exp(-(self._x - xmid)**2 / dg**2)
        # Sod shock tube
        elif init_rho == 'sod':
            result = np.zeros(len(self._x))
            for i in range(0, len(self._x)):
                if i < int(len(self._x) / 2):
                    result[i] = 1
                else:
                    result[i] = 0.125
        # Atmosphere
        elif init_rho == 'atmosphere':
            # Density of air is ~ 1kg/m^3
            result = np.full(len(self._x), 1)
        elif init_rho == 'zero':
            result = np.zeros(self._nx)
        else:
            raise Exception("Choose a valid density")
        return result

    def advect(self, q, ui, dt, nghost):
        nx = len(self._xi) - 1
        if len(self._x) != nx or len(self._xi) != nx + 1 or len(q) != nx or len(ui) != nx + 1 or nghost < 1:
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
        if self._fluxlim == 'donor-cell':
            phi = np.zeros(nx + 1)
        elif self._fluxlim == 'superbee':
            phi = np.zeros(nx + 1)
            for i in range(1, nx):
                a = min([1., 2. * r[i]])
                b = min([2., r[i]])
                phi[i] = max([0., a, b])
        elif self._fluxlim == 'van Leer':
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
                abs(ui[i]) * (1 - abs(ui[i] * dt / (self._x[i] - self._x[i - 1]))) * \
                phi[i] * (q[i] - q[i - 1])

        # Update the cells, except the ghost cells
        for i in range(nghost, nx - nghost):
            q[i] = q[i] - dt * (flux[i + 1] - flux[i]) / (self._xi[i + 1] - self._xi[i])

    def hydrostep(self):
        # Use 2 ghost cells on each side
        nghost = 2
        # Get the number of grid points including the ghost cells
        nx = len(self._x)
        # Impose boundary conditions
        self.set_boundary()
        # Compute the velocities at the cell interfaces
        ui = np.zeros(nx + 1)
        vi = np.zeros(nx + 1)
        for ix in range(1, nx):
            ui[ix] = 0.5 * (self._qrhou[ix] / self._qrho[ix] +
                            self._qrhou[ix - 1] / self._qrho[ix - 1])
            vi[ix] = 0.5 * (self._qrhov[ix] / self._qrho[ix] +
                            self._qrhov[ix - 1] / self._qrho[ix - 1])
        # Advect rho
        self.advect(self._qrho, ui, self._dt, nghost)
        # Advect rhou
        self.advect(self._qrhou, ui, self._dt, nghost)
        # Advect rhoe
        self.advect(self._qrhoe, ui, self._dt, nghost)
        # Advect rhov
        self.advect(self._qrhov, ui, self._dt, nghost)
        # Re-impose boundary conditions
        self.set_boundary()

        # Compute the pressure for each quantity
        u = self._qrhou / self._qrho
        v = self._qrhov / self._qrho
        e_tot = self._qrhoe / self._qrho
        e_kin = u**2 / 2
        e_th = e_tot - e_kin
        p = (gamma - 1.) * self._qrho * e_th

        # Calculate graviational potiental
        V = G * M_star / (self._x + starDist)
        # Now add all external forces, for all cells except the ghost cells
        for ix in range(2, nx - 2):
            self._qrhou[ix] = self._qrhou[ix] - self._dt * \
                (p[ix + 1] - p[ix - 1] + self._qrho[ix + 1] * V[ix + 1] -
                 self._qrho[ix - 1] * V[ix - 1]) / (2 * (self._x[ix + 1] - self._x[ix - 1]))
            self._qrhoe[ix] = self._qrhoe[ix] - self._dt * (p[ix + 1] * u[ix + 1] -
                                                            p[ix - 1] * u[ix - 1] + self._qrho[ix + 1] * V[ix + 1] * u[ix + 1] -
                                                            self._qrho[ix - 1] * V[ix - 1] * u[ix - 1]) / (2 * (self._x[ix + 1] - self._x[ix - 1]))
            # impliment rotational force
            self._qrhou[ix] = self._qrhou[ix] + self._dt * self._qrho[ix]*v[ix]**2/self._x[ix]
        # Re-impose boundary conditions a last time (not strictly necessary)
        self.set_boundary()

    def run(self, time_interval, verbose=False):
        for it in range(1, self._nt+1):
            self._qrho = self._rho[:, it - 1]
            self._qrhou = self._rhou[:, it - 1]
            self._qrhoe = self._rhoe[:, it - 1]
            self._qrhov = self._rhov[:, it - 1]
            # removes ghost cells in dt calculation (can impliment back?)
            cs = np.sqrt(gamma * (gamma - 1) * abs(self._qrhoe[2:-2]/self._qrho[2:-2]))
            dum = self._dx[2:-2] / (cs + abs(self._qrhou[2:-2] / self._qrho[2:-2]))
            self._dt = self._cfl * min(dum)
            self._time[it] = self._time[it - 1] + self._dt
            if self._debug:
                np.savetxt('debug/log{:03d}.dat'.format(self._debug_counter),
                           np.c_['debug text'], header='debug text')
                print("producing log{:03d}".format(self._debug_counter))
                self._debug_counter += 1
            self._print_counter += 1
            if self._print_counter == time_interval:
                print("Time step: {}, Time = {}, dt = {}".format(
                    it, self._time[it]/3.154e7, self._dt/3.154e7))
                np.savetxt('outputdisc/output{:05d}.dat'.format(self._counter),
                           np.c_[self._x, self._qrho, self._qrhou, self._qrhoe, self._qrhov], header='{}'.format(self._time[it - 1]))
                self._print_counter = 0
                self._counter += 1
            self.hydrostep()
            self._rho[:, it] = self._qrho
            self._rhou[:, it] = self._qrhou
            self._rhoe[:, it] = self._qrhoe
            self._rhov[:, it] = self._qrhov


# physical constants - mass in solar masses, length in m
G = 6.67e-11
AU = 1.496e11
M_solar = 1.9891e30
M_star = M_solar
starDist = 100 * AU
m_h = 1.673e-27
mu = 1.3 * m_h  # following Facchini et al. 2016
N = 10**9.2
T = 20
gamma = 7. / 5.


if __name__ == '__main__':
    # run simulation
    print('running main')
    simulation = discSimulation(2500, 10000, 0, 1200*AU, 0.3, 'atmosphere',
                                'atmosphere', 'kep', 'disc', 'van Leer')
    simulation.run(10)
