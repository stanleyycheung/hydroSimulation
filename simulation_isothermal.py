'''
1D Hydrodyanmics evolution code
Isothermal edition - no advection of the energy equation

Written by SC
'''

import numpy as np


class discSimulation:
    def __init__(self, nx, nt, x0, x1, cfl, init_rho, temp, init_v, bc, fluxlim, rotation=True):
        '''
        Main attributes:
            Init attributes:
            _x (numpy arr): grid of x values
            _x0 and _x1 (float): start and end of grid
            _nx (int): number of cells
            _nt (int): number of time steps
            _cfl (float): CFL parameter (usually 0.3)
            _init_v/_init_rho (string): specifies initial velocity and density
            _T (numpy arr): temperature across grid
            _bc (string): specifies nature of boundary conditions
            _rotation (bool): if rotation is considered or not

            Created attributes:
            _rho/_rhou/_rhov (numpy arr(_nx, _nt + 1)): stores data of each time step in column
            _xi (numpy arr(_nx + 1)): cell interfaces
            _dx (numpy arr(_nx)): difference in _x
            _time (numpy arr(_nx + 1)): keeps tracks of the time with dt steps
        '''
        self._nx = nx
        self._nt = nt
        self._x0 = x0
        self._x1 = x1
        self._cfl = cfl
        self._init_v = init_v
        self._bc = bc
        self._fluxlim = fluxlim
        self._T = self.set_T(temp)
        self._rotation = rotation

        # init counters
        self._counter = 1
        self._print_counter = 0
        self._init = True

        # generate arrays needed
        self._x = self._x0 + (self._x1 - self._x0) * (np.arange(self._nx) / (self._nx - 1.))
        # need to check x
        self._rho = np.zeros((self._nx, self._nt + 1))
        self._rhou = np.zeros((self._nx, self._nt + 1))
        if self._rotation:
            self._rhov = np.zeros((self._nx, self._nt + 1))
        self._time = np.zeros(self._nt + 1)
        self._xi = np.zeros(self._nx + 1)
        self._xi[1:self._nx] = 0.5 * (self._x[1:self._nx] + self._x[0:self._nx - 1])
        self._xi[0] = 2 * self._xi[1] - self._xi[2]
        self._xi[self._nx] = 2 * self._xi[self._nx - 1] - self._xi[self._nx - 2]
        self._dx = (self._xi[1:self._nx + 1] - self._xi[0:self._nx])

        # stores initial conditions
        self._v = self.v_initial(init_v)
        self._rho[:, 0] = self.rho_initial(init_rho)
        if self._rotation:
            self._rhov[:, 0] = self._rho[:, 0] * self._v

    def set_boundary(self):
        '''
        Sets _qrho/_qrhou/_qrhov corresponding to the boundary conditions
        Options:
            - 'periodic' - left and right: periodic
            - 'mirror' - left and right: reflective
            - 'mirror/free' - left: reflective; right: free outflow, no inflow
            - 'disc' - left: fixed density and Keplarian rotation; right: free outflow/inflow
            - 'free/ free' (untested) - left and right: free outflow/inflow
        '''
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
            self._qrho[0] = 101325 * m_air / (k * T)
            self._qrho[1] = 101325 * m_air / (k * T)
            self._qrhou[0] = 0
            self._qrhou[1] = 0
            # self._qrho[0] = self._qrho[3]
            # self._qrho[1] = self._qrho[2]
            # self._qrhou[0] = -self._qrhou[3]
            # self._qrhou[1] = -self._qrhou[2]
            # Right BC
            self._qrho[nx - 2] = self._qrho[nx - 3]
            self._qrho[nx - 1] = self._qrho[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[nx - 3]
            self._qrhou[nx - 1] = self._qrhou[nx - 3]
            if self._rotation:
                self._qrhov[0] = -self._qrhov[2]
                self._qrhov[1] = -self._qrhov[3]
                self._qrhov[nx - 2] = self._qrhov[nx - 4]
                self._qrhov[nx - 1] = self._qrhov[nx - 3]
        elif self._bc == 'disc':
            # Left BC
            self._qrho[0] = N * mu * m_h
            self._qrho[1] = N * mu * m_h
            self._qrhou[0] = self._qrhou[2]/self._qrho[2] * self._qrho[0]
            self._qrhou[1] = self._qrhou[3]/self._qrho[3] * self._qrho[1]
            # Right BC
            self._qrho[nx - 2] = self._qrho[nx - 4]
            self._qrho[nx - 1] = self._qrho[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[nx - 4]
            self._qrhou[nx - 1] = self._qrhou[nx - 3]
            if self._rotation:
                self._qrhov[0] = np.sqrt(G*M_star/(self._x[0]+starDist)) * self._qrho[0]
                self._qrhov[1] = np.sqrt(G*M_star/(self._x[1]+starDist)) * self._qrho[1]
                self._qrhov[nx - 2] = self._qrhov[nx - 4]
                self._qrhov[nx - 1] = self._qrhov[nx - 3]
        elif self._bc == 'free/free':
            self._qrho[0] = self._qrho[2]
            self._qrho[1] = self._qrho[3]
            self._qrhou[0] = self._qrhou[2]
            self._qrhou[1] = self._qrhou[3]
            self._qrho[nx - 2] = self._qrho[nx - 4]
            self._qrho[nx - 1] = self._qrho[nx - 3]
            self._qrhou[nx - 2] = self._qrhou[nx - 4]
            self._qrhou[nx - 1] = self._qrhou[nx - 3]
        else:
            raise Exception("Choose valid BC")

    def set_T(self, temp):
        '''
        Sets _T to corresponding value
        Options:
            - 'default': all at 100K except first 2 cells at 20K - used in psuedo disc model
            - 'flat': all at 288.16K - used in atmosphere test
        '''
        if temp == 'default':
            result = np.full(self._nx, 100)
            result[0] = 20
            result[1] = 20
        elif temp == 'flat':
            result = np.full(self._nx, 288.16)
        else:
            raise Exception("Choose a valid temperature")
        return result

    def v_initial(self, init_v):
        '''
        Sets _v to corresponding value
        Options:
            - 'kep': Keplarian rotation
            - 'step': half at 10m/s, half at 5m/s - used to test for advection
            - 'zero': set to 10e-5m/s
        '''
        if init_v == 'kep':
            result = np.sqrt(G*M_star/(self._x+starDist))
        elif init_v == 'step':
            result = np.zeros(self._nx)
            for i in range(0, len(result)):
                if i < int(len(result) / 2):
                    result[i] = 10
                else:
                    result[i] = 5
        elif init_v == 'zero':
            result = np.full(self._nx, 10e-5)
        else:
            raise Exception("Choose a valid initial v")
        return result

    def rho_initial(self, init_rho):
        '''
        Sets _rho to corresponding value
        Options:
            - 'gaussian': set to Gaussian with sigma as 10% of range - used in examples
            - 'sod': half at 1kg/m^3, half at 0.125kg/m^3
            - 'atmosphere': all at 0.788kg/m^3 corresponds to
            - 'exp':
            - 'zero':
        '''
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
        elif init_rho == 'atmosphere':
            # Density of air at 300k
            result = np.full(len(self._x), 0.7880072463768116)
        elif init_rho == 'exp':
            result = 1.2261845353019805 * \
                np.exp(- self._x * 9.80665 * 0.02896968 / 288.16 / 8.314462618)
        elif init_rho == 'zero':
            result = np.full(self._nx, N * mu * m_h / 100)
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
            if self._rotation:
                vi[ix] = 0.5 * (self._qrhov[ix] / self._qrho[ix] +
                                self._qrhov[ix - 1] / self._qrho[ix - 1])
        # Advect rho
        self.advect(self._qrho, ui, self._dt, nghost)
        # Advect rhou
        self.advect(self._qrhou, ui, self._dt, nghost)
        # Advect rhov
        if self._rotation:
            self.advect(self._qrhov, vi, self._dt, nghost)
        # Re-impose boundary conditions
        self.set_boundary()
        # Compute the rotational velocity and pressure
        self._p = self._qrho * k * self._T / (mu * m_h)
        if self._rotation:
            v = self._qrhov / self._qrho
        # Calculate graviational potiental (this could be removed out of loop as this is constant)
        F = G * M_star / (self._x + starDist)**2
        # Now add all external forces, for all cells except the ghost cells
        for ix in range(2, nx - 2):
            self._qrhou[ix] = self._qrhou[ix] - self._dt * \
                (self._p[ix + 1] - self._p[ix - 1]) / (2 * (self._x[ix + 1] - self._x[ix - 1]))
            # gravity implimentation
            self._qrhou[ix] = self._qrhou[ix] - self._dt * self._qrho[ix] * F[ix]
            # impliment rotational force
            if self._rotation:
                self._qrhou[ix] = self._qrhou[ix] + self._dt * self._qrhov[ix]*v[ix]/self._x[ix]
        # Re-impose boundary conditions a last time (not strictly necessary)
        self.set_boundary()

    def run(self, time_interval=10, verbose=True, debug=False):
        if debug:
            self._debug_counter = 1
        for it in range(1, self._nt+1):
            self._qrho = self._rho[:, it - 1]
            self._qrhou = self._rhou[:, it - 1]
            if self._rotation:
                self._qrhov = self._rhov[:, it - 1]

            # Calculates time step with CFL condition
            cs = np.sqrt(k * self._T / (mu * m_h))
            dum = self._dx / (cs + abs(self._qrhou / self._qrho))
            self._dt = self._cfl * min(dum)
            self._time[it] = self._time[it - 1] + self._dt

            # Debug log
            if debug:
                np.savetxt('debug/log{:03d}.dat'.format(self._debug_counter),
                           np.c_[self._qrho, self._qrhou], header='{}'.format(self._dt))
                print("producing log{:03d}".format(self._debug_counter))
                self._debug_counter += 1
            # Deals with saving the initial conditions if time interval is not 1
            if (self._init and time_interval != 1):
                if verbose:
                    print("Time step: {}, Time = {}, dt = {}".format(
                        it, self._time[it], self._dt))
                if self._rotation:
                    np.savetxt('outputdisc_iso_grav/output{:05d}.dat'.format(self._counter),
                               np.c_[self._x, self._qrho, self._qrhou, self._qrhov], header='{}'.format(self._time[it - 1]))
                else:
                    np.savetxt('outputdisc_iso_grav/output{:05d}.dat'.format(self._counter),
                               np.c_[self._x, self._qrho, self._qrhou], header='{}'.format(self._time[it - 1]))
                self._counter += 1
                self._init = False
            # deals with printing and saving results per time interval
            self._print_counter += 1
            if self._print_counter == time_interval:
                if verbose:
                    print("Time step: {}, Time = {}, dt = {}".format(
                        it, self._time[it], self._dt))
                if self._rotation:
                    np.savetxt('outputdisc_iso_grav/output{:05d}.dat'.format(self._counter),
                               np.c_[self._x, self._qrho, self._qrhou, self._qrhov], header='{}'.format(self._time[it - 1]))
                else:
                    np.savetxt('outputdisc_iso_grav/output{:05d}.dat'.format(self._counter),
                               np.c_[self._x, self._qrho, self._qrhou], header='{}'.format(self._time[it - 1]))
                self._print_counter = 0
                self._counter += 1
            # proceeds in one step
            self.hydrostep()
            # stores results
            self._rho[:, it] = self._qrho
            self._rhou[:, it] = self._qrhou
            if self._rotation:
                self._rhov[:, it] = self._qrhov


# physical constants - mass in solar masses, length in m
G = 6.67e-11
AU = 1.496e11
M_solar = 1.9891e30
# M_solar = M_star
M_star = 5.972e24
# starDist = 100 * AU
starDist = 6.371e6
k = 1.38e-23
m_h = 1.673e-27
mu = 1.3  # following Facchini et al. 2016
N = 10**9.2
gamma = 7. / 5.
T = 288.16
m_air = 28.97/1000/(6.02e23)

if __name__ == '__main__':
    '''
    # disc simulation
    simulation = discSimulation(2500, 20000, 0, 1200*AU, 0.3, 'zero',
                                'default', 'zero', 'disc', 'van Leer')
    simulation.run(time_interval=20, debug=False, movie=False)
    '''
    # gravity test
    simulation = discSimulation(500, 50000, 0, 2e4, 0.3, 'exp',
                                'flat', 'kep', 'mirror/free', 'van Leer', rotation=False)
    simulation.run(time_interval=50, debug=False)
