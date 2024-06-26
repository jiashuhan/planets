"""
-----------------------------------
Simple Newtonian N-body simulation
-----------------------------------
"""
import sys, numpy as np
from tqdm.auto import trange
from tqdm import tqdm
from scipy.integrate import RK45

class Particles(object):
    """
    Class of objects interating with gravity only

    Parameters
    ----------
    N: int
        number of object slots to be initialized
    n: float, optional
        (Negative) exponent of force law, i.e. F ~ 1/r^n
    epoch: float, optional
        Epoch at start of simulation in Julian date
    """
    G  = 6.6743e-11 # [m^3/kg/s^2]
    AU = 149597870700 # [m]

    def __init__(self, N=10, n=2, epoch=0):
        self.reset(N=N, n=2, epoch=epoch)

    def reset(self, N=10, n=2, epoch=0):
        self.label = 'particles' # name of simulation run
        self.N     = N
        self.n_obj = 0
        self.step  = 0
        self.obj   = {'M':  np.zeros(N),               # masses
                      'kepler': {},                    # orbital parameters of secondary objects
                      'id': np.zeros(N, dtype='<U21'), # identifiers
                      'id_index': {},                  # id to index
                      'parent': {},                    # identifiers of parent objects
                      'xi': np.zeros((N, 3)),          # initial positions [m]; axes: obj, coord
                      'vi': np.zeros((N, 3))           # initial velocities [m/s]; axes: obj, coord
                      }
        self.n     = n # force index
        self.epoch = 0

    # add additional object
    def add(self, mass, x, y, z, vx, vy, vz, nm, parent):
        """
        Parameters
        ----------
        mass: float
            Mass of object [kg]
        x, y, z: float
            Initial position [m] of object, measured in the
            initial rest frame of the central body.
        vx, vy, vz: float
            Initial velocity [m/s] of object, measured in the
            initial rest frame of the central body.
        nm: str
            Identifier of object
        parent: str
            Identifier of parent object
        """
        self.obj['id_index'][nm] = self.n_obj
        self.obj['parent'][nm] = parent

        if self.n_obj >= self.N: # more objects than number of slots
            self.obj['M']  = np.append(self.obj['M'], mass)
            self.obj['id'] = np.append(self.obj['id'], nm)

            x_ = np.zeros((1,3))
            v_ = np.zeros((1,3))
            x_[0] = [x, y, z]
            v_[0] = [vx, vy, vz]
            self.obj['xi'] = np.append(self.obj['xi'], x_, axis=0)
            self.obj['vi'] = np.append(self.obj['vi'], v_, axis=0)

            self.N += 1
        else:
            self.obj['M'][self.n_obj]  = mass
            self.obj['id'][self.n_obj] = nm
            self.obj['xi'][self.n_obj] = [x, y, z]
            self.obj['vi'][self.n_obj] = [vx, vy, vz]

        self.n_obj += 1

    def update(self, dt, method='rkdp'):
        """
        Update position and velocity of objects

        Parameters
        ----------
        dt: float
            Time step [s]
        method: str
            Integration method; default = 'rkdp' (Dormand-Prince RK5(4))
        """
        V_n = self.obj['v'][self.step] # current velocities; axes: obj,coord
        X_n = self.obj['x'][self.step] # current positions; axes: obj,coord

        if method == 'euler': # symplectic Euler (1st order)
            k1_v = self.acc(X_n) # force based on current state; axes: obj,coord
            
            V_np1 = V_n + dt * k1_v # axes: obj,coord
            X_np1 = X_n + dt * V_np1 # use V_np1 instead of V_n (Forward Euler)

        elif method == 'verlet': # velocity Verlet (2nd order)
            k1_v = self.acc(X_n)

            X_np1 = X_n + dt * V_n + dt**2 * 0.5 * k1_v
            V_np1 = V_n + dt * 0.5 * (k1_v + self.acc(X_np1))

        elif method == 'leapfrog': # Leapfrog (2nd order)
            V_n_half = V_n + 0.5 * dt * self.acc(X_n)
            X_np1 = X_n + dt * V_n_half
            V_np1 = V_n_half + 0.5 * dt * self.acc(X_np1)

        elif method == 'rk4': # Runge-Kutta (4th order)
            k1_v = self.acc(X_n) # force based on current state; axes: obj,coord
            k1_x =          V_n

            k2_v = self.acc(X_n + dt * 0.5 * k1_x) # force based on half-step ahead
            k2_x =          V_n + dt * 0.5 * k1_v

            k3_v = self.acc(X_n + dt * 0.5 * k2_x)
            k3_x =          V_n + dt * 0.5 * k2_v

            k4_v = self.acc(X_n + dt * k3_x)
            k4_x =          V_n + dt * k3_v

            V_np1 = V_n + dt * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6
            X_np1 = X_n + dt * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6

        elif method == 'rkdp': # Dormand-Prince RK5(4) (5th/4th order embedded)
            k1_v = self.acc(X_n) # force based on current state; axes: obj,coord
            k1_x =          V_n

            k2_v = self.acc(X_n + dt * 1/5 * k1_x)
            k2_x =          V_n + dt * 1/5 * k1_v

            k3_v = self.acc(X_n + dt * (3/40 * k1_x + 9/40 * k2_x))
            k3_x =          V_n + dt * (3/40 * k1_v + 9/40 * k2_v)

            k4_v = self.acc(X_n + dt * (44/45 * k1_x - 56/15 * k2_x + 32/9 * k3_x))
            k4_x =          V_n + dt * (44/45 * k1_v - 56/15 * k2_v + 32/9 * k3_v)

            k5_v = self.acc(X_n + dt * (19372/6561 * k1_x - 25360/2187 * k2_x + 64448/6561 * k3_x - 212/729 * k4_x))
            k5_x =          V_n + dt * (19372/6561 * k1_v - 25360/2187 * k2_v + 64448/6561 * k3_v - 212/729 * k4_v)

            k6_v = self.acc(X_n + dt * (9017/3168 * k1_x - 355/33 * k2_x + 46732/5247 * k3_x + 49/176 * k4_x - 5103/18656 * k5_x))
            k6_x =          V_n + dt * (9017/3168 * k1_v - 355/33 * k2_v + 46732/5247 * k3_v + 49/176 * k4_v - 5103/18656 * k5_v)

            k7_x = V_n + dt * (35/384 * k1_v + 500/1113 * k3_v + 125/192 * k4_v - 2187/6784 * k5_v + 11/84 * k6_v)

            # 5th order solution
            X_np1 = X_n + dt * (35/384 * k1_x + 500/1113 * k3_x + 125/192 * k4_x - 2187/6784 * k5_x + 11/84 * k6_x)
            V_np1 = V_n + dt * (35/384 * k1_v + 500/1113 * k3_v + 125/192 * k4_v - 2187/6784 * k5_v + 11/84 * k6_v)

            # 4th order solution
            X_np1_ = X_n + dt * (5719/57600 * k1_x + 7571/16695 * k3_x + 393/640 * k4_x - 92097/339200 * k5_x + 187/2100 * k6_x + 1/40 * k7_x)

            self.obj['x_err'][self.step + 1] = X_np1 - X_np1_ # Dormand-Prince local position error (5th - 4th)

        elif method == 'adaptive_rkdp': # adaptive Dormand-Prince from scipy.integrate
            self.integrator.step()
            self.integrator_finished = self.integrator.status == 'finished'
            self.t.append(self.integrator.t)

            X_np1, V_np1 = self.Y_unflatten(self.integrator.y)

        else:
            raise Exception("Unknown method.")

        self.obj['v'][self.step + 1] = V_np1
        self.obj['x'][self.step + 1] = X_np1

        self.step += 1

    def acc(self, X):
        """
        Computes the gravitational acceleration on each object.

        Parameters
        ----------
        X: numpy.ndarray
            Positions of objects; axes: obj, coord

        Returns
        -------
        acc_vec: 2d array, float
            Acceleration of each object [m/s^2]; axes: obj, coord
        """
        dX = X[:,np.newaxis] - X # axes: obj1,obj2, coord
        r = np.sqrt(np.sum(dX**2, axis=-1)) # axes: obj1, obj2
        _ = np.finfo(float).eps # avoid div by 0; self force removed by dX
        
        acc_vec = self.G * np.sum(self.obj['M'][:,np.newaxis,np.newaxis] * dX / (r[...,np.newaxis]**(self.n+1) + _), axis=0) # axes: obj,coord

        return acc_vec

    # computes position and speed of the center of mass of the system for all snapshots
    def COM(self):
        assert self.step > 0, "No results."

        M = self.obj['M'] # axes: objects
        X = self.obj['x'] # axes: snapshots, objects, coord
        V = self.obj['v'] # axes: snapshots, objects, coord

        X_com = np.sum(M[np.newaxis,:,np.newaxis] * X, axis=1) / np.sum(M) # axes: snapshots, coord
        V_com = np.sum(M[np.newaxis,:,np.newaxis] * V, axis=1) / np.sum(M) # axes: snapshots, coord

        return X_com, V_com # axes: snapshots, coord

    def run(self, duration, res, method='rkdp', max_dt=8640, avg_res_guess=500):
        """
        Runs the simulation and updates the objects. First compute
        forces acting on all objects, then update the position and
        velocitiy of each object.

        Parameters
        ----------
        duration: float
            Duration of simulation [d]
        res: float
            Number of steps per day (for non-adaptive methods), or 
            Maximum local position error [m] (for adaptive methods)
        method: str
            Integration method; default = 'rkdp' (Dormand-Prince RK5(4))
        max_dt: float
            Maximum step size [s] (applies to adaptive methods only); 
            default = 8640 (1/10 of a day)
        avg_res_guess: float
            (Over-) estimated time resolution [1/d] (applies to adaptive 
            methods only); default = 500

        Returns
        -------
        obj: dict
            Properties of time evolved objects.

        Notes
        -----
        (1) There are a total of (Nsteps + 1) steps/snapshots, with the first 
            one (0) being the initial state.
        (2) For non-adaptive methods, the steps are equally spaced in time,
            and the step size is set by the user.
        (3) For adaptive methods, the step size is determined by the RK45 solver,
            with a maximum step size imposed (unlikely to have any effect for 
            sufficiently small local position error). The total number of steps 
            is (intentionally) over-estimated beforehand with a guess for the 
            average time resolution (steps per day) so that there is enough 
            space allocated to store the results.
        (4) An average time resolution of 500 steps per day is a good guess for
            max position error >= 0.1 mm.
        """
        self.step      = 0 # step 0 is the initial state
        self.obj['M']  = self.obj['M'][:self.n_obj]
        self.obj['id'] = self.obj['id'][:self.n_obj]
        self.obj['xi'] = self.obj['xi'][:self.n_obj]
        self.obj['vi'] = self.obj['vi'][:self.n_obj]

        if method != 'adaptive_rkdp': # non-adaptive methods
            Nsteps = int(duration * res) # Total number of steps in the simulation
            dt     = 86400 / res # timestep [s]
            self.t = np.arange(Nsteps + 1) * dt # simulation time [s] for each snapshot; assuming equally spaced in time

            self.obj['x'] = np.zeros((Nsteps + 1, self.n_obj, 3)) # snapshots of positions
            self.obj['v'] = np.zeros((Nsteps + 1, self.n_obj, 3))

            self.obj['x'][self.step] = self.obj['xi'] # first snapshot is initial state
            self.obj['v'][self.step] = self.obj['vi']

            if method == 'rkdp':
                self.obj['x_err'] = np.zeros_like(self.obj['x']) # Dormand-Prince local position error (5th - 4th)

            for _ in trange(Nsteps):
                self.update(dt, method=method)

            if method == 'rkdp':
                print('Max absolute local position error = %.5f m'%np.abs(self.obj['x_err']).max())
                print('Max relative local position error = %.5f'%np.abs(self.obj['x_err'] / (self.obj['x'] + np.finfo(float).eps)).max())

        else: # adaptive methods
            Y0   = self.Y_flatten(self.obj['xi'], self.obj['vi']) # initial state
            xtol = res # Max local position error [m]
            vtol = 1e-3 * xtol # Max local velocity error [m/s]
            atol = np.zeros_like(Y0); atol[:self.n_obj*3] = xtol; atol[self.n_obj*3:] = vtol

            Nsteps = int(duration * avg_res_guess) # Intentionally overestimated total number of steps

            self.t = [0]
            self.integrator_finished = False

            self.obj['x'] = np.zeros((Nsteps + 1, self.n_obj, 3)) # snapshots of positions
            self.obj['v'] = np.zeros((Nsteps + 1, self.n_obj, 3))

            self.obj['x'][self.step] = self.obj['xi'] # first snapshot is initial state
            self.obj['v'][self.step] = self.obj['vi']
            
            self.integrator = RK45(self.RK_derivative, 0, Y0, duration * 86400, max_step=max_dt, rtol=np.finfo(float).eps*100, atol=atol)

            with tqdm(total=Nsteps) as t:
                while not self.integrator_finished:
                    self.update(max_dt, method=method)
                    t.update()

            # truncate arrays in case integrator finished prematurely (self.step < Nsteps)
            self.t = np.array(self.t)[:self.step + 1]
            self.obj['x'] = self.obj['x'][:self.step + 1]
            self.obj['v'] = self.obj['v'][:self.step + 1]

        X_com, V_com = self.COM()

        return {'t': self.t, 'id': self.obj['id'], 'id_index': self.obj['id_index'], \
                'parent': self.obj['parent'], 'x': self.obj['x'], 'v': self.obj['v'], \
                'M': self.obj['M'], 'x_com': X_com, 'v_com': V_com, 'epoch': self.epoch, \
                'n': self.n, 'kepler': self.obj['kepler']}

    def gen_random(self, N, r0, mmax, vmax):
        """
        Generates a list of N random objects.

        Parameters
        ----------
        N: int
            Number of objects to generate.
        r0: float
            Radius of the sphere [m] within which the objects
            are created.
        mmax: float
            Maximum mass [kg] of the objects.
        vmax: float
            Maximum speed [m/s] in each direction.

        Returns
        -------
        names: list
            Names of random objects
        """
        names = []
        for i in range(N):
            m, x, y, z, vx, vy, vz = randomize(r0, mmax, vmax)
            name = str(i)
            self.add(m, x, y, z, vx, vy, vz, name, 'N/A')
            names.append(name)

        return names

    #-----------------------------------------
    # Methods used by adaptive Dormand-Prince
    def Y_flatten(self, X, V):
        Y = np.zeros(self.n_obj * 6)
        Y[:self.n_obj*3] = X.flatten()
        Y[self.n_obj*3:] = V.flatten()

        return Y

    def Y_unflatten(self, Y):
        # Y is one-dimensional with length = N * 6 where N is the number of objects
        X = Y[:self.n_obj*3].reshape(self.n_obj, 3)
        V = Y[self.n_obj*3:].reshape(self.n_obj, 3)

        return X, V

    def RK_derivative(self, t, Y): # No explicit t dependence
        X, V = self.Y_unflatten(Y)
        
        return self.Y_flatten(V, self.acc(X))
    #-----------------------------------------

def randomize(r0, mmax, vmax):
    """
    Generates randomized parameters for creating an object.
    Values are drawn from a uniform distribution.

    Parameters
    ----------
    r0: float
        Radius of the sphere [m] within which the object 
        is created.
    mmax: float
        Maximum mass of the object [kg].
    vmax: float
        Maximum speed [m/s] in each direction.

    Returns
    -------
    params: tuple
        Tuple containing the mass, initial position and initial 
        velocity of the object.
    """
    r = np.random.uniform(0, r0)
    theta = np.random.uniform(0, np.pi)
    phi = np.random.uniform(0, 2 * np.pi)
    mass = np.random.uniform(0, mmax)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    v = np.random.uniform(-1, 1, size=3) * vmax
    return mass, x, y, z, *v
