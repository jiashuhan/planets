"""
-----------------------------------
Simple Newtonian N-body simulation
-----------------------------------
"""
import sys, numpy as np
from tqdm.auto import trange

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
                      'id': np.zeros(N, dtype='<U21'), # identifiers
                      'id_index':  {},                 # id to index
                      'xi': np.zeros((N, 3)),          # initial positions [m]; axes: obj, coord
                      'vi': np.zeros((N, 3))           # initial velocities [m/s]; axes: obj, coord
                      }
        self.n     = n # force index
        self.epoch = 0

    # add additional object
    def add(self, mass, x, y, z, vx, vy, vz, nm):
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
        """
        self.obj['id_index'][nm] = self.n_obj

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

    # Updates position and velocity of object
    def update(self, dt):
        F = self.gravity() # axes: obj,coord
        M = np.expand_dims(self.obj['M'], axis=-1)

        v = self.obj['v'][self.step] + F / M * dt # axes: obj,coord
        self.obj['v'][self.step+1] = v
        # orbits won't close if order swapped
        x = self.obj['x'][self.step] + v * dt
        self.obj['x'][self.step+1] = x

        self.step += 1

    def gravity(self):
        """
        Computes the gravitational force on each object.

        Returns
        -------
        force: 2d array, float
            Force acting on each object [N]; axes: obj, coord
        """
        M = self.obj['M']
        X = self.obj['x'][self.step] # current positions; axes: obj,coord

        dX = X[:,np.newaxis] - X # axes: obj1,obj2, coord
        r = np.sqrt(np.sum(dX**2, axis=-1)) # axes: obj1, obj2
        _ = np.finfo(float).eps # avoid div by 0; self force removed by dX
        
        Fvec = self.G * np.sum(M[:,np.newaxis,np.newaxis] * M[:,np.newaxis] * dX / (r[...,np.newaxis]**(self.n+1) + _), axis=0) # axes: obj,coord

        return Fvec

    # computes position and speed of the center of mass of the system for all snapshots
    def COM(self):
        assert self.step > 0, "No results."

        M = self.obj['M'] # axes: objects
        X = self.obj['x'] # axes: snapshots, objects, coord
        V = self.obj['v'] # axes: snapshots, objects, coord

        X_com = np.sum(M[np.newaxis,:,np.newaxis] * X, axis=1) / np.sum(M) # axes: snapshots, coord
        V_com = np.sum(M[np.newaxis,:,np.newaxis] * V, axis=1) / np.sum(M) # axes: snapshots, coord

        return X_com, V_com # axes: snapshots, coord

    def run(self, num_steps, sample_rate):
        """
        Runs the simulation and updates the objects. First compute
        forces acting on all objects, then update the position and
        velocitiy of each object.

        Parameters
        ----------
        num_steps: int
            Total number of steps in the simulation.
        sample_rate: float
            Number of steps per day.

        Returns
        -------
        obj: dict
            Properties of time evolved objects.
        """
        self.step      = -1
        self.obj['M']  = self.obj['M'][:self.n_obj]
        self.obj['id'] = self.obj['id'][:self.n_obj]
        self.obj['xi'] = self.obj['xi'][:self.n_obj]
        self.obj['vi'] = self.obj['vi'][:self.n_obj]
        self.obj['x']  = np.zeros((num_steps, self.n_obj, 3)) # snapshots of positions
        self.obj['v']  = np.zeros((num_steps, self.n_obj, 3))

        self.obj['x'][self.step] = self.obj['xi']
        self.obj['v'][self.step] = self.obj['vi']

        dt = 86400 / sample_rate # timestep [s]

        for step in trange(num_steps):
            self.update(dt)

        self.sim_time = np.arange(num_steps) * dt # simulation time [s] for each snapshot
        X_com, V_com = self.COM()

        return {'t': self.sim_time, 'id': self.obj['id'], 'id_index': self.obj['id_index'], \
                'x': self.obj['x'], 'v': self.obj['v'], 'x_com': X_com, 'v_com': V_com,
                'epoch': self.epoch, 'n': self.n}

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
        """
        for i in range(N):
            m, x, y, z, vx, vy, vz = randomize(r0, mmax, vmax)
            self.add(m, x, y, z, vx, vy, vz, str(i))

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
