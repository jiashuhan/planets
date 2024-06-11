"""
--------------------------------------------
(DEPRECATED) Simple solar system simulation
--------------------------------------------
*** This is an earlier version modified to work with the newer UI,
    and has been replaced by particles.py
"""
import sys, numpy as np
from tqdm.auto import trange
import warnings
import functools

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("The function '{}' is deprecated.".format(func.__name__),
                      category=DeprecationWarning,
                      stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    return new_func

class Particles(object):
    """
    Ensemble of objects interating with gravity only

    Parameters
    ----------
    N: int
        Number of object slots (maximum possible number of objects)
    n: float, optional
        (Negative) exponent of force law, i.e. F ~ 1/r^n
    epoch: float, optional
        Epoch at start of simulation in Julian date
    """
    G  = 6.6743e-11 # [m^3/kg/s^2]
    AU = 149597870700 # [m]

    @deprecated
    def __init__(self, N=50, n=2, epoch=0):
        self.reset(N=N, n=2, epoch=epoch)

    def reset(self, N=50, n=2, epoch=0):
        self.label = 'particles' # name of simulation run
        self.n_obj = 0
        self.step  = 0
        self.n     = n # force index
        self.epoch = 0
        self.objects = []
        self.old_objects =[]

        self.obj   = {'M':  np.zeros(N),               # masses
                      'id': np.zeros(N, dtype='<U21'), # identifiers
                      'id_index': {},                  # id to index
                      'xi': np.zeros((N, 3)),          # initial positions [m]; axes: obj, coord
                      'vi': np.zeros((N, 3))           # initial velocities [m/s]; axes: obj, coord
                      }

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
        self.objects.append(Particle(mass, x, y, z, vx, vy, vz, nm))
        self.obj['id'][self.n_obj] = nm
        self.obj['M'][self.n_obj] = mass
        self.obj['xi'][self.n_obj] = np.array([x, y, z])
        self.obj['vi'][self.n_obj] = np.array([vx, vy, vz])
        self.n_obj = len(self.objects)

    # computes position and speed of the center of mass of the system for all snapshots
    def COM(self):
        assert self.step > 0, "No results."

        M = self.obj['M'] # axes: objects
        X = self.obj['x'] # axes: snapshots, objects, coord
        V = self.obj['v'] # axes: snapshots, objects, coord

        X_com = np.sum(M[np.newaxis,:,np.newaxis] * X, axis=1) / np.sum(M) # axes: snapshots, coord
        V_com = np.sum(M[np.newaxis,:,np.newaxis] * V, axis=1) / np.sum(M) # axes: snapshots, coord

        return X_com, V_com # axes: snapshots, coord

    @deprecated
    def run(self, Nsteps, sample_rate, rcrit=3.2e8):
        """
        Runs the simulation and updates the objects. First compute
        forces acting on all objects, then update the position and
        velocitiy of each object.

        Parameters
        ----------
        Nsteps: int
            Total number of steps in the simulation.
        sample_rate: float
            Number of steps per day.
        rcrit: float
            Optional. Distance threshold for merger [m]

        Returns
        -------
        obj: dict
            Properties of time evolved objects.
        """
        self.step = 0
        
        dt = 86400 / sample_rate # timestep [s]

        self.old_objects = []
        
        for _ in trange(Nsteps):
            forces_x = []
            forces_y = []
            forces_z = []

            for obj1 in self.objects:
                total_fx = 0
                total_fy = 0
                total_fz = 0

                for obj2 in self.objects:
                    if not (obj1 is obj2) and obj2.mergeable(obj1): # objects are distinct
                        # Merging works only when sample_rate > 100
                        if sample_rate >= 100 and obj1.distance(obj2) < rcrit:
                            self.objects.remove(obj1)
                            self.objects.remove(obj2)
                            obj1, obj2, obj3 = self.merge(obj1, obj2)
                            self.old_objects.append(obj1)
                            self.old_objects.append(obj2)
                            self.objects.append(obj3)
                        else:
                            force = obj1.gravity(obj2, n=self.n)
                            total_fx += force[0]
                            total_fy += force[1]
                            total_fz += force[2]

                forces_x.append(total_fx)
                forces_y.append(total_fy)
                forces_z.append(total_fz)

            for i in range(len(forces_x)): # update all objects
                self.objects[i].update(forces_x[i], forces_y[i], forces_z[i], dt)

            self.step += 1

        self.sim_time = np.arange(Nsteps + 1) * dt # simulation time [s] for each snapshot

        self.obj['M']  = self.obj['M'][:self.n_obj]
        self.obj['id'] = self.obj['id'][:self.n_obj]
        self.obj['xi'] = self.obj['xi'][:self.n_obj]
        self.obj['vi'] = self.obj['vi'][:self.n_obj]
        self.obj['x']  = np.zeros((Nsteps + 1, self.n_obj, 3)) # snapshots of positions
        self.obj['v']  = np.zeros((Nsteps + 1, self.n_obj, 3))

        for i, o in enumerate(self.objects):
            self.obj['M'][i]  = o.mass
            self.obj['id'][i] = o.nm
            self.obj['id_index'][o.nm] = i
            self.obj['x'][:,i,:] = np.array(o.X) # snapshots of positions
            self.obj['v'][:,i,:] = np.array(o.V)

        X_com, V_com = self.COM()        

        return {'t': self.sim_time, 'id': self.obj['id'], 'id_index': self.obj['id_index'], \
                'x': self.obj['x'], 'v': self.obj['v'], 'x_com': X_com, 'v_com': V_com,
                'epoch': self.epoch, 'n': self.n}

    @deprecated
    def merge(self, obj1, obj2):
        """
        Merge two objects if they get too close together. 
        Momentum is conserved.

        Parameters 
        ----------
        obj1, obj2: <class '__main__.Particle'> object
            Input objects to be merged.

        Returns
        -------
        (obj1, obj2, obj3): tuple, <class '__main__.Particle'> object
            input objects and the new object created from merger.
        """
        name = obj1.nm + '+' + obj2.nm
        x = (obj1.pos[0] + obj2.pos[0]) / 2
        y = (obj1.pos[1] + obj2.pos[1]) / 2
        z = (obj1.pos[2] + obj2.pos[2]) / 2
        mass = obj1.mass + obj2.mass
        vx = (obj1.mass * obj1.vel[0] + obj2.mass * obj2.vel[0]) / mass
        vy = (obj1.mass * obj1.vel[1] + obj2.mass * obj2.vel[1]) / mass
        vz = (obj1.mass * obj1.vel[2] + obj2.mass * obj2.vel[2]) / mass
        obj3 = Particle(mass, x, y, z, vx, vy, vz, name)
        obj3.origin1 = obj1.nm
        obj3.origin2 = obj2.nm
        obj1.merged = True
        obj2.merged = True

        return obj1, obj2, obj3

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
        
class Particle(object):
    """
    Objects interating with gravity only

    Parameters
    ----------
    mass: float
        Mass of object in kg
    x, y, z: float
        Initial position (in m) of object, measured in the
        inertial frame of the central body.
    vx, vy, vz: float
        Initial velocity (in m/s) of object, measured in the
        inertial frame of the central body.
    nm: str
        Name of object

    Note
    ----
    Runs much faster with Python lists than numpy.ndarray
    """
    G  = 6.6743e-11 # [m^3/kg/s^2]

    @deprecated
    def __init__(self, mass, x, y, z, vx, vy, vz, nm):
        self.mass = mass
        self.pos = [x, y, z] # current position 
        self.vel = [vx, vy, vz] # current velocity
        self.X = [(self.pos[0], self.pos[1], self.pos[2])] # history of all positions
        self.V = [(self.vel[0], self.vel[1], self.vel[2])]
        self.nm = nm
        self.merged = False # merge status
        self.origin1 = None # name of progenitor 1
        self.origin2 = None # name of progenitor 2

    # check if mergeable with another object
    # objects are mergeable if they share no common origin
    def mergeable(self, other):
        return self.origin1 != other.nm and self.origin2 != other.nm

    # Updates position and velocity of object
    def update(self, force_x, force_y, force_z, dt):
        self.vel[0] += force_x / self.mass * dt
        self.vel[1] += force_y / self.mass * dt
        self.vel[2] += force_z / self.mass * dt
        # orbits won't close if order swapped
        self.pos[0] += self.vel[0] * dt
        self.pos[1] += self.vel[1] * dt
        self.pos[2] += self.vel[2] * dt

        self.X.append((self.pos[0], self.pos[1], self.pos[2]))
        self.V.append((self.vel[0], self.vel[1], self.vel[2]))

    def distance(self, other):
        dx = other.pos[0] - self.pos[0]
        dy = other.pos[1] - self.pos[1]
        dz = other.pos[2] - self.pos[2]
        r = (dx**2 + dy**2 + dz**2)**0.5
        return r

    def gravity(self, other, n=2):
        """
        Computes force between self and another object.

        Parameters
        ----------
        other: <class '__main__.Particle'> object
            The other object to be considered.
        n: float, optional
            (Negative) exponent of force law, i.e. F ~ 1/r^n

        Returns
        -------
        force: tuple, float
            Force vector for self, pointing towards other.
        """
        dx = other.pos[0] - self.pos[0]
        dy = other.pos[1] - self.pos[1]
        dz = other.pos[2] - self.pos[2]
        r = (dx**2 + dy**2 + dz**2)**0.5

        force = self.G * self.mass * other.mass / r**n
        force_x = force * dx / r
        force_y = force * dy / r
        force_z = force * dz / r

        return force_x, force_y, force_z

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
