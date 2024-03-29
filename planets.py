"""
-----------------------------------------------------------------
Simple solar system simulation, with basic graphic user interface
-----------------------------------------------------------------
"""
import sys, numpy as np, matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from tqdm.auto import trange

class particle(object):
    """
    Class of objects interating with gravity only

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

    n: float, optional
        (Negative) exponent of force law, i.e. F ~ 1/r^n

    Note
    ----
    Runs much faster with Python lists than numpy.ndarray
    """
    def __init__(self, mass, x, y, z, vx, vy, vz, nm, n=2):
        self.mass = mass
        self.pos = [x, y, z] # current position 
        self.vel = [vx, vy, vz] # current velocity
        self.path_x = [x] # history of all positions
        self.path_y = [y]
        self.path_z = [z]
        self.nm = nm
        self.merged = False # merge status
        self.origin1 = None # name of progenitor 1
        self.origin2 = None # name of progenitor 2
        self.n = n # force ~ 1/r^n

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
        self.path_x.append(self.pos[0])
        self.path_y.append(self.pos[1])
        self.path_z.append(self.pos[2])

    def distance(self, other):
        dx = other.pos[0] - self.pos[0]
        dy = other.pos[1] - self.pos[1]
        dz = other.pos[2] - self.pos[2]
        r = (dx**2 + dy**2 + dz**2)**(1/2.)
        return r

    def gravity(self, other):
        """
        Computes force between self and another object.

        Parameters
        ----------
        other: <class '__main__.particle'> object
            The other object to be considered.

        Returns
        -------
        force: tuple, float
            Force vector for self, pointing towards other.
        """
        dx = other.pos[0] - self.pos[0]
        dy = other.pos[1] - self.pos[1]
        dz = other.pos[2] - self.pos[2]
        r = (dx**2 + dy**2 + dz**2)**(1/2)
        force = 6.67e-11 * self.mass * other.mass / r**self.n
        force_x = force * dx / r
        force_y = force * dy / r
        force_z = force * dz / r
        return force_x, force_y, force_z

def merge(obj1, obj2):
    """
    Merge two objects if they get too close together. 
    Momentum is conserved.

    Parameters 
    ----------
    obj1, obj2: <class '__main__.particle'> object
        Input objects to be merged.

    Returns
    -------
    (obj1, obj2, obj3): tuple, <class '__main__.particle'> object
        input objects and the new object created from merger.
    """
    name = obj1.nm + '+' + obj2.nm
    x = (obj1.pos[0] + obj2.pos[0]) / 2
    y = (obj1.pos[1] + obj2.pos[1]) / 2
    z = (obj1.pos[2] + obj2.pos[2]) / 2
    mass = obj1.mass + obj2.mass
    vx = (obj1.mass * obj1.vel[0] + obj2.mass * obj2.vel[0])/mass
    vy = (obj1.mass * obj1.vel[1] + obj2.mass * obj2.vel[1])/mass
    vz = (obj1.mass * obj1.vel[2] + obj2.mass * obj2.vel[2])/mass
    obj3 = particle(mass, x, y, z, vx, vy, vz, name)
    obj3.origin1 = obj1.nm
    obj3.origin2 = obj2.nm
    obj1.merged = True
    obj2.merged = True
    return obj1, obj2, obj3

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
    phi = np.random.uniform(0, 2*np.pi)
    mass = np.random.uniform(0, mmax)
    x = r * np.sin(theta)*np.cos(phi)
    y = r * np.sin(theta)*np.sin(phi)
    z = r * np.cos(theta)
    vel = np.random.uniform(-1, 1, size=3)*vmax
    return mass, x, y, z, vel[0], vel[1], vel[2]

def gen_random(N, r0, mmax, vmax, n=2):
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

    n: float, optional
        (Negative) exponent of force law, i.e. F ~ 1/r^n

    Returns
    -------
    objects: list
        List of randomly generated objects.
    """
    objects = []
    for i in range(N):
        mass, x, y, z, vx, vy, vz = randomize(r0, mmax, vmax)
        random_obj = particle(mass, x, y, z, vx, vy, vz, str(i), n)
        objects.append(random_obj)
    return objects

def simulate(objects, num_steps, sample_rate, rcrit=3.2e8):
    """
    Runs the simulation and updates the objects. First compute
    forces acting on all objects, then update the position and
    velocitiy of each object.

    Parameters
    ----------
    objects: list
        List of all objects in the simulation.

    num_steps: int
        Total number of steps in the simulation.

    sample_rate: float
        Number of steps per day.

    rcrit: float
        Optional. Distance threshold for merger [m]

    Returns
    -------
    objects: list
        List of time evolved objects.

    old_objects: list
        List of merged objects.
    """
    old_objects = []
    for step in trange(num_steps):
        forces_x = []
        forces_y = []
        forces_z = []
        for obj1 in objects:
            total_fx = 0
            total_fy = 0
            total_fz = 0
            for obj2 in objects:
                if obj1 != obj2 and obj2.mergeable(obj1):
                    # Merging works only when sample_rate > 100
                    if sample_rate >= 100 and obj1.distance(obj2) < rcrit:
                        objects.remove(obj1)
                        objects.remove(obj2)
                        obj1, obj2, obj3 = merge(obj1, obj2)
                        old_objects.append(obj1)
                        old_objects.append(obj2)
                        objects.append(obj3)
                    else: # if the two are distinct
                        force = obj1.gravity(obj2)
                        total_fx += force[0]
                        total_fy += force[1]
                        total_fz += force[2]
            forces_x.append(total_fx)
            forces_y.append(total_fy)
            forces_z.append(total_fz)
        for i in range(len(forces_x)): # update all objects
            objects[i].update(forces_x[i], forces_y[i], forces_z[i], 86400/sample_rate)
            
    return objects, old_objects

# Animate the orbits; doesn't work with merge
def animate(objects, num_steps, sample_rate, epoch, set_range,
            hide_trace=False, save_file=False, radius=40):
    fig = plt.figure()
    ax = p3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    positions = [] # list of paths; the indices are: [object][direction (xyz),step]
    for obj in objects:
        positions.append(np.array([obj.path_x, obj.path_y, obj.path_z]))
    if not hide_trace:
        paths = [ax.plot(i[0,0:1], i[1,0:1], i[2,0:1])[0] for i in positions]
    else:
        points = [ax.plot([], [], [], '.')[0] for i in positions]
        #points, = ax.plot([], [], [], '.') # plot the points together
    if set_range:
        r = 149597870700*radius
        ax.set_xlim3d([-r,r])
        ax.set_ylim3d([-r,r])
        ax.set_zlim3d([-r,r])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    title = ax.text2D(0.5, 0.95, '', horizontalalignment='center',
                      fontsize=14, transform=ax.transAxes)
    if not hide_trace:
        anim = animation.FuncAnimation(fig, update_paths, num_steps,
                fargs=(positions, paths, title, sample_rate, epoch), interval=1, blit=False)
    else:
        anim = animation.FuncAnimation(fig, no_trace_update, num_steps,
                fargs=(positions, points, title, sample_rate, epoch), interval=1, blit=False)
    if save_file:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=12, metadata=dict(artist='Me'), bitrate=1800)
        anim.save('results/paths.mp4', writer=writer)#, fps=15, extra_args=['-vcodec', 'libx264'])
    plt.show()
    return anim

# Update paths for animation
def update_paths(num_steps, positions, paths, title, sample_rate, epoch):
    for path, position in zip(paths, positions): # for each object
        # Use these to keep temporary trails
        #if num_steps >= 5:
        #    path.set_data(position[0:2, num_steps-4:num_steps])
        #    path.set_3d_properties(position[2, num_steps-4:num_steps])
        #    title.set_text('JD %.1f'%(epoch+num_steps/sample_rate))
        #else:
        path.set_data(position[0:2, :num_steps])
        path.set_3d_properties(position[2, :num_steps])
        title.set_text('JD %.1f'%(epoch+num_steps/sample_rate))
    return paths, title

# Although in FuncAnimation the variable 'num_steps' is provided, only the
# individual iterations of range(num_steps) are fed into 'step'
def no_trace_update(step, positions, points, title, sample_rate, epoch):
    for position, point in zip(positions, points):
        num_steps = position.shape[1] # needs to be looped back manually
        step1 = step%num_steps
        point.set_data(position[0,step1:step1+1], position[1,step1:step1+1])
        point.set_3d_properties(position[2,step1:step1+1])
    title.set_text('JD %.1f'%(epoch+step%num_steps/sample_rate))
    # Use these if you want to plot the points together
    #positions = np.array(positions)
    #points.set_data(positions[:,0,step], positions[:,1,step])
    #points.set_3d_properties(positions[:,2,step])
    return points

# Plot orbits
def plot(objects, old_objects, num_steps, sample_rate, set_range, legend, epoch, radius=40):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for obj in objects:
        ax.plot(obj.path_x, obj.path_y, obj.path_z, label=obj.nm)
    for obj in old_objects:
        ax.plot(obj.path_x, obj.path_y, obj.path_z, label=obj.nm)
    if set_range: # Set range for plot
        r = 149597870700*radius
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)
        ax.set_zlim(-r, r)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('JD %.1f'%(epoch+num_steps/sample_rate))
    if legend:
        plt.legend()
    plt.savefig('results/orbits.pdf', bbox_inches='tight')
    plt.show()
    table_items(objects)

def kep2cart(e, a, i, Ω, ϖ, L, m1, m2=1.989e30):
    """
    Convert Keplerian orbital elements to state vectors in 
    heliocentric cartesian coordinates in the J2000 ecliptic 
    plane. The x-axis points toward the vernal equinox.

    Parameters
    ----------
    e: float
        Eccentricity [rad]

    a: float
        Semi-major axis [au]

    i: float
        Inclination [deg]

    Ω: float
        Longitude of the ascending node [deg]

    ϖ: float
        Longitude of periapsis [deg]

    L: float
        Mean longitude [deg] at given epoch

    m1: float
        Mass of orbiting object [kg]

    m2: float, optional
        Mass of central body, set to 1 solar mass by default

    Returns
    -------
    x, y, z, vx, vy, vz: tuple
        Position [m] and velocity [m/s] measured in the 
        inertial frame of the central body.

    Note
    ----
    e and a define the shape and size of the ellipse. i and 
    Ω define the orientation of the orbital plane. ϖ defines
    the orientation of the ellipse in the orbital plane and 
    L defines the position of the orbiting body along the
    ellipse at a specific time. Alternatively, one can use q
    (periapsis) in place of a, and M (mean anomaly) or T (time
    of perihelion passage) in place of L.

    For planets, the elements used are e, a, i, Ω, ϖ, L.
    """
    ω = (ϖ - Ω)*np.pi/180 # argument of periapsis [rad]
    i = i*np.pi/180
    Ω = Ω*np.pi/180
    M = L - ϖ # mean anomaly [deg]
    E = ecc_anomaly(e, M)*np.pi/180 # eccentric anomaly
    # heliocentric coordinates in orbital plane, z1=0
    x0 = a*(np.cos(E)-e)
    y0 = a*(1-e**2)**(1/2)*np.sin(E)
    # ecliptic coordinates, in au
    x = (np.cos(ω)*np.cos(Ω)-np.sin(ω)*np.sin(Ω)*np.cos(i))*x0+(-np.sin(ω)*np.cos(Ω)-np.cos(ω)*np.sin(Ω)*np.cos(i))*y0
    y = (np.cos(ω)*np.sin(Ω)+np.sin(ω)*np.cos(Ω)*np.cos(i))*x0+(-np.sin(ω)*np.sin(Ω)+np.cos(ω)*np.cos(Ω)*np.cos(i))*y0
    z = np.sin(ω)*np.sin(i)*x0+np.cos(ω)*np.sin(i)*y0
    # find velocity in m/s
    au = 149597870700
    mu = 6.6743e-11*(m1+m2)
    nu = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # true anomaly
    v = np.sqrt(mu/a/au*(1+e*np.cos(E))/(1-e*np.cos(E))) # orbital speed
    theta = np.arcsin(np.sqrt((1-e**2)/(1-e**2*np.cos(E)**2))) # angle between r and v
    vx0, vy0 = v*np.cos(nu+theta), v*np.sin(nu+theta) # velocity in orbital plane
    vx = (np.cos(ω)*np.cos(Ω)-np.sin(ω)*np.sin(Ω)*np.cos(i))*vx0+(-np.sin(ω)*np.cos(Ω)-np.cos(ω)*np.sin(Ω)*np.cos(i))*vy0
    vy = (np.cos(ω)*np.sin(Ω)+np.sin(ω)*np.cos(Ω)*np.cos(i))*vx0+(-np.sin(ω)*np.sin(Ω)+np.cos(ω)*np.cos(Ω)*np.cos(i))*vy0
    vz = np.sin(ω)*np.sin(i)*vx0+np.cos(ω)*np.sin(i)*vy0
    return x*au, y*au, z*au, vx, vy, vz

# Solve for eccentric anomaly [deg] given e [rad] and M [deg]
def ecc_anomaly(e, M, tolerance=0.01):
    e_deg = e*180/np.pi
    E = M+e_deg*np.sin(M*np.pi/180)
    dE = 1
    while dE > tolerance:
        dM = M-(E-e_deg*np.sin(E*np.pi/180))
        dE = dM/(1-e*np.cos(E*np.pi/180))
        E += dE
    return E # [deg]

# Converts state vectors (r in au and v in m/s) into orbital elements
def cart2kep(x, y, z, vx, vy, vz, m1, m2=1.989e30):
    au = 149597870700
    mu = 6.6743e-11*(m1+m2)
    x *= au; y *= au; z *= au
    r = (x**2+y**2+z**2)**(1/2)
    v = (vx**2+vy**2+vz**2)**(1/2)
    eps = v**2/2-mu/r # specific orbital energy
    a = -mu/2/eps
    theta = np.arccos((x*vx+y*vy+z*vz)/r/v) # angle between r and v
    h = r*v*np.sin(theta) # specific angular momentum
    hx, hy, hz = y*vz-z*vy, z*vx-x*vz, x*vy-y*vx
    e = (1-h**2/mu/a)**(1/2)
    E = np.arccos((1-r/a)/e) # eccentric anomaly
    M = E - e*np.cos(E) # mean anomaly
    i = np.arccos(hz/h)
    Om = (np.arctan2(hx, -hy)+2*np.pi)%(2*np.pi)
    nu = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # true anomaly
    w = np.arctan2(z/np.sin(i), x*np.cos(Om)+y*np.sin(Om))-nu
    wb = (w+Om+2*np.pi)%(2*np.pi)
    L = (M+wb+2*np.pi)%(2*np.pi)
    return e, a/au, i*180/np.pi, Om*180/np.pi, wb*180/np.pi, L*180/np.pi

# Keplerian elements (e, a, i, Ω, ϖ, L) and masses of solar system objects.
# J2000, http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
# recreated the Aug 19, 2021 (JD 2459446) Venus-Mars-Mercury conjuction.
# Halley: Epoch 2449400.5 (1994-Feb-17.0), https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=1P
# Hale-Bopp: Epoch 2454724.5 (2008-Sep-15.0)
# must subtract 2144.5 from the JD number in the simulation for Halley, or add 3179.5 for Hale-Bopp
orbital_elements = {'Mercury':[0.20563069, 0.38709893, 7.00487, 48.33167, 77.45645, 252.25084, 3.301e23],
                    'Venus':[0.00677323, 0.72333199, 3.39471, 76.68069, 131.53298, 181.97973, 4.867e24],
                    'Earth':[0.01671022, 1.00000011, 0.00005, -11.26064, 102.94719, 100.46435, 5.972e24],
                    'Mars':[0.09341233, 1.52366231, 1.85061, 49.57854, 336.04084, 355.45332, 6.417e23],
                    'Jupiter':[0.04839266, 5.20336301, 1.30530, 100.55615, 14.75385, 34.40438, 1.898e27],
                    'Saturn':[0.05415060, 9.53707032, 2.48446, 113.71504, 92.43194, 49.94432, 5.683e26],
                    'Uranus':[0.04716771, 19.19126393, 0.76986, 74.22988, 170.96424, 313.23218, 8.681e25],
                    'Neptune':[0.00858587, 30.06896348, 1.76917, 131.72169, 44.97135, 304.88003, 1.024e26],
                    'Pluto':[0.24880766, 39.48168677, 17.14175, 110.30347, 224.06676, 238.92881, 1.303e22],
                    'Halley':[0.96714, 17.834, 162.26, 58.42, 169.75, 208.13, 2.2e14],
                    'Hale-Bopp':[0.99496070, 182.05197034, 89.21708989, 282.94875394, 53.61077446, 55.29041624, 1.3e16]
}
# comet1 = particle(3e14, -4e9, 5e12, 0, 100, 0, 0) # escape
# comet1 = particle(3e14, -4e9, 5e12, 0, 1000, 1000, 0) # capture
# comet1 = particle(3e14, -4e9, 5e12, 2e12, 500, -1000, 2000, 'Some comet')
# Some more presets: objects for demo
test_particles = {'Test Particle 1':[2e30, 1e12, 0, 0, 0, 1e4, 0],  # 1, 2 - fork merge
                  'Test Particle 2':[2e30, -1e12, 0, 0, 0, 1e4, 0],
                  'Test Particle 3':[2e30, 1e12, 0, 0, 0, 1e3, 0],  # 3, 4 - spiral/solenoid
                  'Test Particle 4':[2e30, 1e12, 0, 1e12, -1e3, -1e3, -1e3] 
}

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

objects = []

# Create application and window
app = QApplication(sys.argv)
w = QWidget()
w.resize(640, 480)
w.setWindowTitle('Settings')

# Checkbox for setting range limits to solar system
checkbox1 = QCheckBox(w); checkbox1.move(250, 290)
checkbox1.setChecked(True); checkbox1.setText('Set range at          AU')
textbox0 = QLineEdit(w); textbox0.move(350, 290)
textbox0.resize(25, 20); textbox0.setText('40')
# Checkbox for showing legend
checkbox2 = QCheckBox(w); checkbox2.move(250, 310)
checkbox2.setChecked(False); checkbox2.setText('Show Legend')
# Checkbox for animation
checkbox3 = QCheckBox(w); checkbox3.move(250, 330)
checkbox3.setChecked(True); checkbox3.setText('Animated')
# Checkbox for showing no trace of orbit
checkbox4 = QCheckBox(w); checkbox4.move(250, 350)
checkbox4.setChecked(False); checkbox4.setText('Hide orbits')
# Checkbox for saving animation to file
checkbox5 = QCheckBox(w); checkbox5.move(250, 370)
checkbox5.setChecked(False); checkbox5.setText('Save animation to file')

# Textboxes for entering parameters
label6 = QLabel(w); label6.setText('X (m)'); label6.move(20, 30)
textbox1 = QLineEdit(w); textbox1.move(100, 20)
textbox1.resize(120, 40); textbox1.setText('0')
label7 = QLabel(w); label7.setText('Y (m)'); label7.move(20, 75)
textbox2 = QLineEdit(w); textbox2.move(100, 65)
textbox2.resize(120, 40); textbox2.setText('0')
label8 = QLabel(w); label8.setText('Z (m)'); label8.move(20, 120)
textbox3 = QLineEdit(w); textbox3.move(100, 110)
textbox3.resize(120, 40); textbox3.setText('0')
label9 = QLabel(w); label9.setText('Vx (m/s)'); label9.move(20, 165)
textbox4 = QLineEdit(w); textbox4.move(100, 155)
textbox4.resize(120, 40); textbox4.setText('1e7')
label10 = QLabel(w); label10.setText('Vy (m/s)'); label10.move(20, 210)
textbox5 = QLineEdit(w); textbox5.move(100, 200)
textbox5.resize(120, 40); textbox5.setText('1e7')
label11 = QLabel(w); label11.setText('Vz (m/s)'); label11.move(20, 255)
textbox6 = QLineEdit(w); textbox6.move(100, 245)
textbox6.resize(120, 40); textbox6.setText('1e7')
label12 = QLabel(w); label12.setText('M (kg)'); label12.move(20, 300)
textbox7 = QLineEdit(w); textbox7.move(100, 290)
textbox7.resize(120, 40); textbox7.setText('1e27')
label13 = QLabel(w); label13.setText('Name'); label13.move(20, 345)
textbox8 = QLineEdit(w); textbox8.move(100, 335)
textbox8.resize(120, 40); textbox8.setText('Particle 1')

# Textboxes for entering parameters for random simulation
label14 = QLabel(w); label14.setText('Number'); label14.move(450, 265)
textbox9 = QLineEdit(w); textbox9.move(530, 260)
textbox9.resize(80, 30); textbox9.setText('12')
label15 = QLabel(w); label15.setText('Range (m)'); label15.move(450, 305)
textbox10 = QLineEdit(w); textbox10.move(530, 300)
textbox10.resize(80, 30); textbox10.setText('1e11')
label16 = QLabel(w); label16.setText('Vmax (m/s)'); label16.move(450, 345)
textbox11 = QLineEdit(w); textbox11.move(530, 340)
textbox11.resize(80, 30); textbox11.setText('1e3')
label17 = QLabel(w); label17.setText('Mmax (kg)'); label17.move(450, 385)
textbox12 = QLineEdit(w); textbox12.move(530, 380)
textbox12.resize(80, 30); textbox12.setText('1e27')

# Entering epoch; default: J2000
label18 = QLabel(w); label18.setText('Epoch: JD'); label18.move(20, 420)
textbox16 = QLineEdit(w); textbox16.move(100, 410)
textbox16.resize(120, 40); textbox16.setText('2451545.0')

# Entering number of steps per day
label4 = QLabel(w); label4.setText('Sample rate'); label4.move(250, 25)
textbox13 = QLineEdit(w); textbox13.move(360, 20)
textbox13.resize(50, 30); textbox13.setText('1')

# Entering number of years
label5 = QLabel(w); label5.setText('Duration (years)'); label5.move(250, 65)
textbox14 = QLineEdit(w); textbox14.move(360, 60)
textbox14.resize(50, 30); textbox14.setText('10')

# Entering force law
label2 = QLabel(w); label2.setText('Force index'); label2.move(250, 105)
textbox15 = QLineEdit(w); textbox15.move(360, 100)
textbox15.resize(50, 30); textbox15.setText('2')

# Drop-down list for selecting solar system objects to fill textboxes 1-8
label1 = QLabel(w); label1.setText('Presets'); label1.move(250, 140)
combobox1 = QComboBox(w); combobox1.move(250, 160)
combobox1.addItem('Sun'); combobox1.addItem('Mercury')
combobox1.addItem('Venus'); combobox1.addItem('Earth')
combobox1.addItem('Mars'); combobox1.addItem('Jupiter')
combobox1.addItem('Saturn'); combobox1.addItem('Uranus')
combobox1.addItem('Neptune'); combobox1.addItem('Pluto')
combobox1.addItem('Halley'); combobox1.addItem('Hale-Bopp')
combobox1.addItem('Random object')
def on_activated1(text):
    text_str = str(text)
    if text_str != 'Random object':
        if text_str != 'Sun':
            e = orbital_elements[text_str][0]
            a = orbital_elements[text_str][1]
            i = orbital_elements[text_str][2]
            Om = orbital_elements[text_str][3]
            wb = orbital_elements[text_str][4]
            L = orbital_elements[text_str][5]
            m = orbital_elements[text_str][6]
            x,y,z,vx,vy,vz = kep2cart(e,a,i,Om,wb,L,m)
        else:
            x = y = z = vx = vy = vz = 0
            m = 1.989e30
        textbox1.setText(str(x))
        textbox2.setText(str(y))
        textbox3.setText(str(z))
        textbox4.setText(str(vx))
        textbox5.setText(str(vy))
        textbox6.setText(str(vz))
        textbox7.setText(str(m))
        textbox8.setText(text)
    else:
        data = randomize(1e11, 1e27, 1e3)
        textbox1.setText(str(data[1]))
        textbox2.setText(str(data[2]))
        textbox3.setText(str(data[3]))
        textbox4.setText(str(data[4]))
        textbox5.setText(str(data[5]))
        textbox6.setText(str(data[6]))
        textbox7.setText(str(data[0]))
        textbox8.setText(text)

# Drop-down list for selecting test particles
combobox2 = QComboBox(w); combobox2.move(250, 185)
combobox2.addItem('Test Particle 1')
combobox2.addItem('Test Particle 2')
combobox2.addItem('Test Particle 3')
combobox2.addItem('Test Particle 4')
def on_activated2(text):
    text_str = str(text)
    textbox1.setText(str(test_particles[text_str][1]))
    textbox2.setText(str(test_particles[text_str][2]))
    textbox3.setText(str(test_particles[text_str][3]))
    textbox4.setText(str(test_particles[text_str][4]))
    textbox5.setText(str(test_particles[text_str][5]))
    textbox6.setText(str(test_particles[text_str][6]))
    textbox7.setText(str(test_particles[text_str][0]))
    textbox8.setText(text)

# Table of all particles created for the simulation
label3 = QLabel(w)
label3.setText('List of Objects')
label3.move(485, 25)
table = QTableWidget(w)
table.resize(160, 200)
table.move(450, 50)
table_item = QTableWidgetItem()
def table_items(objects): # Update table of objects
    num_items = len(objects)
    table.setColumnCount(1)
    table.setRowCount(num_items)
    count = 0
    for obj in objects:
        Name = obj.nm
        table.setItem(count, 0, QTableWidgetItem(Name))
        count += 1

# Button for adding new particle
button1 = QPushButton('Add Particle', w)
button1.setToolTip('Add the particle to list')
button1.resize(button1.sizeHint()); button1.move(250, 220)
def on_click_button1():
    try:
        x = float(textbox1.text())
        y = float(textbox2.text())
        z = float(textbox3.text())
        vx = float(textbox4.text())
        vy = float(textbox5.text())
        vz = float(textbox6.text())
        mass = float(textbox7.text())
        name = str(textbox8.text())
        n = float(textbox15.text())
    except ValueError:
        print("Error: invalid parameters.")
        return 1
    new_obj = particle(mass, x, y, z, vx, vy, vz, name, n)
    objects.append(new_obj)
    table_items(objects) # Updates table

# Button for starting the simulation
button2 = QPushButton('Start', w)
button2.setToolTip('Start simulation')
button2.resize(button2.sizeHint()); button2.move(250, 420)
def on_click_button2():
    try:
        sample_rate = int(textbox13.text())
        num_years = int(textbox14.text())
        epoch = float(textbox16.text())
    except ValueError: # If these not entered, set default values
        sample_rate = 1
        num_years = 50
        epoch = 2451545.0
        print("Warning: number of steps not entered; set to default values.")
    num_steps = int(365.25*sample_rate*num_years)
    global objects
    if len(objects) == 0:
        print("Error: no objects entered.")
        return 1
    objects, old_objects = simulate(objects, num_steps, sample_rate)
    plot_size = int(textbox0.text()) # box size in au
    if checkbox3.isChecked():
        animate(objects, num_steps, sample_rate, epoch, checkbox1.isChecked(), 
                checkbox4.isChecked(), checkbox5.isChecked(), plot_size)
    else:
        plot(objects, old_objects, num_steps, sample_rate, 
             checkbox1.isChecked(), checkbox2.isChecked(), epoch, plot_size)

# Button for running simulation with all random particles
button3 = QPushButton('Start Random', w)
button3.setToolTip('Start simulation with all random particles')
button3.resize(button3.sizeHint()); button3.move(470, 420)
def on_click_button3():
    try:
        N = int(textbox9.text())
        r0 = float(textbox10.text())
        vmax = float(textbox11.text())
        mmax = float(textbox12.text())
        n = float(textbox15.text())
    except ValueError:
        N = 12
        r0 = 6e12
        vmax = 5e4
        mmax = 5e30
        n = 2
        print("Warning: initial parameters not entered; set to default values.")
    try:
        sample_rate = int(textbox13.text())
        num_years = int(textbox14.text())
    except ValueError: # If these not entered, set default values
        sample_rate = 1
        num_years = 50
        print("Warning: number of steps not entered; set to default values.")
    epoch = 0 # No need to use J2000 for random objects
    num_steps = int(365.25*sample_rate*num_years)
    objects = gen_random(N, r0, mmax, vmax)
    table_items(objects)
    objects, old_objects = simulate(objects, num_steps, sample_rate)
    plot_size = int(textbox0.text()) # box size in au
    if checkbox3.isChecked():
        animate(objects, num_steps, sample_rate, epoch, checkbox1.isChecked(), 
                checkbox4.isChecked(), checkbox5.isChecked(), plot_size)
    else:
        plot(objects, old_objects, num_steps, sample_rate, 
             checkbox1.isChecked(), checkbox2.isChecked(), epoch, plot_size)

# Button for clearing old list of objects
button4 = QPushButton('Reset', w)
button4.setToolTip('Clear old list of objects')
button4.resize(button4.sizeHint()); button4.move(320, 420)
def on_click_button4():
    global objects
    objects = []
    textbox1.setText('0')
    textbox2.setText('0')
    textbox3.setText('0')
    textbox4.setText('1e7')
    textbox5.setText('1e7')
    textbox6.setText('1e7')
    textbox7.setText('1e27')
    textbox8.setText('Particle 1')
    textbox9.setText('12')
    textbox10.setText('1e11')
    textbox11.setText('1e3')
    textbox12.setText('1e27')
    textbox13.setText('1')
    textbox14.setText('10')
    textbox15.setText('2')
    textbox16.setText('2451545.0')
    table.setColumnCount(0)
    table.setRowCount(0)
    plt.close()

button5 = QPushButton('Quick Add', w)
button5.setToolTip('Add all solar system objects')
button5.resize(button5.sizeHint()); button5.move(250, 250)
def on_click_button5():
    on_activated1('Sun')
    on_click_button1()
    for i in orbital_elements.keys():
        on_activated1(i)
        on_click_button1()

# Connect the signals to the slots
button1.clicked.connect(on_click_button1)
button2.clicked.connect(on_click_button2)
button3.clicked.connect(on_click_button3)
button4.clicked.connect(on_click_button4)
button5.clicked.connect(on_click_button5)
combobox1.activated[str].connect(on_activated1)
combobox2.activated[str].connect(on_activated2)
w.show()

sys.exit(app.exec_())
