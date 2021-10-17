"""
-----------------------------------------------------------------
Simple solar system simulation, with basic graphic user interface
-----------------------------------------------------------------
It's probably more efficient to run the simulation by storing all
position and velocity data in two arrays, instead of creating a 
separate class for all the objects.
"""
import sys, numpy as np, matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

class particle(object):
    """
    Class of objects interating with gravity only

    Parameters
    ----------
    mass: float
        Mass of object in kg

    pos: numpy.ndarray, float
        Initial position vector of object, in m

    vel: numpy.ndarray, float
        Initial velocity vector of object, in m/s

    nm: str
        Name of object
    """
    def __init__(self, mass, pos, vel, nm):
        self.mass = mass
        self.pos = pos # current position 
        self.vel = vel # current velocity
        self.path = self.pos # history of all positions
        self.nm = nm
        self.merged = False # merge status
        self.origin1 = None # name of progenitor 1
        self.origin2 = None # name of progenitor 2

    # check if mergeable with another object
    # objects are mergeable if they share no common origin
    def mergeable(self, other):
        return self.origin1 != other.nm and self.origin2 != other.nm

    # Updates position and velocity of object
    def update(self, force, dt):
        acc = force/self.mass
        self.vel += acc * dt
        self.pos += self.vel * dt # orbits won't close if order swapped
        self.path = np.vstack((self.path, self.pos))

    def gravity(self, other):
        """
        Computes force between self and another object.

        Parameters
        ----------
        other: <class '__main__.particle'> object
            The other object to be considered.

        Returns
        -------
        force: numpy.ndarray, float
            Force vector for self, pointing towards other.
        """
        r = other.pos - self.pos # displacement
        rnorm = np.sqrt(np.sum(r**2))
        force = 6.67e-11*self.mass*other.mass*r/rnorm**3
        return force

def merge(obj1, obj2):
    """
    Merge two objects if they get too close together. Momentum is conserved.

    Parameters 
    ----------
    obj1, obj2: <class '__main__.particle'> object
        Input objects to be merged.

    Returns
    -------
    obj3: <class '__main__.particle'> object
        New object created from merger.
    """
    name = obj1.nm + '+' + obj2.nm
    pos = (obj1.pos + obj2.pos) / 2
    mass = obj1.mass + obj2.mass
    vel = (obj1.mass * obj1.vel + obj2.mass * obj2.vel)/mass 
    obj3 = particle(mass, pos, vel, name)
    obj3.origin1 = obj1.nm
    obj3.origin2 = obj2.nm
    obj1.merged = True
    return obj1, obj2, obj3

def randomize(r0, mmax, vmax):
    """
    Generates randomized parameters for creating an object.
    Values are drawn from a uniform distribution.

    Parameters
    ----------
    r0: float
        Radius of the sphere within which the object is created.

    mmax: float
        Maximum mass of the object.

    vmax: float
        Maximum speed in each direction.

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
    px = r * np.sin(theta)*np.cos(phi)
    py = r * np.sin(theta)*np.sin(phi)
    pz = r * np.cos(theta)
    pos = np.array([px, py, pz])
    vel = np.random.uniform(-1, 1, size=3)*vmax
    return mass, pos, vel

def gen_random(N, r0, mmax, vmax):
    """
    Generates a list of N random objects.

    Parameters
    ----------
    N: int
        Number of objects to generate.

    r0: float
        Radius of the sphere within which the objects are created.

    mmax: float
        Maximum mass of the objects.

    vmax: float
        Maximum speed in each direction.

    Returns
    -------
    objects: list
        List of randomly generated objects.
    """
    objects = []
    for i in range(N):
        mass, pos, vel = randomize(r0, mmax, vmax)
        random_obj = particle(mass, pos, vel, str(i))
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

    num_steps: float
        Total number of steps in the simulation.

    sample_rate: float
        Number of steps per day.

    rcrit: float
        Optional. Distance threshold for merger, in m.

    Returns
    -------
    objects: list
        List of time evolved objects.

    old_objects: list
        List of merged objects.
    """
    old_objects = []
    for step in range(int(num_steps)):
        forces = np.zeros(3) # stores forces acting on all objects
        for obj1 in objects:
            Fsum = np.zeros(3) # net force acting on object 1
            for obj2 in objects:
                if obj1 != obj2 and obj2.mergeable(obj1):
                    r = np.sqrt(np.sum((obj1.pos-obj2.pos)**2))
                    # Merging works only when sample_rate > 100
                    if sample_rate > 100 and r < rcrit:
                        objects.remove(obj1)
                        objects.remove(obj2)
                        obj1, obj2, obj3 = merge(obj1, obj2)
                        old_objects.append(obj1)
                        old_objects.append(obj2)
                        objects.append(obj3)
                    else: # if the two are distinct
                        Fsum += obj1.gravity(obj2)
            forces = np.vstack((forces, Fsum))
        forces = forces[1:] # remove the extra row
        for i in range(forces.shape[0]): # update all objects
            objects[i].update(forces[i], 86400/sample_rate)
        print('%.2f'%(step/num_steps*100), '%') # print progress
    return objects, old_objects

# Animate the orbits; doesn't work with merge
def animate(objects, num_steps, save_file=False):
    fig = plt.figure()
    ax = p3.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    positions = [] # list of paths
    for obj in objects:
        positions.append(obj.path)
    paths = [ax.plot(i[0:1,0], i[0:1,1], i[0:1,2])[0] for i in positions]
    ax.set_xlim3d([-6e12,6e12]); ax.set_xlabel('X')
    ax.set_ylim3d([-6e12,6e12]); ax.set_ylabel('Y')
    ax.set_zlim3d([-6e12,6e12]); ax.set_zlabel('Z')
    ax.set_title('Animated Simulation')
    anim = animation.FuncAnimation(fig, update_paths, int(num_steps), 
                fargs=(positions, paths), interval=1, blit=False)
    if save_file:
        anim.save('results/paths.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
    plt.show()
    return 0

# Update paths for animation
def update_paths(num_steps, positions, paths):
    for path, position in zip(paths, positions):
        path.set_data(np.transpose(position[:num_steps, 0:2]))
        path.set_3d_properties(position[:num_steps, 2])
    return paths

# Plot orbits
def plot(objects, old_objects, solar_range, legend):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for obj in objects:
        ax.plot(obj.path[:,0], obj.path[:,1], obj.path[:,2], label=obj.nm)
    for obj in old_objects:
        ax.plot(obj.path[:,0], obj.path[:,1], obj.path[:,2], label=obj.nm)
    if checkbox1.isChecked(): # Set range for solar system
        ax.set_xlim(-6e12, 6e12)
        ax.set_ylim(-6e12, 6e12)
        ax.set_zlim(-6e12, 6e12)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    if checkbox2.isChecked():
        plt.legend()
    plt.savefig('results/orbits.pdf', bbox_inches='tight')
    plt.show()
    table_items(objects)

# Data for solar system objects
solar_system = {'Sun':[2e30, 0, 0, 0, 0, 0, 0],
                'Mercury':[3.3e23, 4.6e10, 0, 0, 0, 5.66e4, 0], # At perihelion
                'Venus':[4.87e24, 1.08e11, 0, 0, 0, 3.5e4, 0],
                'Earth':[6e24, 1.5e11, 0, 0, 0, 3e4, 0],
                'Mars':[6.4e23, 2.3e11, 0, 0, 0, 2.4e4, 0],
                'Jupiter':[1.9e27, 7.8e11, 0, 0, 0, 13070, 0],
                'Saturn':[5.7e26, 1.4e12, 0, 0, 0, 9690, 0],
                'Uranus':[8.7e25, 2.9e12, 0, 0, 0, 6800, 0],
                'Neptune':[1e26, 4.5e12, 0, 0, 0, 5430, 0],
                'Pluto':[1.3e22, 4.24e12, 0, 1.3e12, 0, 6100, 0], # At perihelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
                'Halley':[2.2e14, -5e12, 0, 1.6e12, 0, 550, 0], # At aphelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/cometfact.html
                'Hale-Bopp':[1.3e16, 4.1e11, 4.1e11, 5.54e13, -77, 77, 0]
}
# comet1 = particle(3e14, np.array([-4e9, 5e12, 0]), np.array([100, 0, 0])) # escape
# comet1 = particle(3e14, np.array([-4e9, 5e12, 0,]) np.array([1000, 1000, 0])) # capture
# comet1 = particle(3e14, np.array([-4e9, 5e12, 2e12]), np.array([500, -1000, 2000]), 'Some comet')
# objects for demonstration
test_particles = {'Test Particle 1':[2e30, 1e12, 0, 0, 0, 1e4, 0],  # 1, 2 - fork
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
w.setWindowTitle('Particle Control')

# Checkbox for setting range limits to solar system
checkbox1 = QCheckBox(w)
checkbox1.move(250, 155)
checkbox1.setChecked(False)
checkbox1.setText('Set range for Solar System')
# Checkbox for showing legend
checkbox2 = QCheckBox(w)
checkbox2.move(250, 320)
checkbox2.setChecked(False)
checkbox2.setText('Show Legend')
# Checkbox for animation
checkbox3 = QCheckBox(w)
checkbox3.move(250, 340)
checkbox3.setChecked(False)
checkbox3.setText('Animated')
# Checkbox for saving animation to file
checkbox4 = QCheckBox(w)
checkbox4.move(250, 360)
checkbox4.setChecked(False)
checkbox4.setText('Save animation to file')

# Textboxes for entering parameters
textbox1 = QLineEdit(w); textbox1.move(20, 20); textbox1.resize(200, 40)
textbox1.setText('Enter Position (x)')
textbox2 = QLineEdit(w); textbox2.move(20, 65); textbox2.resize(200, 40)
textbox2.setText('Enter Position (y)')
textbox3 = QLineEdit(w); textbox3.move(20, 110); textbox3.resize(200, 40)
textbox3.setText('Enter Position (z)')
textbox4 = QLineEdit(w); textbox4.move(20, 155); textbox4.resize(200, 40)
textbox4.setText('Enter Velocity (x)')
textbox5 = QLineEdit(w); textbox5.move(20, 200); textbox5.resize(200, 40)
textbox5.setText('Enter Velocity (y)')
textbox6 = QLineEdit(w); textbox6.move(20, 245); textbox6.resize(200, 40)
textbox6.setText('Enter Velocity (z)')
textbox7 = QLineEdit(w); textbox7.move(20, 290); textbox7.resize(200, 40)
textbox7.setText('Enter Mass')
textbox8 = QLineEdit(w); textbox8.move(20, 335); textbox8.resize(200, 40)
textbox8.setText('Enter Name')

# Textboxes for entering parameters for random simulation
textbox9 = QLineEdit(w)
textbox9.move(450, 260)
textbox9.resize(160, 30)
textbox9.setText('Number of Particles')
textbox10 = QLineEdit(w)
textbox10.move(450, 300)
textbox10.resize(160, 30)
textbox10.setText('Initial Radius')
textbox11 = QLineEdit(w)
textbox11.move(450, 340)
textbox11.resize(160, 30)
textbox11.setText('Initial Max Velocity')
textbox12 = QLineEdit(w)
textbox12.move(450, 380)
textbox12.resize(160, 30)
textbox12.setText('Initial Max Mass')

# Entering number of steps per day
textbox13 = QLineEdit(w); textbox13.move(250, 20)
textbox13.resize(150, 30)
textbox13.setText('Number of steps per day')

# Entering number of years
textbox14 = QLineEdit(w)
textbox14.move(250, 60)
textbox14.resize(150, 30)
textbox14.setText('Number of years')

# Drop-down list for selecting solar system objects to fill textboxes 1-8
label1 = QLabel(w)
label1.setText('Solar System')
label1.move(250, 105)
combobox1 = QComboBox(w)
combobox1.addItem('Sun')
combobox1.addItem('Mercury')
combobox1.addItem('Venus')
combobox1.addItem('Earth')
combobox1.addItem('Mars')
combobox1.addItem('Jupiter')
combobox1.addItem('Saturn')
combobox1.addItem('Uranus')
combobox1.addItem('Neptune')
combobox1.addItem('Pluto')
combobox1.addItem('Halley')
combobox1.addItem('Hale-Bopp')
combobox1.addItem('Random object')
combobox1.move(250, 125)
def on_activated1(text):
    text_str = str(text)
    if text_str != 'Random object':
        textbox1.setText(str(solar_system[text_str][1]))
        textbox2.setText(str(solar_system[text_str][2]))
        textbox3.setText(str(solar_system[text_str][3]))
        textbox4.setText(str(solar_system[text_str][4]))
        textbox5.setText(str(solar_system[text_str][5]))
        textbox6.setText(str(solar_system[text_str][6]))
        textbox7.setText(str(solar_system[text_str][0]))
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
label2 = QLabel(w)
label2.setText('Test Particles')
label2.move(250, 185)
combobox2 = QComboBox(w)
combobox2.addItem('Test Particle 1')
combobox2.addItem('Test Particle 2')
combobox2.addItem('Test Particle 3')
combobox2.addItem('Test Particle 4')
combobox2.move(250, 205)
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
label3.setText('List of Particles')
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
button1.resize(button1.sizeHint())
button1.move(60, 400)
def on_click_button1():
    try:
        pos = np.array([float(textbox1.text()),float(textbox2.text()),float(textbox3.text())])
        vel = np.array([float(textbox4.text()),float(textbox5.text()),float(textbox6.text())])
        mass = float(textbox7.text())
        name = str(textbox8.text())
    except ValueError:
        print("Error: invalid parameters.")
        return 1
    new_obj = particle(mass, pos, vel, name)
    objects.append(new_obj)
    table_items(objects) # Updates table

# Button for starting the simulation
button2 = QPushButton('Start', w)
button2.setToolTip('Start simulation')
button2.resize(button2.sizeHint())
button2.move(250, 240)
def on_click_button2():
    try:
        sample_rate = int(textbox13.text())
        num_years = int(textbox14.text())
    except ValueError: # If these not entered, set default values
        sample_rate = 1
        num_years = 50
        print("Warning: number of steps not entered; set to default values.")
    num_steps = 365.25*sample_rate*num_years
    global objects
    if len(objects) == 0:
        print("Error: no objects entered.")
        return 1
    objects, old_objects = simulate(objects, num_steps, sample_rate)
    if checkbox3.isChecked():
        animate(objects, num_steps, checkbox4.isChecked())
    else:
        plot(objects, old_objects, checkbox1.isChecked(), checkbox2.isChecked())

# Button for running simulation with all random particles
button3 = QPushButton('Random', w)
button3.setToolTip('Start simulation with all random particles')
button3.resize(button3.sizeHint())
button3.move(480, 420)
def on_click_button3():
    try:
        N = int(textbox9.text())
        r0 = float(textbox10.text())
        vmax = float(textbox11.text())
        mmax = float(textbox12.text())
    except ValueError:
        N = 12
        r0 = 6e12
        vmax = 5e4
        mmax = 5e30
        print("Warning: initial parameters not entered; set to default values.")
    try:
        sample_rate = int(textbox13.text())
        num_years = int(textbox14.text())
    except ValueError: # If these not entered, set default values
        sample_rate = 1
        num_years = 50
        print("Warning: number of steps not entered; set to default values.")
    num_steps = 365.25*sample_rate*num_years
    objects = gen_random(N, r0, mmax, vmax)
    objects, old_objects = simulate(objects, num_steps, sample_rate)
    if checkbox3.isChecked():
        animate(objects, num_steps, checkbox4.isChecked())
    else:
        plot(objects, old_objects, checkbox1.isChecked(), checkbox2.isChecked())

# Button for clearing old list of objects
button4 = QPushButton('Reset', w)
button4.setToolTip('Clear old list of objects')
button4.resize(button4.sizeHint()); button4.move(250, 270)
def on_click_button4():
    global objects
    objects = []
    textbox1.setText('Position X')
    textbox2.setText('Position Y')
    textbox3.setText('Position Z')
    textbox4.setText('Velocity X')
    textbox5.setText('Velocity Y')
    textbox6.setText('Velocity Z')
    textbox7.setText('Mass')
    textbox8.setText('Name')
    textbox9.setText('Number of particles')
    textbox10.setText('Initial radius')
    textbox11.setText('Initial max velocity')
    textbox12.setText('Initial max mass')
    textbox13.setText('Steps per day')
    textbox14.setText('Number of years')
    table.setColumnCount(0)
    table.setRowCount(0)
    plt.close()

button5 = QPushButton('Quick add', w)
button5.setToolTip('Add all solar system objects')
button5.resize(button5.sizeHint()); button5.move(60, 430)
def on_click_button5():
    for i in solar_system.keys():
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
