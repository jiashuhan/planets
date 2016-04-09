import sys
#from PyQt4.QtGui import *
#from PyQt4.QtCore import pyqtSlot, QTimer, SIGNAL, SLOT
from PySide.QtGui import *
from PySide.QtCore import QTimer, SIGNAL, SLOT
from numpy import pi, sin, cos
from random import random, uniform
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#NOTE: when using PyQt4 instead of PySide, uncomment all "@pyqtSlot()" lines
bodies = []
G = 6.67e-11
#gravitational constant, everything is in metric system
#create application and window
app = QApplication(sys.argv)
w = QWidget()
w.resize(640, 480)
w.setWindowTitle('Particle Control')
#textboxes for entering values
textbox1 = QLineEdit(w)
textbox1.move(20, 20)
textbox1.resize(200, 40)
textbox1.setText('Enter Position (x)')
textbox2 = QLineEdit(w)
textbox2.move(20, 65)
textbox2.resize(200, 40)
textbox2.setText('Enter Position (y)')
textbox3 = QLineEdit(w)
textbox3.move(20, 110)
textbox3.resize(200, 40)
textbox3.setText('Enter Position (z)')
textbox4 = QLineEdit(w)
textbox4.move(20, 155)
textbox4.resize(200, 40)
textbox4.setText('Enter Velocity (x)')
textbox5 = QLineEdit(w)
textbox5.move(20, 200)
textbox5.resize(200, 40)
textbox5.setText('Enter Velocity (y)')
textbox6 = QLineEdit(w)
textbox6.move(20, 245)
textbox6.resize(200, 40)
textbox6.setText('Enter Velocity (z)')
textbox7 = QLineEdit(w)
textbox7.move(20, 290)
textbox7.resize(200, 40)
textbox7.setText('Enter Mass')
textbox8 = QLineEdit(w)
textbox8.move(20, 335)
textbox8.resize(200, 40)
textbox8.setText('Enter Name')
#button for adding new particle
button1 = QPushButton('Add Particle', w)
button1.setToolTip('Add the particle to list')
button1.resize(button1.sizeHint())
button1.move(60, 400)
#select objects from the solar system
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
#checkbox for setting range limits
checkbox1 = QCheckBox(w)
checkbox1.move(250, 155)
checkbox1.setChecked(False)
checkbox1.setText('Set range for Solar System')
#select test particles
label2 = QLabel(w)
label2.setText('Test Particles')
label2.move(250, 185)
combobox2 = QComboBox(w)
combobox2.addItem('Test Particle 1')
combobox2.addItem('Test Particle 2')
combobox2.addItem('Test Particle 3')
combobox2.addItem('Test Particle 4')
combobox2.move(250, 205)
#list of all particles created for the simulation
label3 = QLabel(w)
label3.setText('List of Particles')
label3.move(485, 25)
table = QTableWidget(w)
table.resize(160, 200)
table.move(450, 50)
table_item = QTableWidgetItem()
#update table
def table_items(bodies):
	num_items = len(bodies)
	table.setColumnCount(1)
	table.setRowCount(num_items)
	count = 0
	for i in bodies:
		Name = i.nm
		table.setItem(count, 0, QTableWidgetItem(Name))
		count += 1
#button for starting the simulation
button2 = QPushButton('Start', w)
button2.setToolTip('Start simulation')
button2.resize(button2.sizeHint())
button2.move(250, 240)
#parameters for randomized simulation
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
#button for running simulation with all random particles
button3 = QPushButton('Random', w)
button3.setToolTip('Start simulation with all random particels')
button3.resize(button3.sizeHint())
button3.move(480, 420)
#button for clearing old list of bodies
button4 = QPushButton('Clear', w)
button4.setToolTip('Clear old list of bodies')
button4.resize(button4.sizeHint())
button4.move(250, 270)
#number of steps per day
textbox13 = QLineEdit(w)
textbox13.move(250, 20)
textbox13.resize(150, 30)
textbox13.setText('Number of steps per day')
#number of years
textbox14 = QLineEdit(w)
textbox14.move(250, 60)
textbox14.resize(150, 30)
textbox14.setText('Number of years')
#checkbox for showing legend
checkbox2 = QCheckBox(w)
checkbox2.move(250, 320)
checkbox2.setChecked(False)
checkbox2.setText('Show Legend')

#create events
#@pyqtSlot()
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
#@pyqtSlot()
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
#@pyqtSlot()
def on_click_button1():
	px = float(textbox1.text())
	py = float(textbox2.text())
	pz = float(textbox3.text())
	vx = float(textbox4.text())
	vy = float(textbox5.text())
	vz = float(textbox6.text())
	mass = float(textbox7.text())
	name = str(textbox8.text())
	new_body = body(mass, px, py, pz, vx, vy, vz, name)
	bodies.append(new_body)
	table_items(bodies)	#updates table
#@pyqtSlot()
def on_click_button2():
	steps_per_day = int(textbox13.text())
	num_years = int(textbox14.text())
	simulate(bodies, 86400/steps_per_day, 365*steps_per_day*num_years, steps_per_day)
#@pyqtSlot()
def on_click_button3():
	num_bodies = int(textbox9.text())
	init_radius = float(textbox10.text())
	init_max_v = float(textbox11.text())
	init_max_mass = float(textbox12.text())
	bodies = randomized_bodies(num_bodies, init_radius, init_max_mass, init_max_v)
	steps_per_day = int(textbox13.text())
	num_years = int(textbox14.text())
	simulate(bodies, 86400/steps_per_day, 365*steps_per_day*num_years, steps_per_day)
#@pyqtSlot()
def on_click_button4():
	global bodies, progress
	bodies = []
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
	progress = 0
	bar.setValue(0)
	plt.close()
#randomly generates set of values for an object
def randomize(init_radius, init_max_mass, init_max_v):
	r = random() * init_radius
	phi = random() * pi
	theta = random() * 2*pi
	mass = random() * init_max_mass
	px = r * sin(phi) * cos(theta)
	py = r * sin(phi) * sin(theta)
	pz = r * cos(phi)
	vx = uniform(-1, 1)*init_max_v
	vy = uniform(-1, 1)*init_max_v
	vz = uniform(-1, 1)*init_max_v
	return [mass, px, py, pz, vx, vy, vz]
#generates a list of random objects
def randomized_bodies(num_bodies, init_radius, init_max_mass, init_max_v):
	bodies = []
	for i in range(num_bodies):
		data = randomize(init_radius, init_max_mass, init_max_v)
		nm = str(i)
		data.append(nm)	
		random_body = body(data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7])
		bodies.append(random_body)
	return bodies
#simulation and plotting
def simulate(bodies, timestep, num_steps, steps_per_day):
	global progress
	timer.start(1000)	#starts timer
	duration = num_steps / steps_per_day
	old_bodies = []
	for i in range(0, num_steps):
		#repeat for given number of steps
		forces_x = []
		forces_y = []
		forces_z = []
		for j in bodies:
			#find force for each body
			total_fx = 0
			total_fy = 0
			total_fz = 0
			for k in bodies:
				#find the other body
				if j != k and k.merged == False and k.origin1 != j.nm and k.origin2 != j.nm:
#					print ((j.position[0]-k.position[0])**2+(j.position[1]-k.position[1])**2+(j.position[2]-k.position[2])**2), j.nm, k.nm#, j.merged
					if ((j.position[0]-k.position[0])**2+(j.position[1]-k.position[1])**2+(j.position[2]-k.position[2])**2) < 1e17:
						new_body = merge(j, k)
						new_body.origin1 = j.nm; new_body.origin2 = k.nm
						if j.merged == False:
							j.merged = True; bodies.remove(j); old_bodies.append(j)
						if k.merged == False:
							k.merged = True; bodies.remove(k); old_bodies.append(k)
						bodies.append(new_body)
					else:
						#check if the two are different
						forces = j.gravity(k)
						total_fx += forces[0]
						#total force in x for each body
						total_fy += forces[1]
						#total force in y for each body
						total_fz += forces[2]
			forces_x.append(total_fx)
			#append force in x for each body
			forces_y.append(total_fy)
			#append force in y for each body
			forces_z.append(total_fz)
		for a in range(len(forces_x)):
			acc_x = forces_x[a] / bodies[a].mass
			acc_y = forces_y[a] / bodies[a].mass
			acc_z = forces_z[a] / bodies[a].mass
			bodies[a].velocity[0] += acc_x * timestep
			bodies[a].velocity[1] += acc_y * timestep
			bodies[a].velocity[2] += acc_z * timestep
			#update velocity
			bodies[a].position[0] += bodies[a].velocity[0] * timestep
			bodies[a].position[1] += bodies[a].velocity[1] * timestep
			bodies[a].position[2] += bodies[a].velocity[2] * timestep
			#update position
			bodies[a].xpos.append(bodies[a].position[0])
			bodies[a].ypos.append(bodies[a].position[1])
			bodies[a].zpos.append(bodies[a].position[2])
			#store position in list
#		time += timestep
	fig = plt.figure()
	ax = fig.gca(projection='3d')
#	print bodies, old_bodies
	for b in bodies:
		#plot orbit for each body
		ax.plot(b.xpos, b.ypos, b.zpos, label=b.nm)
#		print b.xpos, b.ypos, b.zpos
#		plt.xlim((-5e16, 5e16)); plt.ylim((-5e16, 5e16))
	for c in old_bodies:
		ax.plot(c.xpos, c.ypos, c.zpos, label=c.nm)
#		print c.xpos, c.ypos, c.zpos
	if checkbox1.isChecked() == True:
		ax.set_xlim(-6e12, 6e12)
		ax.set_ylim(-6e12, 6e12)
		ax.set_zlim(-6e12, 6e12)
	ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
	if checkbox2.isChecked() == True:
		plt.legend()
	plt.show()
#merge two objects if too close together, momentum is conserved
def merge(body1, body2):
	new_name = body1.nm + '+' + body2.nm
	new_px = (body1.position[0] + body2.position[0]) / 2
	new_py = (body1.position[1] + body2.position[1]) / 2
	new_pz = (body1.position[2] + body2.position[2]) / 2
	new_mass = body1.mass + body2.mass
	new_vx = (body1.mass * body1.velocity[0] + body2.mass * body2.velocity[0]) / new_mass
	new_vy = (body1.mass * body1.velocity[1] + body2.mass * body2.velocity[1]) / new_mass
	new_vz = (body1.mass * body1.velocity[2] + body2.mass * body2.velocity[2]) / new_mass
	new_body = body(new_mass, new_px, new_py, new_pz, new_vx, new_vy, new_vz, new_name)
	return new_body

class body(object):
	def __init__(self, mass, px, py, pz, vx, vy, vz, nm):
		#initiate class
		self.mass = mass
		self.position = [px, py, pz]
		self.velocity = [vx, vy, vz]
		self.xpos = []
		self.ypos = []
		self.zpos = []
		self.nm = nm
		self.merged = False
		self.origin1 = None
		self.origin2 = None

	def gravity(self, other):
		#finds gratational force between two bodies
		distance_x = other.position[0] - self.position[0]
		distance_y = other.position[1] - self.position[1]
		distance_z = other.position[2] - self.position[2]
		#the distance vector points radially outward
		distance = (distance_x**2 + distance_y**2 + distance_z**2)**(1/2.)
		force = G * self.mass * other.mass / (distance**2)
		force_x = force * (distance_x / distance)
		force_y = force * (distance_y / distance)
		force_z = force * (distance_z / distance)
		return [force_x, force_y, force_z]
#data for solar system objects
solar_system = {'Sun':[2e30, 0, 0, 0, 0, 0, 0],
            	'Mercury':[3.3e23, 4.6e10, 0, 0, 0, 5.66e4, 0],#At perihelion
              	'Venus':[4.87e24, 1.08e11, 0, 0, 0, 3.5e4, 0],
 	            'Earth':[6e24, 1.5e11, 0, 0, 0, 3e4, 0],
        	    'Mars':[6.4e23, 2.3e11, 0, 0, 0, 2.4e4, 0],
              	'Jupiter':[1.9e27, 7.8e11, 0, 0, 0, 13070, 0],
              	'Saturn':[5.7e26, 1.4e12, 0, 0, 0, 9690, 0],
              	'Uranus':[8.7e25, 2.9e12, 0, 0, 0, 6800, 0],
              	'Neptune':[1e26, 4.5e12, 0, 0, 0, 5430, 0],
              	'Pluto':[1.3e22, 4.24e12, 0, 1.3e12, 0, 6100, 0],#At perihelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
              	'Halley':[2.2e14, -5e12, 0, 1.6e12, 0, 550, 0],#At aphelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/cometfact.html
              	'Hale-Bopp':[1.3e16, 4.1e11, 4.1e11, 5.54e13, -77, 77, 0]
}
#comet1 = body(3e14, -4e9, 5e12, 0, 100, 0, 0) #escape
#comet1 = body(3e14, -4e9, 5e12, 0, 1000, 1000, 0) #capture
#comet1 = body(3e14, -4e9, 5e12, 2e12, 500, -1000, 2000, 'Some comet')
#objects for demonstration
test_particles = {'Test Particle 1':[2e30, 1e12, 0, 0, 0, 1e4, 0],	#1, 2 - fork
                  'Test Particle 2':[2e30, -1e12, 0, 0, 0, 1e4, 0],
                  'Test Particle 3':[2e30, 1e12, 0, 0, 0, 1e3, 0],	#3, 4 - spiral/solenoid
                  'Test Particle 4':[2e30, 1e12, 0, 1e12, -1e3, -1e3, -1e3]	
}
#create a progressbar class
class QProgBar(QProgressBar):
	value = 0
	#@pyqtSlot()
	def update_value(progressBar):
		progressBar.setValue(progress)
#		print(progress)

bar = QProgBar(w)
bar.resize(150, 80)
bar.setValue(0)
bar.move(250, 350)
#Timer for progressbar
timer = QTimer()
bar.connect(timer, SIGNAL('timeout()'), bar, SLOT('update_value()'))

button1.clicked.connect(on_click_button1)	#connect the signal to the slot
button2.clicked.connect(on_click_button2)
button3.clicked.connect(on_click_button3)
button4.clicked.connect(on_click_button4)
combobox1.activated[str].connect(on_activated1)
combobox2.activated[str].connect(on_activated2)
w.show()

sys.exit(app.exec_())