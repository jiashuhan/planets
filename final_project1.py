from numpy import *
from random import random
from mayavi import mlab

G = 6.67e-11

class body(object):
	def __init__(self, mass, px, py, pz, vx, vy, vz, nm):
		self.mass = mass
		self.position = [px, py, pz]
		self.velocity = [vx, vy, vz]
		self.xpos = []
		self.ypos = []
		self.zpos = []
		self.nm = nm

	def gravity(self, other):
		distance_x = other.position[0] - self.position[0]
		distance_y = other.position[1] - self.position[1]
		distance_z = other.position[2] - self.position[2]
		distance = (distance_x**2 + distance_y**2 + distance_z**2)**(1/2.)
		force = G * self.mass * other.mass / (distance**2)
		force_x = force * (distance_x / distance)
		force_y = force * (distance_y / distance)
		force_z = force * (distance_z / distance)
		return [force_x, force_y, force_z]

def interact(bodies, timestep, step_num):
	forces_x = []
	forces_y = []
	forces_z = []
	for i in bodies:
		#find force
		total_fx = 0
		total_fy = 0
		total_fz = 0
		for j in bodies:
			#find the other body
			if i != j:
				if ((i.position[0]-j.position[0])**2+(i.position[1]-j.position[1])**2+(i.position[2]-j.position[2])**2) < 1e16:
					new_body = merge(i, j)
					bodies.append(new_body)
					del i, j
				else:
					#check if the two are different
					forces = i.gravity(j)
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
	for k in range(len(forces_x)):
		acc_x = forces_x[k] / bodies[k].mass
		acc_y = forces_y[k] / bodies[k].mass
		acc_z = forces_z[k] / bodies[k].mass
		bodies[k].velocity[0] += acc_x * timestep
		bodies[k].velocity[1] += acc_y * timestep
		bodies[k].velocity[2] += acc_z * timestep
		#update velocity
		bodies[k].position[0] += bodies[k].velocity[0] * timestep
		bodies[k].position[1] += bodies[k].velocity[1] * timestep
		bodies[k].position[2] += bodies[k].velocity[2] * timestep
		#update position
	return bodies

def plotting(bodies):
	for b in bodies:
		mlab.points3d(b.xpos, b.ypos, b.zpos)
		#points3d takes a fourth input (a number) that can determine the size/color of the point

def randomize(num_bodies, init_radius, init_max_mass, init_max_v):
	bodies = []
	for i in bodies:
		r = random() * init_radius
		phi = random() * pi
		theta = random() * 2*pi
		mass = random() * init_max_mass
		x = r * sin(phi) * cos(theta)
		y = r * sin(phi) * sin(theta)
		z = r * cos(phi)
		px = r * sin(phi) * cos(theta)
		py = r * sin(phi) * sin(theta)
		pz = r * cos(phi)
		vx = random()*init_max_v
		vy = random()*init_max_v
		vz = random()*init_max_v
		body = particle(mass, px, py, pz, vx, vy, vz)
		bodies.append(body)
	return bodies

################pseudocode here just for reference!!!!!
def simulate():
	step_num = 0
	previous_bodies = initialize(input1)
	while pause = False:
		mlab.clf()
		interact(previous_bodies, timestep, step_num)
		plotting()
		previous_bodies = bodies
		step_num += 1

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

def pause():
	pass

#################################
#Solar System
sun = body(2e30, 0, 0, 0, 0, 0, 0, 'Sun')
mercury = body(3.3e23, 4.6e10, 0, 0, 0, 5.66e4, 0, 'Mercury') #At perihelion
venus = body(4.87e24, 1.08e11, 0, 0, 0, 3.5e4, 0, 'Venus')
earth = body(6e24, 1.5e11, 0, 0, 0, 3e4, 0, 'Earth')
mars = body(6.4e23, 2.3e11, 0, 0, 0, 2.4e4, 0, 'Mars')
jupiter = body(1.9e27, 7.8e11, 0, 0, 0, 13070, 0, 'Jupiter')
saturn = body(5.7e26, 1.4e12, 0, 0, 0, 9690, 0, 'Saturn')
uranus = body(8.7e25, 2.9e12, 0, 0, 0, 6800, 0, 'Uranus')
neptune = body(1e26, 4.5e12, 0, 0, 0, 5430, 0, 'Neptune')
pluto = body(1.3e22, 4.24e12, 0, 1.3e12, 0, 6100, 0, 'Pluto') #At perihelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
#create sun and all planets
halley = body(2.2e14, -5e12, 0, 1.6e12, 0, 550, 0, 'Halley') #At aphelion; http://nssdc.gsfc.nasa.gov/planetary/factsheet/cometfact.html
#comet1 = body(3e14, -4e9, 5e12, 0, 100, 0, 0) #escape
#comet1 = body(3e14, -4e9, 5e12, 0, 1000, 1000, 0) #capture
hale_bopp = body(1.3e16, 4.1e11, 4.1e11, 5.54e13, -77, 77, 0, 'Hale-Bopp')
comet1 = body(3e14, -4e9, 5e12, 2e12, 500, -1000, 2000, 'Some comet')
##################################
#Binary Stars

##################################
#Randomized Accretion Disk

##################################
#escaping comet

##################################
#Star, planet, and moon

##################################
#Cubic Accretion Disk???