from numpy import *
from mayavi import mlab

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

