# -*- coding: utf-8 -*-
"""
Created on Tue Mar 03 16:09:52 2015

@author: Jan Zelmer
@mail: jzelmer@gmx.net
"""

import random
execfile("ragv4.py")


class ShortDistanceNeuron (object):
		
	def __init__ (self, position, bounding_area, area_type, rotation_angle):
		""" Trivial constructor
		"""
		self.position = position
		self.type = "short"
		self.cellbody_radius = 0.2
		self.dendrite_radius = 2
		self.bounding_area = bounding_area
		self.rotation_angle = rotation_angle
		self.area_type = area_type
		self.incoming_connections = 0
		self.outgoing_connections = 0
		self.direction_x = random.random() - 0.5
		self.direction_y = random.random() - 0.5
		self.direction_z = max(-0.5, min (0.5, random.normalvariate(0, 0.2)))
		self.direction = n3Point(self.direction_x, self.direction_y, self.direction_z)
		self.direction = self.direction / self.direction.length()
		self.direction = self.direction.translate_by_degree(rotation_angle[0], rotation_angle[1], rotation_angle[2])
		self.growth_speed = 0.3
		self.axon_flexibility = 0.02
	

	def grow (self, *i):
		""" Is called everytime this neuron shall grow. The return value
			is a direction vector.
		"""
		self.direction.x = random.normalvariate(self.direction.x, self.axon_flexibility) # flexibility of axon
		self.direction.y = random.normalvariate(self.direction.y, self.axon_flexibility)
		self.direction.z = random.normalvariate(self.direction.z, self.axon_flexibility)
		self.direction = self.direction / self.direction.length() * self.growth_speed # growthspeed		
		return self.direction
		
	def can_put_connection(self):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""
		return self.outgoing_connections < 10 and self.bounding_area.lies_inside(self.axon.head)
	
	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		self.outgoing_connections += 1
		
	def can_receive_connection(self):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return self.incoming_connections < 10
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		self.incoming_connections += 1
	
	
class LongDistanceNeuron (object):
		
	def __init__ (self, position, target_area):
		""" Trivial constructor
		"""
						
		self.type = "long"
		self.position = position
		self.cellbody_radius = 0.2
		self.dendrite_radius = 3		
		self.bounding_area = 0
		self.incoming_connections = 0
		self.outgoing_connections = 0
		self.target_area = target_area		
		self.grow_speed_constant = 7
		self.axon_flexibility = 0.02
		self.growth_speed = 0.3
		self.area_type = 0
		
	def grow (self, gradient_field):
		""" modified glücksrad auswahl
		"""
		inverse_distance_to_area = max(gradient_field)
		direction = n3Point(0, 0, 0)
		if inverse_distance_to_area < 8:
			# "white matter"	
			liste = map(lambda x: (max(0, x-0.5))**3, gradient_field)
			summe = random.random() * sum(liste)
			growth_speed = (1.1 - inverse_distance_to_area) * self.grow_speed_constant
			i = 0
			while summe > 0:
				summe += liste[i]
				i += 1
			z = i / 9 -1
			i = i - z * 9		
			if i == 0 or i == 1 or i == 2:
				x = -1
			elif i == 3 or i == 7 or i == 8:
				x = 0
			else:
				x = 1
			if i == 0 or i ==6 or i == 7:
				y = -1 
			elif i == 1 or i == 5 or i == 8:
				y = 0
			else:
				y = 1
			direction.x = x * growth_speed + self.axon.head.x
			direction.y = y * growth_speed + self.axon.head.y
			direction.z = z * growth_speed + self.axon.head.z
		else:
			# "grey matter"
			direction.x = random.normalvariate(self.axon.direction_vector_norm.x , self.axon_flexibility) # flexibility of axon
			direction.y = random.normalvariate(self.axon.direction_vector_norm.y, self.axon_flexibility)
			direction.z = random.normalvariate(self.axon.direction_vector_norm.z, self.axon_flexibility)
			direction = direction / direction.length() * self.growth_speed 
		return direction
			
		
	def can_put_connection(self, target_neuron):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""		
		return self.outgoing_connections < 20 and target_neuron.area_type == self.target_area
	
	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		self.outgoing_connections += 1
		
	def can_receive_connection(self):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return self.incoming_connections < 20
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		self.incoming_connections += 1
		
	
class LayersNeuronGenerator (object):
	""" Very trivial generator. Every method implemented is used by the simulation.
		Every generator class must either directly or indirectly inherit from this.
	"""
	
	def __init__(self, area, area_type, area_target, rotation_angle):
		""" Trivial constructor. Metric and areas should generally exist.
		"""
		self.metric = ThreeDimEuclidMetric()
		self.minimun_iterations = 2
		self.area_type = area_type
		self.area_target = area_target
		points = area.get_points()
		self.bounding_box = area
		self.rotation_angle = rotation_angle
		bottom = points[:4]
		top = points [4:]
		layer6 = n3Rectangle([bottom[0], bottom[1], bottom[2], bottom[3], \
			bottom[0] + (top[0]-bottom[0])/6.0, bottom[1] + (top[1]-bottom[1])/6.0, \
			bottom[2] + (top[2]-bottom[2])/6.0, bottom[3] + (top[3]-bottom[3])/6.0])
		layer6.id = "Layer6"
		layer5 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/6.0, bottom[1] + (top[1]-bottom[1])/6.0, \
			bottom[2] + (top[2]-bottom[2])/6.0, bottom[3] + (top[3]-bottom[3])/6.0, \
			bottom[0] + (top[0]-bottom[0])/3.0, bottom[1] + (top[1]-bottom[1])/3.0, \
			bottom[2] + (top[2]-bottom[2])/3.0, bottom[3] + (top[3]-bottom[3])/3.0])
		layer5.id = "Layer5"
		layer4 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/3.0, bottom[1] + (top[1]-bottom[1])/3.0, \
			bottom[2] + (top[2]-bottom[2])/3.0, bottom[3] + (top[3]-bottom[3])/3.0, \
			bottom[0] + (top[0]-bottom[0])/2.0, bottom[1] + (top[1]-bottom[1])/2.0, \
			bottom[2] + (top[2]-bottom[2])/2.0, bottom[3] + (top[3]-bottom[3])/2.0])
		layer4.id = "Layer4"
		layer3 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/2.0, bottom[1] + (top[1]-bottom[1])/2.0, \
			bottom[2] + (top[2]-bottom[2])/2.0, bottom[3] + (top[3]-bottom[3])/2.0, \
			bottom[0] + (top[0]-bottom[0])/1.5, bottom[1] + (top[1]-bottom[1])/1.5, \
			bottom[2] + (top[2]-bottom[2])/1.5, bottom[3] + (top[3]-bottom[3])/1.5])
		layer3.id = "Layer3"
		layer2 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/1.5, bottom[1] + (top[1]-bottom[1])/1.5, \
			bottom[2] + (top[2]-bottom[2])/1.5, bottom[3] + (top[3]-bottom[3])/1.5, \
			bottom[0] + (top[0]-bottom[0])/1.2, bottom[1] + (top[1]-bottom[1])/1.2, \
			bottom[2] + (top[2]-bottom[2])/1.2, bottom[3] + (top[3]-bottom[3])/1.2])
		layer2.id = "Layer2"
		layer1 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/1.2, bottom[1] + (top[1]-bottom[1])/1.2, \
			bottom[2] + (top[2]-bottom[2])/1.2, bottom[3] + (top[3]-bottom[3])/1.2, \
			top[0], top[1], top[2], top[3]])
		layer1.id = "Layer1"
		self.areas = [layer1, layer2, layer3, layer4, layer5, layer6]
		
	def new_neurons(self, simulation_step, area):
		""" Tells the simulation how many neurons will be added in this step
			for the given area
		"""
		if simulation_step < 50 :
			if   area.id=="Layer1":
				return 1
			elif area.id == "Layer2":
				return 1
			elif area.id == "Layer3":
				return 1
			elif area.id == "Layer4":
				return 1
			elif area.id == "Layer5":
				return 1
			else : #Layer6
				return 1
		else :
			return 0

	def next_neuron(self, simulation_step, area):
		""" Returning a new instance of neuron. This will be placed on the 
			given area. If it collides, it won't be placed, at all.
		"""
		x = random.random() * abs(area.left_front_lower - area.right_front_lower) + area.left_front_lower		
		y = random.random() * abs(area.left_front_lower - area.left_back_lower) + area.left_front_lower		
		z = random.random() * abs(area.left_front_lower - area.right_front_top) + area.left_front_lower				
		if   area.id=="Layer1":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer2":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer3":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer4":
			return LongDistanceNeuron(n3Point(x, y, z), self.area_target)
		elif area.id == "Layer5":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		else : #Layer6
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
	
	
p1 = n3Point(0, 0, 0)
p2 = n3Point(0, 50, 0)
p3 = n3Point(50, 50, 0)
p4 = n3Point(50, 0, 0)
p5 = n3Point(0, 0, 30)
p6 = n3Point(0, 50, 30)
p7 = n3Point(50, 50, 30)
p8 = n3Point(50, 0, 30)

r1 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r1.translate(r1.middle * -1 + n3Point(0, 0, 75))
lg1 = LayersNeuronGenerator(r1, "A1", "A2", (0, 0, 0))

r2 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r2.translate(r2.middle * -1 + n3Point(75, 0, 0))
r2.rotate_by_self(0, 90, 0)
lg2 = LayersNeuronGenerator(r2, "A2", "A3", (0, 90, 0))

r3 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r3.translate(r3.middle * -1 + n3Point(-75, 0, 0))
r3.rotate_by_self(0, -90, 0)
lg3 = LayersNeuronGenerator(r3, "A3", "A4", (0, -90, 0))

r4 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r4.translate(r4.middle * -1 + n3Point(0, 0, -75))
r4.rotate_by_self(0, 180, 0)
lg4 = LayersNeuronGenerator(r4, "A4", "A5", (0, 180, 0))

r5 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r5.translate(r5.middle * -1 + n3Point(0, 75, 0))
r5.rotate_by_self(-90, 0, 0)
lg5 = LayersNeuronGenerator(r5, "A5", "A6", (-90, 0, 0))

r6 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r6.translate(r2.middle * -1 + n3Point(0, -75, 0))
r6.rotate_by_self(90, 0, 0)
lg6 = LayersNeuronGenerator(r6, "A6", "A1", (90, 0, 0))
	
s = Simulation([lg1, lg2, lg3, lg4, lg5, lg6])
s.set_bounding_area = n3Sphere(n3Point(0, 0, 0), 100)
#s.simulate()
#s.print_simulation_meta_data()
#print s.get_statistics()
#d = s.get_distance_matrix()
