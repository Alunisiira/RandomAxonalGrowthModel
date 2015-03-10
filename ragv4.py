# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 17:28:16 2015
Last Update on Fri Feb 20 16:02:53 2015

This is the fourth version of a python script modeling and simulating the 
random axon growth model.
See following papers:
#TODO: papers

@author: Jan Zelmer
@mail: jzelmer@gmx.net

Future possible improvements:
Add bounding box for geometric classes
Adding new geometric objects and its subroutines
Adding get_middle for geometric objects
Discuss and decide the use of different metrics
"""

# coding: utf8

#zest.releaser
#@watman

from math import pi
execfile("rag_utility.py")


class Simulation ():
	""" This is our "Central processing unit". It manages the whole simulation
		including distance matrix, geotree, neuron generators and neurons.
	"""
	
	def __init__ (self, neuron_generators):
		""" Normal constructor. Just defining some variables we need. And 
			constructing our super area, which contains all other areas.
		"""
		self.generators = neuron_generators		
		self.upper_limit_container_capacity = 0.55
		min_x, min_y, min_z, max_x, max_y, max_z = 100000, 100000, 100000, 0, 0, 0
		id_generator = 1
		for generator in neuron_generators:
			self.areas = generator.get_areas()
			for area in self.areas:
				min_x = min( min_x, area.get_points()[0].x)
				min_y = min( min_y, area.get_points()[0].y)
				min_z = min( min_y, area.get_points()[0].z)
				max_x = max( max_x, area.get_points()[2].x)				
				max_y = max( max_y, area.get_points()[2].y)
				max_z = max( max_z, area.get_points()[2].z)
				area.max_space = (area.get_points()[2].x - area.get_points()[0].x) * (area.get_points()[2].y - area.get_points()[0].y) * (area.get_points()[2].z - area.get_points()[0].z)
				area.ocu_space = 0.0
				area.sim_active = True
				if not hasattr(area, 'id'):
					area.id = 'Area%i' %(id_generator)
					id_generator += 1
		self.simulation_area = n3AxisParallelRectangle(n3Point(min_x, min_y, min_z), n3Point(max_x, max_y, max_z))
		self.bounding_box = self.simulation_area.get_bounding_box()
		self.max_distance = self.bounding_box.diagonal_length
		self.tree = 0
		self.verbose = 0
		self.neuron_path = []

	def set_container_capacity_heuristic_limit (self, float_number):
		""" Sets the percentage when a container will be considered full and
			therefore no new neurons will be placed in there.
		"""
		self.upper_limit_container_capacity = float_number

	def set_verbosity (self, integer):
		self.verbose = integer

	def set_bounding_area (self, area):
		self.simulation_area = area
		self.bounding_box = area.get_bounding_box()
		self.max_distance = self.bounding_box.diagonal_length

	def use_geotree (self, boolean, capacity = 50, search_radius = 5):
		""" Tells the simulation whether a indexing structure should be used,
			to speed the computation up. With 300 neurons the factor is 3
			(theoretically and practically about 2).
			And with more neurons the factor increases.
		"""		
		if boolean:
			self.tree = n3GeoTree(capacity, self.simulation_area)
			self.tree_search_radius = search_radius
		else :
			self.tree = 0

	def simulation_step (self):
		""" This method simulates a simulation step. DO NOT CALL DIRECTLY.
		"""		
		self.simulation_step_counter += 1
		added_neurons = 0
		for generator in self.generators:		
			for area in generator.get_areas():
				if area.sim_active:
					neurons_to_add = generator.new_neurons(self.simulation_step_counter, area)					
					degree_of_capacity = area.ocu_space / area.max_space					
					if degree_of_capacity > self.upper_limit_container_capacity:
						if self.verbose:
							print "Error: %s has reached quite its capacity. Skipping this container for the rest of the simulation." %(area.id)
						area.sim_active = False
						neurons_to_add = 0					
					elif degree_of_capacity > 0.35 and self.verbose:
						print "Warning: %s is going into capacity problems." %(area.id) 					
					added_neurons += neurons_to_add
					for i in xrange(neurons_to_add):				
						nneuron = generator.next_neuron(self.simulation_step_counter, area)
						while not self.is_free(nneuron, nneuron.cellbody_radius):
							nneuron = generator.next_neuron(self.simulation_step_counter, area)
						nneuron.axon = n3Line(nneuron.position)
						nneuron.active = True
						nneuron.dist_index = self.dist_counter
						area.ocu_space +=  pi * nneuron.cellbody_radius * nneuron.cellbody_radius
						self.dist_counter += 1
						if self.tree:
							self.tree.add_element(nneuron)
						self.neurons.append(nneuron)
						self.neuron_path.append([nneuron.position])
		for neuron in self.neurons:			
			if neuron.active:
				# grow 
				if neuron.type == "long" :
					x = round(neuron.axon.head.x)
					y = round(neuron.axon.head.y)
					z = round(neuron.axon.head.z)
					target_area = neuron.target_area
					gradient_field = [self.chemical_gradients[x - 1][y - 1][z - 1][target_area], self.chemical_gradients[x - 1][y][z - 1][target_area], \
					self.chemical_gradients[x - 1][y + 1][z - 1][target_area], self.chemical_gradients[x][y + 1][z - 1][target_area], \
					self.chemical_gradients[x + 1][y + 1][z - 1][target_area], self.chemical_gradients[x + 1][y][z - 1][target_area], \
					self.chemical_gradients[x + 1][y - 1][z - 1][target_area], self.chemical_gradients[x][y - 1][z - 1][target_area], \
					self.chemical_gradients[x][y][z - 1][target_area], \
					self.chemical_gradients[x - 1][y - 1][z][target_area], self.chemical_gradients[x - 1][y][z][target_area], \
					self.chemical_gradients[x - 1][y + 1][z][target_area], self.chemical_gradients[x][y + 1][z][target_area], \
					self.chemical_gradients[x + 1][y + 1][z][target_area], self.chemical_gradients[x + 1][y][z][target_area], \
					self.chemical_gradients[x + 1][y - 1][z][target_area], self.chemical_gradients[x][y - 1][z][target_area], \
					self.chemical_gradients[x][y][z][target_area], \
					self.chemical_gradients[x - 1][y - 1][z + 1][target_area], self.chemical_gradients[x - 1][y][z + 1][target_area], \
					self.chemical_gradients[x - 1][y + 1][z + 1][target_area], self.chemical_gradients[x][y + 1][z + 1][target_area], \
					self.chemical_gradients[x + 1][y + 1][z + 1][target_area], self.chemical_gradients[x + 1][y][z + 1][target_area], \
					self.chemical_gradients[x + 1][y - 1][z + 1][target_area], self.chemical_gradients[x][y - 1][z + 1][target_area], \
					self.chemical_gradients[x][y][z + 1][target_area]]						
					neuron.axon.push(neuron.grow(gradient_field))
				elif neuron.type == "short":					
					neuron.axon.push_delta(neuron.grow())
				self.neuron_path[neuron.dist_index].append(neuron.axon.head)
				# make connections
				if neuron.can_put_connection():
					if self.tree :						
						center = neuron.axon.middle
						radius = neuron.axon.get_length() + self.tree_search_radius
						neurons = self.tree.get_points_within_radius(center, radius)
					else :
						neurons = self.neurons
					for ntest in neurons:
						if ntest is not neuron and ntest.can_receive_connection(neuron) and \
							neuron.axon.compute_distance(ntest.position) < ntest.dendrite_radius and \
							self.dmatrix.add(neuron, ntest, self.generators[0].metric.compute_distance):
							neuron.put_connection(ntest)
							ntest.receive_connection(neuron, neuron.axon)
				neuron.active = neuron.can_put_connection() and self.simulation_area.lies_inside(neuron.axon.head)
		if self.verbose:
			print "Step %i: Added %i new neurons." %(self.simulation_step_counter, added_neurons)

	def simulate (self):
		""" Main method to start the simulation.
		"""
		self.neurons = []
		self.dmatrix = Distances(100)
		self.simulation_step_counter = 0
		self.dist_counter = 0
		self.create_chemical_gradients_field()
		if self.tree :
			self.tree = n3GeoTree(self.tree.limit, self.super_area)
		if self.verbose:
			print "Adding Neurons and growing Axons"
		min_iterations = 0
		for gen in self.generators:
			min_iterations = max(min_iterations, gen.get_minimum_iterations())
		for i in xrange(min_iterations):
			self.simulation_step()
		if self.verbose:
			print "Finishing Growth"
		while not self.finished():			
			self.simulation_step()	
			
	def create_chemical_gradients_field (self):
		chemical_gradients = np.zeros(((self.bounding_box.left_front_lower - self.bounding_box.right_front_lower).length(), \
		(self.bounding_box.left_front_lower - self.bounding_box.left_back_lower).length(), \
		(self.bounding_box.left_front_lower - self.bounding_box.left_front_top).length(), len(self.generators) ))
		diagonal_length = self.bounding_box.get_diagonal_length()
		p = n3Point(0, 0, 0)
		for gen in xrange(len(self.generators)):
				gx = gen.chemical_gradient_source.x
				gy = gen.chemical_gradient_source.y
				gz = gen.chemical_gradient_source.z
				for x in xrange(chemical_gradients.shape[0]):
					for y in xrange(chemical_gradients.shape[1]):
						for z in xrange(chemical_gradients.shape[2]):
							p.x, p.y, p.z = gx, gy, gz
							if self.simulation_area.lies_inside(p):
								for inhibitor in self.chemical_inhibitors:
									chemical_gradients[x][y][z][gen] -= max (0, ( 1 - inhibitor.distance_to(p) / inhibitor.radius) * 5  )
								chemical_gradients[x][y][z][gen] += sqrt((gx - x) * (gx -x) + (gy - y) * (gy - y) + (gz -z) * (gz -z)) / diagonal_length * 10	
								chemical_gradients[x][y][z][gen] = max (0, chemical_gradients[x][y][z][gen])
		self.chemical_gradients = chemical_gradients
		
	def print_simulation_meta_data (self):
		""" Prints some statistics onto the console.
		"""
		print "\nNeurons: ---------------------------------------------------\n"
		print "%i Neurons were added and they made %i connections.\n" %(len(self.neurons), self.dmatrix.compute_statistics(self.max_distance, 1)[0])		
		print "Connections: -----------------------------------------------\n"		
		print "The arithmetric mean of all connection lengths is %f." %(self.dmatrix.compute_arithmetic_mean_length())
		print "The median connection has a length of %f."%(self.dmatrix.get_median_length())
		print "The shortest connection has a length of %f and the longest %f.\n" %(self.dmatrix.get_min_length(), self.dmatrix.get_max_length())		
		print "Containers: ------------------------------------------------\n"		
		for gen in self.generators:
			for area in gen.areas:
				if area.sim_active :
					print 'Container "%s" reached %i%s of its capacity during the simulation.' %(area.id, int (area.ocu_space / area.max_space * 100), "%")
				else :
					print 'Container "%s" hit the given limit and was, at some point, excluded from further neuronal placement. Please check the settings for this area.' %(area.id)
				
		print "A capacity of 70 % means the container is full.\n"
	
	def finished (self):
		""" Determines whether our simulation has finished.
		"""
		for neuron in self.neurons:
			if neuron.active:
				return False
		return True
				
	def is_free(self, point, radius):
		""" Computes whether the given area is already (partly) occupied.
		"""
		if self.tree:
			neurons = self.tree.get_elements_within_radius(point, radius)
		else :
			neurons = self.neurons
		for neuron in neurons:			
			if self.generators[0].metric.compute_distance(point.position, neuron.position) < (radius + neuron.cellbody_radius):
				return False		
		return True

	def get_neurons (self):
		""" Returns the list containg all neurons. Please use with care.
		"""
		return self.neurons

	def get_statistics (self, partition_bins = 10):
		""" Computes the statistics for our model. Means make a bar diagramm
			for the distribution of the distances.
		"""		
		return self.dmatrix.compute_statistics(self.max_distance, partition_bins)

	def get_distance_matrix (self):
		""" Returns the distance matrix for this model. The indexes are in
			order in which neurons was placed. Means: First added neuron is 
			index 0, second added neuron is index 1 and so on.
		"""
		return self.dmatrix

	def saveModelNeuronsAsFile (self, path):
		""" Saves our model (neurons) as coordinates in
			a text file
		"""
		with open(path,'w') as f: 		
			for n in self.neurons:
				s = n.get_position().__str__().strip('()') + "\n" 
				f.write(s)

	def saveModelNeuronsWithAxonsAsFile (self, path):
		""" Saves our model (neurons and respecting axons) as coordinates in
			a text file
		"""
		with open(path,'w') as f:
			for n in self.neurons:
				s = n.get_position().__str__().strip('()') + ", " + n.axon.get_head().__str__().strip('()') +  "\n"
				f.write(s)
