# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 20:40:49 2015

@author: Jan Zelmer
@mail: jzelmer@gmx.net
"""

from copy import deepcopy
from math import pow
from math import e
import random
import numpy as np
execfile("rag_geometry.py")

random.seed()

class n3GeoTree(object):
	
	def __init__ (self, limit, cube):
		self.limit = limit
		self.container = []
		self.cube = cube
		self.is_leaf = True
		self.middle_x = cube.get_middle().x
		self.middle_y = cube.get_middle().y
		self.middle_z = cube.get_middle().z

	def add_element (self, element):				
		trie = self
		while not trie.is_leaf:
			if element.x <= trie.middle_x:
				#left
				if element.y <= trie.middle_y:
					#front
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if element.y <= trie.middle_y:
					#front
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		trie.container.append(element)
		if len(trie.container) > trie.limit:
				#split				
				p = trie.cube.get_points()				
				#bottom
				p.insert(1, p[0] + (p[1] - p[0])/2)
				p.insert(3, p[2] + (p[3] - p[2])/2)
				p.insert(5, p[4] + (p[5] - p[4])/2)
				p.insert(7, p[6] + (p[0] - p[6])/2)
				p.insert(8, p[0] + (p[4] - p[0])/2)
				#top				
				p.insert(10, p[9] + (p[10] - p[9])/2)
				p.insert(12, p[11] + (p[12] - p[11])/2)
				p.insert(14, p[13] + (p[14] - p[13])/2)
				p.append(p[15] + (p[9] - p[15])/2)
				p.append(p[9] + (p[13] - p[9])/2)				
				#middle				
				p.insert(9, p[9] + (p[0] - p[9])/2)
				p.insert(10, p[11] + (p[1] - p[11])/2)
				p.insert(11, p[13] + (p[2] - p[13])/2)
				p.insert(12, p[15] + (p[3] - p[15])/2)
				p.insert(13, p[17] + (p[4] - p[17])/2)
				p.insert(14, p[19] + (p[5] - p[19])/2)
				p.insert(15, p[21] + (p[6] - p[21])/2)
				p.insert(16, p[23] + (p[7] - p[23])/2)
				p.insert(17, p[25] + (p[8] - p[25])/2)
				
				trie.cube_left_front_lower = n3GeoTree(self.limit, n3AxisParallelRectangle(p[0], p[17]))
				trie.cube_left_back_lower  = n3GeoTree(self.limit, n3AxisParallelRectangle(p[1], p[12]))
				trie.cube_right_back_lower = n3GeoTree(self.limit, n3AxisParallelRectangle(p[8], p[13]))
				trie.cube_right_front_lower= n3GeoTree(self.limit, n3AxisParallelRectangle(p[7], p[14]))
				trie.cube_left_front_top    = n3GeoTree(self.limit, n3AxisParallelRectangle(p[9], p[26]))
				trie.cube_left_back_top     = n3GeoTree(self.limit, n3AxisParallelRectangle(p[10], p[21]))
				trie.cube_right_back_top    = n3GeoTree(self.limit, n3AxisParallelRectangle(p[17], p[22]))
				trie.cube_right_front_top   = n3GeoTree(self.limit, n3AxisParallelRectangle(p[16], p[23]))				
				
				"""
				trie.cube_left_front_lower = n3GeoTree(self.limit, n3AxisParallelRectangle([p[0], p[1], p[8], p[7], p[9], p[10], p[17], p[16]]))
				trie.cube_left_back_lower  = n3GeoTree(self.limit, n3AxisParallelRectangle([p[1], p[2], p[3], p[8], p[10], p[11], p[12], p[17]]))
				trie.cube_right_back_lower = n3GeoTree(self.limit, n3AxisParallelRectangle([p[8], p[3], p[4], p[5], p[17], p[12], p[13], p[14]]))
				trie.cube_right_front_lower= n3GeoTree(self.limit, n3AxisParallelRectangle([p[7], p[8], p[5], p[6], p[16], p[17], p[14], p[15]]))
				trie.cube_left_front_top    = n3GeoTree(self.limit, n3AxisParallelRectangle([p[9], p[10], p[17], p[16], p[18], p[19], p[26], p[25]]))
				trie.cube_left_back_top     = n3GeoTree(self.limit, n3AxisParallelRectangle([p[10], p[11], p[12], p[17], p[19], p[20], p[21], p[26]]))
				trie.cube_right_back_top    = n3GeoTree(self.limit, n3AxisParallelRectangle([p[17], p[12], p[13], p[14], p[26], p[21], p[22], p[23]]))
				trie.cube_right_front_top   = n3GeoTree(self.limit, n3AxisParallelRectangle([p[16], p[17], p[14], p[15], p[25], p[26], p[23], p[24]]))								
				"""				
				
				# assume the data is nearly linearly distributed, so we don't have to check whether our newly created container gets full
				for el in trie.container:					
					if el.x <= trie.middle_x:
						#left
						if el.y <= trie.middle_y:
							#front
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_left_front_lower.container.append(el)
							else:
								#top
								trie.cube_left_front_top.container.append(el)								
						else:
							#back
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_left_back_lower.container.append(el)
							else:
								#top
								trie.cube_left_back_top.container.append(el)
					else:
						#right
						if el.y <= trie.middle_y:
							#front
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_right_front_lower.container.append(el)
							else:
								#top
								trie.cube_right_front_top.container.append(el)								
						else:
							#back
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_right_back_lower.container.append(el)
							else:
								#top
								trie.cube_right_back_top.container.append(el)
				trie.is_leaf = False
				trie.container = []	
			
	def get_sourrounding_points (self, coords):
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		#return deepcopy(trie.container)
		return trie.container

	def get_points_within_radius (self, coord, radius):
		result = []		
		
		coords = n3Point(coord.x - radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x + radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y - radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y + radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z - radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z + radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		resultset = set(result)
		result = []
		for t in resultset:
			result[len(result):] = t.container
		return result

class n2GeoTree (object):
	"""Index Structure for accessing spatial information in log(n) time.
		It uses 4 (and splitting) axis parallel rectangles as way to access
		spatial data. It works like a tree on those rectangles
	"""
	
	def __init__ (self, limit, rectangle):
		self.limit = limit
		self.container = []
		self.rectangle = rectangle
		self.is_leaf = True
		self.middle_x = rectangle.get_middle().x
		self.middle_y = rectangle.get_middle().y

	def add_element (self, element):				
		trie = self
		while not trie.is_leaf:
			if element.x <= trie.middle_x:
				#left
				if element.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if element.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		trie.container.append(element)
		if len(trie.container) > trie.limit:
				#split				
				p = trie.rectangle.get_points()				
				
				p.insert(1, p[0] + (p[1] - p[0])/2)
				p.insert(3, p[2] + (p[3] - p[2])/2)
				p.insert(5, p[4] + (p[5] - p[4])/2)
				p.insert(7, p[6] + (p[0] - p[6])/2)
				p.insert(8, p[0] + (p[4] - p[0])/2)

				trie.cube_left_front = n2GeoTree(self.limit, n2AxisParallelRectangle(p[0], p[8]))
				trie.cube_left_back  = n2GeoTree(self.limit, n2AxisParallelRectangle(p[1], p[3]))
				trie.cube_right_back = n2GeoTree(self.limit, n2AxisParallelRectangle(p[8], p[4]))
				trie.cube_right_front= n2GeoTree(self.limit, n2AxisParallelRectangle(p[7], p[5]))				
				
				# assume the data is nearly linearly distributed, so we don't have to check whether our newly created container gets full
				for el in trie.container:					
					if element.x <= trie.middle_x:
						#left
						if element.y <= trie.middle_y:
							#front					
							trie = trie.rect_left_front.container.append(el)
						else:
							#back
							trie = trie.rect_left_back.container.append(el)
					else:
						#right
						if element.y <= trie.middle_y:
							#front
							trie = trie.rect_right_front.container.append(el)
						else:
							#back
							trie = trie.rect_right_back.container.append(el)
				trie.is_leaf = False
				trie.container = []	
			
	def get_sourrounding_points (self, coords):
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		#return deepcopy(trie.container)
		return trie.container

	def get_points_within_radius (self, coord, radius):
		result = []		
		
		coords = n3Point(coord.x - radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x + radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y - radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y + radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z - radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z + radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		resultset = set(result)
		result = []
		for t in resultset:
			result[len(result):] = t.container
		return result


class TwoDimEuclidMetric (object):
	""" Simple 2 dimensional euclidean metric to compute the distance between 
	to points
	"""
	
	def compute_distance (self, p1, p2):
		a = p1.x - p2.x
		b = p1.y - p2.y
		return sqrt (a*a + b*b)
		
		
class ThreeDimEuclidMetric (object):
	""" Simple 3 dimensional euclidean metric to compute the distance between 
	to points
	"""
	
	def compute_distance (self, p1, p2):
		a = p1.x - p2.x
		b = p1.y - p2.y
		c = p1.z - p2.z 		
		return sqrt (a*a + b*b + c*c)
		
	
class Noise ():
	""" "Static class" with can be used to noise data
	"""	
	
	@staticmethod
	def add_gauss_noise (mean, variance):
		""" Returning point is gaussion distributed i.r.t. mean and variance
			If result is below zero, zero will be returned instead
		"""
		return max (0, int (random.gauss(mean, variance)))
		
	@staticmethod
	def add_linear_noise(lower_limit, upper_limit):
		""" Returns a point linearly distributed between both limits
			If result is below zero, zero will be returned instead
		"""
		return max (0, int (random.randint(lower_limit, upper_limit)))

"""
	The following lines define functions, which can be used for time distributions
	for the generators. Every class (should) define a get_value method, which
	computes a y (function value) for the given x.
"""

class GaussAlikeFunction ():
	""" This is originally a gauss function. But the area sums not to 1 anymore.
		But is scaled so, that the mean has the function value "highest_point".
		Means a function with very high variance behaves like a constant 
		function: for every x the function value is highest_point
	"""
	
	def __init__ (self, mean, standard_deviation, highest_point, lower_bound, upper_bound ):
		self.mean = mean
		self.highest_point = highest_point
		self.variance = 	standard_deviation * standard_deviation		
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound	
	
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		else:	
			# proudly programmed by Katja			
			#gauss_param = 1 / sqrt(pi * self.variance*self.variance * 2)  * pow (e, - 1 * (pow((x - self.mean) / sqrt(2) / self.variance, 2)))
			#gauss_param = 1 / pow (e, pow (x - self.mean, 2) / 2 / self.variance) / sqrt (2 * pi * self.variance)			
			return self.highest_point / pow (e, pow (x - self.mean, 2) / 2 / self.variance)


class FuzzyFunction ():
	""" This function looks like a pyramid. The mean is the top of the pyramid,
		with lower and upper bound defining the base.
	
	"""	
	
	def __init__ (self, mean, highest_point, lower_bound, upper_bound):
		self.mean = mean
		self.highest_point = highest_point
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound	
	
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		elif x <= self.mean :
			result = self.highest_point * (x - self.lower_bound) / (self.mean - self.lower_bound)
		else :
			result = self.highest_point * (1 - float(x - self.mean) / (self.upper_bound - self.mean))
		return result


class QuadraticFunction ():
	""" It needs 3 Points and then will compute the underlying square function
		of it. 	
	"""
		
	def __init__ (self, p1, p2, p3, lower_bound, upper_bound):		
		""" The computation of the function is according to this website (20.02.2015)
			http://www.arndt-bruenner.de/mathe/10/parabeldurchdreipunkte.htm
		"""
		liste = [p1, p2, p3]		
		liste.sort()
		x1 = liste[0].x
		x2 = liste[1].x
		x3 = liste[2].x
		y1 = liste[0].y
		y2 = liste[1].y
		y3 = liste[2].y
		self.a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)
		self.a /= (x1 - x2) * (x1 - x3) * (x3 - x2)
		self.b = x1 * x1 * (y2 - y3) + x2 * x2 * (y3 - y1)	 + x3 * x3 * (y1 - y2)
		self.b /= (x1 - x2) * (x1 - x3) * (x2 - x3)
		self.c = x1 * x1 * (x2 * y3 - x3 * y2) + x1 * (x3 * x3 * y2 - x2 * x2 * y3) + x2 * x3 * y1 * (x2 - x3)
		self.c /= (x1 - x2) * (x1 - x3) * (x2 - x3)
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound
		
	def get_value (self, x):
		""" Returns the function value for a given x. It is always positiv.
			When below zero, zero is instead returned
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		else:	
			return max(0, self.a * x * x + self.b * x + self.c)

class Scipy1DInterpolation ():
	""" This function uses a Scipy interpolation function to define our function.
		Attention: the range is strictly limited to the given points.
		If you receive:
		ValueError: A value in x_new is above/below the interpolation range.
		it is because of this.
	"""
	
	def __init__ (self, *points):
		from scipy.interpolate import interp1d
		x = []
		y = []
		for p in points:
			x.append(p.x)
			y.append(p.y)
		self.f = interp1d(x, y, kind='cubic')
		
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		return self.f(x) * 1
		
	
class Distances ():
	""" Class which holds the distance matrix and the according functions like
		adding new distances, getting distances and doing the statistics.
		The elements given in add and get have to have the property dist_index
		which donates the index in this matrix.
	"""
	
	def __init__ (self, init_container_size = 10, growth_factor = 1.8):
		""" Simple constructor, just defining variables.
		"""
		self.growth_factor = growth_factor		
		self.container = np.zeros((init_container_size,init_container_size), dtype = float)
			
	def add (self, element1, element2, distance_function):
		""" Adds a distance to the matrix. If the distance is already there,
			no computation is done. If not the distance function is called
			with the 2 given arguments and then written to the matrix.
			If the container can't hold the new elements, the container will
			be extended.
			Finally True will be returned if the connection was computed and
			added and False otherwise. (For example, if distance was already
			computed False will be returned)
		"""
		container_size = len(self.container)
		
		if element1.dist_index > container_size or element2.dist_index > container_size:
			# grow			
			new_container_size = len(self.container) * self.growth_factor
			new_container = np.zeros((new_container_size,new_container_size), dtype = float)			
			for i in xrange(container_size):
				new_container[i][:container_size] = self.container[i]			
			self.container = new_container
			return self.add(element1, element2, distance_function)
		elif	self.container[element1.dist_index-1][element2.dist_index-1] == 0:
			self.container[element1.dist_index-1][element2.dist_index-1] = distance_function(element1, element2)
			return True
		else:
			return False

	def get (self, element1, element2):
		""" Returns the distance between those 2 elements
		"""
		return self.container[element1.dist_index-1][element2.dist_index-1]

	def compute_statistics (self, upper_limit, partition_bins = 10):
		""" Computes the statistics for our model. Means make a bar diagramm
			for the distribution of the distances.
			The maximal distance for the partitioning is also needed.		
		"""
		result = np.zeros((partition_bins), dtype = int)
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					result[int(self.container[i][j]/upper_limit * partition_bins)] += 1
		return result
		
	def compute_arithmetic_mean_length (self):
		""" Computes the arithmetric mean of non-zero entries and returns it.
		"""
		result = 0
		number_of_elements = 0
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					number_of_elements += 1
					result += self.container[i][j]
		return result / float(number_of_elements)
		
	def get_min_length (self):	
		""" Returns the minimal length in this matrix. Zero entries don't count in. 
		"""
		smallest_element = 10000
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0  and self.container[i][j] < smallest_element:
					smallest_element = self.container[i][j]
		return smallest_element
		
	def get_max_length (self):
		""" Returns the maximal length in this matrix
		"""
		return np.amax(self.container)

	def get_median_length (self):
		""" Returns the median length of the distance matrix. Zero entries 
			don't count in.
		"""
		result =[]
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					result.append(self.container[i][j])	
		return np.median(result)

	def get_matrix (self):
		""" Returns the distance matrix
		"""
		return deepcopy(self.container)
		
