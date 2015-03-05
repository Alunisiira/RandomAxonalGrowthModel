# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 17:33:22 2015

@author: Jan Zelmer
@mail: jzelmer@gmx.net

"""
from math import sqrt
from math import radians

class n3Point(object):
	
	def __init__ (self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def __add__ (self, p):
		return self.__class__(self.x + p.x, self.y + p.y, self.z + p.z)

	def __sub__ (self, p):
		return self.__class__(self.x - p.x, self.y - p.y, self.z - p.z)
		
	def __mul__ (self, k):
		return self.__class__(self.x * k, self.y * k, self.z * k)
		
	def __div__ (self, k):
         i = float(k)
         return self.__class__(self.x / i, self.y / i, self.z / i)
		
	def __truediv__ (self, k):
         i = float(k)
         return self.__class__(self.x / i, self.y / i, self.z / i)
		
	def scalar_product (self, v):
		return self.x * v.x + self.y * v.y + self.z * v.z	
	
	def vector_product (self, v):
		return self.__class__(self.y * v.z - self.z * v.y, self.z * v.x - self.x * v.z, self.x * v.y - self.y * v.x)
	
	def length (self):
		return sqrt (self.x * self.x + self.y * self.y + self.z * self.z)
	
	def distance_to (self, p):
		return (self - p).length()
		
	def __str__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")")

	def __repr__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")")
		
	def translate_by_degree (self, alpha, beta, gamma):
		""" Drehung im Uhrzeigersinn
		"""
		alpha = radians (alpha)
		beta = radians (beta)
		gamma = radians (gamma)		
		y1 = self.y * cos (alpha) - self.z * sin (alpha)
		z1 = self.y * sin (alpha)	+ self.z * cos (alpha)		
		x1 = self.x * cos (beta) + z1 * sin (beta)		
		z2 = z1 * cos (beta) - self.x * sin (beta)
		x2 = round(x1 * cos(gamma) - y1 * sin (gamma), 4)
		y2 = round(x1 * sin(gamma) + y1 * cos (gamma), 4)
		return self.__class__(x2, y2, z2)
		
		

class n3Line(object):
	
	def __init__(self, *points):
		""" Simple constructor. Pushes the line by iterating over all points
		provided. So the last point is the head and the one before last is 
		the tail.
		"""		
		self.head = points[0]
		for p in points[1:]:
			self.push(p)
		
	def push (self, p): 
		""" Pushes the line by adding the point p. The new head is point p 
			and the new tail is the old head.
		"""
		self.tail = self.head
		self.head = p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()
		self.middle = self.tail + (self.head - self.tail) / 2 

	def push_delta (self, p):
		""" Pushes the line by adding the coordinates form the point p to the
			actual head as the new head. Finally setting the tail to the old
			head.
		"""
		self.tail = self.head
		self.head = self.head + p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()
		self.middle = self.tail + (self.head - self.tail) / 2 

	def get_middle (self):
		""" Returns the point representing the middle.
		"""
		return self.middle

	def get_head (self):
		""" Returns the point representing the head.
		"""
		return self.head
		
	def get_tail (self):
		""" Returns the point representing the tail.
		"""
		return self.tail

	def compute_distance (self, p):		
		""" Computes the distance between a point p and this line. It is 
			basically the minimum distance over all the points from the 
			line and the point p. This is done by computing the perpendicular
			and then using euclid.		
		"""		
		v = (self.tail - p) * -1 
		r = self.direction_vector.scalar_product(v) / self.direction_vector.length()
		if r < 0:
			return self.tail.distance_to(p)
		elif r > self.direction_vector.length():
			return self.head.distance_to(p)
		else:
			return p.distance_to(self.direction_vector_norm * r + self.tail)

	def get_length (self):
		""" Returns the length of this line
		"""
		return (self.head - self.tail).length()
 
 
class n3Sphere (object):
	
	def __init__ (self, point, radius):
		self.center = point
		self.radius = radius
		
	def lies_in (self, p):
		return self.center.distance_to(p) <= self.radius
    
				
class n3Ellipsoid (object):
	
	def __init__ (self, center, x, y, z):
		self.center = center
		self.x = x
		self.y = y
		self.z = z
		
	def lies_in (self, p):
		return (p.x * p.x / self.x / self.x + p.y * p.y / self.y / self.y + p.z * p.z / self.z / self.z ) <= 1
	

class n3Rectangle(object):
    
	def __init__ (self, points): 
		""" Simple constructor with following ordering of the points.
		"""		
		"""Orders points according to: Most left Bottom Point is first, then
			left upper, then right upper and finally right bottom. 	The 
			points are determined via the distance to coordinate origin
			and a point on the x axis lying to the right of the rectangle.
		"""
		
		max_y = 0
		for i in points:
			if i.y > max_y:
				max_y = i.y
		max_x = 0
		for i in points:
			if i.x > max_x:
				max_x = i.x

		d = []
		for i in points:
			d.append(i.distance_to (n3Point(0,0,0)))
		d.sort(lambda x,y: -1 if x < y else 1)			
		if d[0]==d[1]:
			#idiotic special case, when rectangle was axis parallel and got rotated by 45 degrees
			self.left_front_lower, self.right_back_top = self.init_get_extreme_points(points, n3Point(0, 0.1, 0))
			self.left_back_lower, self.right_front_top = self.init_get_extreme_points(points, n3Point(0, max_y + 0,1, 0))
			self.right_back_lower, self.left_front_top = self.init_get_extreme_points(points, n3Point(max_x, max_y + 0.1, 0))
			self.right_front_lower, self.left_back_top = self.init_get_extreme_points(points, n3Point(max_x, 0.1, 0))
		else:
			self.left_front_lower, self.right_back_top = self.init_get_extreme_points(points, n3Point(0, 0, 0))
			self.left_back_lower, self.right_front_top = self.init_get_extreme_points(points, n3Point(0, max_y, 0))
			self.right_back_lower, self.left_front_top = self.init_get_extreme_points(points, n3Point(max_x, max_y, 0))
			self.right_front_lower, self.left_back_top = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))
		self.inner_radius = ((self.right_front_lower - self.left_front_lower)/2).length()
		self.middle = self.left_front_lower + (self.right_back_top - self.left_front_lower) / 2
		self.outer_radius = (self.middle - self.left_front_lower).length()	
		self.update_coords()			
			
	def update_coords (self):
		self.middle = self.left_front_lower + (self.right_back_top - self.left_front_lower) / 2
		self.norm_vector_side_bottom = (self.right_front_lower - self.left_front_lower).vector_product(self.left_back_lower - self.left_front_lower)
		self.norm_vector_side_bottom = self.norm_vector_side_bottom / self.norm_vector_side_bottom.length()
		self.norm_vector_side_left = (self.left_front_lower - self.left_front_top).vector_product(self.left_back_lower - self.left_front_lower)		
		self.norm_vector_side_left = self.norm_vector_side_left / self.norm_vector_side_left.length()
		self.norm_vector_side_front = (self.left_front_lower - self.right_front_lower).vector_product(self.left_front_top - self.left_front_lower)				
		self.norm_vector_side_front = self.norm_vector_side_front / self.norm_vector_side_front.length()
		self.norm_vector_side_top = self.norm_vector_side_bottom * -1
		self.norm_vector_side_right = self.norm_vector_side_left * -1
		self.norm_vector_side_back = self.norm_vector_side_front * -1
		self.hesse_bottom_d = self.norm_vector_side_bottom.scalar_product(self.left_front_lower)
		self.hesse_left_d = self.norm_vector_side_left.scalar_product(self.left_front_lower)
		self.hesse_front_d = self.norm_vector_side_front.scalar_product(self.left_front_lower)
		self.hesse_back_d = self.norm_vector_side_back.scalar_product(self.right_back_top)
		self.hesse_right_d = self.norm_vector_side_right.scalar_product(self.right_back_top)
		self.hesse_top_d = self.norm_vector_side_top.scalar_product(self.right_back_top)
		
	def init_get_extreme_points(self, liste, p):
         l = []
         for i in liste:
             l.append((i.distance_to (p), i))
         l.sort(lambda x,y: -1 if x < y else 1)
         return (l[0][1], l[7][1])
		
	def get_middle (self):
		""" Returns the point in the middle of the rectangle
		"""
		return self.middle
  
	def get_points(self):
		return [self.left_front_lower, self.left_back_lower, self.right_back_lower, self.right_front_lower, self.left_front_top, self.left_back_top, self.right_back_top, self.right_front_top]
		
	def lies_inside (self, p):
		d = self.middle.distance_to(p)
		if d > self.outer_radius:
			return False
		elif d <= self.inner_radius:
			return True
		else :
			return (p.scalar_product(self.norm_vector_side_back) - self.hesse_back_d) >= 0 and (p.scalar_product(self.norm_vector_side_front) - self.hesse_front_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_left) - self.hesse_left_d) >= 0 and (p.scalar_product(self.norm_vector_side_right) - self.hesse_right_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_bottom) - self.hesse_bottom_d) >= 0 and (p.scalar_product(self.norm_vector_side_top) - self.hesse_top_d) >= 0
				
	def rotate_by_root (self, alpha, beta, gamma):
		self.left_front_lower = self.left_front_lower.translate_by_degree(alpha, beta, gamma)
		self.left_back_lower = self.left_back_lower.translate_by_degree(alpha, beta, gamma)
		self.right_back_lower = self.right_back_lower.translate_by_degree(alpha, beta, gamma)
		self.right_front_lower = self.right_front_lower.translate_by_degree(alpha, beta, gamma)
		self.left_front_top = self.left_front_top.translate_by_degree(alpha, beta, gamma)
		self.left_back_top = self.left_back_top.translate_by_degree(alpha, beta, gamma)
		self.right_back_top = self.right_back_top.translate_by_degree(alpha, beta, gamma)
		self.right_front_top = self.right_front_top.translate_by_degree(alpha, beta, gamma)
		self.update_coords()
		
	def translate (self, p):
		self.left_front_lower = self.left_front_lower + p
		self.left_back_lower = self.left_back_lower + p
		self.right_back_lower = self.right_back_lower + p
		self.right_front_lower = self.right_front_lower + p
		self.left_front_top = self.left_front_top + p
		self.left_back_top = self.left_back_top + p
		self.right_back_top = self.right_back_top + p
		self.right_front_top = self.right_front_top + p
		self.update_coords()
			
	def rotate_by_self (self, alpha, beta, gamma):
		middle = self.middle
		self.translate(middle * -1)
		self.rotate_by_root(alpha, beta, gamma)
		self.translate(middle)
  
		
class n3AxisParallelRectangle (n3Rectangle):
	""" The lines of the rectangle are parallel to the axis.
	"""	
	
	def __init__(self, p1, p7):
		p2 = n3Point(p1.x, p7.y, p1.z)
		p3 = n3Point(p7.x, p7.y, p1.z)
		p4 = n3Point(p7.x, p1.y, p1.z)
		p5 = n3Point(p1.x, p1.y, p7.z)
		p6 = n3Point(p1.x, p7.y, p7.z)
		p8 = n3Point(p7.x, p1.y, p7.z)
		
		points = [p1, p2, p3, p4, p5, p6, p7, p8]		
		
		max_y = 0
		for i in points:
			if i.y > max_y:
				max_y = i.y
		max_x = 0
		for i in points:
			if i.x > max_x:
				max_x = i.x
	
		self.left_front_lower, self.right_back_top = self.init_get_extreme_points(points, n3Point(0, 0, 0))
		self.left_back_lower, self.right_front_top = self.init_get_extreme_points(points, n3Point(0, max_y, 0))
		self.right_back_lower, self.left_front_top = self.init_get_extreme_points(points, n3Point(max_x, max_y, 0))
		self.right_front_lower, self.left_back_top = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))	
		self.middle = self.left_front_lower + (self.right_back_top - self.left_front_lower) / 2
	
	def lies_inside (self, p1):
		return self.left_front_lower.x <= p1.x <= self.right_front_lower.x and  \
			self.left_front_lower.y <= p1.y <= self.left_back_lower.y and \
			self.left_front_lower.z <= p1.z <= self.left_front_top.z


class n2Point(object):
	"""Actually this is a location vector. And therefore some vector operations are supported.
	"""		

	def __init__ (self, x, y):
		""" Simple constructor.
		"""
		self.x = x
		self.y = y
		
		
	def __add__ (self, p):
		""" Simple vector addition. Returns a new vector lengthend by the vector.
		"""
		return self.addition_with_vector(p)

	def addition_with_vector (self, p):
		""" Simple vector addition. Returns a new vector lengthend by the vector.
		"""
		return self.__class__(self.x + p.x, self.y + p.y)
	
	def __sub__(self,p):
		""" Simple vector subraction. Returns a new vector shortened by the vector.
		"""
		return self.subtraction_from_vector(p)
	
	def subtraction_from_vector (self, p):
		""" Simple vector subraction. Returns a new vector shortened by the vector.
		"""
		return self.__class__(self.x - p.x, self.y - p.y)
		
	def scalar_product (self, p):
		""" Vector with vector multiplication. Returns the scalar product.
		"""
		return self.x * p.x + self.y * p.y

	def __mul__ (self, scalar):
		""" Simple vector multiplication. Returns a new vector lengthend by the scalar.
		"""
		return self.__class__(self.x * scalar, self.y * scalar)

	def multiplikation_with_scalar (self, scalar):
		""" Simple vector multiplication. Returns a new vector lengthend by the scalar.
		"""
		return self.__class__(self.x * scalar, self.y * scalar)

	def __div__ (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""	
		return  self.division_from_scalar(scalar)
		
	def __truediv__ (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""
		return  self.division_from_scalar(scalar)

	def division_from_scalar (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""
		return self.__class__(self.x / float(scalar), self.y / float(scalar))
		
	def __str__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + ")")

	def __repr__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + ")")

	def length (self):
		""" Returns the length of this vector.
		"""
		return sqrt ( self.x * self.x + self.y * self.y)

	def distance_to (self, p):
		""" Returns the distance between this and the given point by computing
			the length of the resulting vector.			
		"""
		return (self - p).length()

	def get_position(self):
		""" Returns its position which is basically itself.
		"""
		return self
		
	def get_right_normal (self):
		return self.__class__(self.y, self.x * -1)
		
	def translate_by_degree (self, alpha):
		""" Im Uhrzeigersinn
		"""
		alpha *= -1
		alpha = radians (alpha)
		x = round( self.x * cos (alpha) - self.y * sin (alpha), 4)
		y = round( self.x * sin (alpha) + self.y * cos (alpha), 4)
		return self.__class__(x, y)		


class n2Line (object):
	""" This class defines a line. It is defined by 2 points: head and tail.
		For every push head and tail are newly set and the old ones are 
		forgotten. Finally it implements also a function computing the distance
		to a given point.
	"""		
	
	def __init__(self, *points):
		""" Simple constructor. Pushes the line by iterating over all points
		provided. So the last point is the head and the one before last is 
		the tail.
		"""		
		self.head = points[0]
		for p in points[1:]:
			self.push(p)
		
	def push (self, p): 
		""" Pushes the line by adding the point p. The new head is point p 
			and the new tail is the old head.
		"""
		self.tail = self.head
		self.head = p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()
		self.middle = self.tail + (self.head - self.tail) / 2

	def push_delta (self, p):
		""" Pushes the line by adding the coordinates form the point p to the
			actual head as the new head. Finally setting the tail to the old
			head.
		"""
		self.tail = self.head
		self.head = self.head + p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()
		self.middle = self.tail + (self.head - self.tail) / 2

	def get_middle (self):
		""" Returns the point representing the middle.
		"""
		return self.middle 

	def get_head (self):
		""" Returns the point representing the head.
		"""
		return self.head
		
	def get_tail (self):
		""" Returns the point representing the tail.
		"""
		return self.tail

	def compute_distance (self, p):		
		""" Computes the distance between a point p and this line. It is 
			basically the minimum distance over all the points from the 
			line and the point p. This is done by computing the perpendicular
			and then using euclid.		
		"""		
		v = (self.tail - p) * -1 
		r = self.direction_vector.scalar_product(v) / self.direction_vector.length()
		if r < 0:
			return self.tail.distance_to(p)
		elif r > self.direction_vector.length():
			return self.head.distance_to(p)
		else:
			return p.distance_to(self.direction_vector_norm * r + self.tail)

	def get_length (self):
		""" Returns the length of this line
		"""
		return (self.head - self.tail).length()


class n2Sphere (object):
	
	def __init__ (self, point, radius):
		self.center = point
		self.radius = radius
		
	def lies_in (self, p):
		return self.center.distance_to(p) <= self.radius
    
				
class n2Ellipsoid (object):
	
	def __init__ (self, center, x, y):
		self.center = center
		self.x = x
		self.y = y		
		
	def lies_in (self, p):
		return (p.x * p.x / self.x / self.x + p.y * p.y / self.y / self.y) <= 1


class n2Rectangle (object):
	""" This class implements a rectangle and is implemented based on n2Point"""	
	
	def __init__ (self, p1, p2, p3, p4): 
		""" Simple constructor with following ordering of the points.
		"""		
			
		points = [p1, p2, p3, p4]
		
		max_x = 0
		for i in points:
			if i.y > max_x:
				max_x = i.y

		d = []
		for i in [p1, p2, p3, p4]:
			d.append(i.distance_to (n2Point(0,0)))
		d.sort(lambda x,y: -1 if x < y else 1)			
		if d[0]==d[1]:
			#idiotic special case, when rectangle was axis parallel and got rotated by 45 degrees
			self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0.1, 0))
			self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0.1, 0))
		else:
			self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0, 0))
			self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))
		self.inner_radius = ((self.right_front - self.left_front)/2).length()
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
		self.outer_radius = (self.middle - self.left_front).length()
		self.update_coords()
		
		
	def update_coords (self):
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
		self.norm_vector_side_left = (self.left_back - self.left_front).get_right_normal()
		self.norm_vector_side_left = self.norm_vector_side_left / self.norm_vector_side_left.length()
		self.norm_vector_side_right = self.norm_vector_side_left * -1
		self.norm_vector_side_front = (self.left_front - self.right_front).get_right_normal()
		self.norm_vector_side_front = self.norm_vector_side_front / self.norm_vector_side_front.length()
		self.norm_vector_side_back = self.norm_vector_side_front * -1		
		self.hesse_left_d = self.norm_vector_side_left.scalar_product(self.left_front)
		self.hesse_front_d = self.norm_vector_side_front.scalar_product(self.left_front)
		self.hesse_back_d = self.norm_vector_side_back.scalar_product(self.right_back)
		self.hesse_right_d = self.norm_vector_side_right.scalar_product(self.right_back)
		
	def init_get_extreme_points(self, liste, p):
         l = []
         for i in liste:
             l.append((i.distance_to (p), i))
         l.sort(lambda x,y: -1 if x < y else 1)
         return (l[0][1], l[3][1])
				
	def get_points (self): 
		""" Simple getter method. Returns the points of the rectangle in a list.
			The first element is the point in the left bottom, the second
			is in the left upper corner, the next in the right upper corner
			and finally the last is in the right bottom.		
		"""
		return [self.left_front, self.left_back, self.right_back, self.right_front,]
	
	def get_middle (self):
		""" Returns the point in the middle of the rectangle
		"""
		return self.left_front + (self.right_back - self.left_front) / 2

	def lies_inside (self, p):
		d = self.middle.distance_to(p)
		if d > self.outer_radius:
			return False
		elif d <= self.inner_radius:
			return True
		else :
			return (p.scalar_product(self.norm_vector_side_back) - self.hesse_back_d) >= 0 and (p.scalar_product(self.norm_vector_side_front) - self.hesse_front_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_left) - self.hesse_left_d) >= 0 and (p.scalar_product(self.norm_vector_side_right) - self.hesse_right_d) >= 0
				
	def rotate_by_root (self, alpha):
		self.left_front = self.left_front.translate_by_degree(alpha)
		self.left_back = self.left_back.translate_by_degree(alpha)
		self.right_back = self.right_back.translate_by_degree(alpha)
		self.right_front = self.right_front.translate_by_degree(alpha)		
		self.update_coords()
		
	def translate (self, p):
		self.left_front = self.left_front + p
		self.left_back = self.left_back + p
		self.right_back = self.right_back + p
		self.right_front = self.right_front + p		
		self.update_coords()
			
	def rotate_by_self (self, alpha):
		middle = self.middle
		self.translate(middle * -1)
		self.rotate_by_root(alpha)
		self.translate(middle)
  

class n2AxisParallelRectangle (n2Rectangle):
	""" The lines of the rectangle are parallel to the axis.
	"""	
	
	def __init__(self, p1, p2):
		p3 = n2Point(p1.x, p2.y)
		p4 = n2Point(p2.x, p1.y)
		points = [p1, p2, p3, p4]
		
		max_x = 0
		for i in points:
			if i.y > max_x:
				max_x = i.y
			
		self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0, 0))
		self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
	
	def is_inside (self, p1):
		return self.lbPoint.x <= p1.x <= self.ruPoint.x and self.lbPoint.y <= p1.y <= self.ruPoint.y
