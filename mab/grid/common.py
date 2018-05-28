

class Domain(object):
	def __init__(self, domain=None):
		self.domain = domain
		
	def scale_uniform_point(self, point):
		if self.domain is None:
			return point
		else:
			newpoint = []
			for domain1d, coordinate in zip(self.domain, point):
				x1, x2 = domain1d
				new_coordinate = coordinate * (x2-x1) + x1
				newpoint.append(new_coordinate)
			return tuple(newpoint)
	
	def scale_uniform_length(self, length, dimension):
		if self.domain is None:
			return length
		else:
			domain1d = self.domain[dimension]
			x1, x2 = domain1d
			new_length = length * (x2-x1)
			return new_length
	
	def scale_to_uniform_length(self, length, dimension):
		if self.domain is None:
			return length
		else:
			domain1d = self.domain[dimension]
			x1, x2 = domain1d
			new_length = length / (x2-x1)
			return new_length
			
	def scale_to_uniform_point(self, point):
		if self.domain is None:
			return point
		else:
			newpoint = []
			for domain1d, coordinate in zip(self.domain, point):
				x1, x2 = domain1d
				new_coordinate = (coordinate - x1)/(x2-x1)
				newpoint.append(new_coordinate)
			return tuple(newpoint)			
			