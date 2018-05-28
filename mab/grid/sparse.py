import itertools
from numpy import *
from mab.grid.common import Domain

slow = False
use_numpy = False
class SparseGrid(object):
	def __init__(self, bases, dimension, domain=None, ft=lambda x: x*1., finv=lambda x: x*1.):
		self.ft = ft
		self.finv = finv
		self.dimension = dimension
		#if domain is None:
		#	self.domain = [(0., 1.) for d in range(self.dimension)]
		#else:
		self.domain = Domain(domain)
		self.bases = bases
		self.maxlevel = len(bases)#-1
		levels1d = range(self.maxlevel)
		self.bases = bases
		self.points = []
		pointmap = {}
		self.test = []
		self.indices = [list() for level in range(self.maxlevel)]
		self.indices_levels = [list() for level in range(self.maxlevel)]
		self.indices_values = [list() for level in range(self.maxlevel)]
		self.indices_values_t = [list() for level in range(self.maxlevel)]
		self.indices_points  = [list() for level in range(self.maxlevel)]
		self.indices_weights = [list() for level in range(self.maxlevel)]
		self.level_points = [list() for level in range(self.maxlevel)]
		#for level in levels1d:
		if 1:
			for levels in itertools.product(levels1d, repeat=dimension):
				norm1_level = sum(levels)
				#print "levels", levels
				#print levels
				if norm1_level < self.maxlevel:
					#lowertriangle = len([None for i in range(len(t)-1) if t[i] <= t[i+1]]) == len(t)-1
					#if 1: #lowertriangle:
					if 1:
						#print "norm level", norm1_level, 
						#, lowertriangle
						#for level in t:
						level_indices = list(itertools.product(*[range(len(bases[level])) for level in levels]))
						levels2 = list(itertools.product(*[[level for k in range(len(bases[level]))] for level in levels]))
						#print list(itertools.product(*[[level] for level in levels]))
						#print list(itertools.product(levels))
						points = tuple(itertools.product(*[[basis.x for basis in bases[level]] for level in levels]))
						#print "uniform points", points
						points = [self.domain.scale_uniform_point(point) for point in points]
						#print "scaled points", points
						#weights = tuple(itertools.product(*[[basis.w for basis in bases[level]] for level in levels]))
						weights = tuple(itertools.product(*[[self.domain.scale_uniform_length(basis.w, dimension) for basis in bases[level]] for dimension, level in zip(range(self.dimension), levels)]))
						weights = product(weights, axis=1)
						#print "weights", weights
						#weights = zeros(len(level_indices)) * 0 + 1.
						if slow:
							self.indices[norm1_level].append(level_indices)
							self.indices_levels[norm1_level].append(levels2)
							self.indices_points[norm1_level].append(points)
							self.indices_weights[norm1_level].append(weights)
							values = weights * 0 #zeros(len(weights))
							self.indices_values[norm1_level].append(values)
							self.indices_values_t[norm1_level].append(values * 0)
						else:
							self.indices[norm1_level].extend(list(level_indices))
							self.indices_levels[norm1_level].extend(list(levels2))
							self.indices_points[norm1_level].extend(points)
							self.indices_weights[norm1_level].extend(weights)
							values = weights * 0 #zeros(len(weights))
							self.indices_values[norm1_level].extend(values)
							self.indices_values_t[norm1_level].extend(values * 0)
							
						
						#print "indices, levels, points", level_indices, levels, points
						#print ">>",levels2
						self.test.append((norm1_level, points, levels))
						#print points
						self.points.extend(points)
						self.level_points[norm1_level].extend(points)
						for point in points:
							assert point not in pointmap
							pointmap[point] = 1
							
		if use_numpy:
			names = "level,index,point".split(",")
			names = [k.strip() for k in names]
			types = ["i4", "i4", "f4"]
			dtype_level = zip(names, types)
			names = "value,weight,value_t".split(",")
			names = [k.strip() for k in names]
			types = ["f4", "f4", "f4"]
			dtype_values = zip(names, types)
			self.levelslist = []
			#self.pointslist = []
			self.valueslist = []
			self.purevalueslist = []
			for level in range(self.maxlevel):
				n = len(self.indices[level])
				print len(self.indices[level]), len(self.indices[level][0])
				print "n", n, self.indices[level]
				print self.indices_levels[level]
				levelsarray = zeros((n,self.dimension), dtype=dtype_level).view(recarray)
				valuesarray = zeros((n), dtype=dtype_values).view(recarray)
				purevaluesarray = zeros(n)
				self.levelslist.append(levelsarray)
				self.valueslist.append(valuesarray)
				self.purevalueslist.append(purevaluesarray)
				#for levels, index, values, weights in zip(self.indices_levels[level], self.indices[level], self.indices_values[level], self.indices_weights[level]):
				if 1:
					print "##"
					#print self.indices_levels[level]
					#print array(self.indices_levels[level]).shape
					#print levelsarray.level.shape
					print "p", self.indices_points[level]
					print "p", array(self.indices_points[level])
					levelsarray.level = array(self.indices_levels[level])
					levelsarray.index = array(self.indices[level])
					levelsarray.point = array(self.indices_points[level])
					#valuearray.value = array(self.indices_values[level])
					print valuesarray.shape
					print array(self.indices_weights[level]).shape
					valuesarray.weight = array(self.indices_weights[level])
					#purevaluesarray = 
					#print levelsarray.level[0]
					print levelsarray
					print "test", levelsarray.level
					print "test", levelsarray.index
					print "test", levelsarray.point
					print valuesarray
					#dsa
					
				#self.levellist.append(zeros((dim,n)))
				 
	def ___depr_get_points(self, level, i):
		basis = self.bases[i]
		points = list(itertools.product(*[[basis.x for basis in bases[i]] for i in t]))
		
		
	def feval(self, point, maxlevel=None):
		if maxlevel is None:
			maxlevel = self.maxlevel
		y = 0
		#print "feval", point
		uniform_point = self.domain.scale_to_uniform_point(point)
		for level in range(maxlevel):
			#print "level", level
			if slow:
				for indices, indices_levels, values in zip(self.indices[level], self.indices_levels[level], self.indices_values_t[level]):
					#print "indices, levels, values", indices, indices_levels, values
					for index in range(len(indices)):
						#print "point", point
						#print len(self.bases[level]), indices[index], self.bases[level][0], indices[index][0]
						#print self.bases[level][indices[index][0]]
						#print [self.bases[level][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						#print [self.bases[indices_levels[index][dim]][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						dys = [self.bases[indices_levels[index][dim]][indices[index][dim]](uniform_point[dim]) for dim in range(self.dimension)]
						dy = values[index] * product(dys)
						#print dys, dy 
						y += dy
			else:
				n = len(self.indices[level])
				for i in range(n):
					#overlap_indices = self.bases[self.indices_levels[level][i]]
					#print "point", point
					#print len(self.bases[level]), indices[index], self.bases[level][0], indices[index][0]
					#print self.bases[level][indices[index][0]]
					#print [self.bases[level][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
					#print [self.bases[indices_levels[index][dim]][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
					dys = [self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]](uniform_point[dim]) for dim in range(self.dimension)]
					dy = self.indices_values_t[level][i] * product(dys)
					#print dys, dy 
					y += dy
				
		return self.finv(y)
		
	def feval2(self, point, maxlevel=None):
		if maxlevel is None:
			maxlevel = self.maxlevel
		y = 0
		uniform_point = self.domain.scale_to_uniform_point(point)
		#print "feval", point, uniform_point
		for level in range(maxlevel):
			#print "level", level
			if slow:
				for indices, indices_levels, values in zip(self.indices[level], self.indices_levels[level], self.indices_values[level]):
					#print "indices, levels, values", indices, indices_levels, values
					for index in range(len(indices)):
						#print "point", point
						#print len(self.bases[level]), indices[index], self.bases[level][0], indices[index][0]
						#print self.bases[level][indices[index][0]]
						#print [self.bases[level][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						#print [self.bases[indices_levels[index][dim]][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						dys = [self.bases[indices_levels[index][dim]][indices[index][dim]](uniform_point[dim]) for dim in range(self.dimension)]
						dy = values[index] * product(dys)
						#print dys, dy 
						y += dy
			else:
				n = len(self.indices[level])
				if use_numpy:
					test = [[[self.bases[l][i](uniform_point[dim]) for dim in range(self.dimension)] for i in range(len(self.bases[l]))] for l in range(level+1)]
					for i in range(n):
						p = 1
						for dim in range(self.dimension):
							#p *= test[self.indices_levels[level][i][dim]][self.indices[level][i][dim]][dim]
							l = self.levelslist[level].level[i][dim]
							index = self.levelslist[level].index[i][dim]
							p *= test[l][index][dim]
							# for dim in range(self.dimension)
						dy = self.valueslist[level].value[i] * p #product(dys)
						#print dys, dy 
						y += dy
				else:
					#@dys = [self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]](uniform_point[dim]) for dim in range(self.dimension)]
					#overlap = self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]].overlap(uniform_point[dim]):
					#overlap = [[[self.bases[l][i].overlap(uniform_point[dim]) for dim in range(self.dimension)] for i in range(len(self.bases[l]))] for l in range(level+1)]
					#test = [[[self.bases[l][i](uniform_point[dim]) for dim in range(self.dimension)] for i in range(len(self.bases[l]))] for l in range(level+1)]
					for i in range(n):
						#print n
						overlap_global = True
						if 0:
							for dim in range(self.dimension):
								if not overlap[self.indices_levels[level][i][dim]][self.indices[level][i][dim]][dim]:
									overlap_global = False
									break
							
						if 1:
							for dim in range(self.dimension):
								if not self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]].overlap(uniform_point[dim]):
									overlap_global = False
									break
						if not overlap_global:
							import pdb;
							#pdb.set_trace()
							break
						#print "point", point
						#print len(self.bases[level]), indices[index], self.bases[level][0], indices[index][0]
						#print self.bases[level][indices[index][0]]
						#print [self.bases[level][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						#print [self.bases[indices_levels[index][dim]][indices[index][dim]](point[dim]) for dim in range(self.dimension)]
						#dys = [self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]](uniform_point[dim]) for dim in range(self.dimension)]
						p = 1
						for dim in range(self.dimension):
							#p *= test[self.indices_levels[level][i][dim]][self.indices[level][i][dim]][dim]
							p *= self.bases[self.indices_levels[level][i][dim]][self.indices[level][i][dim]](uniform_point[dim])
							# for dim in range(self.dimension)
						#dys = [test[self.indices_levels[level][i][dim]][self.indices[level][i][dim]][dim] for dim in range(self.dimension)]
						#[self.indices[level][i][dim]](uniform_point[dim]) for dim in range(self.dimension)]
						dy = self.indices_values[level][i] * p #product(dys)
						#dy = self.purevalueslist[level][i] * p
						#print dys, dy 
						y += dy
				
		#print "y =", y
		return (y)				
		#return y
			
						
	def evaluate(self, f):
		#assert self.dimension == 2
		#for t in itertools.product(levels, repeat=dimension):
		#	norm1_level = sum(t)
		#def getter(l):
		#	return reduce(lambda prev, next: prev[next])
		if 0:
			print "evaluate level", 0
			level = 0
			for i, (indices, points) in enumerate(zip(self.indices[level], self.indices_points[level])):
				length = len(indices)
				print "indices,points", indices, points
				#print "i", indices, points
				values = [f(*point) for point in points]
				print values
				#values = [numpy.product([self.bases[0][i](points[i][j])
				#self.indices_values[0].append(values)
				self.indices_values[level][i] = values
				self.indices_values_t[level][i] = self.ft(values)
		 
		#print self.maxlevel
		for level in range(0, self.maxlevel):
			#print "evaluate level", level
			if slow:
				for i, (indices, points) in enumerate(zip(self.indices[level], self.indices_points[level])):
					values = [f(*point) for point in points]
					values = array(values)
					#print "values", values
					#values -= array([self.feval(point, maxlevel=level) for point in points]) 
					#values = [numpy.product([self.bases[0][i](points[i][j])
					#print len(self.indices_values[level]), i
					self.indices_values[level][i] = values
					self.indices_values_t[level][i] = self.ft(values * 1)
			else:
				#length = len(indices)
				#print "indices,points", indices, points
				values = [f(*point) for point in self.indices_points[level]]
				values = array(values)
				#print level, values
				#print "values", values
				#values -= array([self.feval(point, maxlevel=level) for point in points]) 
				#values = [numpy.product([self.bases[0][i](points[i][j])
				#print len(self.indices_values[level]), i
				self.indices_values[level] = values
				self.indices_values_t[level] = self.ft(values * 1)
			if level > 0:
				if slow:
					for i, (indices, points) in enumerate(zip(self.indices[level], self.indices_points[level])):
						x = array([self.feval(point, maxlevel=level) for point in points])
						x2 = array([self.feval2(point, maxlevel=level) for point in points])
						#print "correction", x2, "to", self.indices_values[level][i] 
						self.indices_values[level][i] -= x2
						self.indices_values_t[level][i] -= self.ft(x)
				else:
					#print "values", self.indices_values[level][i]
					x = array([self.feval(point, maxlevel=level) for point in self.indices_points[level]])
					x2 = array([self.feval2(point, maxlevel=level) for point in self.indices_points[level]])
					#print "correction", x2, "to", self.indices_values[level][i] 
					print level, "surplus", list(self.indices_values[level])
					self.indices_values[level] -= x2
					self.indices_values_t[level] -= self.ft(x)
				#self.indices_values_t[level][i] = self.ft(self.indices_values[level][i])
					#print "surplus", self.indices_values[level][i]
					#print "surplus", self.indices_values_t[level][i]
				#for i, (indices, points) in enumerate(zip(self.indices[level], self.indices_points[level])):
					#print "points", [self.feval2(point) for point in points]
				#	pass
		#print self.indices_values
		#self._calcsurplus()
				
	def _calcsurplus(self):
		for level in range(1, self.maxlevel):
			print "converting to surplus", level
			for i, (indices, points) in enumerate(zip(self.indices[level], self.indices_points[level])):
				#print "values", self.indices_values[level][i]
				x = array([self.feval(point, maxlevel=level) for point in points])
				x2 = array([self.feval2(point, maxlevel=level) for point in points])
				#print "correction", x
				self.indices_values[level][i] -= x2  
				self.indices_values_t[level][i] -= self.ft(x)
				#self.indices_values_t[level][i] = self.ft(self.indices_values[level][i])
				#print "surplus", self.indices_values[level][i]
			
		#for t in zip(self.indices, self.indices_levels, self.indices_points, self.indices_values):
		#	print t
		#print 
	
	def intergral(self, maxlevel=None):
		if maxlevel is None:
			maxlevel = self.maxlevel
		I = 0
		for level in range(0, maxlevel):
			for values, weights in zip(self.indices_values[level], self.indices_weights[level]):
				#print values
				#print weights
				I += sum((values) * weights)
				#I += sum(exp(values) + exp(weights))
		return (I)
			
	def marginalize(self, axis=0):
		print "marginalize"
		grid1d = SparseGrid(self.bases, 1, [self.domain.domain[axis]], ft=self.ft, finv=self.finv)
		#grid1dc = SparseGrid(self.bases, 1, ft=self.ft, finv=self.finv)
		print self.points
		print self.indices
		print grid1d.points
		
		print grid1d.indices
		print
		#for indices, levels, points, values in zip(self.indices, self.indices_levels, self.indices_points, self.indices_values):
		for level in range(self.maxlevel):
			for indices, indices_levels, values in zip(self.indices[level], self.indices_levels[level], self.indices_values[level]):
				print "indices, levels, values", indices, indices_levels, values
		for level in range(grid1d.maxlevel):
			for i, (indices, indices_levels) in enumerate(zip(grid1d.indices[level], grid1d.indices_levels[level])):
				#print "indices, levels", indices, indices_levels
				for index, levels in zip(indices, indices_levels):
					print index, level, levels
					for level2 in range(self.maxlevel):
						for indices2, indices_levels2, values2, weights2, points2 in zip(self.indices[level2], self.indices_levels[level2], self.indices_values[level2], self.indices_weights[level2], self.indices_points[level2]):
							for index2, level2, value2, weight2, point2 in zip(indices2, indices_levels2, values2, weights2, points2):
								if (index2[axis] == index[0]) and (level2[axis] == level):
									print "indices, levels, values", index2, level2, value2
									w = weight2 / grid1d.indices_weights[level][i][index]
									#value2 = self.feval2(point2)
									grid1d.indices_values[level][i][index] += (value2) * w
		for level in range(grid1d.maxlevel):
			for i in range(len(grid1d.indices_values[level])):
				grid1d.indices_values_t[level][i] = self.ft(grid1d.indices_values[level][i])
		#grid1d._calcsurplus()
					
		return grid1d
					 
				#for index in range(len(indices)):
			#print indices, levels, points, values
			#print self.indices[level], self.indices_points[level]
			#print " ", grid1d.indices[level]#, grid1d.indices_levels[level]
			
		
		
		
		
		