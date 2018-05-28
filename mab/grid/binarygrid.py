import itertools
import numpy
import mab.parallelize
from mab.grid.sparse import SparseGrid
def scale_domain_1d(domain1d, scale):
	return domain1d[0] + (domain1d[1] - domain1d[0])*scale

def inside_domain1d(domain1d, coordinate):
	return (coordinate >= domain1d[0]) & (coordinate <= domain1d[1])


use_quad = False

class BinaryNode(object):
	def __init__(self, grid, domain, dimension, level=1):
		self.grid = grid
		self.level = level
		self.dimension = dimension
		self.domain = domain
		self.childNodes = []
		self.subgrid, self.fine_subgrid = grid.create_grids(domain, 0, 1)
		self.is_refined = False
		self.refinement_level = 0
		for point in self.subgrid.points:
			self.grid.add_point(point)
		for point in self.fine_subgrid.points:
			self.grid.add_point(point)
		for point in self.subgrid.level_points[0]:
			self.grid.add_corner(point)
		self.midpoints = [domain1d[0] + 0.5 * (domain1d[1] - domain1d[0]) for domain1d in self.domain]
		if 0:
			if use_quad:
				def f1(x):
					return (x-0.5)*(x-1.0)/(0.0-0.5)/(0.0-1.0)
				def f2(x):
					return (x-0.0)*(x-1.0)/(0.5-0.0)/(0.5-1.0)
				def f3(x):
					return (x-0.0)*(x-0.5)/(1.0-0.0)/(1.0-0.5)
				scales_and_basis = [(0., lambda x: f1(x)), (0.5, lambda x: f2(x)/f2(0.5)), (1.0, lambda x: f3(x))]
			else:
				scales_and_basis = [(0., lambda x: 1-x), (1.0, lambda x: x)]
			#fbasis = [, lambda x: x]
			self.corners = []
			self.fs = []
			for scales_and_basis in itertools.product(scales_and_basis, repeat=self.dimension):
				#print scales_and_basis
				scales = [k[0] for k in scales_and_basis]
				basis = [k[1] for k in scales_and_basis]
				def f(point, scales=scales, basis=basis):
					#return numpy.prod([scale*fbasis[i]((domain[i][0]-point[i])/(domain[i][1]-domain[i][0])) for i in range(self.dimension)])
					#print point
					i = 0 
					#print ([((point[i]-self.domain[i][0])/(self.domain[i][1]-self.domain[i][0])) for i in range(self.dimension)])
					#print self.domain
					#print fbasis[i](0.1)((domain[i][0]-point[i])/(domain[i][1]-domain[i][0]))
					x =  [basis[i]((point[i]-self.domain[i][0])/(self.domain[i][1]-self.domain[i][0])) for i in range(self.dimension)]
					if numpy.any(numpy.isinf(x)):
						import pdb; pdb.set_trace()
					return numpy.prod([basis[i]((point[i]-self.domain[i][0])/(self.domain[i][1]-self.domain[i][0])) for i in range(self.dimension)])
				point = tuple([scale_domain_1d(domain, scale) for domain, scale in zip(self.domain, scales)])
				#3 "point", point
				grid.add_point(point)
				self.corners.append(point)
			self.fs.append(f)
			
	def calculate_volumes(self, parentvolume):
		# volume per points in this grid 
		volume = parentvolume / len(self.fine_subgrid.points)
		if self.has_children():
			for node in self.childNodes:
				node.calculate_volumes(volume)
		else:
			for points in self.fine_subgrid.points:
				self.grid.add_volume(points, volume)
		#for corner in self.corners:
		#	self.grid.add_volume(volume)
			
	def has_children(self):
		return len(self.childNodes) != 0
	
	def has_grand_children(self):
		return self.has_children() and self.childNodes[0].has_children()
			
	def inside(self, point):
		insides = [None for coordinate, domain1d in zip(point, self.domain) if inside_domain1d(domain1d, coordinate)]
		return len(insides) == self.dimension
	
	def split(self, levels=1):
		if not self.childNodes: # do not split again!
			scales = [(0., 0.5), (0.5, 1.0)]
			for scales in itertools.product(scales, repeat=self.dimension):
				domain = [(self.domain[k][0] + (self.domain[k][1]-self.domain[k][0])*scale[0], self.domain[k][0] + (self.domain[k][1]-self.domain[k][0])*scale[1])  for scale, k in zip(scales, range(self.dimension))]
				#print scales, domain
				self.childNodes.append(BinaryNode(self.grid, domain, self.dimension, self.level+1))
		if levels > 1:
			for childNode in self.childNodes:
				childNode.split(levels=levels-1)
				
	def refine(self):
		#self.childNodes
		self.is_refined = True
		self.refinement_level += 1
		self.subgrid, self.fine_subgrid = self.grid.create_grids(self.domain, self.refinement_level, self.refinement_level+1)
		for point in self.subgrid.points:
			self.grid.add_point(point)
		for point in self.fine_subgrid.points:
			self.grid.add_point(point)
		
	def interpolate(self, point, fine=True):
		index = 0
		#print "in", self.domain, point, self.midpoints
		if self.childNodes:
			for i, childNode in enumerate(self.childNodes):
				inside = childNode.inside(point)
				#print "inside?", i, inside
				if inside:
					return childNode.interpolate(point, fine=fine)
		else:
			return self.interpolate_self(point, fine=fine)
		
	def interpolate_self(self, point, fine=True):
		if fine:
			return self.fine_subgrid.feval2(point)
		else:
			return self.subgrid.feval2(point)
		#values = [self.grid.values[p] for p in self.corners]
		#return sum([value * f(point) for value, f in zip(values, self.fs)])
		
	def update_grids(self):
		def eval_wrapper(*point): # acts like a function, but takes the values from the dictionary
			return self.grid.values[point]
		#self.subgrid.evaluate(eval_wrapper)
		self.fine_subgrid.evaluate(eval_wrapper)
		for child in self.childNodes:
			child.update_grids()

	def calculate_splitlist(self, splitlist, parent_volume):
		volume = parent_volume / 2**self.grid.dimension 
		#len(self.corners)
		if self.has_children() and not self.has_grand_children():
			if 0:
				N = 5
				scales = numpy.arange(N)/(N-1.0)
				scales = [0.12, 0.41, 1-0.41, 1-0.12] 
				testpoints = []
				for scales in itertools.product(scales, repeat=self.dimension):
					#print scales, self.domain
					point = tuple([scale_domain_1d(domain, scale) for domain, scale in zip(self.domain, scales)])
					#print "(% 7f % 7f)" % point, "% 10f % 10f % 10f" % (self.interpolate(point), self.interpolate_self(point), self.interpolate(point)-self.interpolate_self(point))
					testpoints.append(point)
				values1 = [self.interpolate(p) for p in testpoints]
				values2 = [self.interpolate_self(p) for p in testpoints]
				splitlist.append((self, values1, values2, volume))
			else:
				for child in self.childNodes:
					testpoints = child.fine_subgrid.points
					values1 = [self.interpolate(p) for p in testpoints]
					values2 = [self.interpolate_self(p) for p in testpoints]
					values3 = [self.interpolate(p) for p in testpoints]
					values4 = [self.interpolate(p,fine=False) for p in testpoints]
					#print "=========="
					#print self.domain
					#print testpoints
					#print "values1", values1
					#print "values2", values2
					splitlist.append((child, values1, values2, values3, values4, volume))
					#for point in points:
					#	testpoints.append(point)
					
		else:
			for child in self.childNodes:
				child.calculate_splitlist(splitlist, volume)
	
	def get_point_levels(self):
		for point in self.subgrid.points:
			self.grid.add_point_level(point, self.level+self.refinement_level)
		for point in self.fine_subgrid.points:
			self.grid.add_point_level(point, self.level+self.refinement_level)
		for child in self.childNodes:
			child.get_point_levels()
			
	def check_level_difference(self, max_level_difference=1):
		if self.has_children():
			for child in self.childNodes:
				child.check_level_difference()
		else:
			maxlevel = self.level+self.refinement_level
			for point in self.fine_subgrid.points:
				maxlevel = max(maxlevel, self.grid.get_point_level(point))
			for point in self.subgrid.points:
				maxlevel = max(maxlevel, self.grid.get_point_level(point))
			if maxlevel - max_level_difference > self.level+self.refinement_level: # too large difference
				#print "forced split"
				self.split()
			
			
		
class BinaryGrid(object):
	def __init__(self, dimension, baseslist, domain=None):
		self.dimension = dimension
		if domain is None:
			domain = [(0., 1.)] * dimension
		self.domain = domain
		self.baseslist = baseslist
		
		self.points = list()
		self.corner_points = list()
		self.newpoints = list()
		
		self.rootNode = BinaryNode(self, self.domain, dimension, level=1)
		self.dimension = dimension
		self.volumes = {}
		self.values = {}
		self.maxlevels = {}
		
	def create_grids(self, domain, *refinements):
		subgrids = []
		#for bases in self.baseslist[refinements]:
		for index in refinements:
			bases = self.baseslist[index]
			subgrid = SparseGrid(bases, self.dimension, domain)
			subgrids.append(subgrid)
		return subgrids
		
		
	def __len__(self):
		return len(self.points)
	
	def integrate(self, f=numpy.exp):
		return sum([self.volumes[p]*f(self.values[p]) for p in self.points])
		
	def calculate_volumes(self):
		self.volumes = {}
		for point in self.points:
			self.volumes[point] = 0
		self.rootNode.calculate_volumes(1.0)
	
	def add_volume(self, point, volume):
		self.volumes[point] += volume
		
	def add_point(self, point):
		if point not in self.points:
			self.points.append(point)
			self.newpoints.append(point)
	
	def add_corner(self, point):
		if point not in self.corner_points:
			self.corner_points.append(point)
			
	def add_point_level(self, point, level):
		if point in self.maxlevels:
			self.maxlevels[point] = max(self.maxlevels[point], level)
		else:
			self.maxlevels[point] = level
	def get_point_level(self, point):
		return self.maxlevels[point]
		
	def evaluate(self, f):
		evals = 0
		for point in self.points:
			if point not in self.values:
				self.values[point] = f(point)
				evals += 1
		self.rootNode.update_grids()
		return evals
			
	def interpolate(self, point):
		return self.rootNode.interpolate(point)
	
	def optimize(self, n, max_level, f=numpy.exp, cores=1, info=True):
		self.newpoints = list()
		splitlist = []
		cores = 1
		if cores > 1:
			nodes = self.rootNode.childNodes
			@mab.parallelize.parallelize(cores=cores, info=info, flatten=True)
			def do(i):
				node = nodes[i]
				splitlist = []
				node.calculate_splitlist(splitlist)
				return splitlist
			splitlist = do(range(len(nodes)))
		else: 
			self.rootNode.calculate_splitlist(splitlist, 1.0)
		offset = max(self.values.values())
		def calcerror(node, values1, values2, values3, values4, volume):
			x1 = numpy.array(values1) - offset
			x2 = numpy.array(values2) - offset
			spliterror = max(abs(numpy.exp(x1)-numpy.exp(x2))) * volume # * numpy.log(volume)
			#spliterror = max(abs(numpy.exp(x1)*x1/x2)) * volume # * numpy.log(volume)
			if node.is_refined:
				return spliterror, True
			x1 = numpy.array(values3) - offset
			x2 = numpy.array(values4) - offset
			refineerror = max(abs(numpy.exp(x1)-numpy.exp(x2))) * volume # * numpy.log(volume)
			#refineerror = max(abs(numpy.exp(x1)*x1/x2)) * volume # * numpy.log(volume)
			if refineerror > spliterror:
				return refineerror, False
			else:
				return spliterror, True
		# only split level < max_level
		errorlist = [(node, calcerror(node, values1, values2, values3, values4, volume)) for node, values1, values2, values3, values4, volume in splitlist if node.level < max_level]
		errorlist.sort(key= lambda x: x[1][0]) # sort by error
		errorlist = errorlist[::-1] # increasing order
		i = 0
		while (len(self.newpoints) < n) and (i < len(errorlist)):
			split, node = errorlist[i][1][1], errorlist[i][0]
			node.split()
			if 0:
				if split: 
					print "split"
					node.split()
				else:
					print "refine"
					if node.is_refined:
						print "nice, already refined"
					#	node.split()
					#else:
					node.refine()
				
			self.maxlevels = {}
			self.rootNode.get_point_levels()
			self.rootNode.check_level_difference()
			i += 1
			#errorlist[-2][0].split(levels=2)
		#print "err", errorlist[-1][1]
		#print "err", errorlist[-2][1]
		#for split in splitlist:
		#	print split 
import basis	
class BinaryGridTri(BinaryGrid):
	def __init__(self, dimension, domain=None):
		baseslist = [basis.BasisTriLevels(0., 1., 0) for k in range(0,6)]
		BinaryGrid.__init__(self, dimension=dimension, baseslist=baseslist, domain=domain)
		 
	