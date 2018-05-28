# -*- coding: utf-8 -*-
import itertools
from numpy import *
import os
import mab.gd.logging as logging
import mab.gd.gdfast
import string
import mab.parallelize
#from numpy import nan
nan = numpy.nan
inf = numpy.inf
logger = logging.getLogger("gd.utils")
from numpy import *

def multisum(a, axes):
	correction = 0
	for axis in axes:
		a = numpy.sum(a, axis=axis-correction)
		correction += 1
	return a

class ParameterEstimates(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		
	def run(self, args, opts, scope):
		self.parameterset.load()
		p = numpy.exp(self.parameterset.probability_grid)
		for i, param in enumerate(self.parameterset.parameter_range_list):
			print param.name
			axes = range(self.parameterset.dimension)
			#print axes
			axes.remove(i)
			#print axes
			p1d = multisum(p, axes=axes)
			
			values = numpy.mgrid[param.min:param.max:1j*len(p)]
			#print values
			p1d /= sum(p1d)
			mean = sum(p1d*values)
			variance = sum((values-mean)**2 * p1d)
			sigma = variance**0.5
			print "%g" % mean, "+/-", sigma, "=", (param.f(sigma)-1)*100, "%"
			print "%g" % param.f(mean), "+", "%g" % (param.f(mean+sigma)-param.f(mean)), "-", "%g" % (-(param.f(mean-sigma)-param.f(mean)))

class ParameterValue(object):
	def __init__(self, values, parameter_range_list):
		self.values = values
		self.parameter_range_list = parameter_range_list
		self.values_uniform = [p.scale_to_uniform(v) for v, p in zip(self.values, parameter_range_list)]
		self.values_org = [p.finv(v) for v, p in zip(self.values, parameter_range_list)]
		name_parts = []
		param_names = []
		#self.values = []
		for value, parameter_range in zip(self.values, parameter_range_list):
			#print parameter_range.name, value
			name_parts.append(parameter_range.format_value(value))
			param_names.append(parameter_range.name)
			#self.values.append(value)
		self.name = "/".join(name_parts)
		self.dirname = self.name
		self.items = zip(param_names, self.values)
		self.dict = dict(self.items)
		self.p_prior = 1.

class ParameterRange(object):
	def __init__(self, min, max, name, format="%r", f=lambda x: x, finv=lambda x: x, prior=lambda x: 1, label=None, units=""):
		self.min = float(min)
		self.max = float(max)
		self.name = name
		self.format = format
		self.f = f
		self.finv = finv
		self.label = label
		self.units = units
		self.prior = prior
		#self.values = arange(n) / (n-1.) * (max-min) + min
		
	def scale_from_uniform(self, u):
		return self.f(u * (self.max - self.min) + self.min)
		
	def scale_to_uniform(self, v):
		return (self.finv(v) - self.min) / (self.max - self.min)
		
	def scale_from_uniform_to_range(self, u):
		return u * (self.max - self.min) + self.min
		
	def format_value(self, value):
		return self.format % value
		
class ParameterSetBase(object):
	def __init__(self, modelpath, type, name, parameter_range_list, use_value_directories=True, gridsize=30, postfix=""):
		self.modelpath = modelpath
		self.type = type
		self.name = name
		self.parameter_range_list = parameter_range_list
		self.parameter_ranges = [(p.min, p.max) for p in self.parameter_range_list]
		self.parameter_labels = [p.label for p in self.parameter_range_list] 
		self.use_value_directories = use_value_directories
		self.dirname = os.path.join(self.modelpath, self.type, self.name)
		self.parameter_values = []
		self.parameter_values_new = [] # only use when is_iterative == True
		self.dimension = len(self.parameter_range_list)
		self.gridsize = gridsize
		self.postfix = postfix
		
	def range1d(self, index):
		return self.parameter_range_list[index].min, self.parameter_range_list[index].max
			
	def range2d(self, index1, index2):
		a = self.parameter_range_list[index1].min, self.parameter_range_list[index2].min
		b = self.parameter_range_list[index1].max, self.parameter_range_list[index2].max
		return a,b
			
	def run(self, args, opts, scope):
		self.init()
		
	def init(self):
		self.create_dir()
		if self.use_value_directories:
			self.create_value_dirs()
		
	def create_dir(self):
		if not os.path.exists(self.dirname):
			logger.info("creating directory %s" % self.dirname)
			os.makedirs(self.dirname)
		else:
			logger.info("directory %s already exists" % self.dirname)
	
	def is_iterative(self):
		return False
	
	def get_new_parameter_values(self):
		if self.is_iterative():
			logger.info("iterative method")
			parameter_values = self.parameter_values_new
		else:
			logger.info("non-iterative method")
			parameter_values = self.parameter_values
		return parameter_values
	
	def create_true_parameter(self, scope):
		dirname = os.path.join(self.dirname, "true")
		values = {}
		for pr in self.parameter_range_list:
			name = pr.name
			value = scope[name]
			values[name] = value
			print name, value
		self._create_parameter(dirname, values, "true")
		
	def create_value_dirs(self):
		parameter_values = self.get_new_parameter_values()
		logger.info("creating parameters.ini files for %d models" % len(parameter_values))
		for parameter_value in parameter_values:
			dirname = os.path.join(self.dirname, parameter_value.name)
			values = {}
			for name, value in parameter_value.items:
				values[name] = value
			self._create_parameter(dirname, values, parameter_value.name)
			if 0:
				dirname = os.path.join(self.dirname, parameter_value.name)
				if not os.path.exists(dirname):
					#logger.debug("creating directory %s" % dirname)
					os.makedirs(dirname)
				else:
					#logger.debug("directory %s already exists" % dirname)
					pass
				filename = os.path.join(dirname, "parameters.ini")
				f = open(filename, "w")
				print >>f, "[parameters]"
				for name, value in parameter_value.items:
					print >>f, "%s=%r" % (name, value)
				print >>f, "\n[globals]"
				print >>f, "model_id=%r" % parameter_value.name
				print >>f, "model_id_path=os.path.join(model_set_path, model_id)" 
				f.close()
				#logger.debug("file %s written" % filename)
				
				for name in ["log"]:
					dirname = os.path.join(self.dirname, parameter_value.name, name)
					if not os.path.exists(dirname):
						#logger.debug("creating directory %s" % dirname)
						os.makedirs(dirname)
					else:
						#logger.debug("directory %s already exists" % dirname)
						pass
					
	def _create_parameter(self, dirname, values, name):
		if not os.path.exists(dirname):
			logger.debug("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			#logger.debug("directory %s already exists" % dirname)
			pass
		filename = os.path.join(dirname, "parameters.ini")
		f = open(filename, "w")
		print >>f, "[parameters]"
		for _name, value in values.items(): #parameter_value.items:
			print >>f, "%s=%r" % (_name, value)
		print >>f, "\n[globals]"
		print >>f, "model_id=%r" % name
		print >>f, "model_id_path=os.path.join(model_set_path, model_id)" 
		f.close()
		#logger.info("file %s written" % filename)
		
		for subdirname in ["log"]:
			dirname = os.path.join(self.dirname, name, subdirname)
			if not os.path.exists(dirname):
				logger.debug("creating directory %s" % dirname)
				os.makedirs(dirname)
			else:
				#logger.debug("directory %s already exists" % dirname)
				pass
				
		#logger.info("%d models in total" % len(self.parameter_values))
				
	def collect(self):
		self.valuemap = {}
		n = len(self.parameter_values)
		progressbar = mab.utils.progressbar.ProgressBar(0, n)
		progressbar.update(0)
		i = 0
		for parametervalue in self.parameter_values:
			modelname = parametervalue.name
			dirname = os.path.join(self.dirname, modelname)
			filename = os.path.join(dirname, "results", "logprobability_seperate" +self.postfix +".txt")
			if os.path.exists(filename):
				data = file(filename).read()
				s = str(data)
				#print "s = %r" % s
				#print eval(s)
				values = dict(eval(data))
				logp = values["kin"]# + values["m4"]
				#print logp
				if numpy.isnan(logp):
					logp = -1e10
				if numpy.isinf(logp):
					logp = -1e10
				logger.debug("result: %s = %f" % (filename, logp))
				self.valuemap[tuple(parametervalue.values)] = logp
			else:
				logp = -1e9
				logger.info("not found: %s, using = %f" % (filename, logp))
				self.valuemap[tuple(parametervalue.values)] = logp
					
			#print dirname, logp
			i += 1
			if (i % 10) == 0:
				progressbar.update(i)
		progressbar.update(n)
		values = self.valuemap.values()
		values = numpy.array(values)
		indices = numpy.argsort(values)
		#print values[indices]
		i = indices[-1]
		#print self.parameter_values[i].dirname
		#self.grid()
			#logger.info("fininshed %d out of %d" % (i, len(self.parameterset.parameter_values)))
			#progressbar.update(i)
			
			
class ParameterSetRegular(ParameterSetBase):
	def __init__(self, modelpath, type, name, parameter_range_list, n=2, use_value_directories=True, gridsize=30):
		ParameterSetBase.__init__(self, modelpath, type, name, parameter_range_list, use_value_directories, gridsize=gridsize)
		self.n = n
		uniform_range = numpy.arange(n)/(n-1.)
		ranges = [parameter_range.scale_from_uniform(uniform_range) for parameter_range in parameter_range_list]
		for values in itertools.product(*ranges):
			parameter_value = ParameterValue(values, parameter_range_list)
			self.parameter_values.append(parameter_value)
			#print parameter_value.name, parameter_value.values
			
	def makegrid(self, ps):
		indices = range(self.n)
		grid = numpy.zeros((self.n,) * self.dimension)
		for p, multi_index in zip(ps, itertools.product(*[indices for i in range(self.dimension)])):
			#print multi_index, p
			grid.__setitem__(multi_index, p)
		return grid
		
		
		
		
class CommandCreateTrue(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		
	def run(self, args, opts, scope):
		logger.info("making 'true' model")
		self.parameterset.create_true_parameter(scope)
		
class CommandParameterSetTreeInit(object):
	def __init__(self, tree):
		self.tree = tree
		
	def run(self, args, opts, scope):
		logger.info("running tree initialization")
		self.tree.init()
		
class CommandParameterSetTreeOptimize(object):
	def __init__(self, tree, iteration):
		self.tree = tree
		self.iteration = iteration
		
	def run(self, args, opts, scope):
		logger.info("running tree optimization")
		self.tree.optimize(self.iteration)
		
class ParameterSetTree(ParameterSetBase):
	def __init__(self, modelpath, type, name, parameter_range_list, tree, iteration=0, use_value_directories=True, gridsize=30, postfix=""):
		ParameterSetBase.__init__(self, modelpath, type, name, parameter_range_list, use_value_directories, gridsize=gridsize, postfix=postfix)
		self.tree = tree
		self.valuemap = {}
		self.iteration = iteration
		self.tree.read(iteration=self.iteration)
		for point in self.tree.points:
			p = self._parametervalue_from_uniform(*point)
			self.parameter_values.append(p)
		for point in self.tree.newpoints:
			p = self._parametervalue_from_uniform(*point)
			self.parameter_values_new.append(p)
		logger.info("total models: %d" % len(self.parameter_values))
		logger.info("new models  : %d" % len(self.parameter_values_new))
		self.filename_grid  = os.path.join(self.dirname, "probability_grid" + self.postfix + "-it%03d" % self.iteration +".npy")

		#self.sparsegrid = mab.gd.gdfast.SparseGridD(self.dimension, self.level)
		
		#def f(*uniform_values):
		#	values = [parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)]
		#	#values = tuple([parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)])
		#	parameter_value = ParameterValue(values, self.parameter_range_list)
		#	self.parameter_values.append(parameter_value)
		#	#print parameter_value.name, parameter_value.values
		#	return 0
		#self.tree.eval(f, iteration=self.iteration)
		
	def get_probability(self, parameter_value):
		#return 1
		#print self.valuemap.keys()
		#print self.valuemap[tuple(parameter_value.values)]
		x = self.valuemap.values()
		#print max(x), x
		return numpy.exp(self.valuemap[tuple(parameter_value.values)]-max(self.valuemap.values()))
		
	def get_volume(self, parameter_value):
		return 1; #parameter_value.valuemap[parameter_value.values]
		return parameter_value.volumes[parameter_value.values]
		
	def _parametervalue_from_uniform(self, *uniform_values):
		values = [parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)]
		return ParameterValue(values, self.parameter_range_list)
		
	#def init(self):
	#	self.tree.init()
	def is_iterative(self):
		return True
		
	def optimize(self):
		#print self.parameter_values
		self.load()
		iteration = scope["iteration"]
		self.tree.update(self.f_uniform, iteration=iteration)
		self.tree.optimize(iteration)
		
	def load(self):
		self.tree.read()
		filename = os.path.join(self.dirname, "probability_grid" + self.postfix +".npy")
		logger.info("loading probability_grid from filename %s" % self.filename_grid)
		self.probability_grid = numpy.load(self.filename_grid)
		#print self.probability_grid.shape
		
		filename = os.path.join(self.dirname, "valuegrid" + self.postfix +".npy")
		logger.info("loading valuemap from filename %s" % filename)
		self.valuegrid = numpy.load(filename)
		#print list(self.valuegrid[:,0:self.dimension])
		keys = [tuple(k) for k in self.valuegrid[:,0:self.dimension]]
		values = self.valuegrid[:,self.dimension]
		self.valuemap = dict(zip(keys, values))
		#print numpy.array(keys).shape
		d1 = [k[1] for k in keys]
		d1.sort()
		#print d1
		 
		logger.info("sparse map contains %d points" % len(keys))
		
	def f_uniform(self, *uniform_values):
		values = tuple([parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)])
		values_range = tuple([parameter_range.scale_from_uniform_to_range(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)])
		#print values, self.valuemap[values]
		rlist = self.parameter_range_list
		prior = numpy.prod([rlist[k].prior(values_range[k]) for k in range(len(rlist))])
		#print self.valuemap[values]

		return self.valuemap[values] + numpy.log(prior)

		
	def grid(self):
		#print self.valuemap
		#self.tree.write_solved(self.f_uniform, iteration=self.iteration)
		self.tree.update(self.f_uniform, iteration=self.iteration)
		self.probability_grid = self.tree.grid(iteration=self.iteration)
		#filename = os.path.join(self.dirname, "probability_grid" + self.postfix +".npy")
		logger.info("saving probability_grid to filename %s" % self.filename_grid)
		numpy.save(self.filename_grid, self.probability_grid)
		
		filename = os.path.join(self.dirname, "valuegrid" + self.postfix +".npy")
		N = len(self.valuemap)
		valuegrid = numpy.zeros((N, self.dimension+1))
		valuegrid[:,0:self.dimension] = self.valuemap.keys()
		valuegrid[:,self.dimension] = self.valuemap.values()
		
		logger.info("saving valuegrid to filename %s" % filename)
		numpy.save(filename, valuegrid)
		
	def evidence(self, axis=None, value=None):
		logger.info("loading probability_grid to filename %s" % self.filename_grid)
		grid = numpy.load(self.filename_grid)
		shape = grid.shape
		N = len(self.parameter_range_list)
		deltaXs = [prl.max-prl.min for prl in self.parameter_range_list]
		#dVs = [deltaXs[k]/float(deltaXs[k]) for k in range(N)]
		dVs = [deltaXs[k]/float(shape[k]) for k in range(N)]
		if axis is not None:
			print dVs
			dVs.pop(axis)
			print dVs
			dV = numpy.product(dVs)
			prl = self.parameter_range_list[axis]
			index = (value - prl.min) / deltaXs[axis] * shape[axis]
			print index, shape
			grid = grid[:,:,index] # TODO this is hardcoded to be 3d and use the last axis
			evidence = self.tree.normalization(self.f_uniform, iteration=self.iteration) * numpy.exp(grid).sum() * dV
			if 0:
				from kaplot import box, draw, indexedimage
				grid = numpy.exp(grid)
				print grid
				box()
				print grid.shape
				indexedimage(grid.T)
				draw()
		else:
			dV = numpy.product(dVs)
			print "max", grid.max(),self.tree.normalization(self.f_uniform, iteration=self.iteration), numpy.exp(grid).sum()
			logevidence = self.tree.lognormalization(self.f_uniform, iteration=self.iteration) + numpy.log(numpy.exp(grid).sum() * dV)
			evidence = numpy.exp(self.tree.lognormalization(self.f_uniform, iteration=self.iteration) + numpy.log(numpy.exp(grid).sum())) * dV
		print "N", self.tree.lognormalization(self.f_uniform, iteration=self.iteration), numpy.exp(grid).sum(), dV
		print "-->", self.tree.lognormalization(self.f_uniform, iteration=self.iteration) + numpy.log(numpy.exp(grid).sum())
		print dVs, dV
		print "evidence", evidence, dV, dVs
		print "logevidence", logevidence
		filename = os.path.join(self.dirname, "evidence.txt")
		f = file(filename, "wb")
		f.write(repr(evidence))
		f.close()
		
class Evidence(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		
	def run(self, args, opts, scope):
		print args
		if len(args) == 3:
			self.parameterset.load()
			self.parameterset.evidence(int(args[1]), float(args[2]))
			return 2
		else:
			self.parameterset.load()
			self.parameterset.evidence()
		

class ParameterSetSparse(ParameterSetBase):
	def __init__(self, modelpath, type, name, parameter_range_list, level=2, gridsize=30, use_value_directories=True, postfix=""):
		ParameterSetBase.__init__(self, modelpath, type, name, parameter_range_list, use_value_directories, gridsize, postfix=postfix)
		self.level = level
		self.sparsegrid = mab.gd.gdfast.SparseGridD(self.dimension, self.level)
		
		def f(*uniform_values):
			values = [parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)]
			parameter_value = ParameterValue(values, self.parameter_range_list)
			self.parameter_values.append(parameter_value)
			#print parameter_value.name, parameter_value.values
			return 0
		self.sparsegrid.eval(f)
		
	def compute_weights(self):
		for pv in self.parameter_values:
			pv.weight_dim = [None] * self.dimension
			pv.left = [None] * self.dimension
			pv.right = [None] * self.dimension
			#print pv.weight_dim
		for dim in range(self.dimension):
			groups = dict()
			for pv in self.parameter_values:
				key = pv.values_org[0:dim] + pv.values_org[dim+1:]
				print key
				key = tuple(key)
				#print len(key)
				assert len(key) == self.dimension-1
				if key not in groups:
					print "new group", key
					groups[key] = []
				groups[key].append(pv)
			print "grouping done"
			total = 0
			for key, group in groups.items():
				print len(group), "+",
				total += len(group)
			print "=", total
			print groups.keys()
			for key, group in groups.items():
				print ">>>", dim, key, len(group)
				print "   values ", [pv.values_org for pv in group]
				#print group[0].parameter_range_list
				#if len(group) == 1:
				#	group[0].weight_dim[dim] = group[0].parameter_range_list[dim].max - group[0].parameter_range_list[dim].min
				#	group[0].left[dim] = self.parameter_range_list[dim].min
				#	group[0].right[dim] = self.parameter_range_list[dim].max
				#else:
				if 1:
					xs = numpy.array([pv.values_org[dim] for pv in group])
					indices = numpy.argsort(xs)
					#print xs.shape
					#print indices
					#print len(
					xs = xs[indices]
					print xs
					group = [group[i] for i in indices]
					
					#dx = (xs[1:] - xs[0:-1])/2
					#dx = dx + [dx[-1]] # repeat last value
					
					#scale = 0.8
					scale = 1.
					dx = []
					#mid = (xs[1] + xs[0])/2
					#dx_left = mid - self.parameter_range_list[dim].min
					dx_right = dx_left = (self.parameter_range_list[dim].max-self.parameter_range_list[dim].min)/9
					#dx.append(dx_left)
					#dx = 
					#print "left dx", dx[-1]
					#group[0].left[dim] = self.parameter_range_list[dim].min + scale*dx_left
					#group[0].right[dim] = mid - scale*dx_left
					
					#mids = (xs[0:-1]+xs[1:])/2
					#dxmids = mids[1:] - mids[0:-1]
					#for i in range(1,len(group)-1):
					#	group[i].left[dim] = mids[i-1] + scale * dxmids[i-1]
					#	group[i].right[dim] = mids[i] - scale * dxmids[i-1]
					
					#print "mids", mids, dxmids
					#dx.extend(dxmids)
					#dx.extend([dx_left] * (len(group)-2))
					
					#mid = (xs[-2] + xs[-1])/2
					#dx_right = self.parameter_range_list[dim].max - mid
					#dx.append(dx_right)
					#group[-1].left[dim] = mid + scale * dx_right
					#group[-1].right[dim] = self.parameter_range_list[dim].max - scale*dx_right
					#group[-1].left[dim] = xs[-1] - dx_left/2
					#group[-1].left[dim] = xs[-1] - dx_left/2
					
					
					#Dx = self.parameter_range_list[dim].max - self.parameter_range_list[dim].min
					#print "right dx", dx[-1] 
					#print "!!!!!!!!", sum(dx), Dx, (sum(dx)-Dx)
					#assert abs(sum(dx)-Dx) < 1e-10
					#print dx
					scale = 0.45
					for i in range(len(group)): #pv in group:
						#group[i].weight_dim[dim] = dx[i]
						group[i].left[dim] = max(self.parameter_range_list[dim].min, xs[i]-dx_left*scale)
						group[i].right[dim] = min(self.parameter_range_list[dim].max, xs[i]+dx_left*scale)
					for i in range(len(dx)): #pv in group:
						group[i].weight_dim[dim] = dx[i]
						print group[i].left, group[i].right
		#dsa
		for pv in self.parameter_values:
			#print "values", `pv.values_org`
			#print "values", tuple(pv.values_org)
			#print "weights", pv.weight_dim
			#pv.weight = numpy.prod(pv.weight_dim)
			#print "..", pv.weight
			print pv.left, pv.right
				#for pv in group:
				#	xs = pv
					
			
			
			#group = dict([(pv.values_org[0:i+1,i+1:], pv) for pv in self.parameter_values])
		
	def grid(self):
		def f(*uniform_values):
			values = tuple([parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)])
			return self.valuemap[values]
		self.sparsegrid.eval(f)
		print "eval done"
		w = self.gridsize
		
		#self.parameter_points = transpose(array(self.parameterset.points))
		
		#mab.gd.gdfast.ProfilerStart("sparse")
		self.probability_grid = mab.utils.numpy.mmapzeros((w,) * self.dimension)
		make_linear = False
		if make_linear:
			self.parameterset.make_linear()
		if self.dimension == 1:
			for i in range(w):
				self.probability_grid[i] = self.sparsegrid(i/(w-1.0)) 
		if self.dimension == 2:
			for i in range(w):
				#x1, x2 = self.parameterset.ranges_org[0]
				#x = x1 + (x2-x1) * i / (w-1.)
				#x = self.parameter_range_list[0].scale_from_uniform(i/(w-1.0)) 
				for j in range(w):
					#y1, y2 = self.parameterset.ranges_org[1]
					#y = y1 + (y2-y1) * j / (w-1.)
					#y = self.parameter_range_list[1].scale_from_uniform(j/(w-1.0))
					#print self.parameterset((x, y))
					self.probability_grid[i,j] = self.sparsegrid(*(i/(w-1.0), j/(w-1.0))) 
		if self.dimension == 3:
			x1, x2 = 0., 1.#self.parameterset.ranges_org[0]
			y1, y2 = 0., 1.#self.parameterset.ranges_org[1]
			z1, z2 = 0., 1.#self.parameterset.ranges_org[2]
			logger.info("interpolating grid")
			@mab.parallelize.parallelize(cores=10, info=True)
			def do(i):
				#x = x1 + (x2-x1) * i / (w-1.) 
				for j in range(w):
					#y = y1 + (y2-y1) * j / (w-1.)
					for k in range(w):
						#z = z1 + (z2-z1) * k / (w-1.)
						self.probability_grid[i,j,k] = self.sparsegrid(*(i/(w-1.0), j/(w-1.0), k/(w-1.)))
			xindices = range(w)
			do(xindices)
		if self.dimension == 4:
			x1, x2 = 0., 1.#self.parameterset.ranges_org[0]
			y1, y2 = 0., 1.#self.parameterset.ranges_org[1]
			z1, z2 = 0., 1.#self.parameterset.ranges_org[2]
			logger.info("interpolating grid")
			@mab.parallelize.parallelize(cores=48, info=True)
			def do(i):
				#x = x1 + (x2-x1) * i / (w-1.) 
				for j in range(w):
					#y = y1 + (y2-y1) * j / (w-1.)
					for k in range(w):
						#z = z1 + (z2-z1) * k / (w-1.)
						for l in range(w):
							self.probability_grid[i,j,k,l] = self.sparsegrid(*(i/(w-1.0), j/(w-1.0), k/(w-1.), l/(w-1.)))
			xindices = range(w)
			do(xindices)
		if self.dimension == 5:
			x1, x2 = 0., 1.#self.parameterset.ranges_org[0]
			y1, y2 = 0., 1.#self.parameterset.ranges_org[1]
			z1, z2 = 0., 1.#self.parameterset.ranges_org[2]
			logger.info("interpolating grid")
			@mab.parallelize.parallelize(cores=48, info=True)
			def do(i):
				#x = x1 + (x2-x1) * i / (w-1.) 
				for j in range(w):
					#y = y1 + (y2-y1) * j / (w-1.)
					for k in range(w):
						#z = z1 + (z2-z1) * k / (w-1.)
						for l in range(w):
							for m in range(w):
								self.probability_grid[i,j,k,l,m] = self.sparsegrid(*(i/(w-1.0), j/(w-1.0), k/(w-1.), l/(w-1.), m/(w-1.)))
			xindices = range(w)
			do(xindices)
		#print self.probability_grid
		logger.info("probability grid (min: %f max: %f avg: %f)" % (self.probability_grid.min(), self.probability_grid.max(), self.probability_grid.mean())) 
		self.probability_grid -= self.probability_grid.max()
		#self.probability_grid = numpy.exp(self.probability_grid)
		filename = self.filename_grid
		logger.info("saving probability_grid to filename %s" % filename)
		numpy.save(filename, self.probability_grid)
		
		filename = os.path.join(self.dirname, "valuegrid" + self.postfix +".npy")
		N = len(self.valuemap)
		valuegrid = numpy.zeros((N, self.dimension+1))
		valuegrid[:,0:self.dimension] = self.valuemap.keys()
		valuegrid[:,self.dimension] = self.valuemap.values()
		
		logger.info("saving valuegrid to filename %s" % filename)
		numpy.save(filename, valuegrid)
		
	def load(self):
		filename = os.path.join(self.dirname, "probability_grid" + self.postfix +".npy")
		logger.info("loading probability_grid from filename %s" % filename)
		self.probability_grid = numpy.load(filename)
		
		filename = os.path.join(self.dirname, "valuegrid" + self.postfix +".npy")
		logger.info("loading valuemap from filename %s" % filename)
		self.valuegrid = numpy.load(filename)
		#print list(self.valuegrid[:,0:self.dimension])
		keys = [tuple(k) for k in self.valuegrid[:,0:self.dimension]]
		values = self.valuegrid[:,self.dimension]
		self.valuemap = dict(zip(keys, values))
		logger.info("sparse map contains %d points" % len(keys))
		#print self.valuemap
		#import pdb; pdb.set_trace()
		
		def f(*uniform_values):
			values = tuple([parameter_range.scale_from_uniform(uniform_value) for uniform_value, parameter_range in zip(uniform_values, self.parameter_range_list)])
			#print values
			return self.valuemap[values]
		self.sparsegrid.eval(f)
				
	def logp_uniform(self, *values):
		return self.sparsegrid(*values)
	
	def export_mathematica_marginalized(self, filename, axes):
		def multisum(a, axes):
			correction = 0
			for axis in axes:
				a = numpy.sum(a, axis=axis-correction)
				correction += 1
			return a
		allaxes = range(len(self.probability_grid.shape))
		remove_axes = [] 
		for i in allaxes:
			if i not in axes:
				remove_axes.append(i)
		p = numpy.exp(self.probability_grid)
		p = multisum(p, remove_axes)
		dimension = len(axes)
		w = p.shape[0]
		f = open(filename, "w")
		logger.info("filename: %s" % os.path.abspath(filename))
		if dimension == 3:
			x1, x2 = 0., 1.#self.parameterset.ranges_org[0]
			y1, y2 = 0., 1.#self.parameterset.ranges_org[1]
			z1, z2 = 0., 1.#self.parameterset.ranges_org[2]
			#f.write("{")
			for i in range(w):
				x = self.parameter_range_list[axes[0]].min + (self.parameter_range_list[axes[0]].max-self.parameter_range_list[axes[0]].min) * i/(w-1.0)
				#print x
				#f.write("{")
				for j in range(w):
					y = self.parameter_range_list[axes[1]].min + (self.parameter_range_list[axes[1]].max-self.parameter_range_list[axes[2]].min) * j/(w-1.0)
					#f.write("{")
					for k in range(w):
						z = self.parameter_range_list[axes[2]].min + (self.parameter_range_list[axes[2]].max-self.parameter_range_list[axes[2]].min) * k/(w-1.0)
						f.write("%g %g %g %g\n" % (x, y, z, numpy.log(p[i,j,k])))
						#if (k != (w-1)):
						#	f.write(",")
						#if (i != (w-1)) and (j != (w-1)) and (k != (w-1)):
						#	f.write(",")
					#f.write("}")
					#if (j != (w-1)):
					#	f.write(",")
				#f.write("\n")
				#f.write("}")
				#if (i != (w-1)):
				#	f.write(",")
			#f.write("}")
			f.close()
		else:
			raise "unsupported dimension: %d" % dimension
		
	def _run(self, args, opts, scope):
		filename = args[1]
		axes = eval(args[2])
		self.load()
		self.export_mathematica_marginalized(filename, axes)
		#for axis in axes:

class RunAll(object):
	def __init__(self, parameterset, command, only_new=True):
		self.parameterset = parameterset
		self.command = command
		self.dirname = self.parameterset.dirname
		self.results = {}
		self.only_new = only_new
		
	def run(self, args, opts, scope):
		#values = self.parameterset.parameter_values
		#options = scope.options
		#arguments = scope.arguments
		#progressbar = mab.utils.progressbar.ProgressBar(0, len(self.parameterset.parameter_values))
		#progressbar.update(0)
		#i = 0
		if self.only_new:
			parameter_values = self.parameterset.get_new_parameter_values()
		else:
			parameter_values = self.parameterset.parameter_values
		#@mab.parallelize.parallelize(cores=scope["cores"], info=True)
		#@mab.parallelize.parallelize(cores=4, info=True)
		@mab.parallelize.parallelize(cores=scope["cores"], info=True)
		def do(i):
		#for parametervalue in self.parameterset.parameter_values:
			parametervalue = parameter_values[i]
			modelname = parametervalue.name
			if self.parameterset.use_value_directories:
				dirname = scope["model_set_path"]
				filename = os.path.join(dirname, modelname, "parameters.ini")
				#print filename
				#print filename
				assert os.path.exists(filename), "missing file: %s" % filename
				scope.reset()
				scope.re_readfiles(filename)
				scope.init()
				scope["cores"] = 2
				cmd = scope[self.command]
				scope["cores"] = 2
				scope["progressbar"] = False
				try:
					cmd.run(args, opts, scope)
				except:
					logger.error("error, modelpath %s model_id_path %s" % (opts.modelpath, scope["model_id_path"]))
					raise
				
				#dirname = os.path.join(self.dirname, modelname)
				#model_filename = os.path.join(dirname, opts.type+".ini")
				#model_filename = scope["model_set_path"]
				#assert os.path.exists(model_filename)
				#if 0:
				#	if not os.path.exists(model_filename):
				#		print dirname
				#	else:
				#		continue
				#print model_filename
				#filenames.append(filename)
				#logger.info("modelname: %s" % modelname)
				#scope["options"].modelname = modelname
			if 0:
				scope.reset()
				scope.readfiles(model_filename)
				scope.init()
				key = ""
				if not self.parameterset.use_value_directories:
					for name, value in parametervalue.items:
						#print name, "=", value
						scope[name] = value
						key += "%s=%r" % (name, value)
				command = scope[self.command]
				command.run(args, opts, scope)
				result = scope["result"]
				self.results[key] = results
				i += 1
				#logger.info("fininshed %d out of %d" % (i, len(self.parameterset.parameter_values)))
				#progressbar.update(i)
			#progressbar.update(i+1)
		do(range(len(parameter_values)))
		
class ParameterGridIterator(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		self.dirname = self.parameterset.dirname
		
	def scopelist(self, scope, all=False):
		cleanscope = scope.clone()
		if all:
			parameter_values = self.parameterset.parameter_values
		else:
			parameter_values  = self.parameterset.get_new_parameter_values()
		logger.info("iterating over %d models" % len(parameter_values))
		scopes = []
		for i in range(len(parameter_values)):
			scope = cleanscope.clone()
			#print scope.dict
			parametervalue = self.parameterset.parameter_values[i]
			modelname = parametervalue.name
			#scope.reset()
			#print self.parameterset.use_value_directories
			if self.parameterset.use_value_directories:
				dirname = scope["model_set_path"]
				filename = os.path.join(dirname, modelname, "parameters.ini")
				scope.allfiles.append(filename)
				#filenames = [filename]
				#if opts.include:
				#	filenames.extend(opts.include.split(";"))
				#print filenames
				assert os.path.exists(filename)
				#scope.reset()
				#scope.re_readfiles(*filenames)
				#scope.init()
				#cmd = scope[self.command]
				#scope["progressbar"] = False
				#cmd.run(args, opts, scope)
			else:
				scope.re_readfiles()
				scope.init()
				for name, value in parametervalue.items:
					scope[name] = value
			scope.defaults["probability"] = self.parameterset.get_probability(parametervalue)
			scope.defaults["volume"] = self.parameterset.get_volume(parametervalue)
			scopes.append(scope)
		return scopes
			
			
	def iter(self, scope, all=False):
		cleanscope = scope.clone()
		#print len(self.parameterset.parameter_values)
		#print scope["model_set_path"]
		if all:
			parameter_values = self.parameterset.parameter_values
		else:
			parameter_values  = self.parameterset.get_new_parameter_values()
		logger.info("iterating over %d models" % len(parameter_values))
		for i in range(len(parameter_values)):
			#print i
			#print cleanscope.dict
			#print cleanscope.clone().dict
			#dsa
			scope = cleanscope.clone()
			#print scope.dict
			parametervalue = self.parameterset.parameter_values[i]
			modelname = parametervalue.name
			#scope.reset()
			#print self.parameterset.use_value_directories
			if self.parameterset.use_value_directories:
				dirname = scope["model_set_path"]
				filename = os.path.join(dirname, modelname, "parameters.ini")
				filenames = [filename]
				#if opts.include:
				#	filenames.extend(opts.include.split(";"))
				#print filenames
				assert os.path.exists(filename)
				#scope.reset()
				scope.re_readfiles(*filenames)
				scope.init()
				#cmd = scope[self.command]
				#scope["progressbar"] = False
				#cmd.run(args, opts, scope)
			else:
				scope.re_readfiles()
				scope.init()
				for name, value in parametervalue.items:
					scope[name] = value
			scope["probability"] = self.parameterset.get_probability(parametervalue)
			scope["volume"] = self.parameterset.get_volume(parametervalue)
				
			yield scope
				
		
class RunAllCached(RunAll):
	def __init__(self, parameterset, command, filename_cache="cache"):
		RunAll.__init__(self, parameterset, command)
		self.filename_cache = filename_cache
		self.results = {}
		
	def run(self, args, opts, scope):
		self.dirname
		#filename_cache = 
		#values = self.parameterset.parameter_values
		#options = scope.options
		#arguments = scope.arguments
		parameter_values = self.get_new_parameter_values()
		progressbar = mab.utils.progressbar.ProgressBar(0, len(parameter_values))
		#progressbar.update(0)
		i = 0
		assert not self.parameterset.use_value_directories
		for parametervalue in parameter_values:
			modelname = parametervalue.name
			dirname = os.path.join(self.dirname, modelname)
			scope.reset()
			scope.readfiles()
			scope.init()
			key = ""
			if not self.parameterset.use_value_directories:
				for name, value in parametervalue.items:
					#print name, "=", value
					scope[name] = value
					key += "%s=%r" % (name, value)
			if key not in self.results:
				command = scope[self.command]
				command.run(args, opts, scope)
				result = scope["result"]
				logger.info("running:  %s" % key)
				self.results[key] = results
			else:
				logger.info("skipping: %s" % key)
			i += 1
			logger.info("fininshed %d out of %d" % (i, len(self.parameterset.parameter_values)))
			progressbar.update(i)
		#progressbar.update(i+1)		

class CreateCondorBatchFileCommand(object):
	def __init__(self, parameterset, command, condor_head_template, condor_item_template, rootdir, condor_requirements, filename="condor.batch", use_dagman=False, dagman_filename=None):
		self.parameterset = parameterset
		self.command = command
		self.condor_head_template = condor_head_template
		self.condor_item_template = condor_item_template
		self.rootdir = rootdir
		self.condor_requirements = condor_requirements
		self.dirname = self.parameterset.dirname
		self.filename = filename
		self.use_dagman = use_dagman
		self.dagman_filename = dagman_filename
		assert self.parameterset.use_value_directories
		
	def run(self, args, opts, scope):
		#values = self.parameterset.parameter_values
		#options = scope.options
		#arguments = scope.arguments
		condor_head_template = open(self.condor_head_template).read()
		condor_item_template = open(self.condor_item_template).read()
		model_set_path = scope["model_set_path"]
		filename = os.path.join(model_set_path, "condor", self.filename)
		parameter_values = self.parameterset.get_new_parameter_values()
		if not self.use_dagman:
			logger.info("creating condor file: %s with %d jobs" % (filename, len(parameter_values)))
		print "to submit:"
		print "condor_submit %s" % filename
		#progressbar.update(0)
		params = {}
		params["rootdir"] = self.rootdir
		params["configurations"] = ";".join(scope.configurations)
		params["globals"] = opts.globals
		params["modelpath"] = self.parameterset.modelpath
		params["model_set"] = self.parameterset.name
		params["model_type"] = scope["model_type"]
		#params["modelset"] = self.parameterset.name
		#params["type"] = self.parameterset.type
		params["condor_requirements"] = self.condor_requirements
		params["command"] = self.command
		header = string.Template(condor_head_template).substitute(**params)
		if self.use_dagman:
			pass
		else:
			f = open(filename, "w")
			f.write(header)
		i = 0
		progressbar = mab.utils.progressbar.ProgressBar(0, len(parameter_values))
		filenames = []
		for parametervalue in parameter_values:#[-2:-1]:
			modelname = parametervalue.name
			params["model_id"] = modelname
			dirname = os.path.join(self.dirname, modelname)
			model_filename = os.path.join(dirname, "parameters.ini")
			params["include"] = opts.include +";" + model_filename
			assert os.path.exists(model_filename)
			#data = condor_item_template % params
			if self.use_dagman:
				#model_set_path = scope["model_set_path"]
				#dirname = os.path.join(model_set_path, 
				filename = os.path.join(dirname, self.filename)
				filenames.append(filename)
				f = open(filename, "w")
				f.write(header)
				#f.write("\nLog = /tmp/breddels-log-condor-" +modelname.replace("/", "_") +".txt\n")
			else:
				pass
			data = string.Template(condor_item_template).substitute(**params)
			f.write(data) 
			i += 1
			#logger.info("fininshed %d out of %d" % (i, len(self.parameterset.parameter_values)))
			progressbar.update(i)
		if self.use_dagman:
			filename = os.path.join(self.dirname, "condor", self.dagman_filename)
			logger.info("dagman: %s" % filename)
			f = file(filename, "w")
#Job schw0 models/car/spherical-norot/schw/binned_2p/dm_1.258925e+07/rs_1.000000e-01/condor-single.batch
#Job schw1 models/car/spherical-norot/schw/binned_2p/dm_1.258925e+07/rs_1.333521e-01/condor-single.batch
			i = 0
			jobs = []
			for filename in filenames:
				print >>f, "Job schw%d %s" % (i, filename)
				jobs.append("schw%d" % i)
				i += 1
			print >>f, "Job dummy models/dummy.condor"
			print >>f, "Parent", " ".join(jobs), " CHILD dummy"
		#progressbar.update(i+1)		
class CollectAll(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		assert self.parameterset.parameter_values
		
	def run(self, args, opts, scope):
		#values = self.parameterset.parameter_values
		#options = scope.options
		#arguments = scope.arguments
		self.parameterset.collect()
		self.parameterset.grid()
		
		
		
class RunIterations(object):
	def __init__(self, iterations, command_first, command_post, command, command_post_last, start=0):
		self.iterations = iterations
		self.command_first = command_first
		self.command_post = command_post
		self.command = command
		self.command_post_last = command_post_last
		self.start = start
		
		
	def run(self, args, opts, scope):
		firstscope = scope.clone()
		def run(name, iteration):
			scope = firstscope.clone()
			scope["iteration"] = iteration
			print "RUNNING: %s iteration: %d" % (name, iteration)
			cmd = scope[name]
			cmd.run(args, opts, scope)
			
		for i in range(self.start, self.iterations+1):
			if i == 0:
				run(self.command_first, 0)
			else:
				run(self.command, i)
			
			if i == self.iterations:
				run(self.command_post_last, i)
			else:
				run(self.command_post, i)
		
class RunSchwIterationsDagman(object):
	def __init__(self, mabroot, dirname, iterations, start=0, filename="schw.dag"):
		self.mabroot = mabroot
		self.dirname = dirname
		self.iterations = iterations
		#self.command_first = command_first
		#self.command_post = command_post
		#self.command = command
		#self.command_post_last = command_post_last
		self.start = start
		self.filename = filename

	def format_dagname(self, iteration):
		return "condor.%03d.dag" % (iteration)

	def format_dagname_postprocess(self, iteration):
		return "condor_postprocess.%03d.dag" % (iteration)

	def run(self, args, opts, scope):
		filename = os.path.join(self.dirname, self.filename)
		logger.info("dag filename: %s" % filename)
		file_dag = open(filename, "w")
		firstscope = scope.clone()
		for iteration in range(self.start, self.iterations+1):
			first = iteration == 0
			last = iteration == self.iterations
			filename = self.format_dagname(iteration)
			dagman_X = os.path.join(self.dirname, filename)
			filename = self.format_dagname_postprocess(iteration)
			dagman_postprocess_X = os.path.join(self.dirname, filename)
			# todo --modelpath etc
			
			post = os.path.join(self.mabroot, "bin", "mab-cmd") + " command_schw_dag_post --modelpath=%s --include %s --log-level=warning --configurations=%s --iteration=%d" % (opts.modelpath, opts.include, opts.configurations, iteration)
			pre  = os.path.join(self.mabroot, "bin", "mab-cmd") + " command_schw_dag_pre  --modelpath=%s --include %s --log-level=warning --configurations=%s --iteration=%d" % (opts.modelpath, opts.include, opts.configurations, iteration)
			
			print >>file_dag, "SUBDAG EXTERNAL schwset_%03d %s" % (iteration, dagman_X)
			#if iteration > -1:
			print >>file_dag, "Script Pre schwset_%03d %s" % (iteration, pre)
			#print >>file_dag, "SUBDAG EXTERNAL schwset_postprocess_%03d %s" % (iteration, dagman_postprocess_X)
			print >>file_dag, "Script Post schwset_%03d %s\n" % (iteration, post)
			print >>file_dag
		#for iteration in range(self.start, self.iterations+1):
		#	print >>file_dag, "Parent schwset_%03d Child schwset_postprocess_%03d" % (iteration, iteration)
		for iteration in range(self.start, self.iterations):
			print >>file_dag, "Parent schwset_%03d Child schwset_%03d" % (iteration, iteration+1)
			
		print >>file_dag
		print >>file_dag, "CONFIG /Users/users/breddels/mab/models/dagman.config"
			
#Job schw0 models/car/spherical-norot/schw/binned_2p/dm_1.258925e+07/rs_1.000000e-01/condor-single.batch
#Job schw1 models/car/spherical-norot/schw/binned_2p/dm_1.258925e+07/rs_1.333521e-01/condor-single.batch
		
		#for iteration in range(self.start, self.iterations+1):
		#	"Job schw0 models/car/spherical-norot/schw/binned_2p/dm_1.258925e+07/rs_1.000000e-01/condor-single.batch

			
		file_dag.close()
		

		
				
class CondorDagCommand(object):
	def __init__(self, parameterset, command, condor_head_template, condor_item_template, rootdir, condor_requirements, filename="schw.dag"):
		self.parameterset = parameterset
		self.command = command
		self.condor_head_template = condor_head_template
		self.condor_item_template = condor_item_template
		self.rootdir = rootdir
		self.condor_requirements = condor_requirements
		self.dirname = self.parameterset.dirname
		self.filename = filename
		assert self.parameterset.use_value_directories
		"""
		SUBDAG EXTERNAL schwset0 /Users/users/breddels/phd/models/car/spherical-norot/schw/binned_2p/condor.000.dag
Script Post schwset0 /Users/users/breddels/phd/bin/schw_collect_iterate_prepare_submit.sh 0 1 models/car/spherical-norot /Users/users/breddels/phd/models/car/spherical-norot/schw/binned_2p/results.000.eps /Users/users/breddels/phd/models/car/spherical-norot/schw/binned_2p/condor.001.dag
SUBDAG EXTERNAL schwset1 /Users/users/breddels/phd/models/car/spherical-norot/schw/binned_2p/condor.001.dag
"""

	def format_dagname(self, iteration):
		return "condor.%03d.dag" % (iteration)
		
		
	def run(self, args, opts, scope):
		dagmanscript = ""
		for iteration in range(self.parameterset.max_iterations+1):
			relpath = "/Users/users/breddels/mab"
			dagmanfilename = os.path.join(relpath, self.modelpath, "schw", self.schwsetname, self.format_dagname(iteration))
			if iteration < self.parameterset.skip:
				done = " DONE"
			else:
				done = ""
			dagmanscript += "SUBDAG EXTERNAL schwset%d %s%s\n" % (iteration, dagmanfilename, done)
			resultplot = os.path.join(relpath, self.modelpath, "schw", self.schwsetname, self.format_plotname(iteration))
			
			if iteration < self.parameterset.max_iterations:
				nextdagmanfilename = os.path.join(relpath, self.modelpath, "schw", self.schwsetname, self.format_dagname(iteration+1))
				dagmanscript += "Script Post schwset%d /Users/users/breddels/mab/bin/mab-schw-parameterset-collect-iterate-prepare-submit  %d %d %s %s %s\n" % (iteration, iteration, iteration+1, self.modelpath, resultplot, nextdagmanfilename)
			if iteration == self.parameterset.max_iterations:
				dagmanscript += "Script Post schwset%d /Users/users/breddels/mab/bin/mab-schw-parameterset-collect --iteration=%d --modelpath=%s\n" % (iteration, iteration, self.modelpath)
		for iteration in range(self.parameterset.max_iterations):
			dagmanscript += "Parent schwset%d Child schwset%d\n" % (iteration, iteration+1)
		dagmanscript += "CONFIG /Users/users/breddels/mab/models/dagman.config\n"
		filename = os.path.join(self.modelpath, "schw", self.schwsetname, "schw.dag")
		logger.info("writing dagman file (%s)" % filename)
		logger.info("to submit:\ncondor_submit_dag -f -no_recurse %s" % filename)
		file(filename, "w").write(dagmanscript)
	
