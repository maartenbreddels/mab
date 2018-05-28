# -*- coding: utf-8 -*-
import mab.gd.logging as logging
import os
import string
from numpy import *
import numpy
from kaplot import *
usevtk = False
if usevtk:
	from vtk import *
import scipy.ndimage
import scipy.interpolate as interpolate
import mab.gd.schw.configure
import pickle
from mab.grid.basis import *
from mab import fastsg

logger = logging.getLogger("gd.schw.parametersweep")

def logschwdirname(modelpath, schwsetname):
	dirname = os.path.join(modelpath, schwsetname)
	logger.info("schw path: %s" % dirname)

"""class ParameterSetSingle(object):
	def __init__(self, modelpath, schwsetname, postfix, **params):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.params = params
		self.postfix = postfix
		dirname = os.path.join(modelpath, "schw", schwsetname, "single")
		self.paramnames = [name for name, value in params.items()]
		values = [value for name, value in params.items()]
		self.paramlist = [("single", dirname, values, values)]
		self.max_iterations = 1
		self.skip = 0
		
	def load(self, iteration):
		pass"""

class ParameterSetTrue(object):
	def __init__(self, modelpath, runtemplate, postfix="", schwsetname="true"):
		self.modelpath = modelpath
		self.runtemplate = runtemplate
		self.schwsetname = schwsetname
		self.postfix = postfix
		self.dirname = os.path.join(modelpath, "schw", self.schwsetname)
		self.modeldirname = os.path.join(modelpath, "schw", self.schwsetname, "true")
		if 0:
			self.paramnames = [name for name, value in params.items()]
			values = [value for name, value in params.items()]
			self.paramlist = [("single", dirname, values, values)]
			self.max_iterations = 1
			self.skip = 0
			
	def load(self, iteration=0):
		logger.info("does not need load")
		
	def init(self, iteration=0):
		logger.info("does not need init")
		
	def prepare(self, iteration=0):
		if not os.path.exists(self.dirname):
			logger.info("creating directory: %s" % self.dirname)
			os.makedirs(self.dirname)
		if not os.path.exists(self.modeldirname):
			logger.info("creating directory: %s" % self.modeldirname)
			os.makedirs(self.modeldirname)
			
		dirname = self.modeldirname
		
		name = "true"
		for subdirname in ["log", "intermediate", "results"]:
			dirnamesub = os.path.join(dirname, subdirname)
			if not os.path.exists(dirnamesub):
				logger.info("creating directory %s" % dirnamesub)
				os.makedirs(dirnamesub)
			else:
				logger.info("directory %s already exists" % dirnamesub)
		
		script = file(self.runtemplate).read()
		vars = {}
		vars["modelpath"] = os.path.abspath(self.modelpath)
		vars["schwsetname"] = self.schwsetname
		vars["schwmodelname"] = name
		script = script % vars
		filename = os.path.join(self.modelpath, "schw", self.schwsetname, name, "run.sh")
		logger.info("writing runscript to filename: %s" % filename)
		file(filename, "w").write(script)
		os.chmod(filename, 0700)			
		

class ParameterValues(object):
	def __init__(self, modelpath,  schwsetname, schwmodelname, values_org, values_trans, p_prior):
		self.modelpath  = modelpath
		self.schwsetname = schwsetname
		self.modelname  = self.schwmodelname = schwmodelname
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.values_org = values_org
		self.values_trans = values_trans
		self.p_prior = p_prior
	

		
class ParameterSetSweep(object):
	isregular = True
	def __init__(self,  modelpath, schwsetname, paramformats, postfix="", **extra_params):
		self.modelpath = modelpath
		#self.runtemplate = runtemplate
		self.schwsetname = schwsetname
		self.paramformats = paramformats
		#self.execute = execute
		self.extra_params = extra_params
		#self.condor_params = {"rootdir":rootdir, "condor_requirements":condor_requirements, "schwsetname":setname, "modelpath":modelpath}
		#self.condor_head_template = condor_head_template
		#self.condor_item_template = condor_item_template
		self.postfix = postfix
		
		self.paramlist = self.createparamlist(self.paramformats[0], self.paramformats[1:], "", [], [], [])
		self.paramnames = [param_name for param_name, param_format, param_transformation, param_values in self.paramformats]
		self.ranges_org = [k[-1] for k in paramformats]
		#self.ranges_new = [k[-1] for k in paramformats]
		self.ranges_new = [k[-2](k[-1]) for k in paramformats]
		
		
		shape = []
		for paramformat in self.paramformats:
			name = paramformat[0]
			values = paramformat[-1]
			shape.append(len(values))
		self.shape = tuple(shape)
		
		#print self.paramlist
		#for p in self.paramlist:
		#	print p
		#	
		#print len(self.paramlist)
		
	def createparamlist(self, paramformat_head, paramformat_tail, name, indices, valueslist_org, valueslist_new):
		paramlist = []
		#print paramformat_head
		param_name, param_format, param_transformation, param_values = paramformat_head
		for i, value_org in enumerate(param_values):
			value_new = param_transformation(value_org)
			indices_new = list(indices) + [i]
			valueslist_new_sub = list(valueslist_new) + [value_new]
			valueslist_org_sub = list(valueslist_org) + [value_org]
			if name:
				name_new = name + "/" + (param_format % value_new)
			else:
				name_new = (param_format % value_new)
			#print name_new
			if len(paramformat_tail) >= 1:
				paramlist.extend(self.createparamlist(paramformat_tail[0], paramformat_tail[1:], name_new, indices_new, valueslist_org_sub, valueslist_new_sub))
			else:
				dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name_new)
				#paramlist.append((name_new, dirname, indices_new, values_new))
				paramlist.append((name_new, dirname, valueslist_org_sub, valueslist_new_sub, indices_new))
		return paramlist

class ParameterSetHyperOctree(object):
	isregular = False
	def __init__(self, modelpath, schwsetname, basename, paramformats, max_newpoints, max_level, binary, doall=False, init=False, initial_subdivides=2, postfix="", skip=0, **extra_params):
		
		self.modelpath = modelpath
		#self.runtemplate = runtemplate
		self.schwsetname = schwsetname
		self.paramformats = paramformats
		self.doall = doall
		self.basename = basename
		self.postfix = postfix
		self.skip = skip
		self.dimension = len(self.paramformats)
		
		self.initial_subdivides = initial_subdivides
		self.max_newpoints = max_newpoints
		self.max_level = max_level
		
		self.paramnames = [k[0] for k in self.paramformats]
		self.functions_transform = [k[4] for k in self.paramformats]
		self.paramlabels = [k[5] for k in self.paramformats]
		#self.paramlabels = [k.split(".")[-1] for k in self.paramnames]
		self.ranges_org = [k[3] for k in self.paramformats]
		self.dirname = filename = os.path.join(self.modelpath, "schw", self.schwsetname)
		self.binary = binary
		self.max_iterations = len(max_newpoints)
		
	def clean(self):
		def _clean(pattern):
			dirname = self.dirname
			cmd = "find %s -iname '%s' | xargs rm" % (dirname, pattern)
			print cmd
			os.system(cmd)
		_clean("logprob*")
		_clean("finished.orblib")
		_clean("finished.initial")
		_clean("schw.dag.*")
		_clean("condor.???.dag.dagman.*")
		_clean("condor.???.dag.lib.*")
		#_clean("condor.??[1-9].dag.*")
		
	def format_name(self, iteration):
		return "%s.%03d" % (self.basename, iteration)
	
	def init(self):
		#dirname = self.parameterset.dirname
		
		dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		if not os.path.exists(dirname):
			logger.info("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.info("directory %s already exists" % dirname)
		
		output_filename = os.path.join(self.dirname, self.format_name(0))
		cmd = "%s -d %d -i -o %s -s %d " % (self.binary, self.dimension, output_filename, self.initial_subdivides)
		cmd += " -- "
		for a, b in self.ranges_org:
			cmd += str(a) + " "
		for a, b in self.ranges_org:
			cmd += str(b) + " "
		print cmd
		if os.system(cmd) != 0:
			return
		
	def iterate(self, iteration, optimize=True):
		if optimize:
			inputname = os.path.join(self.dirname, self.format_name(iteration))
			outputname = os.path.join(self.dirname, self.format_name(iteration+1))
			rest = "-- 0 0 0 %d %d" % (self.max_newpoints[iteration], self.max_level)
			template = "%s -d %d --output=%s --input=%s --limit-level-difference=1 --optimize " + rest
			cmd = template % (self.binary, self.dimension, outputname, inputname)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
		else:
			outputname_known = os.path.join(self.dirname, self.format_name(iteration+1))+ ".known"
			
			inputname_known = os.path.join(self.dirname, self.format_name(iteration))+ ".known"
			cmd = "cp %s %s" % (inputname_known, outputname_known)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
			
			inputname_solved = os.path.join(self.dirname, self.format_name(iteration))+ ".solved"
			cmd = "cat %s >> %s" % (inputname_solved, outputname_known)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
			
		
	def load(self, iteration, load_unknown=True, load_known=False, load_solved=False):
		self.paramlist = []
		self.paramvalues = []
		self.volumes = []
		filename = os.path.join(self.dirname, self.format_name(iteration))
		if not os.path.exists(self.dirname):
			logger.info("creating directory: %s" % self.dirname)
			os.makedirs(self.dirname)
		self.dimension = len(self.ranges_org)
		logger.info("dimension = %d" % self.dimension)
		if 0:
			output_filename = os.path.join(dirname, inputname)
			cmd = "%s -d %d -i -o %s -s %d " % (self.binary, self.dimension, output_filename, intial_subdivides)
			cmd += " -- "
			for a, b in self.ranges_org:
				cmd += str(a) + " "
			for a, b in self.ranges_org:
				cmd += str(b) + " "
			print cmd
			if os.system(cmd) != 0:
				return
			
			
		#self.
		#self.execute = execute
		#self.extra_params = extra_params
		#self.condor_params = {"rootdir":rootdir, "condor_requirements":condor_requirements, "schwsetname":setname, "modelpath":modelpath}
		#self.condor_head_template = condor_head_template
		#self.condor_item_template = condor_item_template
		#self.postfix = postfix
		#self.paramlist = self.createparamlist(self.paramformats[0], self.paramformats[1:], "", [], [])
		#self.paramnames = [param_name for param_name, param_format, param_values in self.paramformats]

		def do(name, skiplast):
			inputname = self.format_name(iteration)
			filename = os.path.join(self.modelpath, "schw", self.schwsetname, inputname + name)
			lines = file(filename).readlines()
			for line in lines:
				values = [float(k) for k in line.split()]
				values, volume = values[:-1-skiplast], values[-1-skiplast]
				newvalues = []
				dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
				names = []
				p_prior = 1.
				for value, paramformat in zip(values, self.paramformats):
					param_name, param_format, param_transformation, param_range = paramformat[:4]
					newvalue = param_transformation(value)
					name = param_format % newvalue
					newvalues.append(newvalue)
					names.append(name)
					param_prior = paramformat[6]
					p_prior *= param_prior(value)
					
				name = "/".join(names)
				dirname = os.path.join(dirname, name)
				self.paramlist.append((name, dirname, values, newvalues))
				self.paramvalues.append(ParameterValues(self.modelpath, self.schwsetname, name, values, newvalues, p_prior))
				self.volumes.append(volume)
		if load_unknown:
			do(".unknown", 0)
		if load_known:
			do(".known", 1)
		if load_solved:
			do(".solved", 1)
		#if doall:
		#	do(".known", 1)
		
from mab.grid.binarygrid import *

class ParameterSetHyperOctreePython(ParameterSetHyperOctree):
	def init(self):
		dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		if not os.path.exists(dirname):
			logger.info("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.info("directory %s already exists" % dirname)
		
		domain = self.ranges_org
		baseslist = [BasisTriLevels(0., 1., 1+k) for k in range(10)]
		self.grid = BinaryGrid(self.dimension, domain=domain, baseslist=baseslist)
		self.grid.rootNode.split(self.initial_subdivides)
		logger.info("initial nodes: %d" % len(self.grid))
		self.refresh()
		
	def refresh(self):
		self.paramlist = []
		self.paramvalues = []
		self.volumes = []
		for i in range(len(self.grid)):
			values = self.grid.points[i]
			volume = 1. # FIXME: implement volumes
			newvalues = []
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
			names = []
			p_prior = 1.
			for value, paramformat in zip(values, self.paramformats):
				param_name, param_format, param_transformation, param_range = paramformat[:4]
				newvalue = param_transformation(value)
				name = param_format % newvalue
				newvalues.append(newvalue)
				names.append(name)
				param_prior = paramformat[5]
				p_prior *= param_prior(value)
			name = "/".join(names)
			dirname = os.path.join(dirname, name)
			self.paramlist.append((name, dirname, values, newvalues))
			self.paramvalues.append(ParameterValues(self.modelpath, self.schwsetname, name, values, newvalues, p_prior))
			self.volumes.append(volume)
		
		
	def iterate(self, iteration, optimize=True, cores=1):
		self.grid.optimize(self.max_newpoints[iteration], max_level=self.max_level, cores=cores)
		logger.info("new number of nodes: %d" % len(self.grid.newpoints))
		logger.info("total number of nodes: %d" % len(self.grid))
		self.refresh()
		
	def feed(self, values):
		assert len(self.grid.points) == len(values)
		for point, value in zip(self.grid.points, values):
			self.grid.values[point] = value
	def load(self, iteration, load_unknown=True, load_known=False, load_solved=False):
		pass
				
class ParameterSetTri(ParameterSetHyperOctree):
	isregular = False
	def __init__(self, modelpath, schwsetname, basename, paramformats, max_level, initial_level=2, postfix="", skip=0, **extra_params):
		
		self.modelpath = modelpath
		#self.runtemplate = runtemplate
		self.schwsetname = schwsetname
		self.paramformats = paramformats
		#self.doall = doall
		self.basename = basename
		self.postfix = postfix
		self.skip = skip
		self.dimension = len(self.paramformats)
		
		self.initial_level = initial_level
		self.max_level = max_level
		
		self.paramnames = [k[0] for k in self.paramformats]
		self.functions_transform = [k[4] for k in self.paramformats]
		self.paramlabels = [k[5] for k in self.paramformats]
		self.parampriors = [k[6] for k in self.paramformats]
		#self.paramlabels = [k.split(".")[-1] for k in self.paramnames]
		self.ranges_org = [k[3] for k in self.paramformats]
		self.dirname = filename = os.path.join(self.modelpath, "schw", self.schwsetname)
		#self.binary = binary
		self.max_iterations = 0
		#self.max_iterations = len(max_newpoints)
		
	def __call__(self, point):
		point_normalized = [(point[k] - self.ranges_org[k][0])/(self.ranges_org[k][1]-self.ranges_org[k][0]) for k in range(self.dimension)]
		#print "get", point_normalized, self.grid(*point)
		return self.grid(*point_normalized)
	
	def eval_normalized(self, point_normalized):
		#point_normalized = [(point[k] - self.ranges_org[k][0])/(self.ranges_org[k][1]-self.ranges_org[k][0]) for k in range(self.dimension)]
		#print "get", point_normalized, self.grid(*point)
		return self.grid(*point_normalized)
	
	def create_grid(self):
		return mab.gd.gdfast.RegularGridD(self.dimension, self.initial_level)
	def init(self):
		dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		if not os.path.exists(dirname):
			logger.info("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.info("directory %s already exists" % dirname)
		
		domain = self.ranges_org
		baseslist = [BasisTriLevels(0., 1., 1+k) for k in range(10)]
		self.grid = self.create_grid()
		
		#self.grid.rootNode.split(self.initial_subdivides)
		#logger.info("initial nodes: %d" % len(self.grid))
		self.refresh()
		
	def refresh(self):
		self.paramlist = []
		self.paramvalues = []
		self.volumes = []
		self.points = []
		def f(*point_normalized):
			point = [point_normalized[k]*(self.ranges_org[k][1]-self.ranges_org[k][0]) + self.ranges_org[k][0] for k in range(self.dimension)]
			point = tuple(point)
			#print "p", point
			self.points.append(point)
			return -1000
		self.grid.eval(f)
		logger.info("number of points: %d" % len(self.points))
		for i in range(len(self.points)):
			values = self.points[i]
			volume = 1. # FIXME: implement volumes
			newvalues = []
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
			names = []
			p_prior = 1.
			for value, paramformat in zip(values, self.paramformats):
				param_name, param_format, param_transformation, param_range = paramformat[:4]
				newvalue = param_transformation(value)
				name = param_format % newvalue
				newvalues.append(newvalue)
				names.append(name)
				param_prior = paramformat[6]
				p_prior *= param_prior(value)
			name = "/".join(names)
			dirname = os.path.join(dirname, name)
			self.paramlist.append((name, dirname, values, newvalues))
			self.paramvalues.append(ParameterValues(self.modelpath, self.schwsetname, name, values, newvalues, p_prior))
			self.volumes.append(volume)
		
		
	def iterate(self, iteration, optimize=True, cores=1):
		logger.info("ignoring iteration: %d" % iteration)
		#self.grid.optimize(self.max_newpoints[iteration], max_level=self.max_level, cores=cores)
		#logger.info("new number of nodes: %d" % len(self.grid.newpoints))
		#logger.info("total number of nodes: %d" % len(self.grid))
		#self.refresh()
		
	def feed(self, values):
		print "feed", len(self.points), len(values)
		self.maxlog = max(values)
		print "maxlog", self.maxlog
		assert len(self.points) == len(values)
		#print zip(self.points, values)
		#for point, value in zip(self.points, values):
		#	self.grid.values[point] = value
		mapping = dict(zip(self.points, values))
		def f(*point_normalized):
			point = [point_normalized[k]*(self.ranges_org[k][1]-self.ranges_org[k][0]) + self.ranges_org[k][0] for k in range(self.dimension)]
			point = tuple(point)
			#print point, "=", mapping[point]
			return mapping[point]
		self.grid.eval(f)
		#print "DDDDDDDDDDDDDDDDDDDDDDDDDDDD"
		#print self.grid(0., 0.)
		#print self.grid(0., 0.)
		#print self.grid(1., 1.)
		
	def load(self, iteration, load_unknown=True, load_known=False, load_solved=False):
		pass

class ParameterSetSparse(ParameterSetTri):
	isregular = False
	
	def create_grid(self):
		#basis = mab.gd.gdfast.BasisSetHierTriangular()
		return mab.gd.gdfast.SparseGridD(self.dimension, self.initial_level)

	def make_linear(self):
		offset = self.maxlog
		#offset = self.grid(*((0.5,) * self.dimension))
		#offset = self.grid(*((0.5,) * self.dimension))
		def f(*point_normalized):
			return exp(self.grid(*point_normalized) - offset)
		N = self.initial_level+9
		gridlinear = mab.gd.gdfast.SparseGridD(self.dimension, N)
		gridlinear.eval(f)
		for i in range(N):
			print i, gridlinear.integrate(i)
		self.gridold = self.grid
		self.grid = gridlinear

class ParameterSetSparseFast(ParameterSetTri):
	isregular = False
	def __init__(self, modelpath, schwsetname, basename, paramformats, max_level, initial_level=2, postfix="", skip=0, **extra_params):
		
		self.modelpath = modelpath
		#self.runtemplate = runtemplate
		self.schwsetname = schwsetname
		self.paramformats = paramformats
		#self.doall = doall
		self.basename = basename
		self.postfix = postfix
		self.skip = skip
		self.dimension = len(self.paramformats)
		
		self.initial_level = initial_level
		self.max_level = max_level
		
		self.paramnames = [k[0] for k in self.paramformats]
		self.functions_transform = [k[4] for k in self.paramformats]
		self.paramlabels = [k[5] for k in self.paramformats]
		#self.paramlabels = [k.split(".")[-1] for k in self.paramnames]
		self.ranges_org = [k[3] for k in self.paramformats]
		self.dirname = filename = os.path.join(self.modelpath, "schw", self.schwsetname)
		#self.binary = binary
		self.max_iterations = 0
		#self.max_iterations = len(max_newpoints)
		
	def __call__(self, point):
		point_normalized = [(point[k] - self.ranges_org[k][0])/(self.ranges_org[k][1]-self.ranges_org[k][0]) for k in range(self.dimension)]
		#print "get", point_normalized, self.grid(*point)
		return self.grid(*point_normalized)
	
	
	def init(self):
		dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		if not os.path.exists(dirname):
			logger.info("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.info("directory %s already exists" % dirname)
		
		domain = self.ranges_org
		baseslist = [BasisTriLevels(0., 1., 1+k) for k in range(10)]
		#self.grid = self.create_grid()
		
		#self.grid.rootNode.split(self.initial_subdivides)
		#logger.info("initial nodes: %d" % len(self.grid))
		self.refresh()
		
	def refresh(self):
		self.paramlist = []
		self.paramvalues = []
		self.volumes = []
		self.points = []
		def f(*point_normalized):
			point = [point_normalized[k]*(self.ranges_org[k][1]-self.ranges_org[k][0]) + self.ranges_org[k][0] for k in range(self.dimension)]
			point = tuple(point)
			#print "p", point
			self.points.append(point)
			return -1000
		self.func = fastsg.PyFunction(self.dimension, f)
		self.grid = fastsg.SparseGrid(self.initial_level, self.func)
		#self.grid.eval(f)
		logger.info("number of points: %d" % len(self.points))
		for i in range(len(self.points)):
			values = self.points[i]
			volume = 1. # FIXME: implement volumes
			newvalues = []
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
			names = []
			p_prior = 1.
			for value, paramformat in zip(values, self.paramformats):
				param_name, param_format, param_transformation, param_range = paramformat[:4]
				newvalue = param_transformation(value)
				name = param_format % newvalue
				newvalues.append(newvalue)
				names.append(name)
				param_prior = paramformat[5]
				p_prior *= param_prior(value)
			name = "/".join(names)
			dirname = os.path.join(dirname, name)
			self.paramlist.append((name, dirname, values, newvalues))
			self.paramvalues.append(ParameterValues(self.modelpath, self.schwsetname, name, values, newvalues, p_prior))
			self.volumes.append(volume)
		
		
	def iterate(self, iteration, optimize=True, cores=1):
		logger.info("ignoring iteration: %d" % iteration)
		#self.grid.optimize(self.max_newpoints[iteration], max_level=self.max_level, cores=cores)
		#logger.info("new number of nodes: %d" % len(self.grid.newpoints))
		#logger.info("total number of nodes: %d" % len(self.grid))
		#self.refresh()
		
	def feed(self, values):
		print "feed", len(self.points), len(values)
		self.maxlog = max(values)
		print "maxlog", self.maxlog
		assert len(self.points) == len(values)
		#print values, self.points
		#print zip(self.points, values)
		#for point, value in zip(self.points, values):
		#	self.grid.values[point] = value
		mapping = dict(zip(self.points, values))
		def f(*point_normalized):
			point = [point_normalized[k]*(self.ranges_org[k][1]-self.ranges_org[k][0]) + self.ranges_org[k][0] for k in range(self.dimension)]
			point = tuple(point)
			#print point, "=", mapping[point]
			return mapping[point]
		self.func = fastsg.PyFunction(self.dimension, f)
		self.grid = fastsg.SparseGrid(self.initial_level, self.func)
		self.grid.hierarchize()
		#self.grid.eval(f)
		#print "DDDDDDDDDDDDDDDDDDDDDDDDDDDD"
		#print self.grid(0., 0.)
		#print self.grid(0., 0.)
		#print self.grid(1., 1.)
		
	def load(self, iteration, load_unknown=True, load_known=False, load_solved=False):
		pass
	isregular = False
	
	def make_linear(self):
		offset = self.maxlog
		#offset = self.grid(*((0.5,) * self.dimension))
		#offset = self.grid(*((0.5,) * self.dimension))
		def f(*point_normalized):
			return exp(self.grid(*point_normalized) - offset)
		gridlinear = mab.gd.gdfast.SparseGridD(self.dimension, self.initial_level+5)
		gridlinear.eval(f)
		self.gridold = self.grid
		self.grid = gridlinear
		

class ParameterSetExecutorCondor(object):
	def __init__(self, modelpath, runtemplate, rootdir, condor_requirements, condor_head_template, condor_item_template, parameterset, execute=True, output_true=False, **extra_params):
		self.modelpath = modelpath
		self.runtemplate = runtemplate
		self.schwsetname = parameterset.schwsetname
		self.execute = execute
		self.extra_params = extra_params
		self.condor_params = {"rootdir":rootdir, "condor_requirements":condor_requirements, "schwsetname":self.schwsetname, "modelpath":modelpath}
		self.condor_head_template = condor_head_template
		self.condor_item_template = condor_item_template
		self.output_true = output_true
		#self.postfix = postfix
		self.parameterset = parameterset
		
		#self.paramlist = self.createparamlist(self.paramformats[0], self.paramformats[1:], "", [], [])
	def format_name(self, iteration):
		return "condor.%03d.batch" % (iteration)
	
	def format_plotname(self, iteration):
		return "results.%03d.eps" % (iteration)
	
	def format_dagname(self, iteration):
		return "condor.%03d.dag" % (iteration)
	
	def prepare(self, iteration=0):
		self.parameterset.load(iteration)
		execute = self.execute
		
		dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		condor_arguments_list = []
		if not os.path.exists(dirname):
			logger.debug("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.debug("directory %s already exists" % dirname)
			
		condortxt = ""
		txt = file(self.condor_head_template).read()
		condortxt += string.Template(txt).substitute(**self.condor_params)
		condortxt += "\n"
		
		
		filename = os.path.join(dirname, "runall.%03d.sh" % iteration)
		logger.info("all scripts collected in: %s" % filename)
		file_runall = open(filename, "w") 
		
		if self.output_true:
			name = "true"
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name)
			
			for subdirname in ["log", "intermediate", "results"]:
				dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name, subdirname)
				if not os.path.exists(dirname):
					logger.debug("creating directory %s" % dirname)
					if execute:
						os.makedirs(dirname)
				else:
					logger.debug("directory %s already exists" % dirname)
			
			script = file(self.runtemplate).read()
			vars = {}
			vars["modelpath"] = os.path.abspath(self.modelpath)
			vars["schwsetname"] = self.schwsetname
			vars["schwmodelname"] = name
			script = script % vars
			filename = os.path.join(self.modelpath, "schw", self.schwsetname, name, "run.sh")
			print >>file_runall, filename
			logger.debug("writing runscript to filename: %s" % filename)
			if execute:
				file(filename, "w").write(script)
				os.chmod(filename, 0700)
			
			
			
		single_condor_filenames = []
		
		
		for paramlist in self.parameterset.paramlist:
			name, dirname, orig_values, true_values = paramlist[:4]
			#dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name_new)
			#print dirname, self.modelpath, self.schwsetname, name_new
			if not os.path.exists(dirname):
				logger.debug("creating directory %s" % dirname)
				if execute:
					os.makedirs(dirname)
			else:
				logger.debug("directory %s already exists" % dirname)
			
			for subdirname in ["log", "intermediate", "results"]:
				dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name, subdirname)
				if not os.path.exists(dirname):
					logger.debug("creating directory %s" % dirname)
					if execute:
						os.makedirs(dirname)
				else:
					logger.debug("directory %s already exists" % dirname)
			
			filename = os.path.join(self.modelpath, "schw", self.schwsetname, name, "galaxy.ini")
			logger.debug("writing parameters to filename: %s" % filename)
			if execute:
				f = open(filename, "w")
				print >>f, "[parameters]"
				#for key, value in params_new.items():
				for key, value in zip(self.parameterset.paramnames, true_values):
					print >>f, "%s=%r" % (key, value)
				
			script = file(self.runtemplate).read()
			vars = {}
			vars["modelpath"] = os.path.abspath(self.modelpath)
			vars["schwsetname"] = self.schwsetname
			vars["schwmodelname"] = name
			script = script % vars
			filename = os.path.join(self.modelpath, "schw", self.schwsetname, name, "run.sh")
			print >>file_runall, filename
			logger.debug("writing runscript to filename: %s" % filename)
			if execute:
				file(filename, "w").write(script)
				os.chmod(filename, 0700)
				
			condor_arguments = dict(self.condor_params)
			condor_arguments["schwmodelname"] = name
			#condor_arguments_list.append(condor_arguments)
			
			filename = os.path.join(self.modelpath, "schw", self.schwsetname, name, "condor-single.batch")
			single_condor_filenames.append(filename)
			condor_singletxt = ""
			txt = file(self.condor_head_template).read()
			condor_singletxt += string.Template(txt).substitute(**self.condor_params)
			condor_singletxt += "\n"
			
			txt = file(self.condor_item_template).read()
			condor_singletxt_item = string.Template(txt).substitute(**condor_arguments) + "\n"
			condor_singletxt += condor_singletxt_item
			#condor_singletxt += "\n"
			if execute:
				file(filename, "w").write(condor_singletxt)
			
			condortxt += condor_singletxt_item
		
			
		filename = os.path.join(self.modelpath, "schw", self.schwsetname, self.format_name(iteration))
		logger.info("writing condor batch (%s)" % filename)
		file(filename, "w").write(condortxt)
		
		
		#if iteration == 0:
			#for iteration in range(self.parameterset.max_iterations):
		if 1:
			#jobindex = 0
				#single_condor_filenames = single_condor_filenames[:10]
			#for iteration in range(self.parameterset.max_iterations):
				nextcondorfilename = os.path.join(self.modelpath, "schw", self.schwsetname, self.format_name(iteration+1))
				nextdagmanfilename = os.path.join(self.modelpath, "schw", self.schwsetname, self.format_dagname(iteration+1))
				resultplot = os.path.join(self.modelpath, "schw", self.schwsetname, self.format_plotname(iteration))
				dummycondor = os.path.join("models", "dummy.condor")
				jobnames = []
				dagmanscript = "#generated file\n"
				for i, filename in enumerate(single_condor_filenames):
					jobname = "schw%d" % i
					jobnames.append(jobname)
					dagmanscript += "Job %s %s\n" % (jobname, filename)
				dagmanscript += "Job dummy %s\n" % (dummycondor)
				dagmanscript += "Parent %s Child dummy\n" % " ".join(jobnames)
				#dagmanscript += "Parent dummy Child %s\n" % nextdagmanfilename`
				dagmanscript += "CONFIG /Users/users/breddels/mab/models/dagman.config\n"
				#if iteration > 0:
				#	dagmanscript += "Script Pre schwjob /Users/users/breddels/phd/schw_iterate_prepare.sh %d %d\n" % (iteration-1, iteration)
				#if iteration == 0:
				#	dagmanscript += "Script Pre schwjob ~/phd/bin/schw_init_prepare.sh\n"
				#if iteration < self.parameterset.max_iterations-1:
				#	dagmanscript += "Script Post dummy /Users/users/breddels/phd/bin/schw_collect_iterate_prepare_submit.sh %d %d %s %s\n" % (iteration, iteration+1, resultplot, nextdagmanfilename)
				#if iteration == self.parameterset.max_iterations-1:
				#	dagmanscript += "Script Post dummy /Users/users/breddels/phd/bin/schw_collect.sh %d %s\n" % (iteration, resultplot)
					#dagmanscript += "Script Post dummy python /Users/users/breddels/phd/bin/schw_parameterset_collect.py --hardcopy=%s --iteration=%d\n" % (resultplot, iteration)
		
				filename = os.path.join(self.modelpath, "schw", self.schwsetname, self.format_dagname(iteration))
				logger.info("writing dagman file (%s)" % filename)
				file(filename, "w").write(dagmanscript)
				os.system("condor_submit_dag -f -no_submit %s" % filename)
		if iteration == 0:
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
				
		
		total = 1
		#for paramformat in self.paramformats:
		#	name = paramformat[0]
		#	values = paramformat[-1]
		#	logger.info("parameter: %s, [%e, ..., %e] (#=%d)" % (name, min(values), max(values), len(values)))
		#	total *= len(values)
		total = len(self.parameterset.paramlist)
		logger.info("%d models" % total)
		
		 

"""	def collect(self):
		filename = os.path.join(self.modelpath, "schw", self.setname, "quadtree.b.0.solved")
		f = file(filename, "w")
		for name, dirname, orig_values, true_values in self.parameterset.paramlist:
			dirname = os.path.join(self.modelpath, "schw", self.setname, name)
			filename = os.path.join(dirname, "results/logprobability" +self.postfix +".txt")
			if os.path.exists(filename):
				#print name, file(filename).read().strip()
				logp = eval(file(filename).read().strip()) #* 1e100
				p = exp(logp)
				print name, p, logp
				line = " ".join(["%e" % k for k in list(orig_values) + [p]])
				#print line
				print >>f, line
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
			#print indices, p"""
			

#def create_anisotropy_grid(paramlist, plist, volumelist, postfix="", betamin=-3, betamax=1, rmin=0, rmax=1.5):
def create_anisotropy_grid(paramlist, plist, volumelist, nr, postfix="", betamin=-3, betamax=1):
	Nbeta = 500
	#nr = 100
	#Nbeta = 10
	#Nr = 10
	#xmin=0
	#xmax=1.5 # * 2610.9660574412533 
	#rs = (arange(nr) + 0.5 )/ (nr) * (rmax - rmin) + rmin
	betagrid = zeros((nr, Nbeta))
	betagrid_prior = zeros((nr, Nbeta))
	#print betagrid
	#for paramlist in self.parameterset.paramlist:
	#print ">>>>>>>>>>>>>", paramlist, plist, volumelist, "<<<<<<<<<<<<<"
	#print ">>>>>>>>>>>>>", len(paramlist), len(plist), len(volumelist), "<<<<<<<<<<<<<"
	logger.info("constructing anisotropy grid")
	for param, p, volume in zip(paramlist, plist, volumelist):
		#print p
		#print param
		dirname = param[1]
		name, dirname, orig_values, true_values, p_prior = param[:5]
#	for name, dirname, indices, values in self.paramlist[1:]:
		filename = os.path.join(dirname, "results", "solution_moments3d" +postfix +".npy")
		#print dirname, filename, os.path.exists(filename)
		if os.path.exists(filename):
			solution_moments3d = numpy.load(filename)
			varvr = solution_moments3d[4]
			varvphi = solution_moments3d[5]
			varvtheta = solution_moments3d[6]
			betas = 1 - (varvphi + varvtheta)/(2*varvr)
			if 0:
				print betas
				print varvr
				print varvphi
				print varvtheta
			betaindices = ((betas-betamin)/(betamax-betamin) * Nbeta).astype(int)
			#print ".", ~isnan(varvtheta)
			good_betaindices = (betaindices >= 0) & (betaindices < Nbeta) & (~isnan(varvtheta))
			betaindices = betaindices[good_betaindices]
			#p = pgrid[indices[0], indices[1]]
			#print good_betaindices,betaindices
			#print "lalala", p, volume
			#print solution_moments3d.shape
			#print good_betaindices,betaindices, len(good_betaindices)
			#print good_betaindices, good_betaindices.shape
			betagrid[good_betaindices,betaindices] += p * volume
			betagrid_prior[good_betaindices,betaindices] += p * volume * p_prior
			#print betagrid
			#print len(solution_moments3d)
			#f = interpolate.interp1d(oldrs, betagrid[:,i])
		else:
			#logger.error("file '%s' missing (model not finished?)" % filename)
			pass	
	print sort(plist)[-20:]
	print len(paramlist), len(plist), len(volumelist)
	#dsajdksla
	return betagrid, betagrid_prior
			
def sigma_grid(paramlist, plist, volumelist, sigma_max, NR, Nsigma=100, postfix=""):
	sigmagrid = zeros((NR, Nsigma))
	sigmagrid_prior = zeros((NR, Nsigma))
	logger.info("constructing sigma grid")
	for param, p, volume in zip(paramlist, plist, volumelist):
		#print p
		#print param
		dirname = param.dirname #[1]
		#name, dirname, orig_values, true_values, p_prior = param[:5]
		filename = os.path.join(dirname, "results", "solution_projectedmoments" +postfix +".npy")
		if os.path.exists(filename):
			moments2d = numpy.load(filename)
			sigmas = moments2d[2]**0.5
			#print sigmas
			sigma_indices = (sigmas/sigma_max * Nsigma).astype(int)
			good_sigma_indices = (sigma_indices >= 0) & (sigma_indices < Nsigma)# & (~isnan(varvtheta))
			sigma_indices = sigma_indices[good_sigma_indices]
			#print sigmagrid.shape, good_sigma_indices.shape, sigma_indices.shape
			sigmagrid[good_sigma_indices,sigma_indices] += p * volume
			sigmagrid_prior[good_sigma_indices,sigma_indices] += p * volume * param.p_prior
		else:
			#logger.error("file '%s' missing (model not finished?)" % filename)
			pass	
	print sort(plist)[-20:]
	print len(paramlist), len(plist), len(volumelist)
	#dsajdksla
	return sigmagrid, sigmagrid_prior

def create_rho_grid(paramvalues, plist, volumelist, logrmin=-3, logrmax=3, logrhomin=-3, logrhomax=3):
	Nr = 200
	Nrho = 200
	logrs = (arange(Nr) + 0.5 )/ (Nr) * (logrmax - logrmin) + logrmin
	#rs = (arange(Nr) + 0.5 )/ (Nr) * (rmax - rmin) + rmin
	rs = 10**logrs
	logrhogrid = zeros((len(logrs), Nrho))
	#for paramlist in self.parameterset.paramlist:
	Ntot = len(paramvalues)
	i = 0
	for parametervalue, p, volume in zip(paramvalues, plist, volumelist):
		ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
		dm_density = scope["dm_density"]
		print dm_density.alpha
		print i, "out of", Ntot
		i += 1
		#if os.path.exists(filename):
		if 1:
			logrhos = log10(dm_density.densityr(rs))
			rho_indices = ((logrhos-logrhomin)/(logrhomax-logrhomin) * Nrho).astype(int)
			#print ".", ~isnan(varvtheta)
			good_rho_indices = (rho_indices >= 0) & (rho_indices < Nrho)# & (~isnan(varvtheta))
			rho_indices = rho_indices[good_rho_indices]
			logrhogrid[good_rho_indices, rho_indices] += p * volume
	return logrhogrid
		

def create_Menc_and_slope_grid(paramvalues, plist, volumelist, logrmin=-3, logrmax=3, logMmin=0, logMmax=10, slopemin=-1, slopemax=3, rmin=0, rmax=1.5, logarithmic=True, configuration=None):
	Nr = 200
	NM = 500
	Nslope = 500
	Nslope2 = 500*5
	Ntanr = 200
	if logarithmic:
		logrs = (arange(Nr) + 0.5 )/ (Nr) * (logrmax - logrmin) + logrmin
		#rs = (arange(Nr) + 0.5 )/ (Nr) * (rmax - rmin) + rmin
		rs = 10**logrs
	else:
		rs = (arange(Nr) + 0.5 )/ (Nr) * (rmax - rmin) + rmin
	slopes = (arange(Nslope2)+0.5)/Nslope2 * (slopemax - slopemin) + slopemin
		
	slopegrid = zeros((len(rs), Nslope))
	slopegrid_prior = zeros((len(rs), Nslope))
	logMencgrid = zeros((len(rs), NM))
	logMencgrid_prior = zeros((len(rs), NM))
	r_at_slopegrid = zeros((Nslope2, Ntanr))
	r_at_slopegrid_prior = zeros((Nslope2, Ntanr))
	#for paramlist in self.parameterset.paramlist:
	Ntot = len(paramvalues)
	i = 0
	logger.info("constructing enclosed mass and slope grid")
	print len(rs), NM, Nslope
	print `configuration`
	for parametervalue, p, volume in zip(paramvalues, plist, volumelist):
		#ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
		#dm_profile = scope["dm_density"]
		filename = os.path.join(parametervalue.modelpath, "schw", parametervalue.schwsetname, parametervalue.schwmodelname, "galaxy.ini")
		assert os.path.exists(filename), "filename %s missing" % filename
		#print "1"
		configuration.reset()
		#print "2"
		configuration.readfiles(filename)
		#print "3"
		configuration.init()
		#print "4"
		dm_density = configuration["dm_density"]
		#print filename
		#print dm_density.enclosed_mass(1.)
		#print "4b"
		#dm_profile = configuration["dm_profile"]
		#print "4c"
		if (i % 100) == 0:
			print i, "out of", Ntot
			#print dm_density.rs
		i += 1
		#if os.path.exists(filename):
		if 1:
			logMenc = log10(dm_density.enclosed_mass(rs))
			#print "5"
			Menc_indices = ((logMenc-logMmin)/(logMmax-logMmin) * NM).astype(int)
			#print ".", ~isnan(varvtheta)
			good_Menc_indices = (Menc_indices >= 0) & (Menc_indices < NM)# & (~isnan(varvtheta))
			Menc_indices = Menc_indices[good_Menc_indices]
			#print list(good_Menc_indices), good_Menc_indices.shape
			#volume = 1.
			logMencgrid[good_Menc_indices, Menc_indices] += p * volume
			logMencgrid_prior[good_Menc_indices, Menc_indices] += p * volume * parametervalue.p_prior
			
			slope = dm_density.logslope(rs)
			#print "6"
			#print slope
			slope_indices = ((slope-slopemin)/(slopemax-slopemin) * Nslope).astype(int)
			good_slope_indices = (slope_indices >= 0) & (slope_indices < Nslope)# & (~isnan(varvtheta))
			slope_indices = slope_indices[good_slope_indices]
			slopegrid[good_slope_indices, slope_indices] += p * volume
			slopegrid_prior[good_slope_indices, slope_indices] += p * volume * parametervalue.p_prior
			
			rs_at_slope = dm_density.r_at_slope(slopes)
			mask = isfinite(rs_at_slope) & (rs_at_slope >= 0)
			tanrs = arctan(rs_at_slope[mask])*2/pi
			rs_index = (tanrs*(Ntanr-1)+0.5).astype(int)
			r_at_slopegrid[mask, rs_index] +=  p * volume
			r_at_slopegrid_prior[mask, rs_index] +=  p * volume * parametervalue.p_prior
			
			
	print sort(plist)[-20:]
	volumelist = array(volumelist)
	plist = array(plist)
	print sort(plist*volumelist)[-20:]
	print len(paramvalues), len(plist), len(volumelist)
	return logMencgrid, slopegrid, r_at_slopegrid, logMencgrid_prior, slopegrid_prior, r_at_slopegrid_prior

"""def create_rhoslope_grid(paramvalues, plist, volumelist, logrmin=-3, logrmax=3, slopemin=-1, slopemin=3):
	Nr = 200
	Nslope = 500
	logrs = (arange(Nr) + 0.5 )/ (Nr) * (logrmax - logrmin) + logrmin
	#rs = (arange(Nr) + 0.5 )/ (Nr) * (rmax - rmin) + rmin
	rs = 10**logrs
	slopegrid = zeros((len(logrs), Nslope))
	#for paramlist in self.parameterset.paramlist:
	Ntot = len(paramvalues)
	i = 0
	for parametervalue, p, volume in zip(paramvalues, plist, volumelist):
		ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_profile")
		dm_profile = scope["dm_profile"]
		#print dm_density.alpha
		print i, "out of", Ntot
		i += 1
		#if os.path.exists(filename):
		if 1:
			logMenc = log10(dm_profile.enclosed_mass(rs))
			Menc_indices = ((logMenc-logMmin)/(logMmax-logMmin) * NM).astype(int)
			#print ".", ~isnan(varvtheta)
			good_Menc_indices = (Menc_indices >= 0) & (Menc_indices < NM)# & (~isnan(varvtheta))
			Menc_indices = Menc_indices[good_Menc_indices]
			logMencgrid[good_Menc_indices, Menc_indices] += p * volume
	print sort(plist)[-20:]
	volumelist = array(volumelist)
	plist = array(plist)
	print sort(plist*volumelist)[-20:]
	print len(paramvalues), len(plist), len(volumelist)
	return logMencgrid"""

class ParameterSetCollectHyperOctree(object):
	def __init__(self, modelpath, parameterset, binned_data, betamin, betamax, logrmin, logrmax, logmassmin, logmassmax, logslopemin, logslopemax, nr, rmin, rmax, configuration, sigma_max, **extra_params):
		self.modelpath = modelpath
		self.schwsetname = parameterset.schwsetname
		logschwdirname(self.modelpath, self.schwsetname) 
		self.parameterset = parameterset
		#self.postfix = postfix
		self.binned_data = binned_data
		self.betamin, self.betamax = betamin, betamax
		self.logrmin, self.logrmax = logrmin, logrmax
		self.rmin, self.rmax = 0., 1.5
		self.logmassmin, self.logmassmax = logmassmin, logmassmax
		self.logslopemin, self.logslopemax = logslopemin, logslopemax
		self.dirname = os.path.join(self.modelpath, "schw", self.parameterset.schwsetname)
		self.configuration = configuration
		self.sigma_max = sigma_max
		
		self.paramlabels = self.parameterset.paramlabels
		self.dimension = self.parameterset.dimension
		self.probability_range = self.parameterset.ranges_org
		self.parameterset.load(0)
		
		parametervalue = self.parameterset.paramvalues[0]
		self.nr = nr
		#ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="storage_3d")
		#storage_3d = scope["storage_3d"] 
		#self.rs = storage_3d.x
		
		self.anisotropy_range = (rmin, self.betamin), (rmax, self.betamax)
		self.mass_enclosed_range = (self.logrmin, self.logmassmin), (self.logrmax, self.logmassmax)
		self.logslope_range = (self.logrmin, self.logslopemin), (self.logrmax, self.logslopemax)
		
		self.mass_enclosed_range_linear = (self.rmin, self.logmassmin), (self.rmax, self.logmassmax)
		self.logslope_range_linear = (self.rmin, self.logslopemin), (self.rmax, self.logslopemax)
		
		assert not self.parameterset.isregular, "parameter set should not be a regular set" 
	
	def getpoints(self, filename):
		lines = file(filename).readlines()
		if len(lines) == 0:
			return None
		values = array([[float(k) for k in line.split()] for line in lines])
		print filename, values.shape
		return transpose(values)
	
	def save(self, iteration=0):
		resultdir = os.path.join(self.dirname, "results")
		logger.info("ensuring directory exists: %s" % resultdir)
		if not os.path.exists(resultdir):
			os.makedirs(resultdir)
		
		filename = os.path.join(resultdir, "probability_grid_%03d.npy" % iteration)
		logger.info("storing probability grid: %s" % filename)
		numpy.save(filename, self.probability_grid)
		
		filename = os.path.join(resultdir, "probability_grid_prior_%03d.npy" % iteration)
		logger.info("storing probability grid: %s" % filename)
		numpy.save(filename, self.probability_grid_prior)
		
		filename = os.path.join(resultdir, "anisotropy_grid_%03d.npy" % iteration)
		logger.info("storing anisotropy grid: %s" % filename)
		numpy.save(filename, self.anisotropy_grid)
		
		filename = os.path.join(resultdir, "anisotropy_grid_prior_%03d.npy" % iteration)
		logger.info("storing anisotropy grid: %s" % filename)
		numpy.save(filename, self.anisotropy_grid_prior)
		
		filename = os.path.join(resultdir, "mass_enclosed_grid_%03d.npy" % iteration)
		logger.info("storing mass_enclosed grid: %s" % filename)
		numpy.save(filename, self.mass_enclosed_grid)
	
		filename = os.path.join(resultdir, "mass_enclosed_grid_prior_%03d.npy" % iteration)
		logger.info("storing mass_enclosed grid: %s" % filename)
		numpy.save(filename, self.mass_enclosed_grid_prior)
		
		filename = os.path.join(resultdir, "mass_enclosed_linear_grid_%03d.npy" % iteration)
		logger.info("storing mass_enclosed grid: %s" % filename)
		numpy.save(filename, self.mass_enclosed_grid_linear)
	
		filename = os.path.join(resultdir, "mass_enclosed_linear_grid_prior_%03d.npy" % iteration)
		logger.info("storing mass_enclosed grid: %s" % filename)
		numpy.save(filename, self.mass_enclosed_grid_linear_prior)
	
		filename = os.path.join(resultdir, "parameter_points_%03d.npy" % iteration)
		logger.info("storing parameter_points: %s" % filename)
		numpy.save(filename, self.parameter_points)
		
		filename = os.path.join(resultdir, "logslope_grid_%03d.npy" % iteration)
		logger.info("storing logslope grid: %s" % filename)
		numpy.save(filename, self.logslope_grid)
		
		filename = os.path.join(resultdir, "logslope_grid_prior_%03d.npy" % iteration)
		logger.info("storing logslope grid: %s" % filename)
		numpy.save(filename, self.logslope_grid_prior)
		
		filename = os.path.join(resultdir, "logslope_linear_grid_%03d.npy" % iteration)
		logger.info("storing logslope grid: %s" % filename)
		numpy.save(filename, self.logslope_grid_linear)
		
		filename = os.path.join(resultdir, "logslope_linear_grid_prior_%03d.npy" % iteration)
		logger.info("storing logslope grid: %s" % filename)
		numpy.save(filename, self.logslope_grid_linear_prior)
		
		filename = os.path.join(resultdir, "sigma_grid_prior_%03d.npy" % iteration)
		logger.info("storing sigma grid (prior): %s" % filename)
		numpy.save(filename, self.sigma_grid_prior)
		
		filename = os.path.join(resultdir, "sigma_grid_%03d.npy" % iteration)
		logger.info("storing sigma grid: %s" % filename)
		numpy.save(filename, self.sigma_grid)
		
		filename = os.path.join(resultdir, "r_at_slope_grid_%03d.npy" % iteration)
		logger.info("storing r_at_slope grid: %s" % filename)
		numpy.save(filename, self.r_at_slope_grid)
		
		filename = os.path.join(resultdir, "r_at_slope_grid_prior_%03d.npy" % iteration)
		logger.info("storing r_at_slope (prior) grid: %s" % filename)
		numpy.save(filename, self.r_at_slope_grid_prior)
		
		
	def load(self, iteration=0):
		resultdir = os.path.join(self.dirname, "results")
		
		filename = os.path.join(resultdir, "probability_grid_%03d.npy" % iteration)
		logger.info("loading probability grid: %s" % filename)
		self.probability_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "anisotropy_grid_%03d.npy" % iteration)
		logger.info("loading anisotropy grid: %s" % filename)
		self.anisotropy_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "mass_enclosed_grid_%03d.npy" % iteration)
		logger.info("loading mass_enclosed grid: %s" % filename)
		self.mass_enclosed_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "mass_enclosed_linear_grid_%03d.npy" % iteration)
		logger.info("loading mass_enclosed grid: %s" % filename)
		self.mass_enclosed_grid_linear = numpy.load(filename)
		
		filename = os.path.join(resultdir, "parameter_points_%03d.npy" % iteration)
		logger.info("loading parameter_points: %s" % filename)
		self.parameter_points = numpy.load(filename)
		
		filename = os.path.join(resultdir, "logslope_grid_%03d.npy" % iteration)
		logger.info("loading logslope grid: %s" % filename)
		self.logslope_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "logslope_linear_grid_%03d.npy" % iteration)
		logger.info("loading logslope grid: %s" % filename)
		self.logslope_grid_linear = numpy.load(filename)
		
		# with prior
		filename = os.path.join(resultdir, "probability_grid_prior_%03d.npy" % iteration)
		if not os.path.exists(filename):
			return
		logger.info("loading probability grid: %s" % filename)
		self.probability_grid_prior = numpy.load(filename)
		
		filename = os.path.join(resultdir, "mass_enclosed_grid_prior_%03d.npy" % iteration)
		if not os.path.exists(filename):
			return
		logger.info("loading mass_enclosed grid: %s" % filename)
		self.mass_enclosed_grid_prior = numpy.load(filename)
		
		filename = os.path.join(resultdir, "mass_enclosed_linear_grid_prior_%03d.npy" % iteration)
		logger.info("loading mass_enclosed grid: %s" % filename)
		self.mass_enclosed_grid_linear_prior = numpy.load(filename)
		
		filename = os.path.join(resultdir, "anisotropy_grid_prior_%03d.npy" % iteration)
		logger.info("loading anisotropy grid: %s" % filename)
		self.anisotropy_grid_prior = numpy.load(filename)
		

		filename = os.path.join(resultdir, "logslope_grid_prior_%03d.npy" % iteration)
		logger.info("loading logslope grid: %s" % filename)
		self.logslope_grid_prior = numpy.load(filename)
		
		
		filename = os.path.join(resultdir, "logslope_linear_grid_prior_%03d.npy" % iteration)
		logger.info("loading logslope grid: %s" % filename)
		self.logslope_grid_linear_prior = numpy.load(filename)
		
		filename = os.path.join(resultdir, "sigma_grid_%03d.npy" % iteration)
		if not os.path.exists(filename):
			return
		logger.info("storing sigma grid: %s" % filename)
		self.sigma_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "sigma_grid_prior_%03d.npy" % iteration)
		logger.info("storing sigma grid (prior): %s" % filename)
		self.sigma_grid_prior = numpy.load(filename)
	
		filename = os.path.join(resultdir, "r_at_slope_grid_%03d.npy" % iteration)
		if not os.path.exists(filename):
			return
		logger.info("loading r_at_slope grid: %s" % filename)
		self.r_at_slope_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "r_at_slope_grid_prior_%03d.npy" % iteration)
		logger.info("loading r_at_slope (prior) grid: %s" % filename)
		self.r_at_slope_grid_prior = numpy.load(filename)
	
	def collect(self, iteration=0):
		self.parameterset.load(iteration)
		self.binned_data.load()
		
		inputname = self.parameterset.format_name(iteration)
		setdirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		filename = os.path.join(setdirname, inputname+".solved")
		f = file(filename, "w")
		logger.info("writing file: %s" % filename)
		ps = []
		print len(self.parameterset.paramlist), len(self.parameterset.volumes)
		self.probability_set = "m2"
		#self.probability_set = None
		#self.probability_set = "m2"
		logps = []
		for (name, dirname, orig_values, true_values), volume in zip(self.parameterset.paramlist, self.parameterset.volumes):
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name)
			if self.probability_set is None:
				filename = os.path.join(dirname, "results/logprobability" +self.parameterset.postfix +".txt")
			else:
				filename = os.path.join(dirname, "results/logprobability_seperate" +self.parameterset.postfix +".txt")
			if os.path.exists(filename):
				#print name, file(filename).read().strip()
				if self.probability_set is None:
					logp = eval(file(filename).read().strip())# + 100 #* 1e100
				else:
					values = dict(eval(file(filename).read()))
					logp = values[self.probability_set]
				#print logp
				p = exp(logp)
				if isnan(p):
					print filename
				#ps.append(0)
				#print name, p, logp
				#line = " ".join(["%1.30e" % k for k in list(orig_values) + [volume, logp]])
				#print line
				#print >>f, line
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
				logp = -100000
				#sys.exit(-1)
			logp += 0 #4000 -650+ 1.354397999999999956344254314899e+04-9000
			if isinf(logp):
				logp = -10000
			if isnan(logp):
				logp = -10000
			logps.append(logp)
			line = " ".join(["%1.30e" % k for k in list(orig_values) + [volume, logp]])
			print >>f, line
		f.close()	
		#width, height = 300, 310
		w = 200
		#w = 17
		logps = array(logps)
		print logps
		ps = exp(logps)
		#dsada
		print where(isnan(ps))
		width, height = w, w
		inputname = os.path.join(self.modelpath, "schw", self.schwsetname, self.parameterset.format_name(iteration))
		cmd = "~/src/gdfast/build/src/hyperquadtree -d %d --grid --input=%s " % (self.dimension, inputname)
		for i in range(self.dimension):
			cmd += " %d " % w
		print cmd
		if os.system(cmd) != 0:
			return
		inputname = os.path.join(self.modelpath, "schw", self.schwsetname, self.parameterset.format_name(iteration)+".grid")
		#lines = file(inputname).readlines()
		#if dimension is 
		#values = array([[float(k) for k in line.split()] for line in lines])
		print "reading in data..."
		data = file(inputname).read()
		#print data
		gridlist = eval(data)
		grid = array(gridlist)
		#print grid
		#print gridlist
		
		print grid.shape, width, height
		print "min/max", grid.min(), grid.max(), grid.mean()
		grid[isnan(grid)] = 0
		#grid **= 0.02
		if 0:
			grid = log(grid)
			grid[isinf(grid)] = min(grid[~isinf(grid)])
			print grid.mean()
			print grid
			grid[grid>-510] = -600
		#fsdklfsd
		#grid = scipy.ndimage.gaussian_filter(grid, 4.5)
		#grid = exp(grid)
		self.parameterset.load(iteration, load_unknown=False, load_known=True, load_solved=True)
		paramlist = []
		volumes = []
		logps = []
		#paramnames = [k[0] for k in self.parameterset.paramformats]
		#self.ranges_org = [k[-1] for k in self.paramformats]
		lines = []
		filebase = os.path.join(self.modelpath, "schw", self.schwsetname, self.parameterset.format_name(iteration))
		#if not self.parameterset.doall:
		filename = filebase+".known"
		lines.extend(file(filename).readlines())
		filename = filebase+".solved"
		lines.extend(file(filename).readlines())
		points = []
		pointlist = []
		for line in lines:
			values = [float(k) for k in line.split()]
			values, volume, logp = values[:-2], values[-2], values[-1]
			#p = exp(0)
			points.append(values)
			pointlist.append(values)
			newvalues = []
			dirname = os.path.join(self.modelpath, "schw", self.parameterset.schwsetname)
			names = []
			#print "lenghs", len(values), len(self.parameterset.paramformats)
			assert len(values) == len(self.parameterset.paramformats)
			p_prior = 1.
			for value, paramformat in zip(values, self.parameterset.paramformats):
				param_name, param_format, param_transformation, param_range = paramformat[:4]
				param_prior = paramformat[6]
				p_prior *= param_prior(value)
				newvalue = param_transformation(value)
				name = param_format % newvalue
				newvalues.append(newvalue)
				names.append(name)
			name = "/".join(names)
			dirname = os.path.join(dirname, name)
			paramlist.append((name, dirname, values, newvalues, p_prior))
			volumes.append(volume)
			#xs.append(values[0])
			#ys.append(values[1])
			logps.append(logp)
			#print line
			#sys.exit()
		logps = array(logps)
		logps -= logps.max()
		#print logps
		ps = exp(logps)
		points = transpose(array(points))
		
		self.parameter_points = points
		#grid = scipy.ndimage.gaussian_filter(grid, 1.5)
		
		
		self.probability_grid = grid
		self.probability_grid_prior = grid * 1.0
		
		import itertools
		for index in itertools.product(range(w), repeat=self.dimension):
			p_prior = 1.
			for i, paramformat in zip(index, self.parameterset.paramformats):
				param_name, param_format, param_transformation, param_range = paramformat[:4]
				value = param_range[0] + (param_range[1] - param_range[0]) * i / (w-1.)
				param_prior = paramformat[6]
				p_prior *= param_prior(value)
			#numpy.put(self.probability_grid_prior, [index], [
			#print index, p_prior
			self.probability_grid_prior.__setitem__(index, p_prior*self.probability_grid_prior.__getitem__(index))
				
		
			
		#self.parameterset.ranges_org[i][0], self.parameterset.ranges_org[j][0]), (self.parameterset.ranges_org[i][1]
		#
		#print dir(self.binned_data)
		self.sigma_grid, self.sigma_grid_prior = sigma_grid(self.parameterset.paramvalues, ps, volumes, self.sigma_max, NR=self.binned_data.n_constraints, postfix=self.parameterset.postfix)
		#return
		#nr = len(self.rs)
		self.anisotropy_grid, self.anisotropy_grid_prior = create_anisotropy_grid(paramlist, ps, volumes, self.nr, postfix=self.parameterset.postfix, betamin=self.betamin, betamax=self.betamax)
		
		
		#self.mass_enclosed_grid = create_Menc_grid(self.parameterset.paramvalues, ps, volumes, -2, 1, 2.5, 9)
		self.mass_enclosed_grid, self.logslope_grid, self.r_at_slope_grid, self.mass_enclosed_grid_prior, self.logslope_grid_prior, self.r_at_slope_grid_prior = \
			create_Menc_and_slope_grid(self.parameterset.paramvalues, ps, volumes, self.logrmin, self.logrmax, self.logmassmin, self.logmassmax, self.logslopemin, self.logslopemax, configuration=self.configuration)

		self.mass_enclosed_grid_linear, self.logslope_grid_linear, _, self.mass_enclosed_grid_linear_prior, self.logslope_grid_linear_prior, _ = \
			create_Menc_and_slope_grid(self.parameterset.paramvalues, ps, volumes, self.logrmin, self.logrmax, self.logmassmin, self.logmassmax, self.logslopemin, self.logslopemax, logarithmic=False, rmin=self.rmin, rmax=self.rmax, configuration=self.configuration)
			

		
		print "iteration", iteration
		return
		
		
		#sys.exit(0)
		
		
		#f = open(filename, "wb")
		#pickle.dump((
		
		#box()
		#vsplit(box)
		print dimension
		if 0:
			document(size="45cm,25cm")
			page(fontsize="10pt")
			mozaic(3,3,box)
			
			filename = os.path.join(self.dirname, "selection.txt")
			print filename
			
			if os.path.exists(filename):
				#print filename, "exists"
				selections = file(filename).readlines()
				selections = [k.strip() for k in selections if k.strip()]
				indices = []
				print len(selections)
				#print selections
				for i in range(len(pointlist)):
					#print paramlist[i][0]
					if paramlist[i][0] in selections:
						x, y = pointlist[i]
						#print logps[i]
						indices.append(i)
						#print x, y
						#symbol(x, y, symbolname="circle", color="red")
				print len(indices)
				#for param, p, volume in zip(paramlist, plist, volumelist):
				#for name, dirname, indices, values in self.paramlist[1:]:
				parametervalue = self.parameterset.paramvalues[0]
				ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="storage_3d")
				storage_3d = scope["storage_3d"]
				
				
				r = storage_3d.x
				rborders = storage_3d.xborders
				ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="aperture")
				aperture = scope["aperture"]
				aperture.load()
				
				ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="light_model")
				light_model= scope["light_model"]
				R = light_model.arcsec_to_kpc(aperture.aperture_rcenters)
				
				ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="binned_data")
				binned_data = scope["binned_data"]
				binned_data.load()
				
				mass_r = array([light_model.cumdensityr(r1, r2) for r1, r2 in zip(rborders[:-1], rborders[1:])])
				select(0, 0)
				if 0:
					graph(log10(r), log10(mass_r), color="blue")
					vfill(-2, 0, color="blue", alpha=0.1)
					vfill(-1, 0, color="blue", alpha=0.1)
				
				
				select(0, 0)
				vfill(-2, 0, color="blue", alpha=0.1)
				vfill(-1, 0, color="blue", alpha=0.1)
				select(0, 1)
				vfill(-2, 0, color="blue", alpha=0.1)
				vfill(-1, 0, color="blue", alpha=0.1)
				
				chisqsm2 = []
				chisqsm4 = []
				for i, index in enumerate(indices):
					
					#print p
					color = Color(float(i)/len(indices), 0, 0)
					name, dirname, orig_values, true_values = paramlist[index][:4]
					#print name, dirname
					
					schwmodelname = name
					ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, schwmodelname, objectname="dm_density_twoslope")
					dm_density = scope["dm_density_twoslope"]
					
					if 1:
						select(0,0)
						logrs = arange(-2, 1., 0.05)
						rs = 10**logrs
						Ms = dm_density.enclosed_mass(rs)
						graph(logrs, log10(Ms), color=color)
						
						select(0,1)
						logrs = arange(-2, 1., 0.05)
						rs = 10**logrs
						rho = dm_density.densityr(rs)
						graph(logrs, log10(rho), color=color)
					
					if 1:
						filename = os.path.join(dirname, "results", "solution_moments3d" +self.parameterset.postfix +".npy")
						#print filename
						if os.path.exists(filename):
							#print "\texists"
							solution_moments3d = numpy.load(filename)
							varvr = solution_moments3d[4]
							varvphi = solution_moments3d[5]
							varvtheta = solution_moments3d[6]
							betas = 1 - (varvphi + varvtheta)/(2*varvr)
							mask = ~isnan(betas)
							#graph(arange(len(betas))[mask], betas[mask])
							select(1,0)
							graph((r[mask]), betas[mask], color=color)
							#xlim(0, 1.5)
							xlim(0, 1.5)
							ylim(-2,1)
							#select(i % 6, i/6)
							#graph((r[mask]), solution_moments3d[0][mask], color=color)
							select(0, 1)
							#graph((r[mask]), log10(solution_moments3d[0][mask]), color=color)
							
							#print betas
							if 0:
								print betas
								print varvr
								print varvphi
								print varvtheta
								betaindices = ((betas-betamin)/(betamax-betamin) * Nbeta).astype(int)
								#print ".", ~isnan(varvtheta)
								good_betaindices = (betaindices >= 0) & (betaindices < Nbeta) & (~isnan(varvtheta))
								betaindices = betaindices[good_betaindices]
								#p = pgrid[indices[0], indices[1]]
								#print good_betaindices,betaindices
								#print p, volume
								betagrid[good_betaindices,betaindices] += p * volume
					if 1:
						filename = os.path.join(dirname, "results", "solution_projectedmoments" +self.parameterset.postfix +".npy")
						print filename
						assert os.path.exists(filename)
						if os.path.exists(filename):
							#print "\texists"
							solution_projectedmoments = numpy.load(filename)
							#import pdb
							#pdb.set_trace()
							sigmalos = solution_projectedmoments[2]**0.5
							m2 = solution_projectedmoments[2]
							m4 = solution_projectedmoments[4]
							kappa = m4/m2**2
							#varvphi = solution_moments3d[5]
							#varvtheta = solution_moments3d[6]
							#betas = 1 - (varvphi + varvtheta)/(2*varvr)
							#mask = ~isnan(betas)
							#graph(arange(len(betas))[mask], betas[mask])
							#print R
							#print sigmalos
							#select(i % 6, i/6)
							#graph(R, sigmalos)
							#print sigmalos.shape, binned_data.moments.shape, R.shape, kappa.shape
							#print R
							if 0:
								filename = os.path.join(dirname, "results/logprobability_seperate" +self.parameterset.postfix +".txt")
								logpvalues = dict(eval(file(filename).read()))
								#logp = values[self.probability_set]
								#chis
								
								chisq = sum(((solution_projectedmoments[2] - binned_data.moments[2])/binned_data.e_moments[2])**2)
								print chisq/len(solution_projectedmoments[2]),
								#chisqsm2.append(chisq/len(solution_projectedmoments[2]))
								chisqsm2.append(logpvalues["m2"])
								chisq = sum(((solution_projectedmoments[4] - binned_data.moments[4])/binned_data.e_moments[4])**2)
								print chisq/len(solution_projectedmoments[2])
								#chisqsm4.append(chisq/len(solution_projectedmoments[2]))
								chisqsm4.append(logpvalues["m4"])
								print logpvalues
							#print chisq, logps[index]
							#print sum((solution_projectedmoments[0] - binned_data.moments[0])**2)
							#xlim(0, 1.5)
							#ylim(0, 25)
							select(1,1)
							graph(R, sigmalos, color=color)
							scatter(R, binned_data.moments[2]**0.5, symbolName="circlesolid", color="blue")
							xlim(0, 1.5)
							ylim(0, 25)
							select(1,2)
							graph(R, kappa, color=color)
							print kappa
							
							m2 = binned_data.moments[2]
							m4 = binned_data.moments[4]
							e_m2 = binned_data.e_moments[2]
							e_m4 = binned_data.e_moments[4]
							kappa = m4/m2**2
							e_kappa = sqrt(1/m2**4 * e_m4**2 + (2 * m4/m2**4)**2*e_m2**2)
							
							scatter(R, m4/m2**2, symbolName="circlesolid", color="blue")
							errorbars(R, kappa, yerr=e_kappa)

							ylim(-1, 5)
							xlim(0,1.5)
							
							select(2,2)
							m2 = solution_projectedmoments[2]
							m4 = solution_projectedmoments[4]
							graph(R, m4, color=color)
							scatter(R, binned_data.moments[4], symbolName="circlesolid", color="blue")
							errorbars(R, binned_data.moments[4], yerr=binned_data.e_moments[4])
							#print betas
					
					
				select(2,0)
				#chisqsm2 = array(chisqsm2)
				#chisqsm4 = array(chisqsm4)
				#graph(chisqsm2, color="red")
				#graph(chisqsm4, color="green")
				#graph(chisqsm2+chisqsm4/5, color="blue")
				#ylim(0, max(chisqsm2+chisqsm4))
			draw()
			sys.exit(0) 
		if self.dimension == 2:
			document(size="25cm,20cm")
			#document(size="15cm,15cm")
			page(fontsize="16pt")
			mozaic(3,2,box)
			def makeresize(i, j):
				return (self.parameterset.ranges_org[i][0], self.parameterset.ranges_org[j][0]), (self.parameterset.ranges_org[i][1], self.parameterset.ranges_org[j][1])
			if 1:
				select(0, 0)
				probimage2d(grid, 0, 1, resize=makeresize(0, 1), colormap="whiteblue", drawcontourlines=True)
				scatter(points[0], points[1], color="black")#, symbolsize="1pt")
				filename = os.path.join(self.dirname, "selection.txt")
				#print filename
				selections = []
				if os.path.exists(filename):
					#print filename, "exists"
					selections = file(filename).readlines()
					selections = [k.strip() for k in selections if k.strip()]
					#print selections
					for i in range(len(pointlist)):
						#print paramlist[i][0]
						if paramlist[i][0] in selections:
							selections.append(paramlist[i][0])
							x, y = pointlist[i]
							print logps[i]
							#print x, y
							symbol(x, y, symbolname="circle", color="red") 

				labels(self.parameterset.paramlabels[0], self.parameterset.paramlabels[1]) 
				#labels("M<sub>1kpc</sub> - solar masses", "log rs - kpc")
				def onmouse(x, y, options, window):
					current.container = curr
					vx, vy = current.container.getDocument().windowToViewport(x, y, current.container.borderViewport)
					#px = kaplot.utils.convertPixelsTo(x, "cm")
					#py = kaplot.utils.convertPixelsTo(y, "cm")
					wx, wy = current.container.windowToWorld(x, y)
					#print wx, wy
					
					if options["leftdown"]:
						nearest_index = 0
						rmin = 1e6
						for i, point in enumerate(pointlist):
							#print i, point
							x, y = point
							r = sqrt((x-wx)**2+(y-wy)**2)
							if r < rmin:
								nearest_index = i
								rmin = r
						print "->", paramlist[nearest_index], pointlist[nearest_index], wx, wy, vx, vy
						selections.append(paramlist[nearest_index][0])
						# xs[nearest_index], ys[nearest_index], ps[nearest_index]
					
					return False # dont block
				curr = current.container
				current.container.onmouse = onmouse
				def onkey(x, y, keycode, character, options, window):
					if keycode == 25:
						filename = os.path.join(self.dirname, "selection.txt")
						print "writing to", filename
						f = open(filename, "w")
						for selection in selections:
							print >>f, selection
					return False
				current.container.handleKeyboardEvent = onkey 
				
				if 1:
					select(2, 0)
					probgraph(grid, 0, resize=self.parameterset.ranges_org[0])
					labels(self.parameterset.paramlabels[0], "p")
					
					select(2, 1)
					probgraph(grid, 1, resize=self.parameterset.ranges_org[1])
					labels(self.parameterset.paramlabels[1], "p")
					
					#select(2, 0)
					#probgraph(grid, 2, resize=self.parameterset.ranges_org[2])
					#labels(self.parameterset.paramlabels[2], "p")
				
				
				if 0:
					select(1, 0)
					probimage2d(grid, 1, 2, resize=makeresize(1, 2))
					scatter(points[1], points[2])
					labels(self.parameterset.paramlabels[1], self.parameterset.paramlabels[2])
					
					select(0, 1)
					probimage2d(grid, 0, 1, resize=makeresize(0, 1))
					scatter(points[0], points[1])
					labels(self.parameterset.paramlabels[0], self.parameterset.paramlabels[1]) 
					
			if 0:
				select(0,1)
				title("Density profile", fontsize="15pt") 
				labels("log r - kpc", "&rho;(r)", fontsize="15pt")
				
				#logr = arange(-3, 3, 0.01)
				#r = 10**logr
				#for i, parametervalue in enumerate(self.parameterset.paramvalues):
				#	ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
				#	dm_density = scope["dm_density"]
				#	print dm_density.alpha
				#	print i, "out of", len(self.parameterset.paramvalues)
				rhogrid = create_rho_grid(self.parameterset.paramvalues, ps, volumes, -3, 3, -7, 12)
				s = 2.
				rhogrid = scipy.ndimage.gaussian_filter(rhogrid, [s, s])
				if 0:
					for i in range(rhogrid.shape[1]):
						if sum(rhogrid[:,i]) > 0:
							rhogrid[:,i] /= rhogrid[:,i].sum()
				if 1:
					for i in range(rhogrid.shape[0]):
						if sum(rhogrid[i]) > 0:
							rhogrid[i] /= rhogrid[i].max()
				
				import kaplot
				probimage2d(rhogrid, 0, 1, resize=((-3, -7), (3, 12)), colormap="whiterainbow")
				#kaplot.grid(xinterval=1., yinterval=1., color="lightgrey")
				kaplot.grid(xinterval=1., yinterval=1., alpha=0.2)
				
			if 0:
				select(1,1)
				title("Enclosed mass", fontsize="15pt") 
				labels("log r - kpc", "M(&lt;r)", fontsize="15pt")
				
				#logr = arange(-3, 3, 0.01)
				#r = 10**logr
				#for i, parametervalue in enumerate(self.parameterset.paramvalues):
				#	ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
				#	dm_density = scope["dm_density"]
				#	print dm_density.alpha
				#	print i, "out of", len(self.parameterset.paramvalues)
				Menc_grid = create_Menc_grid(self.parameterset.paramvalues, ps, volumes, -2, 1, 2.5, 9)
				s = 1.
				Menc_grid = scipy.ndimage.gaussian_filter(Menc_grid, [s, s])
				if 0:
					for i in range(Menc_grid.shape[1]):
						if sum(Menc_grid[:,i]) > 0:
							Menc_grid[:,i] /= Menc_grid[:,i].sum()
				if 1:
					for i in range(Menc_grid.shape[0]):
						if sum(Menc_grid[i]) > 0:
							Menc_grid[i] /= Menc_grid[i].max()
				
				import kaplot
				probimage2d(Menc_grid, 0, 1, resize=((-2, 2.5), (1, 9)), colormap="whiterainbow")
				#kaplot.grid(xinterval=1., yinterval=1., color="lightgrey")
				kaplot.grid(xinterval=1., yinterval=1., alpha=0.2)				
				
			if 0:
				select(1,1)
				title("Enclosed mass", fontsize="15pt") 
				labels("log r - kpc", "M(&lt;r)", fontsize="15pt")
				
				#logr = arange(-3, 3, 0.01)
				#r = 10**logr
				#for i, parametervalue in enumerate(self.parameterset.paramvalues):
				#	ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
				#	dm_density = scope["dm_density"]
				#	print dm_density.alpha
				#	print i, "out of", len(self.parameterset.paramvalues)
				sorted = argsort(ps)[::-1]
				N = 10
				logrs = arange(-3, 3, 0.1)
				rs = 10**logrs
				#rs = arange(1e-4, 1.5, 0.01)
				#logrs = log10(rs)
				indices = sorted[:N]
				print indices, len(ps), len(self.parameterset.paramvalues)
				for i in indices:
					parametervalue = self.parameterset.paramvalues[i]
					#p
					#volume in zip(paramvalues, plist, volumelist):
					ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_profile")
					dm_density = scope["dm_profile"]
					print dm_density.alpha
					print i, "out of", N
					i += 1
					#if os.path.exists(filename):
					if 1:
						logMenc = log10(dm_density.enclosed_mass(rs))
						#Menc_indices = ((logMenc-logMmin)/(logMmax-logMmin) * NM).astype(int)
						#print ".", ~isnan(varvtheta)
						#good_Menc_indices = (Menc_indices >= 0) & (Menc_indices < NM)# & (~isnan(varvtheta))
						#Menc_indices = Menc_indices[good_Menc_indices]
						#logMencgrid[good_Menc_indices, Menc_indices] += p * volume
						graph(logrs,logMenc)
					
					
				if 0:
					Menc_grid = create_Menc_grid(self.parameterset.paramvalues, ps, volumes, -2, 1, 2.5, 9)
					s = 1.
					Menc_grid = scipy.ndimage.gaussian_filter(Menc_grid, [s, s])
					if 0:
						for i in range(Menc_grid.shape[1]):
							if sum(Menc_grid[:,i]) > 0:
								Menc_grid[:,i] /= Menc_grid[:,i].sum()
					if 1:
						for i in range(Menc_grid.shape[0]):
							if sum(Menc_grid[i]) > 0:
								Menc_grid[i] /= Menc_grid[i].max()
					
					import kaplot
					probimage2d(Menc_grid, 0, 1, resize=((-2, 2.5), (1, 9)), colormap="whiterainbow")
					#kaplot.grid(xinterval=1., yinterval=1., color="lightgrey")
					kaplot.grid(xinterval=1., yinterval=1., alpha=0.2)				
				
				#print ps[indices]
				
			if 1:
				select(1,0)
				print len(paramlist), len(ps), len(volumes)
				#jdklsa
				betagrid = create_anisotropy_grid(paramlist, ps, volumes, nr, postfix=self.parameterset.postfix)
				s = 3.01
				betagrid = scipy.ndimage.gaussian_filter(betagrid, [s*0.2, s])
				if 0:
					for i in range(betagrid.shape[1]):
						if sum(betagrid[:,i]) > 0:
							betagrid[:,i] /= betagrid[:,i].sum()
				if 1:
					for i in range(betagrid.shape[0]):
						if sum(betagrid[i]) > 0:
							betagrid[i] /= betagrid[i].max()
				print betagrid.min(), betagrid.max()
				#probimage2d(betagrid, 0, 1, resize=((0, -3), (1.5, 1)), colormap="whiteblue", drawcontourlines=True)
				probimage2d(betagrid, 0, 1, resize=((-3, -3), (3, 1)), colormap="whiteblue", drawcontourlines=False)
				#print betagrid[0][0]
				#indexedimage(betagrid, colormap="whiteblack")
				hline(-0.5)
				#title("Anisotropy parameter as\na function or radius.", fontsize="15pt")
				vline(-1.5)
				vline(0.0) 
				labels("r - kpc", "anisotropy", fontsize="15pt")

			draw()
		if self.dimension == 3:
			document(size="30cm,20cm")
			if usevtk:
				image = vtkImageData()
				image.SetDimensions(w, w, w)
				image.SetScalarTypeToDouble ()
				image.SetNumberOfScalarComponents(1)
				for i in range(w):
					for j in range(w):
						for k in range(w):
							image.SetScalarComponentFromDouble(i, j, k, 0, grid[i,j,k])
				if 0:
					image.SetSpacing(
							self.parameterset.ranges_org[0][1]-self.parameterset.ranges_org[0][0],
							self.parameterset.ranges_org[1][1]-self.parameterset.ranges_org[1][0],
							self.parameterset.ranges_org[2][1]-self.parameterset.ranges_org[2][0],
							)
					image.SetOrigin(
							self.parameterset.ranges_org[0][0],
							self.parameterset.ranges_org[1][0],
							self.parameterset.ranges_org[2][0],
							)
					
				
				writer = vtkXMLImageDataWriter()
				writer.SetFileName("vtk-image-test.vti");
				writer.SetInput(image);
				writer.Write();
				#sys.exit(0)
				#"""
		
				
				
			mozaic(3,3,box)
			def makeresize(i, j):
				return (self.parameterset.ranges_org[i][0], self.parameterset.ranges_org[j][0]), (self.parameterset.ranges_org[i][1], self.parameterset.ranges_org[j][1])
			if 1:
				select(0, 0)
				color = None
				colormap = "whiteblue"
				drawcontourlines = True
				probimage2d(grid, 0, 2, resize=makeresize(0, 2), color=color, colormap=colormap, drawcontourlines=drawcontourlines)
				scatter(points[0], points[2])
				labels(self.parameterset.paramlabels[0], self.parameterset.paramlabels[2]) 
				
				select(1, 0)
				probimage2d(grid, 1, 2, resize=makeresize(1, 2), color=color, colormap=colormap, drawcontourlines=drawcontourlines)
				scatter(points[1], points[2])
				labels(self.parameterset.paramlabels[1], self.parameterset.paramlabels[2])
				
				select(0, 1)
				probimage2d(grid, 0, 1, resize=makeresize(0, 1), color=color, colormap=colormap, drawcontourlines=drawcontourlines)
				scatter(points[0], points[1])
				labels(self.parameterset.paramlabels[0], self.parameterset.paramlabels[1]) 
				
				select(0, 2)
				probgraph(grid, 0, resize=self.parameterset.ranges_org[0])
				labels(self.parameterset.paramlabels[0], "p")
				
				select(1, 1)
				probgraph(grid, 1, resize=self.parameterset.ranges_org[1])
				labels(self.parameterset.paramlabels[1], "p")
				
				select(2, 0)
				probgraph(grid, 2, resize=self.parameterset.ranges_org[2])
				labels(self.parameterset.paramlabels[2], "p")
				
			if 0:
				select(1,2)
				#title("Density profile", fontsize="15pt") 
				labels("log r - kpc", "&rho;(r)", fontsize="15pt")
				
				#logr = arange(-3, 3, 0.01)
				#r = 10**logr
				#for i, parametervalue in enumerate(self.parameterset.paramvalues):
				#	ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
				#	dm_density = scope["dm_density"]
				#	print dm_density.alpha
				#	print i, "out of", len(self.parameterset.paramvalues)
				rhogrid = create_rho_grid(self.parameterset.paramvalues, ps, volumes, -3, 3, -7, 12)
				s = 2.
				rhogrid = scipy.ndimage.gaussian_filter(rhogrid, [s, s])
				if 0:
					for i in range(rhogrid.shape[1]):
						if sum(rhogrid[:,i]) > 0:
							rhogrid[:,i] /= rhogrid[:,i].sum()
				if 1:
					for i in range(rhogrid.shape[0]):
						if sum(rhogrid[i]) > 0:
							rhogrid[i] /= rhogrid[i].max()
				
				import kaplot
				probimage2d(rhogrid, 0, 1, resize=((-3, -7), (3, 12)), colormap="whiterainbow")
				#kaplot.grid(xinterval=1., yinterval=1., color="lightgrey")
				kaplot.grid(xinterval=1., yinterval=1., alpha=0.2)
				
			if 0:
				select(2,1)
				title("Enclosed mass", fontsize="15pt") 
				labels("log r - kpc", "M(&lt;r)", fontsize="15pt")
				
				#logr = arange(-3, 3, 0.01)
				#r = 10**logr
				#for i, parametervalue in enumerate(self.parameterset.paramvalues):
				#	ini, scope = mab.gd.schw.configure.loadini(parametervalue.modelpath, parametervalue.schwsetname, parametervalue.schwmodelname, objectname="dm_density")
				#	dm_density = scope["dm_density"]
				#	print dm_density.alpha
				#	print i, "out of", len(self.parameterset.paramvalues)
				Menc_grid = create_Menc_grid(self.parameterset.paramvalues, ps, volumes, -3, 3, 0, 13)
				s = 2.
				Menc_grid = scipy.ndimage.gaussian_filter(Menc_grid, [s, s])
				if 0:
					for i in range(Menc_grid.shape[1]):
						if sum(Menc_grid[:,i]) > 0:
							Menc_grid[:,i] /= Menc_grid[:,i].sum()
				if 1:
					for i in range(Menc_grid.shape[0]):
						if sum(Menc_grid[i]) > 0:
							Menc_grid[i] /= Menc_grid[i].sum()
				
				import kaplot
				probimage2d(Menc_grid, 0, 1, resize=((-3, 0), (3, 13)), colormap="whiterainbow")
				#kaplot.grid(xinterval=1., yinterval=1., color="lightgrey")
				kaplot.grid(xinterval=1., yinterval=1., alpha=0.2)				
				
				
				
			if 0:
				select(2,2)
				betagrid = create_anisotropy_grid(paramlist, ps, volumes, postfix=self.parameterset.postfix)
				s = 5.1
				betagrid = scipy.ndimage.gaussian_filter(betagrid, [s*0.2, s])
				if 0:
					for i in range(betagrid.shape[1]):
						if sum(betagrid[:,i]) > 0:
							betagrid[:,i] /= betagrid[:,i].sum()
				if 0:
					for i in range(betagrid.shape[0]):
						if sum(betagrid[i]) > 0:
							betagrid[i] /= betagrid[i].sum()
				print betagrid.min(), betagrid.max()
				probimage2d(betagrid, 0, 1, resize=((0, -3), (1.5, 1)), colormap="whiterainbow")
				#print betagrid[0][0]
				#indexedimage(betagrid, colormap="whiteblack")
				hline(-0.5)
				#title("Anisotropy parameter as\na function or radius.", fontsize="15pt") 
				labels("r - kpc", "&beta; - anisotropy", fontsize="15pt")
				
				
				
				
			if 0:
				select(2, 2)
				probimage2d(grid[:,:,-int(6*70/17)], 0, 1, resize=makeresize(0, 1))
				scatter(points[0], points[1])
				labels(self.parameterset.paramnames[0], self.parameterset.paramnames[1])
				select(1, 2)
				probimage2d(grid[:,:,-1], 0, 1, resize=makeresize(0, 1))
				scatter(points[0], points[1])
				labels(self.parameterset.paramnames[0], self.parameterset.paramnames[1])
				select(2, 1)
				probimage2d(grid[:,:,0], 0, 1, resize=makeresize(0, 1))
				scatter(points[0], points[1])
				labels(self.parameterset.paramnames[0], self.parameterset.paramnames[1])
			
			
			if 0:
				for i in range(3): 
					select(i, 2)
					probimage2d(grid[:,:,8+i], 0, 1, resize=makeresize(0, 1))
					scatter(points[0], points[1])
					labels(self.parameterset.paramnames[0], self.parameterset.paramnames[1]) 
			draw()
		if 0: #dimension == 2:
			mozaic(2,2,box)
			x1, x2 = self.parameterset.ranges_org[0]
			y1, y2 = self.parameterset.ranges_org[1]
			#probimage2d(values, 0, 1, arange(width), arange(height)) #, colormap="whiteblack")
			resize = (x1, y1), (x2, y2)
			#values = scipy.ndimage.gaussian_filter(values, 4.0)
			#indexedimage(transpose(values), colormap="whiterainbow", resize=resize)
			#indexedimage(transpose(values), colormap="whiterainbow", resize=resize)
			#values *= 1e55
			values = grid
			print values.min(), values.max()
			im = probimage2d(values, 0, 1, resize=resize, color="blue", drawcontourlines=True) #, colormap="whiteblack")
			filename = os.path.join(setdirname, self.parameterset.inputname+".solved")
			x, y, volumes, values = self.getpoints(filename)
			n_models = len(x)
			scatter(x, y)
			
			# this file is empty the first time
			filename = os.path.join(setdirname, self.parameterset.inputname+".known")
			res = self.getpoints(filename)
			if res != None:
				x, y, volumes, values = res
				n_models += len(x)
				scatter(x, y)
			print "# models", n_models
			title("Joined pdf of NFW mass vs scale radius.", fontsize="15pt") 
			labels("log(M<sub>DM</sub>) - solar masses", "log(r<sub>s</sub>) - kpc", fontsize="15pt")
			innercolorbar(im, colormap="whiterainbow")
			
			#paramlist = []
			#volumes = []
			#ps = []
			#paramnames = [k[0] for k in self.parameterset.paramformats]
			#self.ranges_org = [k[-1] for k in self.paramformats]
			#lines = []
			#filebase = os.path.join(self.modelpath, "schw", self.schwsetname, self.parameterset.inputname)
			#filename = filebase+".known"
			#lines.extend(file(filename).readlines())
			#filename = filebase+".solved"
			#lines.extend(file(filename).readlines())
			if 0:
				xs = []
				ys = []
				for line in lines:
					values = [float(k) for k in line.split()]
					values, volume, p = values[:-2], values[-2], values[-1]
					newvalues = []
					dirname = os.path.join(self.modelpath, "schw", self.parameterset.schwsetname)
					names = []
					for value, paramformat in zip(values, self.parameterset.paramformats):
						param_name, param_format, param_transformation, param_range = paramformat
						newvalue = param_transformation(value)
						name = param_format % newvalue
						newvalues.append(newvalue)
						names.append(name)
					name = "/".join(names)
					dirname = os.path.join(dirname, name)
					paramlist.append((name, dirname, values, newvalues))
					volumes.append(volume)
					xs.append(values[0])
					ys.append(values[1])
					ps.append(p)
				
			#scales = array(ps)
			#scales /= scales.max()
			self.binned_data.load()
			self.binned_data.aperture.load()
			if 0:
				for scale, line in zip(scales, lines):
					values = [float(k) for k in line.split()]
					values, volume, p = values[:-2], values[-2], values[-1]
					newvalues = []
					dirname = os.path.join(self.modelpath, "schw", self.parameterset.schwsetname)
					names = []
					for value, paramformat in zip(values, self.parameterset.paramformats):
						param_name, param_format, param_transformation, param_range = paramformat
						newvalue = param_transformation(value)
						name = param_format % newvalue
						newvalues.append(newvalue)
						names.append(name)
					name = "/".join(names)
					dirname = os.path.join(dirname, name)
					select(0,1)
					
					filename = os.path.join(dirname, "results", "solution_projectedmoments"+self.parameterset.postfix+".npy")
					moments = load(filename)
					#print moments.shape
					sigma_los = sqrt(moments[2])
					index = arange(len(sigma_los))
					indices = self.binned_data.aperture.aperture_rborders
					#scatter(indices, sigma_los, symbolName="squaresolid", alpha=scale/2,color="red")
					
				sigma_los = sqrt(self.binned_data.moments[2])
				index = arange(len(sigma_los))
				scatter(indices, sigma_los, symbolName="squaresolid", alpha=0.8)
				sigma_max = max(sigma_los) * 1.4
				ylim(0, sigma_max)
				
			#print ps, len(ps), len(volumes), len(paramlist)
			def onmouse(x, y, options, window):
				vx, vy = current.container.getDocument().windowToViewport(x, y, current.container.borderViewport)
				#px = kaplot.utils.convertPixelsTo(x, "cm")
				#py = kaplot.utils.convertPixelsTo(y, "cm")
				wx, wy = current.container.windowToWorld(x, y)
				#print wx, wy
				if options["leftdouble"]:
					nearest_index = 0
					rmin = 1e6
					for i, (x, y) in enumerate(zip(xs, ys)):
						r = sqrt((x-wx)**2+(y-wy)**2)
						if r < rmin:
							nearest_index = i
							rmin = r
							
					print "->", xs[nearest_index], ys[nearest_index], ps[nearest_index]
				
				return False # dont block
			current.container.onmouse = onmouse
			select(1,0)
			print ps
			betagrid = create_anisotropy_grid(paramlist, ps, volumes, postfix=self.parameterset.postfix)
			s = 7
			print betagrid.min(), betagrid.max()
			print betagrid[betagrid>betagrid.min()]
			betagrid = scipy.ndimage.gaussian_filter(betagrid, [s*0.1, s])
			if 0:
				for i in range(betagrid.shape[1]):
					if sum(betagrid[:,i]) > 0:
						betagrid[:,i] /= betagrid[:,i].sum()
			if 1:
				for i in range(betagrid.shape[0]):
					if sum(betagrid[i]) > 0:
						betagrid[i] /= betagrid[i].sum()
			
			probimage2d(betagrid, 0, 1, resize=((0, -3), (1.5, 1)), colormap="whiterainbow", drawcontourlines=False)
			hline(-0.5)
			title("Anisotropy parameter as a function or radius.", fontsize="15pt") 
			labels("r - kpc", "&beta; - anisotropy", fontsize="15pt")
			select(0)
			draw()
		#select(1)
		#print len(lines)
			#print indices, p"""			
			
class ParameterSetCollectRegular(object):
	def __init__(self, modelpath, parameterset, **extra_params):
		self.modelpath = modelpath
		self.schwsetname = parameterset.schwsetname
		self.parameterset = parameterset
		logschwdirname(self.modelpath, self.schwsetname) 
		#self.postfix = postfix
		assert self.parameterset.isregular, "parameter set should be a regular set" 
		
		
	def collect(self):
		shape = self.parameterset.shape
		print shape
		pgrid = numpy.zeros(tuple(shape))
		logpgrid = numpy.zeros(tuple(shape)) - 1e500
		#betas =
			
		for paramlist in self.parameterset.paramlist:
			name, dirname, orig_values, true_values, indices = paramlist
			print indices[0:3]
			#dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name_new)
			#print dirname, self.modelpath, self.schwsetname, name_new
			
			filename = os.path.join(dirname, "results/probability" +self.parameterset.postfix +".txt")
			if os.path.exists(filename):
				p = eval(file(filename).read().strip())
				#pgrid[indices[0], indices[1], indices[2]] = p
				pgrid.__setitem__(tuple(indices), p)
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
			#print indices, p
			
			filename = os.path.join(dirname, "results/logprobability" +self.parameterset.postfix +".txt")
			if os.path.exists(filename):
				logp = eval(file(filename).read().strip())
				#logpgrid[indices[0], indices[1], indices[2]] = logp
				logpgrid.__setitem__(tuple(indices), logp)
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
			#print indices, p
		pgrid2 = exp(logpgrid)
		#print pgrid.max(), pgrid.min()
		#print pgrid2.max(), pgrid2.min()
		#print pgrid
		pgrid /= pgrid.max()
		pgrid2 = exp(logpgrid-logpgrid.max())
		#pgrid = pgrid2
		if 0: 
			Nbeta = 500
			betamin, betamax = -3, 1
			#nr = 700
			#rmax = 1.0
			#xmin=(1.556303+0.45)
			#xmax=(5.181438-0.75)
			#oldrs = arange(0, nr)/(nr-1.)*rmax + 0.03
			#rs = arange(0, nr)/(nr-1.)*rmax + 0.03
			#linbetagrid = zeros((len(rs),Nbeta))
			#oldrs = self.storage3d
			nr = 50
			xmin=0
			xmax=1.5 # * 2610.9660574412533 
			rs = (arange(nr) + 0.5 )/ (nr) * (xmax - xmin) + xmin
			betagrid = zeros((len(rs), Nbeta))
			print "# models", len(self.parameterset.paramlist)
			for paramlist in self.parameterset.paramlist:
				name, dirname, orig_values, true_values, indices = paramlist
		#	for name, dirname, indices, values in self.paramlist[1:]:
				filename = os.path.join(dirname, "results", "solution_moments3d" +self.parameterset.postfix +".npy")
				if os.path.exists(filename):
					solution_moments3d = numpy.load(filename)
					varvr = solution_moments3d[4]
					varvphi = solution_moments3d[5]
					varvtheta = solution_moments3d[6]
					betas = 1 - (varvphi + varvtheta)/(2*varvr)
					if 0:
						print betas
						print varvr
						print varvphi
						print varvtheta
					betaindices = ((betas-betamin)/(betamax-betamin) * Nbeta).astype(int)
					#print ".", ~isnan(varvtheta)
					good_betaindices = (betaindices >= 0) & (betaindices < Nbeta) & (~isnan(varvtheta))
					betaindices = betaindices[good_betaindices]
					p = pgrid[indices[0], indices[1]]
					#print good_betaindices,betaindices
					betagrid[good_betaindices,betaindices] += p
					#print len(solution_moments3d)
					#f = interpolate.interp1d(oldrs, betagrid[:,i])
				else:
					#logger.error("file '%s' missing (model not finished?)" % filename)
					pass
		
		#print pgrid.max(), pgrid.min()
		#print pgrid2.max(), pgrid2.min()
		#print pgrid2
		from kaplot import *
		import scipy.misc.pilutil
		im = pgrid2
		import scipy.ndimage
		#im = scipy.ndimage.gaussian_filter(im, 1.0)
		if 0:
			im = scipy.misc.pilutil.toimage(im)
			width, height = 500, 500
			import Image
			im = im.resize((width, height), Image.BILINEAR)
			im = scipy.misc.pilutil.fromimage(im)
			#im = pgrid2
		#im = scipy.ndimage.gaussian_filter(im, 0.5)
		#im = scipy.misc.pilutil.imresize(pgrid2, (200,200))
		#box()
		plot_anisotropy = False
		if plot_anisotropy:
			vsplit(box)
		else:
			#box()
			mozaic(3,3,box)
		#indexedimage(transpose(im), colormap="whiterainbow")
		#print log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])
		#print (self.paramformats[0][-1]), (self.paramformats[1][-1])
		#probimage2d(im, 0, 1, log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])) #, colormap="whiteblack")
		print im.min(), im.max()
		select(0, 0)
		probimage2d(im, 0, 1, self.parameterset.ranges_org[0], self.parameterset.ranges_org[1]) #, colormap="whiteblack")
		select(1,0)
		probimage2d(im, 1, 2, self.parameterset.ranges_org[1], self.parameterset.ranges_org[2]) #, colormap="whiteblack")
		select(0,1)
		probimage2d(im, 0, 2, self.parameterset.ranges_org[0], self.parameterset.ranges_org[2]) #, colormap="whiteblack")
		#ylim(-1, 1)
		if plot_anisotropy:
			select(1)
			#probimage2d(betagrid, 0, 1, log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])) #, colormap="whiteblack")
			#betagrid = scipy.ndimage.gaussian_filter(betagrid, [2, 3])
			betagrid = scipy.ndimage.gaussian_filter(betagrid, [2, 2])
			if 0:
				for i in range(betagrid.shape[1]):
					if sum(betagrid[:,i]) > 0:
						betagrid[:,i] /= betagrid[:,i].sum()
			if 0:
				for i in range(betagrid.shape[0]):
					if sum(betagrid[i]) > 0:
						betagrid[i] /= betagrid[i].sum()
			resize = (xmin, betamin), (xmax, betamax)
			indexedimage(transpose(betagrid), colormap="whiterainbow", resize=resize)
			betas = (arange(Nbeta) + 0.5 )/ (Nbeta) * (betamax- betamin) + betamin
			#probimage2d(betagrid, 0, 1, rs, betas) #, colormap="whiteblack")
		draw()
		

class ParameterSetIteratorCondor(object):
	def __init__(self, modelpath, parametersweep, basename, iterations, initial_subdivides, max_newpoints, binary):
		self.modelpath = modelpath
		self.parametersweep = parametersweep
		self.parameterset = self.parametersweep.parameterset
		self.basename = basename
		self.iterations = iterations
		self.initial_subdivides = initial_subdivides
		self.max_newpoints = max_newpoints
		self.binary = binary
		
	def format_name(self, iteration):
		return "%s.%03d" % (self.basename, iteration)
	
	def prepare(self):
		dirname = self.parameterset.dirname
		
		output_filename = os.path.join(dirname, self.format_name(0))
		cmd = "%s -d %d -i -o %s -s %d " % (self.binary, self.parameterset.dimension, output_filename, self.initial_subdivides)
		cmd += " -- "
		for a, b in self.parameterset.ranges_org:
			cmd += str(a) + " "
		for a, b in self.parameterset.ranges_org:
			cmd += str(b) + " "
		print cmd
		if os.system(cmd) != 0:
			return
		
	def iterate(self, iteration):
		pass
		
		#--limit-level-difference=1
		
		
