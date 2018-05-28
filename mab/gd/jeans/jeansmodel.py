from numpy import *
import numpy
import mab.gd.logging as logging
import mab.parallelize
from scipy import special
import pickle
import mab.utils.numpy
import os
import sys
from kaplot import *

store_cache = True
logger = logging.getLogger("gd.jeans.jeansmodel")

def h3(x):
	return x**3-3*x

def h4(x):
	return x**4-6*x**2+3
def h6(x):
	return x**6-15*x**4+45*x**2-15

fac3 = 1*2*3
fac4 = fac3*4
fac6 = fac3 * 4 * 5 * 6


def makeresize(parameterset, i, j):
	return (parameterset.ranges_org[i][0], parameterset.ranges_org[j][0]), (parameterset.ranges_org[i][1], parameterset.ranges_org[j][1])

class JeansModel(object):
	def __init__(self, modelpath, modelsetname, parameterset, projection_matrix, observation, use_fourth_moment=False, filters=[], postfix=""):
		self.modelpath = modelpath
		self.parameterset = parameterset
		self.projection_matrix = projection_matrix
		self.observation = observation
		self.schwsetname = modelsetname #"biased_10k_jeans2"
		self.postfix = postfix
		
		u1=-3.
		u2=3.
		length=250
		self.u = (arange(length)+0.5)/(length+0.) * (u2-u1) + u1
		self.r = 10**self.u
		self.du = (u2-u1)/length
		
		self.dimension = self.parameterset.dimension
		
		#self.schwsetname = "biased_jeans2"
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname)
		
		self.probability_range = parameterset.ranges_org 
		self.labels = parameterset.paramnames
		self.paramlabels = parameterset.paramlabels
		
		self.use_m4 = use_fourth_moment
		
		self.modelcache = {}
		resultdir = os.path.join(self.dirname, "results")
		if self.use_m4:
			self.modelcachefilename = os.path.join(resultdir, "modelcache_m4_d%d%s_mod.pickle" % (self.parameterset.dimension, self.postfix))
		else:
			self.modelcachefilename = os.path.join(resultdir, "modelcache_d%d%s.pickle" % (self.parameterset.dimension, self.postfix))
		
		
	def save(self, iteration=0):
		resultdir = os.path.join(self.dirname, "results")
		logger.info("ensuring directory exists: %s" % resultdir)
		if not os.path.exists(resultdir):
			os.makedirs(resultdir)
		
		filename = os.path.join(resultdir, "probability_grid_d%d%s%03d.npy" % (self.parameterset.dimension, self.postfix, iteration))
		logger.info("storing probability grid: %s" % filename)
		numpy.save(filename, self.probability_grid)
		if store_cache:
			logger.info("storing modelcache: %s" % self.modelcachefilename)
			fp = open(self.modelcachefilename, "wb")
			pickle.dump(self.modelcache, fp)
			fp.close()
		
		filename = os.path.join(resultdir, "parameter_points_d%d%s%03d.npy" % (self.parameterset.dimension, self.postfix, iteration))
		logger.info("storing parameter_points: %s" % filename)
		numpy.save(filename, self.parameter_points)
		
		#numpy.save(filename, self.modelcache)
		
#self.paramvalues		
		
		return
		
		filename = os.path.join(resultdir, "anisotropy_grid_%03d.npy" % iteration)
		logger.info("storing anisotropy grid: %s" % filename)
		numpy.save(filename, self.anisotropy_grid)
		
		filename = os.path.join(resultdir, "mass_enclosed_grid_%03d.npy" % iteration)
		logger.info("storing mass_enclosed grid: %s" % filename)
		numpy.save(filename, self.mass_enclosed_grid)
	
		filename = os.path.join(resultdir, "logslope_grid_%03d.npy" % iteration)
		logger.info("storing logslope grid: %s" % filename)
		numpy.save(filename, self.logslope_grid)
		
	def load(self, iteration=0):
		resultdir = os.path.join(self.dirname, "results")
		
		#filename = os.path.join(resultdir, "probability_grid_%03d.npy" % iteration)
		filename = os.path.join(resultdir, "probability_grid_d%d%s%03d.npy" % (self.parameterset.dimension, self.postfix, iteration))
		logger.info("loading probability grid: %s" % filename)
		self.probability_grid = numpy.load(filename)
		
		if os.path.exists(self.modelcachefilename):
			logger.info("loading modelcache: %s" % self.modelcachefilename)
			#self.modelcache = numpy.load(filename)
			fp = open(self.modelcachefilename, "rb")
			self.modelcache = pickle.load(fp)
			fp.close()
		
		#filename = os.path.join(resultdir, "parameter_points_%03d.npy" % iteration)
		filename = os.path.join(resultdir, "parameter_points_d%d%s%03d.npy" % (self.parameterset.dimension, self.postfix, iteration))
		logger.info("loading parameter_points: %s" % filename)
		self.parameter_points = numpy.load(filename)
		
		return
		
		filename = os.path.join(resultdir, "anisotropy_grid_%03d.npy" % iteration)
		logger.info("loading anisotropy grid: %s" % filename)
		self.anisotropy_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "mass_enclosed_grid_%03d.npy" % iteration)
		logger.info("loading mass_enclosed grid: %s" % filename)
		self.mass_enclosed_grid = numpy.load(filename)
		
		filename = os.path.join(resultdir, "logslope_grid_%03d.npy" % iteration)
		logger.info("loading logslope grid: %s" % filename)
		self.logslope_grid = numpy.load(filename)		
		
	def __collect(self, iteration=0):
		filename = os.path.join(self.modelpath, "logps.npy")
		logpss = load(filename)
		
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
		for logp, (name, dirname, orig_values, true_values), volume in zip(list(logpss), self.parameterset.paramlist, self.parameterset.volumes):
			dirname = os.path.join(self.modelpath, "schw", self.schwsetname, name)
			logp += 37700 #$7000/2000.*len(self.observation.stars)
			if isnan(logp):
				logp = -1000
			line = " ".join(["%1.30e" % k for k in list(orig_values) + [volume, logp]])
			print >>f, line
		f.close()
		#sys.exit(0)	
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
		errorcode = os.system(cmd)
		if errorcode != 0:
			print "errorcode:", errorcode
			raise Exception("error gridding")
		inputname = os.path.join(self.modelpath, "schw", self.schwsetname, self.parameterset.format_name(iteration)+".grid")
		#lines = file(inputname).readlines()
		#if dimension is 
		#values = array([[float(k) for k in line.split()] for line in lines])
		print "reading in data..."
		data = file(inputname).read()
		#print data
		gridlist = eval(data)
		grid = array(gridlist)
		self.probability_grid = p = grid
		
		#print gridlist
		print grid.shape
		if 0:
			from kaplot import *
			box()
			
			#prob((grid), colormap="whiterainbow")
			resize = makeresize(self.parameterset, 0, 1)
			fill = "blue"
			probimage2d(self.probability_grid, 0, 1, resize=resize, colormap=None, color="blue", drawcontourlines=True, premultiply=True, fill=fill)
			draw()
		
	def _collect(self, iteration=0, cores=4):
		w = 80
		
		self.parameter_points = transpose(array(self.parameterset.points))
		
		#mab.gd.gdfast.ProfilerStart("sparse")
		self.probability_grid = mab.utils.numpy.mmapzeros((w,) * self.parameterset.dimension)
		make_linear = False
		if make_linear:
			self.parameterset.make_linear()
		if self.parameterset.dimension == 2:
			for i in range(w):
				x1, x2 = self.parameterset.ranges_org[0]
				x = x1 + (x2-x1) * i / (w-1.) 
				for j in range(w):
					y1, y2 = self.parameterset.ranges_org[1]
					y = y1 + (y2-y1) * j / (w-1.)
					#print self.parameterset((x, y))
					self.probability_grid[i,j] = self.parameterset((x, y)) 
		if self.parameterset.dimension == 3:
			x1, x2 = 0., 1.#self.parameterset.ranges_org[0]
			y1, y2 = 0., 1.#self.parameterset.ranges_org[1]
			z1, z2 = 0., 1.#self.parameterset.ranges_org[2]
			logger.info("interpolating grid")
			@mab.parallelize.parallelize(cores=cores, info=info)
			def do(i):
				x = x1 + (x2-x1) * i / (w-1.) 
				for j in range(w):
					y = y1 + (y2-y1) * j / (w-1.)
					for k in range(w):
						z = z1 + (z2-z1) * k / (w-1.)
						self.probability_grid[i,j,k] = self.parameterset.eval_normalized((x, y, z))
			xindices = range(w)
			do(xindices)
		if self.parameterset.dimension == 4:
			x1, x2 = self.parameterset.ranges_org[0]
			y1, y2 = self.parameterset.ranges_org[1]
			z1, z2 = self.parameterset.ranges_org[2]
			v1, v2 = self.parameterset.ranges_org[3]
			logger.info("interpolating grid")
			#info = False
			@mab.parallelize.parallelize(cores=cores, info=info)
			def do(ij):
					i = ij%w
					j = ij/w
					x = x1 + (x2-x1) * i / (w-1.) 
					y = y1 + (y2-y1) * j / (w-1.)
					for k in range(w):
						z = z1 + (z2-z1) * k / (w-1.)
						for l in range(w):
							v = v1 + (v2-v1) * l / (w-1.)
							#print (x, y, z, v), self.parameterset((x, y, z, v)) 
							self.probability_grid[i,j,k,l] = self.parameterset((x, y, z, v))
			xindices = range(w**2)
			do(xindices)
		if self.parameterset.dimension == 5:
			x1, x2 = self.parameterset.ranges_org[0]
			y1, y2 = self.parameterset.ranges_org[1]
			z1, z2 = self.parameterset.ranges_org[2]
			a1, a2 = self.parameterset.ranges_org[3]
			b1, b2 = self.parameterset.ranges_org[4]
			logger.info("interpolating grid")
			#info = False
			@mab.parallelize.parallelize(cores=cores, info=info)
			def do(ijk):
					i = ijk%w
					j = (ijk/w)%w
					k = (ijk/w)/w
					x = x1 + (x2-x1) * i / (w-1.) 
					y = y1 + (y2-y1) * j / (w-1.)
					z = z1 + (z2-z1) * k / (w-1.)
					for l in range(w):
						a = a1 + (a2-a1) * l / (w-1.)
						for m in range(w):
							b = b1 + (b2-b1) * m / (w-1.)
							#print (x, y, z, v), self.parameterset((x, y, z, v)) 
							self.probability_grid[i,j,k,l,m] = self.parameterset((x, y, z, a, b))
			xindices = range(w**3)
			do(xindices)
		#mab.gd.gdfast.ProfilerStop()
		if not make_linear:
			#print self.probability_grid 
			self.probability_grid -= self.probability_grid.max()
			#print self.probability_grid 
			self.probability_grid = exp(self.probability_grid)
			mask = numpy.isnan(self.probability_grid)
			self.probability_grid[mask] = 0 
			#print self.probability_grid 
		
	def init(self, iteration=0):
		self.parameterset.init()
		self.parameterset.load(iteration=iteration)
		self.projection_matrix.load()
		self.observation.load()
		self.aperture = self.projection_matrix.gridR
		
		self.M = self.projection_matrix.matrices
		
		#print self.Mr.shape
		#print self.projection_matrix.gridR.radial_surface_densities.shape
		
		light_model = self.projection_matrix.light_model
		self.R = R = light_model.arcsec_to_kpc(self.projection_matrix.gridR.aperture_rcenters)
		dR = light_model.arcsec_to_kpc(self.projection_matrix.gridR.aperture_rcenters[1] - self.projection_matrix.gridR.aperture_rcenters[0])
		
		
		self.massr = light_model.densityr(self.r)/1e6 * self.r**3 * 4 * pi * log(10) * self.du
		#massR = tensordot(self.Ms, self.massr, axes=[(1,), (0,)])
		self.massR = (2*pi*R*dR)*light_model.densityR(R,M=1.)
		print "sum", sum(self.massr)
		
		for i in range(self.M.shape[0]):
			self.M[i] = (self.M[i].T /self.massR).T
			self.M[i] = (self.M[i] *self.massr)
		
		self.Ms = self.M[0]
		self.Mr = self.M[0] - self.M[1]
		self.Mt = self.M[1]
		#self.Mr = (self.Mr.T /self.massR).T
		#self.Mt = (self.Mt.T /self.massR).T
		#self.Mr = (self.Mr *self.massr)
		#self.Mt = (self.Mt *self.massr)
		
		stars = self.observation.stars
		for filter in self.filters:
			stars = filter(stars)
		logger.info("using %d stars/observations" % len(stars))
		#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		#print stars
		stars_inrange = []
		for i in range(len(stars)):
			star = stars[i]
			#print star
			if self.aperture.inrange(stars.xi[i], stars.eta[i]):
				stars_inrange.append(star)
		use_numpy = not isinstance(stars, mab.cvsfile.CsvObject)
		if use_numpy:
			stars_inrange = array(stars_inrange, dtype=stars.dtype).view(recarray)
		else:
			stars_inrange = mab.cvsfile.CsvObject(stars_inrange)
		#import pdb;pdb.set_trace()
		logger.info("stars in aperture range     : %d" % len(stars_inrange))
		logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		stars = stars_inrange
		aperture_indices = []
		for i in range(len(stars)):
			aperture_indices.append(self.aperture.aperture.findindex(stars.xi[i], stars.eta[i]))
		self.aperture_indices = array(aperture_indices)
		#print self.aperture_indices
			
		real_data = True #False
		if not real_data:
			sigma_v = 2.01
			numpy.random.seed(8)
			for i in range(len(stars)):
				stars.vlos[i] = stars.vlos_true[i] + numpy.random.normal(0, sigma_v)
				#star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d.Nv);
				#outlier = True
				#for losvd in losvds:
				#	if losvd[star.v_index, star.aperture_index] != 0:
				#		outlier = False
				#		break
				#star.is_outlier = outlier
				
				#print star.aperture_index
			#indices = array([star.aperture_index for star in stars])
		self.vlos = stars.vlos - self.vmean #array([star.vlos for star in stars])
		self.vlos_sigma = stars.e_vlos #array([star.e_vlos for star in stars])
		print self.aperture_indices.min(), self.aperture_indices.max()
		#print self.storage_2d.aperture
		#moment0 = self.storage_2d.aperture.radial_surface_densities
		#moment0 = self.aperture.radial_surface_densities
		#print vlos_sigma
		#print indices[argsort(indices)]
		#print stars.rc[argsort(indices)]
		
		if os.path.exists(self.modelcachefilename):
			logger.info("loading modelcache: %s" % self.modelcachefilename)
			fp = open(self.modelcachefilename, "rb")
			self.modelcache = pickle.load(fp)
			fp.close()
			#self.modelcache = numpy.load(filename)
			logger.info("model cache contains %d models" % len(self.modelcache))
		
				
		
	def test(self, init, scope):
		self.scope = scope
		mozaic(2,1,box)
		
		scope.flush()
		scope.init()
		scope["dm_density_twoslope.rs"] = 0.5
		scope["dm_density_twoslope.alpha"] = 1
		#scope["jeans.beta"] = -0.5 
		jeans = scope["jeans"]
		self.drawsigmar(jeans, color="black")
		
		scope.flush()
		scope.init()
		scope["dm_density_twoslope.rs"] = 10**-0.85
		scope["dm_density_twoslope.alpha"] = -0.56
		#scope["jeans.beta"] = -0.5 
		jeans = scope["jeans"]
		self.drawsigmar(jeans, color="green", linestyle="dot")
		
		scope.flush()
		scope.init()
		scope["dm_density_twoslope.rs"] = 10**-0.57
		scope["dm_density_twoslope.alpha"] = 0.37
		#scope["jeans.beta"] = -0.5 
		jeans = scope["jeans"]
		self.drawsigmar(jeans, color="red", linestyle="dash")
		
		scope.flush()
		scope.init()
		scope["dm_density_twoslope.rs"] = 10**-0.392
		scope["dm_density_twoslope.alpha"] = 0.860
		#scope["jeans.beta"] = -0.5 
		jeans = scope["jeans"]
		self.drawsigmar(jeans, color="blue", linestyle="dot")
		
		select(0, 0)
		ylim(0, 15)
		select(1, 0)
		ylim(0, 15)
		draw()
			
	def drawsigmar(self, jeans, **kwargs):
		m2 = array([jeans.m2(k) for k in self.r])
		sigmar = sqrt(m2)
		beta = jeans.beta
		m2_los = tensordot(self.M[0]-beta*self.M[1], m2, axes=[(1,), (0,)])
		sigma_los = m2_los**0.5
		select(0, 0)
		graph(self.u, sigmar, **kwargs)
		select(1, 0)
		graph(self.R, sigma_los, **kwargs)

	def run(self, ini, scope, cores=1, info=True):
		self.scope = scope
		#if 1:
		logm = 7.91
		if 0:
			scope.flush()
			scope.init()
			logprior = 0
			paramnames = ['dm_density_twoslope.rs', 'dm_density_twoslope.alpha']#, 'jeans.beta']
			values = [10**-0.5, 0.54, -0.5]
			values = [3.16227766017, -0.25]

			for name, value in zip(paramnames, values):
				scope[name] = value
			#logprior = log(paramvalue.p_prior)
			#scope["jeans.beta"] = -0.5 
			jeans = scope["jeans"]
			print self.run_one(0, test=True)
			sys.exit(0)
			#print "beta = ", jeans.beta 
			beta = jeans.beta
			logr = arange(-2, 2, 0.1)
			r = 10**logr
			m2 = numpy.array([jeans.m2(k) for k in r])
			m4 = numpy.array([jeans.m4(k) for k in r])
			gamma = m4/m2**2
			graph(logr, gamma)
			draw()
			sys.exit(0)
			
		if 0:
			if 1:
				i = 0
				logm += 0.01
				logrs = -0.4796
				self.parameterset.paramvalues[i].values_org = list(self.parameterset.paramvalues[i].values_org)
				self.parameterset.paramvalues[i].values_trans = list(self.parameterset.paramvalues[i].values_trans)
				self.parameterset.paramvalues[i].values_org[0] = logm
				self.parameterset.paramvalues[i].values_trans[0] = 10**logm
				self.parameterset.paramvalues[i].values_org[1] = logrs
				self.parameterset.paramvalues[i].values_trans[1] = 10**logrs
				print self.run_one(0)
				logrs = -0.4795
				self.parameterset.paramvalues[i].values_org = list(self.parameterset.paramvalues[i].values_org)
				self.parameterset.paramvalues[i].values_trans = list(self.parameterset.paramvalues[i].values_trans)
				self.parameterset.paramvalues[i].values_org[0] = logm
				self.parameterset.paramvalues[i].values_trans[0] = 10**logm
				self.parameterset.paramvalues[i].values_org[1] = logrs
				self.parameterset.paramvalues[i].values_trans[1] = 10**logrs
				print self.run_one(0)
				sys.exit(0)
			while 1:
				logm += 0.01
				print "run, logm=", logm
				if 0:
					for i, p in enumerate(self.parameterset.paramvalues):
						v = p.values_org
						if (v[0] >= 7.97) and (v[0] < 7.98) and (v[1] > -0.52) and (v[1] < -0.25):
							print p.modelname, v
							print self.modelcache[p.modelname]
				
				#logrss = arange(-0.5, -0.45, 0.01/10)
				#logrss = arange(-0.55, -0.45, 0.01/100)
				logrss = arange(-0.4798, -0.4794, 0.01/1000)
				#logrss = arange(-0.5, 0.0, 0.01/4)
				logps = []
				for i, logrs in enumerate(logrss):
					#logm = 8.0 #7.9749999999999996
					self.parameterset.paramvalues[i].values_org = list(self.parameterset.paramvalues[i].values_org)
					self.parameterset.paramvalues[i].values_trans = list(self.parameterset.paramvalues[i].values_trans)
					self.parameterset.paramvalues[i].values_org[0] = logm
					self.parameterset.paramvalues[i].values_trans[0] = 10**logm
					self.parameterset.paramvalues[i].values_org[1] = logrs
					self.parameterset.paramvalues[i].values_trans[1] = 10**logrs
				@mab.parallelize.parallelize(cores=cores, info=info)
				def wrap(index):
					return self.run_one(index)
				logps = wrap(range(len(logrss)))
				print logps
				#logp = self.run_one(0)
				#print logrs, logp
				#	logps.append(logps)
				logps = array(logps)
				logps -= logps.max()
				ps = exp(logps)
				clear()
				box()
				graph(logrss, ps)
				draw()
					
				
			sys.exit(0)
		def do():
			all_indices = range(len(self.parameterset.paramvalues))
			# skip cached models
			known_indices = [i for i in all_indices if self.parameterset.paramvalues[i].modelname in self.modelcache]
			unknown_indices = [i for i in all_indices if self.parameterset.paramvalues[i].modelname not in self.modelcache]
			if 0:
				self.run_one(0)
			else:
				@mab.parallelize.parallelize(cores=cores, info=info)
				def wrap(index):
					return self.run_one(index)
				if cores == 1:
					for i in unknown_indices:
						print i
						print self.run_one(i)
				else:
					logps_solved = wrap(unknown_indices)
					logps_solved = array(logps_solved)
					logps = zeros(len(all_indices))
					logps[unknown_indices] = logps_solved
					logps[known_indices] = [self.modelcache[self.parameterset.paramvalues[i].modelname] for i in known_indices]
					for i in unknown_indices:
						self.modelcache[self.parameterset.paramvalues[i].modelname] = logps[i]
					filename = os.path.join(self.modelpath, "logps.npy")
					save(filename, logps)
					resultdir = os.path.join(self.dirname, "results")
					if store_cache:
						logger.info("storing modelcache: %s (%d models)" % (self.modelcachefilename, len(self.modelcache)))
						fp = open(self.modelcachefilename, "wb")
						pickle.dump(self.modelcache, fp)
						fp.close()
					
					print "logps", logps
				self.parameterset.feed(logps)
			
		print self.parameterset.paramnames
		do()
		for i in range(self.parameterset.max_iterations):
			#print "iterate"
			logger.info("iteration %d (out of %i)" % (i+1,self.parameterset.max_iterations)) 
			self.parameterset.iterate(i,cores=cores)
			do()
			#self.parameterset.refresh()
		
		
		self._collect(iteration=0, cores=cores)
		
	def run_one(self, index, test=False):
		scope = self.scope
		if not test:
			paramnames = self.parameterset.paramnames
			paramvalue = self.parameterset.paramvalues[index]
			#paramprior = self.parameterset.parampriors
			values = paramvalue.values_trans
			scope.flush()
			scope.init()
			logprior = 0
			for name, value in zip(paramnames, values):
				scope[name] = value
			logprior = log(paramvalue.p_prior)
			#scope["jeans.beta"] = -0.5 
		jeans = scope["jeans"]
		#print "beta = ", jeans.beta 
		beta = jeans.beta(self.r)
		
		#jeans = jeans.fast()
		#sigmar = jeans.sigmar(self.r)
		sigmar = array([jeans.m2(k)**0.5 for k in self.r])
		#g = 1 - 2 * self.beta * R**2/r**2 + self.beta*(1+self.beta)/2*R**4/r**4
		
		use_m4 = self.use_m4
		if use_m4:
			m4 = array([jeans.m4(k) for k in self.r])
			kappa4 = m4##/ sigmar**4
		

		
		varr = sigmar**2# * self.massr
		vart =(1-beta)*varr
		
		
		varr_los = tensordot(self.Mr, varr, axes=[(1,), (0,)])
		vart_los = tensordot(self.Mt, vart, axes=[(1,), (0,)])
		#light_model = jeans.light_model
		#varR_los = (varr_los + vart_los)
		m2 = sigmar**2
		m2_los = tensordot(self.M[0]-beta*self.M[1], m2, axes=[(1,), (0,)])
		if use_m4:
			#g = 1 - 2 * self.beta * R**2/r**2 + self.beta*(1+self.beta)/2*R**4/r**4
			m4_los = tensordot(self.M[0]-2*beta*self.M[1]+beta*(1+beta)/2*self.M[3], m4, axes=[(1,), (0,)])
		# m2*(1-beta*
		varR_los = m2_los
		
		m2_los = varR_los[self.aperture_indices] + self.vlos_sigma**2
		if use_m4:
			m4_los = m4_los[self.aperture_indices] + 3*self.vlos_sigma**4
			kappa4 = m4_los
		sigma_los = sqrt(m2_los)
		#print sigmasqs
		#print sigmasqs.shape
		logps = sum(- self.vlos**2/(2*m2_los) - log(sqrt(m2_los*2*pi)))
		if use_m4:
			u = (kappa4/sigma_los**4-3)
			g = m4_los/sigma_los**4
			
			if 0:
				a1, a2, a3, b1, b2,b3 = [-0.10971413, -0.03265063,  0.02697869,  -0.52812814, -0.05336269, -0.0592689] # unnormalized
				gamma = g
				mask = gamma < 3
				#print "%e" % gamma
				sigma = sigma_los
				sigma_mod = sigma * 1.0
				sigma_mod[mask] = sigma[mask] * (1 + (3-gamma[mask])**2*a1 + (3-gamma[mask])**3*a2 + (3-gamma[mask])**4*a3)
				kappa2 = sigma**2
				kappa4 = gamma * kappa2**2
				kappa4_mod = kappa4 * 1.
				kappa4_mod[mask] = kappa4[mask] * (1 + (3-gamma[mask])**2*b1 + (3-gamma[mask])**3*b2 + (3-gamma[mask])**4*b3)
				gamma_mod = gamma * 1.
				gamma_mod[mask] = kappa4_mod[mask]/kappa2[mask]**2
				u = gamma_mod * 1.-3
				sigma_los = sigma_mod
			
			#kappa4 = kappa4_mod
			#print list(g)
			#print 
			#print g.min(), g.max(), mean(g)
			#u = minimum(u, 2)
			#u = maximum(u, -1)
			#u = maximum(u, -0.25)
			#u = minimum(u, 0.25)
			kappa4 = (u+3)*sigma_los**4
			s4 = (kappa4 - 3 * sigma_los**4)/sigma_los**6
			b1 = 0 #self.sigma*(1./fac3 * self.s3 * h3(x/self.sigma))
			b2 = sigma_los**2 * (1./fac4 * s4 * h4(self.vlos/sigma_los))
			#print list(g)
			#print b2.min(), b2.max(), mean(b2)
			extra = b2
			#mask = (1+b2)>0
			mask = u > 0.
			mask = b2 > 0.
			z = b2
			z[mask] = log(1+b2[mask])
			
			if 0:
				u = (kappa4/sigma_los**4-3)
				#print b2[mask].max(), u.min(), u.max()
				assert numpy.all(1+b2[mask]) > 0, zip(b2[mask], u[mask], s4[mask])
				oldb2 = b2 * 1 
				b2[mask] = log(1+b2[mask])
				#print b2.min(), b2.max(), mean(b2)
				nans = numpy.isnan(b2)
				if isnan(sum(b2)) and not test:
					for name, value in zip(paramnames, values):
						print name, value
					
				#assert not isnan(sum(b2)), zip(b2[nans], oldb2[nans], u[nans], s4[nans], (self.vlos/sigma_los)[nans], self.vlos[nans], sigma_los[nans])
				print logps, sum(b2)
			else:
				assert not isnan(sum(z)), zip(b2[nans], oldb2[nans], u[nans], s4[nans], (self.vlos/sigma_los)[nans], self.vlos[nans], sigma_los[nans])
			r = logps + sum(z)
			import pdb;
			#pdb.set_trace()
			return r + logprior
		else:
			return logps + logprior

				
		
		if 0:
			from kaplot import *
			box()
			graph(sqrt(varR_los))
			sigmar = jeans.sigma_los(self.R)
			graph(sigmar, color="green")
			draw()
		if 0:
			R = light_model.arcsec_to_kpc(self.projection_matrix.gridR.aperture_rcenters)
			dR = light_model.arcsec_to_kpc(self.projection_matrix.gridR.aperture_rcenters[1] - self.projection_matrix.gridR.aperture_rcenters[0])
			varR_los = (varr_los + vart_los)#/self.massR # /(2*pi*R*dR)*(jeans.light_model.densityR(R)/1e6)#massR*(2*pi*R*dR)# * 10
			mozaic(2,2,box)
			area = 1 #R / jeans.light_model.densityR(R,M=1.) # * self.projection_matrix.gridR.radial_surface_densities# / 100
			sigmar = jeans.sigma_los(R)
			graph(sigmar, color="green")
			graph(sqrt(varR_los*area))
			print sigmar.shape, varR_los.shape
			graph(sigmar**2/(varR_los*area), color="red")
			select(0,1)
			graph(self.u, jeans.light_model.densityr(self.r))
			select(1,1)
			graph(self.projection_matrix.gridR.radial_surface_densities)
			select(1,0)
			graph(self.massR/(2*pi*R*dR), color="red")
			
			graph(jeans.light_model.densityR(R)/1e6, color="green", linestyle="dot")
			draw()
		#print u
		
	