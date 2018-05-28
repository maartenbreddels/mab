# -*- coding: utf-8 -*-
from numpy import *
import os
import mab.gd.logging as logging
import scipy.interpolate
import mab.utils.numpy
logger = logging.getLogger("gd.schw.derived")

class CommandDerivedProfiles(object):
	def __init__(self, profile_calculators, profile_values):
		self.profile_calculators = profile_calculators
		self.profile_values = profile_values
		
	def run(self, args, opts, scope):
		for profile_calculator, profile_value in zip(self.profile_calculators, self.profile_values):
			values = profile_calculator.calculate()
			profile_value.set_values(values)
			profile_value.save()
	
class ProfileAperture(object):
	def __init__(self, aperture, name="aperture", label="nolabel"):
		aperture.load()
		self.xmin = aperture.aperture_rcenters_kpc[0]
		self.xmax = aperture.aperture_rcenters_kpc[-1]
		self.aperture = aperture
		self.label = label
		self.name = name
		
		self.centers = self.aperture.aperture_rcenters_kpc
		self.u_centers = self.centers
		self.borders = self.aperture.aperture_rborders_kpc
		self.u_borders = self.borders
		
		self.N = len(self.centers)
			
class ProfileLog(object):
	def __init__(self, logrmin, logrmax, N, name="logr", label="nolabel"):
		self.xmin = logrmin
		self.xmax = logrmax
		self.label = label
		
		u = (arange(N)+0.5) / N
		logr = u * (logrmax - logrmin) + logrmin
		self.centers = 10**logr
		self.u_centers = logr
		
		u = (arange(N+1.)) / (N)
		logr = u * (logrmax - logrmin) + logrmin
		self.borders = 10**logr
		self.u_borders = logr
		
		self.name = name
		self.N = N
			
class ProfileLog(object):
	def __init__(self, logrmin, logrmax, N, name="logr", label="nolabel"):
		self.xmin = logrmin
		self.xmax = logrmax
		self.label = label
		
		u = (arange(N)+0.5) / N
		logr = u * (logrmax - logrmin) + logrmin
		self.centers = 10**logr
		self.u_centers = logr
		
		u = (arange(N+1.)) / (N)
		logr = u * (logrmax - logrmin) + logrmin
		self.borders = 10**logr
		self.u_borders = logr
		
		self.name = name
		self.N = N
	
class ProfileLinear(object):
	def __init__(self, rmin, rmax, N, name="r", label="nolabel"):
		self.xmin = rmin
		self.xmax = rmax
		self.label = label
		
		u = (arange(N)+0.5) / N
		self.centers = self.u_centers = self.centers = u * (rmax - rmin) + rmin
		
		u = (arange(N+1.)) / (N)
		self.borders = self.u_borders = self.borders = u * (rmax - rmin) + rmin
		
		self.name = name
		self.N = N
		
class ProfileValues(object):
	def __init__(self, profile, dirname, name, postfix=""):
		self.profile = profile
		self.dirname = dirname
		filename = "_".join(["profile", self.profile.name, name + postfix]) + ".npy"
		self.filename = os.path.join(self.dirname, filename)
		
	def set_values(self, values):
		self.values = values
		
	def load(self):
		self.values = load(self.filename)
	
	def save(self):
		save(self.filename, self.values)
		
		
class ProfileGrid(object):
	def __init__(self, dirname, profile, profile_values_name, ymin, ymax, N, parameterset_iterator, name, postfix="", label="nolabel"):
		self.profile = profile
		self.profile_values_name = profile_values_name
		self.parameterset_iterator = parameterset_iterator
		self.N = N
		self.ymin = ymin
		self.ymax = ymax
		self.levels = [0.682689492137, 0.954499736104, 0.997300203937]
		self.label = label
		
		self.grid = zeros((profile.N, N))
		self.contours = zeros((profile.N, len(self.levels), 2))
		self.contour_mask = zeros((profile.N), dtype=bool)
		self.median = zeros((profile.N))
		basefilename = "_".join(["profile_grid", self.profile.name, name]) + postfix
		self.filename_grid = os.path.join(dirname, basefilename + "_grid.npy")
		self.filename_contours = os.path.join(dirname, basefilename + "_contours.npy")
		self.filename_contour_mask = os.path.join(dirname, basefilename + "_contour_mask.npy")
		self.filename_median = os.path.join(dirname, basefilename + "_median.npy")
		self.resize = [(profile.xmin, ymin), (profile.xmax, ymax)]
		
	def load(self):
		logger.info("loading grid: %s" % self.filename_grid)
		self.grid = load(self.filename_grid)
		logger.info("loading contours: %s" % self.filename_contours)
		self.contours = load(self.filename_contours)
		logger.info("loading contour mask: %s" % self.filename_contour_mask)
		self.contour_mask = load(self.filename_contour_mask)
		logger.info("loading median: %s" % self.filename_median)
		self.median = load(self.filename_median)
		
	def save(self):
		logger.info("saving grid: %s" % self.filename_grid)
		save(self.filename_grid, self.grid)
		logger.info("saving contours: %s" % self.filename_contours)
		save(self.filename_contours, self.contours)
		logger.info("saving contour mask: %s" % self.filename_contour_mask)
		save(self.filename_contour_mask, self.contour_mask)
		logger.info("saving median: %s" % self.filename_median)
		save(self.filename_median, self.median)
		
	def run(self, args, opts, scope):
		self.parameterset_iterator.parameterset.load()
		self.calculate(scope)
		self.save()
		
	def calculate_contours(self):
		Nx, Ny = self.grid.shape
		assert self.profile.N == Nx
		u = (arange(Ny+1.)) / (Ny)
		y_borders = u * (self.ymax - self.ymin) + self.ymin
		for x_index in range(Nx):
			y = concatenate([[0], cumsum(self.grid[x_index,:])])
			if max(y) > 0:
				y /= max(y)
				f = scipy.interpolate.interp1d(y, y_borders)
				ymedian = f(0.5)
				self.median[x_index] = ymedian
				for level_index, level in enumerate(self.levels):
					ylow = f(0.5-level/2)
					yhigh = f(0.5+level/2)
					self.contours[x_index, level_index, 0] = ylow
					self.contours[x_index, level_index, 1] = yhigh
				self.contour_mask[x_index] = True
			else:
				self.contour_mask[x_index] = False
	
	def calculate(self, scope):
		N = self.N
		ymin = self.ymin
		ymax = self.ymax
		scopes = self.parameterset_iterator.scopelist(scope, all=True)
		#for scope in self.parameterset_iterator.iter(scope, all=True):
		grids = mab.utils.numpy.mmapzeros((len(scopes), self.profile.N, N))
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def do(i):
			scope = scopes[i] 
			scope.re_readfiles()
			scope.init()
			profile_values = scope[self.profile_values_name]
			profile_values.load()
			values = profile_values.values 
			#print scope["model_id_path"]
			#print profile_values.values.shape
			#print values
			indices = ((values-ymin)/(ymax-ymin) * N).astype(int)
			good_indices = (indices >= 0) & (indices < N)
			indices = indices[good_indices]
			#print self.grid.shape
			#print good_indices, indices
			#print good_indices.shape
			#print indices.shape
			#self.grid[0,indices]
			#self.grid[good_indices, indices] += scope["probability"] * scope["volume"]
			#print scope["model_id_path"], scope["probability"]
			grids[i, good_indices, indices] += scope["probability"] * scope["volume"]
		do(range(len(scopes)))
		self.grid = sum(grids, axis=0)
			

		self.calculate_contours()
		print self.grid[-1,-1]
			
				#print name 
		
class MassEnclosed(object):
	def __init__(self, density, profile):
		self.density = density
		self.profile = profile
		
	def calculate(self):
		return log10(self.density.enclosed_mass(self.profile.centers))
		
class Logslope(object):
	def __init__(self, density, profile):
		self.density = density
		self.profile = profile
		
	def calculate(self):
		return self.density.logslope(self.profile.centers)
		
		
class Anisotropy(object):
	def __init__(self, profile, storage_3d, solution):
		self.profile = profile
		self.storage_3d = storage_3d
		self.solution = solution
		
	def calculate(self):
		self.storage_3d.load()
		self.solution.load()
		allmoments = self.storage_3d.moments3d
		moments = self.solution.calculate_solution_moments(allmoments)
		mask = moments[0] > 0
		moments[1:,mask] /= moments[0,mask]
		varvr = moments[self.storage_3d.v20index]
		varvphi = moments[self.storage_3d.v02index]/2
		#varvtheta = moments[6]
		betas = varvphi * 0 
		betas[mask] = 1 - (varvphi[mask])/(varvr[mask])
		#print moments.shape
		#print betas.shape
		return betas
		
class VelocityDispersion(object):
	def __init__(self, profile, storage_2d, solution):
		self.profile = profile
		self.storage_2d = storage_2d
		self.solution = solution
		
	def calculate(self):
		self.storage_2d.load()
		self.solution.load()
		allmoments = self.storage_2d.projectedmoments
		moments = self.solution.calculate_solution_moments(allmoments)
		moments[1:] /= moments[0]
		vlosvar = moments[2]
		return vlosvar**0.5
				
class VelocityMoment4(object):
	def __init__(self, profile, storage_2d, solution):
		self.profile = profile
		self.storage_2d = storage_2d
		self.solution = solution
		
	def calculate(self):
		self.storage_2d.load()
		self.solution.load()
		allmoments = self.storage_2d.projectedmoments
		moments = self.solution.calculate_solution_moments(allmoments)
		moments[1:] /= moments[0]
		m4 = moments[4]
		return m4**0.25
				
class Kurtosis(object):
	def __init__(self, profile, storage_2d, solution):
		self.profile = profile
		self.storage_2d = storage_2d
		self.solution = solution
		
	def calculate(self):
		self.storage_2d.load()
		self.solution.load()
		allmoments = self.storage_2d.projectedmoments
		moments = self.solution.calculate_solution_moments(allmoments)
		moments[1:] /= moments[0]
		m4 = moments[4]
		m2 = moments[2]
		return m4/m2**2
		

class MdynoverMstar(object):
	def __init__(self, profile, dm_density, light_profile):
		self.profile = profile
		self.dm_density = dm_density
		self.light_profile = light_profile
		
	def calculate(self):
		r = self.profile.centers
		Menc_dm = self.dm_density.enclosed_mass(r)
		Menc_light = self.light_profile.enclosed_mass(r)
		return log10(Menc_dm/Menc_light)
		
