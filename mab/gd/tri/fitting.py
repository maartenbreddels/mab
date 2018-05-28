# -*- coding: utf-8 -*-
from numpy import *
from kaplot import *
import scipy


#lass Fit


class FitPlummer(object):
	def __init__(self, ellipticity, scale, orblib, aperture):
		self.ellipticity = ellipticity
		self.scale = scale
		self.orblib = orblib
		self.aperture = aperture
		
	def run(self, args, opts, scope):
		self.orblib.load()
		self.aperture.load()
		x, y = self.aperture.xgrid, self.aperture.ygrid
		R = sqrt(x**2+y**2/self.ellipticity)
		density = self.scale**2 / (pi * (self.scale**2+R**2)**2)
		density = density / sum(density)
		
		#box()
		mozaic(3,2,box)
		select(0, 0)
		
		colormap = "blackwhite"
		indexedimage(density, colormap=colormap)
		
		densityflat = density.reshape(-1)
		buildingblocks = array([self.orblib.grid[k].sum(axis=1) for k in range(self.orblib.grid.shape[0])])
		
		weights = self.fit(densityflat, buildingblocks)
		#densityfit = sum([buildingblocks * weight for weight in weights])
		densityfit = dot(weights, buildingblocks)
		densityfit = densityfit/sum(densityfit)
		density = density/sum(density)
		select(0, 1)
		print buildingblocks.shape
		print weights.shape
		print densityfit.shape
		densityfit = densityfit.reshape(density.shape)
		indexedimage(densityfit, colormap=colormap)
		select(1,0)
		diff = densityfit-density
		indexedimage(diff, colormap=colormap)
		print abs(diff).max(), density.max(), densityfit.max()
		select(1,1)
		rdiff = diff/density
		indexedimage(rdiff, colormap=colormap)
		print abs(rdiff).max(), density.max(), densityfit.max()
		
		select(2, 0)
		for i in range(density.shape[0])[::4]:
			graph(density[i,:], color="red")
			graph(densityfit[i,:], color="black")
		select(2, 1)
		for i in range(density.shape[0])[::4]:
			graph(density[:,i], color="red")
			graph(densityfit[:,i], color="black")
		draw()
	
	def fit(self, target, vectors):
		#x, error = scipy.optimize.nnls(vectors.T, target)
		#print x, 
		#return x
		#mask = target > 0
		print target.shape, vectors.shape
		if 0:
			mass_matrix = vectors[:,mask] * 1
			totalmass_matrix = sum(vectors[:,mask], axis=1) * 1
		else:
			mass_matrix = vectors * 1
			totalmass_matrix = sum(vectors, axis=1) * 1
		mass_matrix = array(mass_matrix, copy=True, order='C')
		totalmass_matrix = array(totalmass_matrix, copy=True, order='C')
		target = array(target, copy=True, order='C')
		print totalmass_matrix.shape, mass_matrix.shape
		print totalmass_matrix.dtype, mass_matrix.dtype
		from mab.gd import gdfast_schw
		import pyublas
		#opt_matrix_mass = gdfast_schw.OptimizationMatrix(mass_matrix, totalmass_matrix* 0 + 1)
		opt_matrix_mass = gdfast_schw.OptimizationMatrixN(mass_matrix, target, totalmass_matrix)
		#opt_matrix_mass = gdfast_schw.OptimizationMatrix(pyublas.why_not(mass_matrix), pyublas.why_not(totalmass_matrix))
		weights = ones(vectors.shape[0])
		x = weights
		logL1 = opt_matrix_mass.logp(x)
		
		result = dot(x, vectors)
		logL2 = sum(target* log(result))
		#import pdb
		#pdb.set_trace()
		
		
		def f_(weight):
			weight = 10**weight
			result = dot(weight, vectors)
			chisq = sum(((target-result)**2))
			return chisq
			
		def f(u):
			#x = weight
			global last_fitweights
			x = weight = exp(u)
			#weight = array(weight, copy=True, order='C')
			#result = dot(weight, vectors)
			#logL = sum(target[mask] * log(result[mask]))
			#print logL
			logL = opt_matrix_mass.logp(x)
			grad = weight * 0
			opt_matrix_mass.dlogpdx(x, grad)
			grad = grad * x
			print logL, grad[0:10]
			last_fitweights = u
			return -logL, -grad
			
		
		weights = ones(vectors.shape[0])
		bounds = None
		try: 
			fitweights = scipy.optimize.fmin_l_bfgs_b(f, log(weights), None, bounds=bounds, approx_grad=False, iprint=1,factr=1e-2,maxfun=200000)[0]
		except KeyboardInterrupt:
			fitweights = exp(last_fitweights)
		else:
			fitweights = exp(fitweights)
		
		print fitweights
		return fitweights
		
		
		