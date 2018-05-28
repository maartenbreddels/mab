# -*- coding: utf-8 -*-
import cvxopt.solvers
import os
import numpy
from numpy import *
import mab.gd.logging as logging
import mab.utils.numpy
import sys
import scipy.optimize
from mab.gd import gdfast_schw
import numpy.linalg

logger = logging.getLogger("gd.schw.solution")


class OrbitWeightCalculatorSeperable(object):
	def __init__(self, modelpath, dfgrid, galaxy, name="orbitweights", plot=True):
		self.modelpath = modelpath
		self.dfgrid = dfgrid
		self.plot = plot
		self.galaxy = galaxy
		self.name = name
		
	def calculate(self, cores=1, info=True):
		filename = os.path.join(self.modelpath, "df", "fE.npy")
		logger.info("loading f1(E) as: %s" % filename)
		dff = load(filename)
		
		filename = os.path.join(self.modelpath, "df", "fE-E.npy")
		logger.info("loading corresponding energy as: %s" % filename)
		dfEs = load(filename)
		#orbitweights = mab.utils.numpy.mmapzeros((self.dfgrid.n_I1, self.dfgrid.n_I2))
		orbitweights_inner_products = mab.utils.numpy.mmapzeros((self.dfgrid.dof))
		df_inner_products = mab.utils.numpy.mmapzeros((self.dfgrid.dof))
		@mab.parallelize.parallelize(cores=cores, info=info)
		def orbitweight(I1, I2):
			dither1 = self.dfgrid.dither
			dither2 = self.dfgrid.dither
			for i in range(dither1):
				E = self.dfgrid.subgrid.Es[I1*dither1 + i]
				dE = (self.dfgrid.subgrid.E_borders[I1*dither1 + i + 1] - self.dfgrid.subgrid.E_borders[I1*dither1 + i])
				Lmax = self.dfgrid.profile_model.Lmax_at_E(E)
				#for j in range(len(dfgrid.subgrid.ls)):
				for j in range(dither2):
					L = (self.dfgrid.subgrid.ls[I2*dither2 + j])*Lmax
					dL = (self.dfgrid.subgrid.l_borders[I2*dither2 + j + 1] - self.dfgrid.subgrid.l_borders[I2*dither2 + j])*Lmax
					
					#if opts.true:
					#	g = galaxy.gEL(E, L)
					#	mass = g * galaxy.fEL(E, L) * dE * dL
					#if opts.truenumerical:
					if 1:
						Ei = argmin(abs(dfEs-E))
						fEL = dff[Ei] * self.galaxy.fL(L) #* L ** (-2*self.anisotropy_beta)
						g = self.dfgrid.profile_model.gEL(E, L)
						mass = g * fEL * dE * dL
						#countG = g * dE *  dL
					u = float(i + 0.5) / dither1
					v = float(j + 0.5) / dither2
					#print u, v
					for k in range(self.dfgrid.dof_per_cell):
						orbitweights_inner_products[self.dfgrid.dof_index(I1, I2, k)] += mass * self.dfgrid.basis(k, u, v)
						df_inner_products[self.dfgrid.dof_index(I1, I2, k)] += fEL * dE * dL * self.dfgrid.basis(k, u, v)
				
		nE = self.dfgrid.n_I1
		nL = self.dfgrid.n_I2
		I1s = [i for i in range(nE) for j in range(nL)] #[:1]
		I2s = [j for i in range(nE) for j in range(nL)] #[:1]
		#print len(I1s), len(dfgrid.subgrid.E_borders), len(dfgrid.subgrid.Es)
		orbitweight(I1s, I2s)
		#print orbitweights_inner_products
		self.orbitweights = self.dfgrid.solve_coordinates(orbitweights_inner_products)
		self.df = self.dfgrid.solve_coordinates(df_inner_products)
		#print self.orbitweights
		if 0:
			for I1 in I1s:
				for I2 in I2s:
					orbitweight(I1, I2)
					
	def load(self):
		filename_base = os.path.join(self.modelpath, "df", self.name) 
		filename_npy = filename_base + ".npy"
		self.orbitweights = numpy.load(filename_npy)
		filename = os.path.join(self.modelpath, "df", "df.npy")
		self.df = numpy.load(filename)
		
	def save(self):
		filename = os.path.join(self.modelpath, "df", "df.npy")
		logger.info("storing distribution function in file: %s" % filename)
		numpy.save(filename, self.df)
		filename_base = os.path.join(self.modelpath, "df", self.name) 
		filename_png = filename_base + ".png" 
		filename_npy = filename_base + ".npy"
		logger.info("storing orbitweights in file: %s" % filename_npy) 
		numpy.save(filename_npy, self.orbitweights)
		import kaplot
		kaplot.box()
		n = 200
		x = (arange(n)+0.5)/(n)
		us = self.dfgrid.umin + x * (self.dfgrid.umax - self.dfgrid.umin)
		ls = 1 * x
		#logrs = self.dfgrid.subgrid.logrs
		#logrs = arange(len(logrs)) / (len(logrs)-1.)
		#ls = self.dfgrid.subgrid.ls
		#print ls, logrs
		#print self.orbitweights
		df = array([[self.dfgrid(self.orbitweights, u, l) for u in us] for l in ls])
		print ls
		#print df
		#print df.shape
		kaplot.indexedimage((df), colormap="whiterainbow")
		print self.orbitweights.shape
		#kaplot.indexedimage((self.orbitweights.reshape(10, 30)))
		kaplot.hardcopy(filename_png)


		
class Solution(object):
	
	def run(self, args, opts, scope):
		self.load()
		self.solve()
	
	def store_moments(self, moments2d, moments3d, x):
		moments2d = tensordot(x, moments2d, axes=[(0,), (0,)])
		moments3d = tensordot(x, moments3d, axes=[(0,), (0,)])
		moments2d[1:] /= moments2d[0]
		moments3d[1:] /= moments3d[0]
		
		#filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		#numpy.save(filename, orbitweights)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +".npy")
		numpy.save(filename, moments2d)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		numpy.save(filename, moments3d)
		 
		
		
class Likelihood2(object):
	def __init__(self, modelpath, light_model, aperture, profile_model, schwsetname, schwmodelname, storage_2d, storage_3d, fitdensity2d, fitdensity3d, observation, dfgrid, max_iterations, regularization, regularization_delta=1., use_jeans=False, jeans_fraction=0.05, postfix=""):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.aperture = aperture
		#self.storage_2d_binned = storage_2d_binned
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.observation = observation
		self.dfgrid = dfgrid
		self.max_iterations = max_iterations
		self.regularization = regularization
		self.regularization_delta = regularization_delta
		self.use_jeans = use_jeans
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.likelihood")
		self.postfix = postfix
		
	def run(self, args, opts, scope):
		self.load()
		self.solve()
		
	def load(self):
		self.observation.load()
		self.storage_2d.init()
		self.storage_2d.load()
		#self.storage_2d_binned.init()
		#self.storage_2d_binned.load()
		self.storage_3d.init()
		self.storage_3d.load()
		#self.storage_2d.aperture.load()
		self.aperture.load()
	
	def solve(self):
		return self.findsolution()
	
	def findsolution(self, solution=None):
		stars = self.observation.stars
		
		stars = self.observation.stars
		logger.info("using %d stars/observations" % len(stars))
		#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		#print stars
		stars_inrange = []
		for i in range(len(stars)):
			star = stars[i]
			#print star
			if self.aperture.inrange(stars.xi[i], stars.eta[i]):
				stars_inrange.append(star)
		stars_inrange = array(stars_inrange, dtype=stars.dtype).view(recarray)
		#import pdb;pdb.set_trace()
		logger.info("stars in aperture range     : %d" % len(stars_inrange))
		logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		stars = stars_inrange
		aperture_indices = []
		for i in range(len(stars)):
			aperture_indices.append(self.storage_2d.aperture.aperture.findindex(stars.xi[i], stars.eta[i]))
		self.aperture_indices = array(aperture_indices)
		#print self.aperture_indices
			
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
		self.vlos = stars.vlos#array([star.vlos for star in stars])
		self.vlos_sigma = stars.e_vlos #array([star.e_vlos for star in stars])
		print self.aperture_indices.min(), self.aperture_indices.max()
		
		if 0:
			
			self.logger.info("using %d stars/observations" % len(stars))
			#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
			mask_inrange = [self.storage_2d.aperture.inrange(stars.xi[i], stars.eta[i]) for i in range(len(stars))]
			stars_inrange = stars[mask_inrange]
			self.logger.info("stars in aperture range     : %d" % len(stars_inrange))
			self.logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
			stars = stars_inrange
			for star in stars:
				star.aperture_index = self.storage_2d.aperture.findindex(star.xi, star.eta)
				
			sigma_v = 2.01
			numpy.random.seed(8)
			for star in stars:
				star.vlos = star.vlos_true + numpy.random.normal(0, sigma_v)
				#star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d.Nv);
				#outlier = True
				#for losvd in losvds:
				#	if losvd[star.v_index, star.aperture_index] != 0:
				#		outlier = False
				#		break
				#star.is_outlier = outlier
				
				#print star.aperture_index
		indices = self.aperture_indices #array([star.aperture_index for star in stars])
		vlos = self.vlos #array([star.vlos for star in stars])
		vlos_sigma = self.vlos_sigma #array([star.e_vlos for star in stars])
		print indices.min(), indices.max()
		print self.storage_2d.aperture
		#moment0 = self.storage_2d.aperture.radial_surface_densities
		moment0 = self.aperture.radial_surface_densities
		print vlos_sigma
		print indices[argsort(indices)]
		print stars.rc[argsort(indices)]
		#dsa
		
		# expectation value if gaussian: <-x^2/(2*sigma^2) - log(sqrt(sigma^2*2*pi))>
		sigma_approx = 10.
		scale = (-0.5 + log(1/(sqrt(2*pi)*sigma_approx))) * len(stars)
		exp_value = scale
		scale = 1.
		print "scale =", scale
		Greg, hreg = self.regularization.assemble_system()
		delta = self.regularization.regularization_delta
		Ptreg = Greg/delta
		Preg = tensordot(Ptreg, Ptreg, axes=([1],[1]))
		print Greg.shape, Preg.shape, hreg.shape
		qreg = -tensordot(array(hreg)/delta, Ptreg, axes=([0], [1]))
		#qreg = array(hreg)/delta
		
		Greg/delta
		#moment4 = True
		
		def logp(x, dlogpdx=False, seperate=False):
			moments = tensordot(x, self.storage_2d.projectedmoments, axes=[(0,), (0,)])
			moments[1:,:] /= moment0
			#print "d", sum((moment0-moments[0])**2)
			#print moments[2].shape, indices.shape
			sigmasqs = moments[2,indices] + vlos_sigma**2
			sigmas = sqrt(sigmasqs)
			#print sigmasqs
			#print sigmasqs.shape
			logps = - vlos**2/(2*sigmasqs) - log(sqrt(sigmasqs*2*pi))
			wk = 1./abs(scale)
			w1 = 1e3/abs(scale)
			wlight = 1e2/abs(scale)
			wall = 1e0
			wreg = -1.0*10#/105#/5
			almost_one = 1.000001
			devfromone = w1 * (sum(x) - almost_one)**2
			#double logLnorm = -pow((1.000001-x.sum())/error_x, 2);
			
			light = wlight * sum(((moment0-moments[0])/moment0)**2)
			#print x.shape, Preg.shape, qreg.shape
			reg = (0.5 * dot(x, dot(x, Preg)) + sum(qreg*x)) * wreg
			#light = 1e5 * sum(((moment0-moments[0]))**2)
			#print sum(logps),devfromone,light 
			#return (sum(logps)/abs(scale)  - devfromone - light) * 1e3
			y = (sum(logps)*wk  - devfromone - light + reg) * wall
			#y = reg
			if seperate:
				return sum(logps), -devfromone, -light, reg
			#y = - light
			
			if dlogpdx:
				dmu2dx = (self.storage_2d.projectedmoments[:,2,:] / moment0)[:,indices]
				#print dmu2dx.shape
				dsigmadx = 1/sigmas * dmu2dx/2
				dydx1 = wk * sum((vlos**2/(sigmas**3.) - 1/sigmas) * (dsigmadx), axis=1)
				#print sum(x)-1
				dydx2 = -w1 * 2*(sum(x)-almost_one) * (x*0+1);
				#print moment0.shape
				#print self.storage_2d.projectedmoments[:,0,indices].shape
				dydx3 = wlight * sum(2*(moment0-moments[0])/moment0**2 * self.storage_2d.projectedmoments[:,0,:], axis=1)
				dydx_reg = (dot(x, Preg) + qreg) * wreg
				dydx = (dydx1+ dydx2 + dydx3 + dydx_reg) * wall
				#dydx = dydx_reg
				return y, dydx
			else:
				return y
		
		
			 
		
		x = ones(160)
		if 1:
			x *= 1.0
			dx = 1e-4
			for index in [1,30, 50, 100]:
				dxs = x * 0
				dxs[index] = dx
				#print "y", logp(x,False)
				#print "y(+dx)", logp(x+dx,False)
				print (logp(x+dxs,False) - logp(x,False))/dx,
				y, dydx = logp(x,True)
				print dydx[index]
			print dydx.shape
			print x.shape
			#sys.exit(0)
		
		x /= sum(x)
		print sum(x)
		print "logp", -logp(x)
		#return x
		if 1:
			filename = os.path.join(self.modelpath, "orbitweights.npy")
			orbitweights = load(filename)
			c = orbitweights.flatten() 
			c /= sum(c)
			print sum(c)
			#x = c
			#x = c
			print "logp", -logp(c)
		filename = os.path.join(self.dirname, "results", "orbitweights" +self.postfix +".npy")
		if os.path.exists(filename):
			x = load(filename)
		#else:
		if 0:
			print "sumx", sum(x)	
			dx = 1e-7
			for index in [1,30, 50, 100]:
				dxs = x * 0
				dxs[index] = dx
				#print "y", logp(x,False)
				#print "y(+dx)", logp(x+dx,False)
				#print logp(x+dxs,False), logp(x,False)
				#print (logp(x+dxs,False) - logp(x-dxs,False))/(2*dx),
				print (logp(x+dxs,False) - logp(x,False))/(dx),
				y, dydx = logp(x,True)
				print dydx[index]
				#import pdb
				#pdb.set_trace()
			print dydx.shape
			print x.shape
			#sys.exit(0)
		#x += (random.random(len(x)) - 0.5) * 2 * 1e-3
		if 0:
			print "X", logp(x)
			sys.exit(0)
		#dsa
		if 1:
			bounds = [(1e-10, 1) for i in range(len(x))]
			#bounds = None
			approx_grad = False
			if approx_grad:
				def f(x):
					return -logp(x,False)
			else:
				def f(x):
					y, dy = logp(x,True)
					return -y, -dy

			x = scipy.optimize.fmin_l_bfgs_b(f, x, None, bounds=bounds, approx_grad=approx_grad, factr=1000., pgtol=1e-6, iprint=1,maxfun=self.max_iterations)[0]
		print sum(x)	
		print "saved", filename
		save(filename, x)
		print "res 1)", logp(x)
		print "res 2)", logp(x, seperate=True)
		logL = logp(x)-exp_value
		self.logger.info("logL = %f (=%e)" % (logL, logL))
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (logL)
		f.close()
		
		return x
		sys.exit(0)
			
			
		
		
		
		
		return None
		for star in stars:
			star.aperture_index = self.storage_2d.aperture.findindex(star.xi, star.eta)
		
		vmax = self.storage_2d.vmax
		print "vmax", vmax
		delta_v = 2*vmax/self.storage_2d.Nv
		print "res", delta_v
		#sys.exit(0)
		losvds = self.storage_2d.losvds
		stars_invrange = stars.filter(lambda star: abs(star.vlos) < vmax)
		self.logger.info("stars in velocity range     : %d" % len(stars_invrange))
		self.logger.info("stars outside velocity range: %d" % (len(stars)-len(stars_invrange)))
		stars = stars_invrange
		sigma_v = 2.01
		numpy.random.seed(8)
		for star in stars:
			star.vlos = star.vlos_true + numpy.random.normal(0, sigma_v)
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d.Nv);
			outlier = True
			for losvd in losvds:
				if losvd[star.v_index, star.aperture_index] != 0:
					outlier = False
					break
			star.is_outlier = outlier
		stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
				
		self.logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		self.logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		
		self.obs_losvd = losvds[0] * 0
		for star in stars:
			self.obs_losvd[star.v_index, star.aperture_index] += 1
		
		self.logger.info("non zero losv gridpoints: %d" % sum(self.obs_losvd > 0))
		self.logger.info("(check: %f)" % sum(self.obs_losvd[self.obs_losvd > 0]))
		
		Rborders = arange(self.storage_2d.NR+1) / (0.0+self.storage_2d.NR) * (self.storage_2d.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		#print Rborders
		self.rho2d_target = array([self.light_model.cumdensityR(R1, R2) for R1, R2 in zip(R1s, R2s)])
		
		rho2ds = sum(losvds, axis=1)
		rho2dmatrix = sum(losvds, axis=1)
		
		self.vlos_centers = ((arange(self.storage_2d.Nv)+0.5)/(0.+self.storage_2d.Nv)) * 2 * self.storage_2d.vmax - self.storage_2d.vmax
		#print "vlos centers"
		#print self.vlos_centers
		#sys.exit(0) 
		
		for i in range(losvds.shape[0]):
			#print sum(losvds[i])
			for j in range(losvds.shape[1]):
				#dens = sum(losvds[i,:,j])
				#if dens > 0:
				#print self.rho2d_target.shape
				#print losvds.shape 
				#losvds[i,j,:] /= self.rho2d_target
				pass
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [sigma_v/delta_v, 14])
			losvds[i] /= (delta_v * delta_R)
		#self.storage_2d.losvds = losvds
		#sys.exit(0)
		
		
		v_indices = [star.v_index for star in stars]
		aperture_indices = [star.aperture_index for star in stars]
		
		pmatrix = array(list(losvds[:,v_indices, aperture_indices]))
		
		rho2d_error = self.rho2d_target.max() * 0.01
		error_x = 1e-4
		
		
		if 1:
			filename = os.path.join(self.modelpath, "orbitweights.npy")
			orbitweights = load(filename)
			c = orbitweights.flatten() 
			c /= sum(c)
			x = c
			self.true_losvd = numpy.tensordot(self.storage_2d.losvds, c, axes=[(0,),(0,)])
			self.true_rho2d = numpy.tensordot(rho2ds, c, axes=[(0,),(0,)])
		else:
			x = None
		x[x<=0] = 1e-20
		x = x * 0 + 1
		x /= sum(x)
		u = log(x+1e-9)
		
		print pmatrix.shape, rho2dmatrix.shape, x.shape, self.rho2d_target.shape, (self.rho2d_target*0+rho2d_error).shape, error_x
		opt = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, True, True, True)
		opt1 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, True, False, False)
		opt2 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, False, True, False)
		opt3 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, False, False, True)
		opts = [opt1, opt2, opt3]
		
		self.x = x
		if 1:
			for i1 in [0,1,2]:
				for i2 in [0,1,2]:
					self.u1 = log(self.x)
					self.u2 = log(self.x)
					du = 1e-7
					#i1 = 0
					self.u1[i2] -= du/2
					self.u2[i2] += du/2
					w1 = -opt.likelihood(self.u1)
					w2 = -opt.likelihood(self.u2)
					
					self.grad1 = x * 0
					self.grad2 = x * 0
					
					
					opt.dfdx(self.u1, self.grad1)
					opt.dfdx(self.u2, self.grad2)
					print "grad: %13f %13f %13f" % ((w2-w1)/du, self.grad1[i2], self.grad2[i2])
					#print "hess: %13f %13f " % ((self.grad2[i1] - self.grad1[i1])/du, self.hessian[i1,i2])
					
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		#if solution is not None:
		#	x = solution
		#	u = log(x+1e-10)
		#else:
		if 1:
			if 0:
				filename = os.path.join(self.modelpath, "orbitweights.npy")
				orbitweights = load(filename)
				orbitweights = array(orbitweights.flatten()) 
				orbitweights[orbitweights<0] = 0
				orbitweights /= sum(orbitweights)
				x = orbitweights
				u = log(x+1e-20) 
			elif os.path.exists(filename) and 0:
				self.logger.info("orbitweights found at: %s" % filename)
				x = numpy.load(filename)
				u = log(x)
			else:
			#if 1:
				def f_and_g(x):
					u = log(x)
					#print u
					w = -opt.likelihood(u)
					grad = u * 0
					opt.dfdx(u, grad)
					return w, grad/x
				if 1:
					#u = scipy.optimize.fmin_l_bfgs_b(f_and_g, log(x+1e-9), None, iprint=1,factr=1e3,maxfun=1000)[0]
					bounds = [(1e-10, 1) for i in range(len(x))]
					#bounds = None
					approx_grad = False
					x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=approx_grad, iprint=1,factr=1e2,maxfun=200000)[0]
					#x /= sum(x)
					u = log(x)
				else:
					print u
					opt.optimize(10000, 10, u)
					x = exp(u)
					x /= sum(x)
					u = log(x)
				numpy.save(filename, x)
		print u
		logL = opt.likelihood(u)
		print "log L", logL, -6966.89732413
		
		print "log L k/d/n", opt1.likelihood(u), opt2.likelihood(u), opt3.likelihood(u)
		print "\t", -10815.8321688 -14.8476681244 -0.00208193164913
		
		#filename = os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt")
		#f = file(filename, "w")
		#print >>f, "%e" % logL
		
		logL = opt.likelihood(u)
		self.logger.info("logL = %f (=%e)" % (logL, logL))
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (logL)
		f.close()
		
		for opt_single in opts:
			logL1 = opt_single.likelihood(u)
			self.logger.info("logL_i = %f (=%e)" % (logL1, logL1))

		#sys.exit(0)
		self.x = x
		self.c = c
		self.hessian = zeros((len(x), len(x)))
		opt.hessian(self.hessian, log(x), False)
		print "sum(x) = ", sum(x)
		self.schw_losvd = numpy.tensordot(self.storage_2d.losvds, x, axes=[(0,),(0,)])
		self.schw_rho2d = numpy.tensordot(rho2ds, x, axes=[(0,),(0,)])
		print "sum(losvd) = ", sum(self.schw_losvd)
					
		#sys.exit(0)
		#u = log(x+1e-9)
		#w = -opt.likelihood(u)
		#grad = u * 0
		#opt.dfdx(u, grad)
		#print w
		#print grad
		moments = self.storage_2d_binned.projectedmoments
		#print moments.shape, x.shape
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_los = moments[2]
		mask = moments[0] > 0
		#print var_los 
		#print x
		#print moments
		var_los[mask] /= moments[0][mask]
		self.sigma_los = sqrt(var_los)
		
		moments = self.storage_3d.moments3d
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_r = moments[4]
		var_t = moments[5]
		mask = moments[0] > 0
		#print x
		#print moments.shape
		var_r[mask] /= moments[0][mask]
		var_t[mask] /= moments[0][mask]
		#print var_r.shape
		#efjsdl
		self.sigma_r = sqrt(var_r)
		self.sigma_t = sqrt(var_t)
		self.beta = 1 - var_t/(2*var_r)
		self.store_moments(self.storage_2d_binned.projectedmoments, self.storage_3d.moments3d, x)
		
				
	
		
		
class Likelihood(Solution):
	def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d, storage_2d_binned, storage_3d, fitdensity2d, fitdensity3d, observation, dfgrid, regularization, regularization_delta=1., use_jeans=False, jeans_fraction=0.05, postfix=""):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.storage_2d_binned = storage_2d_binned
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.observation = observation
		self.dfgrid = dfgrid
		self.regularization = regularization
		self.regularization_delta = regularization_delta
		self.use_jeans = use_jeans
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.likelihood")
		self.postfix = postfix
		
	def load(self):
		self.observation.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_2d_binned.init()
		self.storage_2d_binned.load()
		self.storage_3d.init()
		self.storage_3d.load()
	
	def solve(self):
		return self.findsolution()
	
	def findsolution(self, solution=None):
		stars = self.observation.stars
		self.logger.info("using %d stars/observations" % len(stars))
		stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		self.logger.info("stars in aperture range     : %d" % len(stars_inrange))
		self.logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		stars = stars_inrange
		for star in stars:
			star.aperture_index = self.storage_2d.aperture.findindex(star.xi, star.eta)
		
		vmax = self.storage_2d.vmax
		print "vmax", vmax
		delta_v = 2*vmax/self.storage_2d.Nv
		print "res", delta_v
		#sys.exit(0)
		losvds = self.storage_2d.losvds
		stars_invrange = stars.filter(lambda star: abs(star.vlos) < vmax)
		self.logger.info("stars in velocity range     : %d" % len(stars_invrange))
		self.logger.info("stars outside velocity range: %d" % (len(stars)-len(stars_invrange)))
		stars = stars_invrange
		sigma_v = 2.01
		numpy.random.seed(8)
		for star in stars:
			star.vlos = star.vlos_true + numpy.random.normal(0, sigma_v)
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d.Nv);
			outlier = True
			for losvd in losvds:
				if losvd[star.v_index, star.aperture_index] != 0:
					outlier = False
					break
			star.is_outlier = outlier
		stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
				
		self.logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		self.logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		
		self.obs_losvd = losvds[0] * 0
		for star in stars:
			self.obs_losvd[star.v_index, star.aperture_index] += 1
		
		self.logger.info("non zero losv gridpoints: %d" % sum(self.obs_losvd > 0))
		self.logger.info("(check: %f)" % sum(self.obs_losvd[self.obs_losvd > 0]))
		
		Rborders = arange(self.storage_2d.NR+1) / (0.0+self.storage_2d.NR) * (self.storage_2d.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		#print Rborders
		self.rho2d_target = array([self.light_model.cumdensityR(R1, R2) for R1, R2 in zip(R1s, R2s)])
		
		rho2ds = sum(losvds, axis=1)
		rho2dmatrix = sum(losvds, axis=1)
		
		self.vlos_centers = ((arange(self.storage_2d.Nv)+0.5)/(0.+self.storage_2d.Nv)) * 2 * self.storage_2d.vmax - self.storage_2d.vmax
		#print "vlos centers"
		#print self.vlos_centers
		#sys.exit(0) 
		
		for i in range(losvds.shape[0]):
			#print sum(losvds[i])
			for j in range(losvds.shape[1]):
				#dens = sum(losvds[i,:,j])
				#if dens > 0:
				#print self.rho2d_target.shape
				#print losvds.shape 
				#losvds[i,j,:] /= self.rho2d_target
				pass
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [sigma_v/delta_v, 14])
			losvds[i] /= (delta_v * delta_R)
		#self.storage_2d.losvds = losvds
		#sys.exit(0)
		
		
		v_indices = [star.v_index for star in stars]
		aperture_indices = [star.aperture_index for star in stars]
		
		pmatrix = array(list(losvds[:,v_indices, aperture_indices]))
		
		rho2d_error = self.rho2d_target.max() * 0.01
		error_x = 1e-4
		
		
		if 1:
			filename = os.path.join(self.modelpath, "orbitweights.npy")
			orbitweights = load(filename)
			c = orbitweights.flatten() 
			c /= sum(c)
			x = c
			self.true_losvd = numpy.tensordot(self.storage_2d.losvds, c, axes=[(0,),(0,)])
			self.true_rho2d = numpy.tensordot(rho2ds, c, axes=[(0,),(0,)])
		else:
			x = None
		x[x<=0] = 1e-20
		x = x * 0 + 1
		x /= sum(x)
		u = log(x+1e-9)
		
		print pmatrix.shape, rho2dmatrix.shape, x.shape, self.rho2d_target.shape, (self.rho2d_target*0+rho2d_error).shape, error_x
		opt = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, True, True, True)
		opt1 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, True, False, False)
		opt2 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, False, True, False)
		opt3 = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, False, False, True)
		opts = [opt1, opt2, opt3]
		
		self.x = x
		if 1:
			for i1 in [0,1,2]:
				for i2 in [0,1,2]:
					self.u1 = log(self.x)
					self.u2 = log(self.x)
					du = 1e-7
					#i1 = 0
					self.u1[i2] -= du/2
					self.u2[i2] += du/2
					w1 = -opt.likelihood(self.u1)
					w2 = -opt.likelihood(self.u2)
					
					self.grad1 = x * 0
					self.grad2 = x * 0
					
					
					opt.dfdx(self.u1, self.grad1)
					opt.dfdx(self.u2, self.grad2)
					print "grad: %13f %13f %13f" % ((w2-w1)/du, self.grad1[i2], self.grad2[i2])
					#print "hess: %13f %13f " % ((self.grad2[i1] - self.grad1[i1])/du, self.hessian[i1,i2])
					
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		#if solution is not None:
		#	x = solution
		#	u = log(x+1e-10)
		#else:
		if 1:
			if 0:
				filename = os.path.join(self.modelpath, "orbitweights.npy")
				orbitweights = load(filename)
				orbitweights = array(orbitweights.flatten()) 
				orbitweights[orbitweights<0] = 0
				orbitweights /= sum(orbitweights)
				x = orbitweights
				u = log(x+1e-20) 
			elif os.path.exists(filename) and 0:
				self.logger.info("orbitweights found at: %s" % filename)
				x = numpy.load(filename)
				u = log(x)
			else:
			#if 1:
				def f_and_g(x):
					u = log(x)
					#print u
					w = -opt.likelihood(u)
					grad = u * 0
					opt.dfdx(u, grad)
					return w, grad/x
				if 1:
					#u = scipy.optimize.fmin_l_bfgs_b(f_and_g, log(x+1e-9), None, iprint=1,factr=1e3,maxfun=1000)[0]
					bounds = [(1e-10, 1) for i in range(len(x))]
					#bounds = None
					approx_grad = False
					x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=approx_grad, iprint=1,factr=1e2,maxfun=200000)[0]
					#x /= sum(x)
					u = log(x)
				else:
					print u
					opt.optimize(10000, 10, u)
					x = exp(u)
					x /= sum(x)
					u = log(x)
				numpy.save(filename, x)
		print u
		logL = opt.likelihood(u)
		print "log L", logL, -6966.89732413
		
		print "log L k/d/n", opt1.likelihood(u), opt2.likelihood(u), opt3.likelihood(u)
		print "\t", -10815.8321688 -14.8476681244 -0.00208193164913
		
		#filename = os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt")
		#f = file(filename, "w")
		#print >>f, "%e" % logL
		
		logL = opt.likelihood(u)
		self.logger.info("logL = %f (=%e)" % (logL, logL))
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (logL)
		f.close()
		
		for opt_single in opts:
			logL1 = opt_single.likelihood(u)
			self.logger.info("logL_i = %f (=%e)" % (logL1, logL1))

		#sys.exit(0)
		self.x = x
		self.c = c
		self.hessian = zeros((len(x), len(x)))
		opt.hessian(self.hessian, log(x), False)
		print "sum(x) = ", sum(x)
		self.schw_losvd = numpy.tensordot(self.storage_2d.losvds, x, axes=[(0,),(0,)])
		self.schw_rho2d = numpy.tensordot(rho2ds, x, axes=[(0,),(0,)])
		print "sum(losvd) = ", sum(self.schw_losvd)
					
		#sys.exit(0)
		#u = log(x+1e-9)
		#w = -opt.likelihood(u)
		#grad = u * 0
		#opt.dfdx(u, grad)
		#print w
		#print grad
		moments = self.storage_2d_binned.projectedmoments
		#print moments.shape, x.shape
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_los = moments[2]
		mask = moments[0] > 0
		#print var_los 
		#print x
		#print moments
		var_los[mask] /= moments[0][mask]
		self.sigma_los = sqrt(var_los)
		
		moments = self.storage_3d.moments3d
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_r = moments[4]
		var_t = moments[5]
		mask = moments[0] > 0
		#print x
		#print moments.shape
		var_r[mask] /= moments[0][mask]
		var_t[mask] /= moments[0][mask]
		#print var_r.shape
		#efjsdl
		self.sigma_r = sqrt(var_r)
		self.sigma_t = sqrt(var_t)
		self.beta = 1 - var_t/(2*var_r)
		self.store_moments(self.storage_2d_binned.projectedmoments, self.storage_3d.moments3d, x)
		
				
				
	def findsolution_(self):
		self.observation.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_2d_binned.init()
		self.storage_2d_binned.load()
		self.storage_3d.init()
		self.storage_3d.load()
		stars = self.observation.stars
		self.logger.info("using %d stars/observations" % len(stars))
		stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		self.logger.info("stars in aperture range     : %d" % len(stars_inrange))
		self.logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		stars = stars_inrange
		for star in stars:
			star.aperture_index = self.storage_2d.aperture.findindex(star.xi, star.eta)
		
		vmax = self.storage_2d.vmax
		losvds = self.storage_2d.losvds
		stars_invrange = stars.filter(lambda star: abs(star.vlos) < vmax)
		self.logger.info("stars in velocity range     : %d" % len(stars_invrange))
		self.logger.info("stars outside velocity range: %d" % (len(stars)-len(stars_invrange)))
		stars = stars_invrange
		for star in stars:
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d.Nv);
			outlier = True
			for losvd in losvds:
				if losvd[star.v_index, star.aperture_index] != 0:
					outlier = False
					break
			star.is_outlier = outlier
		stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
				
		self.logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		self.logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		
		self.obs_losvd = losvds[0] * 0
		for star in stars:
			self.obs_losvd[star.v_index, star.aperture_index] += 1
		
		self.logger.info("non zero losv gridpoints: %d" % sum(self.obs_losvd > 0))
		self.logger.info("(check: %f)" % sum(self.obs_losvd[self.obs_losvd > 0]))
		
		rho2ds = sum(losvds, axis=1)
		
		Rborders = arange(self.storage_2d.NR+1) / (0.0+self.storage_2d.NR) * (self.storage_2d.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		#print Rborders
		self.rho2d_target = array([self.light_model.cumdensityR(R1, R2) for R1, R2 in zip(R1s, R2s)])
		
		delta_v = 2*vmax/self.storage_2d.Nv
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		for i in range(losvds.shape[0]):
			#print sum(losvds[i])
			for j in range(losvds.shape[1]):
				#dens = sum(losvds[i,:,j])
				#if dens > 0:
				#print self.rho2d_target.shape
				#print losvds.shape 
				#losvds[i,j,:] /= self.rho2d_target
				pass
			#losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [sigma_v/delta_v, 1e-4])
			losvds[i] /= (delta_v * delta_R)
		
		
		v_indices = [star.v_index for star in stars]
		aperture_indices = [star.aperture_index for star in stars]
		
		pmatrix = array(list(losvds[:,v_indices, aperture_indices]))
		rho2dmatrix = sum(losvds, axis=1)
		
		rho2d_error = self.rho2d_target.max() * 0.01
		error_x = 1e-1
		
		
		if 1:
			filename = os.path.join(self.modelpath, "orbitweights.npy")
			orbitweights = load(filename)
			c = orbitweights.flatten() 
			c /= sum(c)
			x = c
			self.true_losvd = numpy.tensordot(self.storage_2d.losvds, c, axes=[(0,),(0,)])
			self.true_rho2d = numpy.tensordot(rho2ds, c, axes=[(0,),(0,)])
		else:
			x = None
		x[x<=0] = 1e-20
		x = x * 0 + 1
		x /= sum(x)
		u = log(x+1e-9)
		
		print pmatrix.shape, rho2dmatrix.shape, x.shape, self.rho2d_target.shape, (self.rho2d_target*0+rho2d_error).shape, error_x
		# kin/light/norm
		opt = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, x, self.rho2d_target, self.rho2d_target*0+rho2d_error, error_x, True, True, True)
		
		def f_and_g(u):
			w = -opt.likelihood(u)
			grad = u * 0
			opt.dfdx(u, grad)
			return w, grad
		#u = log(x+1e-9)
		#w = -opt.likelihood(u)
		#grad = u * 0
		#opt.dfdx(u, grad)
		#print w
		#print grad
		
		if 0:
			filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
			#filename = os.path.join(self.dirname, "../../../orbitweights.npy")
			if os.path.exists(filename):
				self.logger.info("orbitweights found at: %s" % filename)
				x = numpy.load(filename)
				if 0:
					x /= sum(x)
					x *= 1.001
					u = log(x)
			else:
				#u = scipy.optimize.fmin_l_bfgs_b(f_and_g, log(x+1e-9), None, iprint=1,factr=1e3,maxfun=1000)[0]
				print u
		if 1:
			def f_and_g(x):
				u = log(x)
				#print u
				w = -opt.likelihood(u)
				grad = u * 0
				opt.dfdx(u, grad)
				return w, grad/x
			#u = scipy.optimize.fmin_l_bfgs_b(f_and_g, log(x+1e-9), None, iprint=1,factr=1e3,maxfun=1000)[0]
			bounds = [(1e-10, 1) for i in range(len(x))]
			#bounds = None
			approx_grad = False
			x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=approx_grad, iprint=1,factr=1e2,maxfun=10000)[0]
			#x /= sum(x)
			u = log(x)
		if 0:
			print x.shape
			print sum(x)
			#x /= sum(x)
			#x *= 0.99
			#u = log(x)
			print "start:", -opt.likelihood(u)
			#opt.optimize(100, 1, u)
			opt.optimize(10000, 10, u)
			x = exp(u)
			#numpy.save(filename, x)
		numpy.save(filename, x)
		logL = opt.likelihood(u)
		self.logger.info("logL = %f (=%e)" % (logL, logL))
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (logL)
		f.close()
		
		self.x = x
		self.c = c
		self.hessian = zeros((len(x), len(x)))
		opt.hessian(self.hessian, log(x), False)
		print "sum(x) = ", sum(x)
		self.schw_losvd = numpy.tensordot(self.storage_2d.losvds, x, axes=[(0,),(0,)])
		self.schw_rho2d = numpy.tensordot(rho2ds, x, axes=[(0,),(0,)])
		#print numpy.linalg.eigvals(self.hessian)
		#L = numpy.linalg.cholesky(self.hessian)
		#import pdb;
		#pdb.set_trace()
		#sys.exit(0)
		
		
		
		for i1 in [0,1,2]:
			for i2 in [0,1,2]:
				self.u1 = log(self.x)
				self.u2 = log(self.x)
				du = 1e-7
				#i1 = 0
				self.u1[i2] -= du/2
				self.u2[i2] += du/2
				w1 = -opt.likelihood(self.u1)
				w2 = -opt.likelihood(self.u2)
				
				self.grad1 = x * 0
				self.grad2 = x * 0
				
				
				opt.dfdx(self.u1, self.grad1)
				opt.dfdx(self.u2, self.grad2)
				#print "grad: %13f %13f %13f" % ((w2-w1)/du, self.grad1[i2], self.grad2[i2])
				#print "hess: %13f %13f " % ((self.grad2[i1] - self.grad1[i1])/du, self.hessian[i1,i2])
		
		self.cov = numpy.matrix(self.hessian).I
		u, s, vd = numpy.linalg.svd(self.hessian)
		self.u = u
		self.vd = vd
		self.s = s
		#print s
		#print sqrt(s)
		#print sqrt(1/s)
		
		S = zeros((len(s), len(s)))
		for i in range(10):
			S[-1-i,-1-i] = 1/s[-1-i]
			#S[i,i] = 1/s[i]
		#print S
		#self.hessian2 = dot(u, dot(S, vd))
		#self.cov2 = numpy.matrix(self.hessian2).I
		self.cov2 = numpy.matrix(vd).T * S * numpy.matrix(u).T
		u, s, vd = numpy.linalg.svd(self.cov2)
		self.u2 = u
		self.vd2 = vd
		self.s2 = s
		import pdb
		#pdb.set_trace()
		#print numpy.linalg.det(self.cov2)
		#sys.exit(0)
		#print "eigenvalues", sqrt(s)
		#print numpy.linalg.eigvals(self.cov2)
		u, s, vd = numpy.linalg.svd(self.cov)
		self.u3 = u
		self.vd3 = vd
		self.s3 = s
		
		moments = self.storage_2d_binned.projectedmoments
		#print moments.shape, x.shape
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_los = moments[2]
		mask = moments[0] > 0
		#print var_los 
		#print x
		#print moments
		var_los[mask] /= moments[0][mask]
		self.sigma_los = sqrt(var_los)
		
		moments = self.storage_3d.moments3d
		moments = numpy.tensordot(moments, x, axes=[(0,),(0,)])
		var_r = moments[4]
		var_t = moments[5]
		mask = moments[0] > 0
		#print x
		#print moments.shape
		var_r[mask] /= moments[0][mask]
		var_t[mask] /= moments[0][mask]
		#print var_r.shape
		#efjsdl
		self.sigma_r = sqrt(var_r)
		self.sigma_t = sqrt(var_t)
		self.beta = 1 - var_t/var_r
		

		return x
		
		
		"""pmatrix = array(list(losvds[:,starvlosis, starRis]))
		p_rho2dmatrix = sum(losvds[:,:,starRis], axis=1)
		rho2dmatrix = array(list(sum(losvds, axis=1)))
		print pmatrix.shape
		print p_rho2dmatrix.shape
		from mab.gd import gdfast_schw 
		import pyublas
		import pdb;
		#pdb.set_trace()
		l = [pmatrix, rho2dmatrix, orbitweights_target, rho2d_true, rho2d_true*0+rho2d_error]
		[pyublas.why_not(k) for k in l]
		print [k.shape for k in l]
		o = gdfast_schw.OptimizationProblemSchw(pmatrix, rho2dmatrix, orbitweights_target, rho2d_true, rho2d_true*0+rho2d_error, error_x)
				
					
		rho2d_error = rho2d_true.max() * 0.01
		fudgefactor = 1e0
		error_x = 1e-2"""
		
		
		#sys.exit(0)

class SolutionMoments(object):
	
	def load(self):
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		self.orbitweights = numpy.load(filename)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +".npy")
		self.moments_solution = numpy.load(filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		self.moments3d_solution = numpy.load(filename)
		
	def calculate_solution_moments(self, moments):
		return tensordot(self.orbitweights, moments, axes=[(0,), (0,)])
		
	def solve(self):
		x = orbitweights = self.findsolution()
		
		moments = tensordot(orbitweights, self.storage_2d.projectedmoments, axes=[(0,), (0,)])
		moments3d = tensordot(orbitweights, self.storage_3d.moments3d, axes=[(0,), (0,)])
		moments[1:] /= moments[0]
		moments3d[1:] /= moments3d[0]
		
		self.moments2d_solution = moments
		self.moments3d_solution = moments3d
		
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		numpy.save(filename, orbitweights)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +".npy")
		numpy.save(filename, moments)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		numpy.save(filename, moments3d)
		
		
		self.logger.info("sum of orbitweight %f" % sum(orbitweights))
		self.logger.info("total light (2d) %f" % sum(moments[0]))
		
		
		f = file(os.path.join(self.dirname, "results/statistics_qp" +self.postfix +".txt"), "w")
		#print "lenghts", len(kinematics.m2), len(varvlos)
		dof = len(self.binned_data.moments[2])
		chisq_m2_check = sum( ((self.binned_data.moments[2]-moments[2])/self.binned_data.e_moments[2])**2 )
		chisq_m4_check = sum( ((self.binned_data.moments[4]-moments[4])/self.binned_data.e_moments[4])**2 )
		self.logger.info("(check) chisq 2nd moment: %f (reduced: %f, dof: %d)" % (chisq_m2_check, chisq_m2_check/dof, dof))
		self.logger.info("(check) chisq 4th moment: %f (reduced: %f, dof: %d)" % (chisq_m4_check, chisq_m4_check/dof, dof))
		#print chis**2
		#print kinematic_m2
		#print moments[2]
		#print e_kinematic_m2
		chisq_tot, chisqs = self.calc_chisq(x, print_info=True, returnall=True)
		print >>f, "chisq, redchsiq, dof"
		#chisq = sum(chis**2)
		print >>f, chisq_tot, ",", chisq_tot/dof,",", dof
				
		#print "chisq:",chisq_tot, chisq
		
		chisq = chisq_m2_check + chisq_m4_check
		f = file(os.path.join(self.dirname, "results/probability" +self.postfix +".txt"), "w")
		print >>f, "%e" % exp(-0.5*chisq)
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (-0.5*chisq)
		
		logps = [(name, -0.5*chisq) for name, chisq in chisqs] 
		f = file(os.path.join(self.dirname, "results/logprobability_seperate" +self.postfix +".txt"), "w")
		print >>f, repr(logps)
		
		return orbitweights
		
	def calc_chisq(self, x, print_info=False, returnall=False):
		dof = len(self.binned_data.moments[2])
		chisq_tot = 0
		chisqs = []
		
		chisq_m2 = dot(dot(x, self.P_m2),x)
		chisq_tot += chisq_m2
		chisqs.append(("m2", chisq_m2))
		if print_info:
			self.logger.info("chisq 2nd moment: %f (reduced: %f, dof: %d)" % (chisq_m2, chisq_m2/dof, dof))
		 
		chisq_m4 = dot(dot(x, self.P_m4),x)
		chisq_tot += chisq_m4
		chisqs.append(("m4", chisq_m4))
		if print_info: 
			self.logger.info("chisq 4th moment: %f (reduced: %f, dof: %d)" % (chisq_m4, chisq_m4/dof, dof))
		if self.regularization: 
			chisq_reg = dot(dot(x, self.P_reg),x) 
			chisqs.append(("reg", chisq_reg))
			chisq_tot += chisq_reg
			if print_info: 
				self.logger.info("chisq regular.  : %f (reduced: %f)" % (chisq_reg, chisq_reg/dof))
			#self.logger.info("chisq 2nd moment: %f (reduced: %f)" % (chisq_m22, chisq_m22/dof))
		if print_info:
			self.logger.info("chisq total     : %f (reduced: %f)" % (chisq_tot, chisq_tot/dof))
		#filename = os.path.join(self.dirname, ...
		if returnall:
			return chisq_tot, chisqs
		else:
			return chisq_tot
	
	def run(self, args, opts, scope):
		self.init()
		self.solve()

	

class QP(SolutionMoments):
	def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d, storage_3d, fitdensity2d, fitdensity3d, binned_data, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix=""):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.binned_data = binned_data
		self.dfgrid = dfgrid
		self.regularization = regularization
		#self.regularization_delta = regularization_delta
		self.require_stability = require_stability
		self.use_jeans = use_jeans
		self.use_fourth_moment = use_fourth_moment
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.qp")
		self.postfix = postfix

		
	def init(self):
		self.storage_2d.load()
		self.storage_3d.load()
		self.binned_data.load()
		
	
	def findsolution(self):
		logger.info("using quadratic programming to find a solution")
		n_nodes = self.storage_2d.projectedmoments.shape[0]
		nI1nodes = self.dfgrid.n_I1
		nI2nodes = self.dfgrid.n_I2
		
		constraints = len(self.binned_data.moments[0])
		dof_3d = self.storage_3d.moments3d.shape[2]
		logger.info("dof 3d: %d" % dof_3d)
		logger.info("unknown weights: %d" % n_nodes)
		logger.info("constraints: %d" % constraints)
		logger.debug("shape of projectedmoments: %r" % (self.storage_2d.projectedmoments.shape,))
		rescale = 1.
		
		Gpos = cvxopt.spmatrix(-1, range(n_nodes), range(n_nodes))
		#print "Gpos", array(Gpos).shape, Gpos
		hpos = cvxopt.matrix(0, size=(n_nodes, 1))
		assert not (self.fitdensity2d and self.fitdensity3d), "cannot fit both 2d and 3d density"
		
		
		indices1 = []
		indices2 = []
		values = []
		h_stability = []
		#self.require_stability = True
		if self.require_stability: #self.dfgrid.order == 0:
			filename = os.path.join(self.modelpath, "data", "DoS.npy")
			DoS = load(filename)
			
			assert self.dfgrid.order == 0
			n = 0
			for i2 in range(nI2nodes):
				for i1 in range(nI1nodes-1):
						lowE = self.dfgrid.index_to_orbitnr(i1,i2, 0)
						highE = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
						# -lowE + highE <= 0
						indices1.append(lowE)
						indices2.append(n)
						values.append(-1./DoS[lowE])
						
						indices1.append(highE)
						indices2.append(n)
						values.append(1./DoS[highE])
						h_stability.append(0.0)
						n += 1
			G_stability_matrix = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		
		
		
		densitiesprojected = self.storage_2d.projectedmoments[:,0,:]
		densityprojected_target = self.binned_data.moments[0]
		
		moments3d = self.storage_3d.moments3d
		
		rborders = self.storage_3d.xborders
		r = self.storage_3d.x
		
		#light_model
		#density_target = self.storage_3d.moments3d[:,0,:]
		#print density_target.shape
		#density_target
		
		kinematic_m2 = self.binned_data.moments[2]
		schw_m2 = self.storage_2d.projectedmoments[:,2,:]
		e_kinematic_m2 = self.binned_data.e_moments[2]
		kinematic_m4 = self.binned_data.moments[4]
		schw_m4 = self.storage_2d.projectedmoments[:,4,:]
		e_kinematic_m4 = self.binned_data.e_moments[4]
		logger.debug("shape of densitiesprojected: %r" % (densitiesprojected.shape,))
		logger.debug("shape of densityprojected_target: %r" % (densityprojected_target.shape,))
		
		light_model = self.light_model
		Erborders = self.dfgrid.r_borders
		#density_Etarget = array([galaxy.stellar_profile.cumdensityr(galaxy.arcsec_to_kpc(Erborders[i]), galaxy.arcsec_to_kpc(Erborders[i+1])) for i in range(self.dfgrid.n_I1)])
		density_Etarget = array([light_model.light_profile.cumdensityr((Erborders[i]), (Erborders[i+1]), M=1.) for i in range(self.dfgrid.n_I1)])
		if 0:
			from kaplot import *
			box()
			#graph(arange(len(density_Etarget)), density_Etarget)
			l = 0.8
			E = self.dfgrid.E_borders[len(self.dfgrid.E_borders)/2]
			Lmax = self.profile_model.Lmax_at_E(E)
			#Ls = [l * self.profile_model.Lmax_at_E(E) for E in self.dfgrid.Es]
			v = []
			#for i, (L, E) in enumerate(zip(Ls, self.dfgrid.Es)):
			#	deltaE = self.dfgrid.E_borders[i+1] - self.dfgrid.E_borders[i]
			#	print deltaE
			#	v.append(self.profile_model.gEL(E,L) / deltaE)
			for l in arange(0, 1, 0.01):
				L = l * Lmax
				v.append(self.profile_model.gEL(E,L))
				print v
			graph(arange(len(v)), v)
			draw()
	
		# this part constrainst the density to within a certain range
		if self.fitdensity3d: # this part is for testing, it uses the 3d density
			diff = (density_target*0.01).max()
			Gmin = cvxopt.matrix(-transpose(densities))
			hmin = cvxopt.matrix(minimum(-(density_target-diff/2), 0), size=(Nr3d, 1))
			Gmax = cvxopt.matrix(transpose(densities))
			hmax = cvxopt.matrix(density_target+diff/2, size=(Nr3d, 1))
		if self.fitdensity2d:
			diff = max(densityprojected_target*0.01)
			Gmin = cvxopt.matrix(-1/rescale*transpose(densitiesprojected))
			hmin = cvxopt.matrix(minimum(-(densityprojected_target-diff/2), 0), size=(constraints, 1))
			Gmax = cvxopt.matrix(1/rescale*transpose(densitiesprojected))
			hmax = cvxopt.matrix(densityprojected_target+diff/2, size=(constraints, 1))
			
		if self.use_jeans:
			dr = rborders[1:] - rborders[:-1]
			rc = (r[1:] + r[:-1])/2
			drc = r[1:] - r[:-1]
			nu = moments3d[:,0,:] / (4*pi*r**2*dr)
			mask = moments3d[:,0,:] > 0
			
			varvr = moments3d[:,4,:] * 0
			varvr[mask] = moments3d[:,4,:][mask]/moments3d[:,0,:][mask]
			
			varvphi = moments3d[:,5,:] * 0
			varvphi[mask] = moments3d[:,5,:][mask]/moments3d[:,0,:][mask]
			
			varvtheta = moments3d[:,6,:] * 0
			varvtheta[mask] = moments3d[:,6,:][mask]/moments3d[:,0,:][mask]
			
			#varvphi = moments3d[:,5,:]/moments3d[:,0,:]
			#varvtheta = moments3d[:,6,:]/moments3d[:,0,:]
			assert all(~isnan(varvr))
			assert all(~isnan(varvphi))
			assert all(~isnan(varvtheta))
			#print varvr.shape, nu.shape
			rho_varvr = varvr * nu
			#print rho_varvr.shape, drc.shape
			jeans_term1 = (rho_varvr[:,1:] - rho_varvr[:,:-1]) / drc
			
			varvrc = (varvr[:,:-1] + varvr[:,1:])/2
			varvphic = (varvphi[:,:-1] + varvphi[:,1:])/2
			varvthetac = (varvtheta[:,:-1] + varvtheta[:,1:])/2
			
			nuc = (nu[:,:-1] + nu[:,1:])/2
			#jeans_term2 = 2 * -0.5 * nuc * varvrc/rc
			jeans_term2 = nuc * (2*varvrc - varvphic - varvthetac) / rc
			#print sum(moments3d[0,0])
			
			#dphidr = rc * 0.
			#for i in range(len(rc)):
			dphidr = self.profile_model.dphidr(rc)
			#print jeans_term1.shape, jeans_term2.shape
			jeans =  jeans_term1 + jeans_term2
			jeans_term_rhs = -self.light_model.densityr(rc)/self.light_model.light_profile.M * self.profile_model.dphidr(rc)
			#print nu.shape, jeans.shape, jeans_term_rhs.shape
			#print jeans,jeans_term_rhs
			J = (jeans/jeans_term_rhs)
			#print J.max(), J.min()
			
			
			JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs))-1)
			Jhmin = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			
			#JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs)))
			#Jhmin = cvxopt.matrix(1.25, size=(constraints-1, 1))
			
			JGmax = cvxopt.matrix(-(transpose(jeans/jeans_term_rhs)-1))
			Jhmax = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			#JGmax = cvxopt.matrix((transpose(jeans/jeans_term_rhs)-1))
			#Jhmax = cvxopt.matrix(1.51, size=(constraints-1, 1))
			#print Jhmax
			#print jeans
			
			#print array(JGmin).shape
			#print array(JGmax).shape
			#print " "
			#print array(Jhmin).shape
			#print array(Jhmax).shape
			
			
			#target = graph(xc, log10(nuc * dphidr), color="orange", linestyle="dash")
			import sys
			#sys.exit(0)
			
		if self.regularization:
			Greg, hreg = self.regularization.assemble_system()
			delta = self.regularization.regularization_delta
		elif 0:
			# regularisation part
			indices1 = []
			indices2 = []
			values = []
			hreg = []
			#delta = 10000.1
			logger.info("regularization delta: %.5f" % self.regularization_delta)
			delta =  self.regularization_delta #.regdelta #sqrt(delta / (2*orbits))
			#delta = 1.0
			#delta = sqrt(delta / (2*orbits))
			logger.info("internal regularization delta: %.5f" % delta)
			n = 0
			s = 1.
			scale_i1 = 1 * s
			scale_i2 = 4 * s
			rscale = 1e-4
			
			#s = 1
			#scale_i1 = 1 * s
			#scale_i2 = 8 * s
			#rscale = 1e-4
			
			#s = 0.01
			#scale_i1 = 1 * s
			#scale_i2 = 1 * s
			#rscale = 1e-4
			if 0:
				for i in range(4):
					print "=" * 10
					for vi in range(7):
						for ui in range(7):
							u = ui/6.
							v = vi/6.
							print "% 8.6e" % self.dfgrid.basis(i, u, v),
						print
					print "=" * 10
			#sys.exit(0)
			if self.dfgrid.order == 0:
				for i1 in range(1,nI1nodes-1):
					for i2 in range(nI2nodes):
						extras = 1 #/ ((i2 + 0.5) / (nI2nodes))
						o1 = self.dfgrid.index_to_orbitnr(i1-1,i2, 0) 
						o2 = self.dfgrid.index_to_orbitnr(i1+0,i2, 0)
						o3 = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
						# minimize -x_{i-1,j} + 2 x_{i,j} -x_{i+1,j} 
						r1 = 1 + (random.random() - 0.5) * rscale
						r2 = 1 + (random.random() - 0.5) * rscale
						r3 = 1 + (random.random() - 0.5) * rscale
						values.append(-1.*r1/density_Etarget[i1-1]/scale_i1*extras)
						assert o1 < n_nodes
						assert o2 < n_nodes
						assert o3 < n_nodes, "o3 = %d, n_nodes = %d" % (o3, n_nodes)
						indices1.append(o1)
						indices2.append(n)
						
						values.append(2.*r2/density_Etarget[i1]/scale_i1*extras)
						indices1.append(o2)
						indices2.append(n)
						
						values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
						indices1.append(o3)
						indices2.append(n)
						rh = (random.random() - 0.5) * delta * 1e-5
						hreg.append(rh)
						n += 1
				for i1 in range(0,nI1nodes):
					for i2 in range(1,nI2nodes-1):
						extras1 = 1 #/ ((i2- 0.5) / (nI2nodes))
						extras2 = 1 #/ ((i2 + 0.) / (nI2nodes))
						extras3 = 1 #/ ((i2 + 0.5) / (nI2nodes))
						o1 = self.dfgrid.index_to_orbitnr(i1,i2-1, 0) 
						o2 = self.dfgrid.index_to_orbitnr(i1,i2, 0)
						o3 = self.dfgrid.index_to_orbitnr(i1,i2+1, 0)
						# minimize -x_{i,j-1} + 2 x_{i,j} -x_{i,j+1} 
						r1 = 1 + (random.random() - 0.5) * rscale
						r2 = 1 + (random.random() - 0.5) * rscale
						r3 = 1 + (random.random() - 0.5) * rscale
						values.append(-1.*r1/density_Etarget[i1]/scale_i2*extras1)
						indices1.append(o1)
						indices2.append(n)
						
						values.append(2.*r2/density_Etarget[i1]/scale_i2*extras2)
						indices1.append(o2)
						indices2.append(n)
						
						values.append(-1.*r3/density_Etarget[i1]/scale_i2*extras3)
						indices1.append(o3)
						indices2.append(n)
						rh = (random.random() - 0.5) * delta * 1e-5
						hreg.append(rh)
						n += 1
			else:
				s = 14.0
				scale_i1 = 1 * s
				scale_i2 = 0.1 * s
				#s = 2.8
				#scale_i1 = 1.0 * s
				#scale_i2 = 1.05 * s
				#s = 0.40
				#scale_i1 = 1.0 * s
				#scale_i2 = 0.5 * s
				rscale = 1e-4
				assert self.dfgrid.order == 1
				for i1 in range(nI1nodes-1):
					for i2 in range(nI2nodes):
						extras = 1 # / ((i2 + 0.5) / (nI2nodes))
						#i1 = 0
						#i2 = 0
						dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
						dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
						dof_index20 = self.dfgrid.dof_index(i1+1, i2, 1)
						dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
						dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
						dof_index21 = self.dfgrid.dof_index(i1+1, i2, 3)
						if (i1 == 0) and (i2 == 0):
							print dof_index00, dof_index10, dof_index20, dof_index01, dof_index11, dof_index21
						
						if 1:
							#print dof_index00, dof_index10 
							r1 = 1 + (random.random() - 0.5) * rscale
							r2 = 1 + (random.random() - 0.5) * rscale
							r3 = 1 + (random.random() - 0.5) * rscale
							
							values.append(-1.*r1/density_Etarget[i1]/scale_i1*extras)
							indices1.append(dof_index00)
							indices2.append(n)
							
							values.append(2.*r2/((density_Etarget[i1]+density_Etarget[i1+1])/2)/scale_i1*extras)
							#values.append(2.*r2/density_Etarget[i1]/scale_i1)
							indices1.append(dof_index10)
							indices2.append(n)
							
							values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
							indices1.append(dof_index20)
							indices2.append(n)
							
							rh = (random.random() - 0.5) * delta * 1e-5
							hreg.append(rh)
							n += 1
							
							
							values.append(-1.*r1/density_Etarget[i1]/scale_i1*extras)
							indices1.append(dof_index01)
							indices2.append(n)
							
							#values.append(2.*r2/density_Etarget[i1]/scale_i1)
							values.append(2.*r2/((density_Etarget[i1]+density_Etarget[i1+1])/2)/scale_i1*extras)
							indices1.append(dof_index11)
							indices2.append(n)
							
							values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
							indices1.append(dof_index21)
							indices2.append(n)
							
							rh = (random.random() - 0.5) * delta * 1e-5
							hreg.append(rh)
							n += 1
				for i1 in range(nI1nodes):
					for i2 in range(nI2nodes-1):
						extras = 1 # / ((i2 + 0.5) / (nI2nodes-1))
						#i1 = 0
						#i2 = 0
						dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
						dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
						dof_index02 = self.dfgrid.dof_index(i1, i2+1, 2)
						dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
						dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
						dof_index12 = self.dfgrid.dof_index(i1, i2+1, 3)
						if (i1 == 0) and (i2 == 0):
							print dof_index00, dof_index10, dof_index01, dof_index11, dof_index02, dof_index12
						if 1:
							#print dof_index00, dof_index10 
							r1 = 1 + (random.random() - 0.5) * rscale
							r2 = 1 + (random.random() - 0.5) * rscale
							r3 = 1 + (random.random() - 0.5) * rscale
							
							if i1 > 0:
								dl = (density_Etarget[i1] + density_Etarget[i1-1])/2
							else:
								dl = density_Etarget[i1]
							if i1 < nI1nodes-1:
								dr = (density_Etarget[i1] + density_Etarget[i1+1])/2
								#dr = density_Etarget[i1+1]
							else:
								dr = density_Etarget[i1]
							values.append(-1.*r1/dl/scale_i2*extras)
							indices1.append(dof_index00)
							indices2.append(n)
							
							values.append(2.*r2/dl/scale_i2*extras)
							indices1.append(dof_index01)
							indices2.append(n)
							
							values.append(-1.*r3/dl/scale_i2*extras)
							indices1.append(dof_index02)
							indices2.append(n)
							
							rh = (random.random() - 0.5) * delta * 1e-5
							hreg.append(rh)
							n += 1
							
							
							values.append(-1.*r1/dr/scale_i2*extras)
							indices1.append(dof_index10)
							indices2.append(n)
							
							values.append(2.*r2/dr/scale_i2*extras)
							indices1.append(dof_index11)
							indices2.append(n)
							
							values.append(-1.*r3/dr/scale_i2*extras)
							indices1.append(dof_index12)
							indices2.append(n)
							
							rh = (random.random() - 0.5) * delta * 1e-5
							hreg.append(rh)
							n += 1
							
				
			# regularisation matrix
			#print len(values), len(indices1), len(indices2), (n_nodes, n), max(indices1), max(indices2)
			Greg = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		
		#print "jeans?", self.use_jeans
		if 0: #regularisation:
			G = cvxopt.matrix([Gpos, Gmin, Gmax, Greg])
			h = cvxopt.matrix([hpos, hmin, hmax, hreg])
		else:
			Glist = []
			hlist = []
			if self.use_jeans:
				#G = cvxopt.matrix([Gpos, Gmin, Gmax])
				#h = cvxopt.matrix([hpos, hmin, hmax])
				G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin, JGmax])
				h = cvxopt.matrix([hpos, hmin, hmax, Jhmin, Jhmax])
				
				print array(Gpos).shape
				print array(Gmin).shape
				print array(Gmax).shape
				print array(JGmin).shape
				print array(JGmax).shape
				print array(G).shape
				print " "
				print array(hpos).shape
				print array(hmin).shape
				print array(hmax).shape
				print array(Jhmin).shape
				print array(Jhmax).shape
				print array(h).shape
				
				#G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin])
				#h = cvxopt.matrix([hpos, hmin, hmax, Jhmin])
			else:
				Glist.extend([Gpos, Gmin, Gmax])
				hlist.extend([hpos, hmin, hmax])
				if self.require_stability:
					Glist.append(cvxopt.matrix(transpose(G_stability_matrix)))
					hlist.append(cvxopt.matrix(h_stability))
					
				G = cvxopt.matrix(Glist)
				h = cvxopt.matrix(hlist)
		
		P = cvxopt.spmatrix(0., range(n_nodes), range(n_nodes))
		q = cvxopt.matrix(0., size=(n_nodes, 1))
		if 0:
			P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
			x = array(sigmar**2)
			mask = density_target > 0
			x[mask] *= density_target[mask] 
			q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
		else:
			if 0:
				Pt = hstack((varvloss/kinematics.em2, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
				x *= hstack((densityprojected_target2, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			elif 0:
				print varvloss.shape, n, n_nodes, Greg.shape
				Pt = hstack((Greg/delta, varvloss/kinematics.em2))#, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((array(hreg)/delta, kinematics.m2/kinematics.em2))#, kinematics.m4/kinematics.em4))
				x *= hstack((ones(n), densityprojected_target2))#, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			else:
				#print varvloss.shape, n, orbits, Greg.shape
				#print Greg.shape
				#print kinematic_m2.shape
				#print schw_m2.shape
				Pt_list = []
				if self.regularization:
					Pt_list.append(Greg/delta)
				Pt_list.append((kinematic_m2- schw_m2/densityprojected_target)/(e_kinematic_m2))
				#Pt_list.append((kinematic_m2*densityprojected_target- schw_m2)/(e_kinematic_m2))
				if self.use_fourth_moment:
					Pt_list.append((kinematic_m4- schw_m4/densityprojected_target)/(e_kinematic_m4))
					#Pt = 1/rescale*hstack((Greg/delta, kinematic_m2/e_kinematic_m2*densitiesprojected- schw_m2/e_kinematic_m2))
					#Pt = 1/rescale*hstack((Greg/delta,\
					#(kinematic_m2*densitiesprojected- schw_m2)/(densityprojected_target*e_kinematic_m2),\
					#(kinematic_m4*densitiesprojected - schw_m4)/(densityprojected_target*e_kinematic_m4) ))
					#Pt = 1/rescale*hstack((Greg/delta,\
					#(kinematic_m2- schw_m2/densityprojected_target)/(e_kinematic_m2),\
					#(kinematic_m4- schw_m4/densityprojected_target)/(e_kinematic_m4) ))
					
					#, kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
				#else:
					#Pt = hstack(((kinematic_m2/e_kinematic_m2*densityprojected_target- schw_m2/e_kinematic_m2),\
					#	(kinematic_m4/e_kinematic_m4*densityprojected_target- schw_m4/e_kinematic_m4)))
				#	Pt = hstack((	(kinematic_m2 - schw_m2/densityprojected_target)/e_kinematic_m2,\
				#					(kinematic_m4 - schw_m4/densityprojected_target)/e_kinematic_m4 ))
							# kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
				Pt = hstack(tuple(Pt_list))
				test = tensordot(Pt, Pt, axes=([1],[1]))
				P = cvxopt.matrix(test)
				x_list = []
				if self.regularization:
					x_list.append(array(hreg)/delta)
					#x = hstack((array(hreg)/delta, zeros(constraints), zeros(constraints)))#, zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
					#x *= hstack((ones(n), zeros(constraints), zeros(constraints)))
				#else:
					#x = hstack((zeros(constraints), zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
					#x *= hstack((zeros(constraints), zeros(constraints)))
				x_list.append(zeros(constraints))
				if self.use_fourth_moment:
					x_list.append(zeros(constraints))
				x = hstack(tuple(x_list))
				#print Pt.shape, x.shape
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
				#Pt = ((kinematics.m2/kinematics.em2*densitiesprojected - varvloss/kinematics.em2))
				#P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				#x = ((kinematics.m2/kinematics.em2))
				#x *= ((densityprojected_target2))
			
			
		A = cvxopt.matrix(1/rescale, size=(1, n_nodes))
		b = cvxopt.matrix([1.0])
		
		orbitweights, chisq_solution = self.solve_problem(P, q, G, h, A, b)
		#cvxopt.solvers.options["maxiters"] = 500
		#solution = cvxopt.solvers.qp(P, q, G, h, A, b)
		#chisq_solution = solution['primal objective']
		#import pdb
		#pdb.set_trace()
		
		#print solution['x']
		#orbitweights = array(solution['x']).flatten()
		if self.regularization:
			self.P_reg = tensordot(Greg/delta, Greg/delta, axes=[(1,), (1,)])
		Ph = (kinematic_m2 - schw_m2/densityprojected_target)/e_kinematic_m2
		self.P_m2 = tensordot(Ph, Ph, axes=[(1,), (1,)])
		Ph = (kinematic_m4 - schw_m4/densityprojected_target)/e_kinematic_m4
		self.P_m4 = tensordot(Ph, Ph, axes=[(1,), (1,)])
		P = tensordot(Pt, Pt, axes=[(1,), (1,)])
		chisq_tot = tensordot(orbitweights, tensordot(orbitweights, P, axes=[(0,), (0,)]), axes=[(0,),(0,)]) +\
				tensordot(orbitweights, q, axes=[(0,), (0,)])
		logger.info("chisq_tot: %s %s" % (chisq_tot, chisq_solution*2))
		self.hessian = P
		
		
		return orbitweights
		
	def solve_problem(self, P, q, G, h, A, b):
		cvxopt.solvers.options["maxiters"] = 500
		cvxopt.solvers.options['show_progress'] = False
		solution = cvxopt.solvers.qp(P, q, G, h, A, b)
		chisq_solution = solution['primal objective']
		orbitweights = array(solution['x']).flatten()
		return orbitweights, chisq_solution
		
		
class KnownWeights(QP):
	"""def __init__(self, modelpath, schwsetname, schwmodelname, dfgrid, storage_2d, storage_3d, binned_data, postfix=""):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dfgrid = dfgrid
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.binned_data = binned_data
		self.light_model = dfgrid.light_model
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.postfix = postfix
		self.logger = logging.getLogger("gd.schw.solution.known")"""
		
		
	"""def load(self):
		self.binned_data.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_3d.load()"""
	
	def load(self):
		self.storage_2d.load()
		self.storage_3d.load()
		self.binned_data.load()
	
	def solve_problem(self, P, q, G, h, A, b):
		filename = os.path.join(self.modelpath, "df", "orbitweights.npy")
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		return c, "?"

class KnownWeights2(SolutionMoments):
	def __init__(self, modelpath):
		self.modelpath = modelpath
	"""def __init__(self, modelpath, schwsetname, schwmodelname, dfgrid, storage_2d, storage_3d, binned_data, postfix=""):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dfgrid = dfgrid
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.binned_data = binned_data
		self.light_model = dfgrid.light_model
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.postfix = postfix
		self.logger = logging.getLogger("gd.schw.solution.known")"""
		
		
	"""def load(self):
		self.binned_data.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_3d.load()"""
	
	def load(self):
		pass
		#self.storage_2d.load()
		#self.storage_3d.load()
		#self.binned_data.load()
	
	def solve(self):
		filename = os.path.join(self.modelpath, "df", "orbitweights.npy")
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		return c
			
if 0:
			"""if 1:
				for k in range(4):
					print ">> ", self.dfgrid.dof_index(0,0,k), self.dfgrid.dof_index(1,0,k), self.dfgrid.dof_index(0,1,k)
				#for i1 in range(1,nI1nodes-1):
				if 0:
					for i1 in range(nI1nodes):
						for i2 in range(nI2nodes):
							extras = 1 # / ((i2 + 0.5) / (nI2nodes))
							#i1 = 0
							#i2 = 0
							dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
							dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
							dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
							dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
							if 1:
								#print dof_index00, dof_index10 
								r1 = 1 + (random.random() - 0.5) * rscale
								r2 = 1 + (random.random() - 0.5) * rscale
								
								values.append(-1.*r1/density_Etarget[i1]/scale_i1*extras)
								indices1.append(dof_index00)
								indices2.append(n)
								
								values.append(1.*r2/density_Etarget[i1]/scale_i1*extras)
								#values.append(2.*r2/density_Etarget[i1]/scale_i1)
								indices1.append(dof_index10)
								indices2.append(n)
								
								rh = (random.random() - 0.5) * delta * 1e-5
								hreg.append(rh)
								n += 1
								
								r1 = 1 + (random.random() - 0.5) * rscale
								r2 = 1 + (random.random() - 0.5) * rscale
								
								values.append(-1.*r1/density_Etarget[i1]/scale_i2*extras)
								indices1.append(dof_index01)
								indices2.append(n)
								
								values.append(1.*r2/density_Etarget[i1]/scale_i2*extras)
								indices1.append(dof_index11)
								indices2.append(n)
								
								rh = (random.random() - 0.5) * delta * 1e-5
								hreg.append(rh)
								n += 1
								
				if 0:
					for i1 in range(nI1nodes-1):
						for i2 in range(nI2nodes):
							extras = 1 # / ((i2 + 0.5) / (nI2nodes))
							#i1 = 0
							#i2 = 0
							dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
							dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
							dof_index20 = self.dfgrid.dof_index(i1+1, i2, 1)
							dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
							dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
							dof_index21 = self.dfgrid.dof_index(i1+1, i2, 3)
							if (i1 == 0) and (i2 == 0):
								print dof_index00, dof_index10, dof_index20, dof_index01, dof_index11, dof_index21
							
							if 1:
								#print dof_index00, dof_index10 
								r1 = 1 + (random.random() - 0.5) * rscale
								r2 = 1 + (random.random() - 0.5) * rscale
								r3 = 1 + (random.random() - 0.5) * rscale
								
								values.append(-1.*r1/density_Etarget[i1]/scale_i1*extras)
								indices1.append(dof_index00)
								indices2.append(n)
								
								values.append(2.*r2/((density_Etarget[i1]+density_Etarget[i1+1])/2)/scale_i1*extras)
								#values.append(2.*r2/density_Etarget[i1]/scale_i1)
								indices1.append(dof_index10)
								indices2.append(n)
								
								values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
								indices1.append(dof_index20)
								indices2.append(n)
								
								rh = (random.random() - 0.5) * delta * 1e-5
								hreg.append(rh)
								n += 1
								
								
								values.append(-1.*r1/density_Etarget[i1]/scale_i1*extras)
								indices1.append(dof_index01)
								indices2.append(n)
								
								#values.append(2.*r2/density_Etarget[i1]/scale_i1)
								values.append(2.*r2/((density_Etarget[i1]+density_Etarget[i1+1])/2)/scale_i1*extras)
								indices1.append(dof_index11)
								indices2.append(n)
								
								values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
								indices1.append(dof_index21)
								indices2.append(n)
								
								rh = (random.random() - 0.5) * delta * 1e-5
								hreg.append(rh)
								n += 1
								
					if 0:		
						for i1 in range(nI1nodes):
							for i2 in range(nI2nodes-1):
								extras = 1 # / ((i2 + 0.5) / (nI2nodes-1))
								#i1 = 0
								#i2 = 0
								dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
								dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
								dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
								dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
								dof_index02 = self.dfgrid.dof_index(i1, i2+1, 2)
								dof_index12 = self.dfgrid.dof_index(i1, i2+1, 3)
								if (i1 == 0) and (i2 == 0):
									print dof_index00, dof_index10, dof_index01, dof_index11, dof_index02, dof_index12
								if 1:
									#print dof_index00, dof_index10 
									r1 = 1 + (random.random() - 0.5) * rscale
									r2 = 1 + (random.random() - 0.5) * rscale
									r3 = 1 + (random.random() - 0.5) * rscale
									
									if i1 > 0:
										dl = (density_Etarget[i1] + density_Etarget[i1-1])/2
									else:
										dl = density_Etarget[i1]
									if i1 < nI1nodes-1:
										#dr = (density_Etarget[i1] + density_Etarget[i1+1])/2
										dr = density_Etarget[i1+1]
									else:
										dr = density_Etarget[i1]
									values.append(-1.*r1/dl/scale_i2*extras)
									indices1.append(dof_index00)
									indices2.append(n)
									
									values.append(2.*r2/dl/scale_i2*extras)
									indices1.append(dof_index01)
									indices2.append(n)
									
									values.append(-1.*r3/dl/scale_i2*extras)
									indices1.append(dof_index02)
									indices2.append(n)
									
									rh = (random.random() - 0.5) * delta * 1e-5
									hreg.append(rh)
									n += 1
									
									
									values.append(-1.*r1/dr/scale_i2*extras)
									indices1.append(dof_index10)
									indices2.append(n)
									
									values.append(2.*r2/dr/scale_i2*extras)
									indices1.append(dof_index11)
									indices2.append(n)
									
									values.append(-1.*r3/dr/scale_i2*extras)
									indices1.append(dof_index12)
									indices2.append(n)
									
									rh = (random.random() - 0.5) * delta * 1e-5
									hreg.append(rh)
									n += 1
									
						
			elif 0:
				for i1 in range(0,nI1nodes-2):
					for i2 in range(nI2nodes):
						dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
						dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
						dof_index20 = self.dfgrid.dof_index(i1+1, i2, 1)
						dof_index30 = self.dfgrid.dof_index(i1+2, i2, 1)
						dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
						dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
						dof_index21 = self.dfgrid.dof_index(i1+1, i2, 3)
						dof_index31 = self.dfgrid.dof_index(i1+2, i2, 3)
						if (i1 == 0) and (i2 == 0):
							print
							print dof_index00, dof_index10, dof_index20, dof_index30
							print dof_index01, dof_index11, dof_index21, dof_index31
							print
						
						#o1 = self.dfgrid.index_to_orbitnr(i1-1,i2, 0) 
						#o2 = self.dfgrid.index_to_orbitnr(i1+0,i2, 0)
						#o3 = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
						# minimize -x_{i-1,j} + 2 x_{i,j} -x_{i+1,j} 
						r1 = 1 + (random.random() - 0.5) * rscale
						r2 = 1 + (random.random() - 0.5) * rscale
						r3 = 1 + (random.random() - 0.5) * rscale
						values.append(-0.25*r1/density_Etarget[i1]/scale_i1)
						indices1.append(dof_index00)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i1)
						indices1.append(dof_index01)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i1)
						indices1.append(dof_index11)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i1)
						indices1.append(dof_index11)
						indices2.append(n)
						
						values.append(0.5*r2/density_Etarget[i1+1]/scale_i1)
						indices1.append(dof_index10)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1+1]/scale_i1)
						indices1.append(dof_index11)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1+1]/scale_i1)
						indices1.append(dof_index20)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1+1]/scale_i1)
						indices1.append(dof_index21)
						indices2.append(n)
						
						values.append(-0.25*r3/density_Etarget[i1+2]/scale_i1)
						indices1.append(dof_index20)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1+2]/scale_i1)
						indices1.append(dof_index21)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1+2]/scale_i1)
						indices1.append(dof_index30)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1+2]/scale_i1)
						indices1.append(dof_index31)
						indices2.append(n)
						
						rh = (random.random() - 0.5) * delta * 1e-5
						hreg.append(rh)
						n += 1
				for i1 in range(0,nI1nodes):
					for i2 in range(0,nI2nodes-2):
						dof_index00 = self.dfgrid.dof_index(i1, i2, 0)
						dof_index10 = self.dfgrid.dof_index(i1, i2, 1)
						dof_index01 = self.dfgrid.dof_index(i1, i2, 2)
						dof_index11 = self.dfgrid.dof_index(i1, i2, 3)
						dof_index02 = self.dfgrid.dof_index(i1, i2+1, 2)
						dof_index12 = self.dfgrid.dof_index(i1, i2+1, 3)
						dof_index03 = self.dfgrid.dof_index(i1, i2+2, 2)
						dof_index13 = self.dfgrid.dof_index(i1, i2+2, 3)
						if (i1 == 0) and (i2 == 0):
							print
							print dof_index00, dof_index01, dof_index02, dof_index03
							print dof_index10, dof_index11, dof_index12, dof_index13
							print
						#o1 = self.dfgrid.index_to_orbitnr(i1,i2-1, 0) 
						#o2 = self.dfgrid.index_to_orbitnr(i1,i2, 0)
						#o3 = self.dfgrid.index_to_orbitnr(i1,i2+1, 0)
						# minimize -x_{i,j-1} + 2 x_{i,j} -x_{i,j+1} 
						r1 = 1 + (random.random() - 0.5) * rscale
						r2 = 1 + (random.random() - 0.5) * rscale
						r3 = 1 + (random.random() - 0.5) * rscale
						values.append(-0.25*r1/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index00)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index10)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index01)
						indices2.append(n)
						values.append(-0.25*r1/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index11)
						indices2.append(n)
						
						values.append(0.5*r2/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index01)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index11)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index02)
						indices2.append(n)
						values.append(0.5*r2/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index12)
						indices2.append(n)
						
						values.append(-0.25*r3/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index02)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index12)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index03)
						indices2.append(n)
						values.append(-0.25*r3/density_Etarget[i1]/scale_i2)
						indices1.append(dof_index13)
						indices2.append(n)
						
						rh = (random.random() - 0.5) * delta * 1e-5
						hreg.append(rh)
						n += 1
			else:"""
			
