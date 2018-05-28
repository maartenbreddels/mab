# -*- coding: utf-8 -*-
import os
import mab.gd.logging as logging
import emcee
from numpy import *
import numpy
import scipy
from mab.gd import gdfast_schw
from kaplot import *

logger = logging.getLogger("gd.schw.solution2")

class Dummy(object):
	pass

dummy = Dummy()

def lnprob(u):
	x = exp(u)
	x /= sum(x)
	logp =  sum([k.logp(x) for k in dummy.opts])
	#print logp
	return logp

def domcmc(x, opts):
	dummy.opts = opts
	N = len(x)
	x0 = x
	ndim = N
	nwalkers = 2*ndim
	grad = zeros(N)
	for opt in opts:
		opt.dlogpdx(x, grad)
	def gen():
		x = array([x0[i]*(1+1e-8*(random.random()*2-1)) for i in range(ndim)])
		x /= sum(x)
		return x
	p0 = [log(gen()) for i in xrange(nwalkers)]
	dummy.opts = opts
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[], threads=50)
	result = sampler.run_mcmc(p0, 1000)
	print sampler.flatchain.shape
	x = sampler.flatchain
	print "x shape", x.shape
	print "acceptance_fraction", sampler.acceptance_fraction
	logprob = array(sampler.lnprobability.flat)
	mozaic(2,2,box)
	#m = x.mean(axis=0).reshape(20, 8)
	#s = x.std(axis=0).reshape(20, 8)
	m = x.mean(axis=0).reshape(8, 20)
	s = x.std(axis=0).reshape(8, 20)
	#mozaic(2,2,box)
	select(0, 0)
	indexedimage(x0.reshape(8,20))
	select(0, 1)
	indexedimage(m)
	select(1, 0)
	indexedimage(s)
	select(1,1)
	histogram(logprob, bincount=100)
	l = lnprob(log(x0))
	vline(l, color="red")
	xlim(l-40, l+1)
	draw()
	import pdb
	pdb.set_trace()
	dsa
	

class Discrete(object):
	def __init__(self, modelpath, light_model, aperture_light, profile_model, schwsetname, schwmodelname, storage_2d_m0, storage_2d_m2, storage_2d_m4, storage_3d, storage_2d_losvd,  fitdensity2d, fitdensity3d, observation, binned_data_m2, binned_data_m4, dfgrid, max_iterations=1000, regularization=None, postfix=""):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d_m0 = storage_2d_m0
		self.storage_2d_m2 = storage_2d_m2
		self.storage_2d_m4 = storage_2d_m4
		self.storage_3d = storage_3d
		self.storage_2d_losvd = storage_2d_losvd
		self.aperture_light = aperture_light
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		#self.storage_2d_binned = storage_2d_binned
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.observation = observation
		self.dfgrid = dfgrid
		self.max_iterations = max_iterations
		self.regularization = regularization
		#self.regularization_delta = regularization_delta
		#self.use_jeans = use_jeans
		#self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.likelihood")
		self.postfix = postfix
		
	def run(self, args, opts, scope):
		self.init()
		self.solve(scope)
		
	def init(self):
		self.observation.load()
		self.storage_2d_m0.init()
		self.storage_2d_m0.load()
		self.storage_2d_m2.init()
		self.storage_2d_m2.load()
		self.storage_2d_m4.init()
		self.storage_2d_m4.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		#self.storage_2d_binned.init()
		#self.storage_2d_binned.load()
		self.storage_3d.init()
		self.storage_3d.load()
		self.storage_2d_losvd.load()
		#self.storage_2d.aperture.load()
		self.aperture_light.load()
	
	def solve(self, scope):
		stars = self.observation.stars
		self.logger.info("using %d stars/observations" % len(stars))
		#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		stars_inrange = stars.filter(lambda star: self.storage_2d_losvd.aperture.inrange(star.xi, star.eta))
		self.logger.info("stars in aperture range     : %d" % len(stars_inrange))
		self.logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		vmax = self.storage_2d_losvd.vmax
		#print "vmax", vmax
		delta_v = 2*vmax/self.storage_2d_losvd.Nv
		#print "res", delta_v
		
		for star in stars:
			star.aperture_index = self.storage_2d_losvd.aperture.findindex(star.xi, star.eta)
		
	
		losvds = self.storage_2d_losvd.losvds
		stars_invrange = stars.filter(lambda star: abs(star.vlos) < vmax)
		self.logger.info("stars in velocity range     : %d" % len(stars_invrange))
		self.logger.info("stars outside velocity range: %d" % (len(stars)-len(stars_invrange)))
		stars = stars_invrange
		
		sigma_v = 2.01
		numpy.random.seed(8)
		for star in stars:
			star.vlos = star.vlos_true# + numpy.random.normal(0, sigma_v)
			#star.vlos = star.vlos_true + numpy.random.normal(0, sigma_v)
			#print star.vlos, vmax, self.storage_2d_losvd.Nv
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d_losvd.Nv);
			outlier = True
			for losvd in losvds:
				if losvd[star.v_index, star.aperture_index] != 0:
					outlier = False
					break
			star.is_outlier = outlier
		stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
		
		
		self.logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		self.logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		
		
		
		
		Rborders = arange(self.storage_2d_losvd.NR+1) / (0.0+self.storage_2d_losvd.NR) * (self.storage_2d_losvd.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		#print Rborders
		self.rho2d_target = array([self.light_model.cumdensityR(R1, R2, M=1.) for R1, R2 in zip(R1s, R2s)])
		
		
		rho2ds = sum(losvds, axis=1)
		rho2dmatrix = sum(losvds, axis=1)
		rho3dmatrix = self.storage_3d.moments3d[:,0,:]
		#rho2dmatrix = self.storage_2d_m0.moments[:,0,:]
		r1s = self.storage_3d.rborders[:-1]
		r2s = self.storage_3d.rborders[1:]
		delta_r = r2s[0] - r1s[0]
		#R1s = self.storage_2d.rborders[:-1]
		#R2s = self.storage_2d.rborders[1:]
		self.rho3d_target = array([self.light_model.cumdensityr(r1, r2, M=1.) for r1, r2 in zip(r1s, r2s)])
		
		#for i in range(losvds.shape[0]):
			#for j in range(losvds.shape[1]):
		#	print i, sum(losvds[i]),
		for i in range(losvds.shape[0]):
			#print sum(losvds[i])
			for j in range(losvds.shape[2]):
				#dens = sum(losvds[i,:,j])
				#if dens > 0:
				#print self.rho2d_target.shape
				#print losvds.shape 
				#losvds[i,j,:] /= self.rho2d_target
				#losvds[i,:,j] = scipy.ndimage.gaussian_filter(losvds[i,:,j], sigma_v/delta_v, mode='constant')
				pass
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [sigma_v/delta_v, 0.51])
			losvds[i] /= (delta_v * delta_R)
		#print 
		#print losvds.shape, delta_v, delta_R, Rborders[0], Rborders[-1]
		#print Rborders
		#for i in range(losvds.shape[0]):
			#for j in range(losvds.shape[1]):
		#	print i, sum(losvds[i]*delta_v*delta_R),
			
		v_indices = [star.v_index for star in stars]
		aperture_indices = [star.aperture_index for star in stars]
		
		#print losvds.shape
		pmatrix = array(list((losvds/(self.rho2d_target/delta_R))[:,v_indices, aperture_indices]))
		pmatrix = array(list((losvds)[:,v_indices, aperture_indices]))
		pmatrix = pmatrix * 1.
		#print pmatrix.shape
		
		rho2d_error = self.rho2d_target.max() * 0.000001*0.5 * 0.1
		error_x = 1e-3
		if 1:
			filename = os.path.join(self.modelpath, "df/orbitweights_tang.npy")
			orbitweights = load(filename)
			c = orbitweights.flatten() 
			c /= sum(c)
			x = c
			xtrue = x
		if 0:
			self.x0 = x
			self.true_losvd = numpy.tensordot(self.storage_2d_losvd.losvds, c, axes=[(0,),(0,)])
			self.true_rho2d = numpy.tensordot(self.storage_2d_losvd.masses, c, axes=[(0,),(0,)])
			#self.true_rho2d = numpy.tensordot(rho2ds, c, axes=[(0,),(0,)])
			self.true_rho3d = numpy.tensordot(rho3dmatrix, c, axes=[(0,),(0,)])
			filename = os.path.join(self.modelpath, "df/losvd_tang.npy")
			save(filename, self.true_losvd)
			filename = os.path.join(self.modelpath, "df/masses_tang.npy")
			save(filename, self.true_rho2d)
			dsa
			if 0:
				#graph(self.true_rho3d)
				#graph(self.rho3d_target, color="red")
				#avg = 
				graph((self.true_rho3d-self.rho3d_target)/self.rho3d_target.max(), color="red")
				draw()
			#import pdb; pdb.set_trace()
		else:
			filename = os.path.join(self.modelpath, "df/losvd_tang.npy")
			self.true_losvd = load(filename)
			filename = os.path.join(self.modelpath, "df/masses_tang.npy")
			self.true_rho2d = load(filename)
			self.true_losvd = scipy.ndimage.gaussian_filter(self.true_losvd, [sigma_v/delta_v, 0.51])
			self.true_losvd /= (delta_v * delta_R)
		debug = False
		if 0:
			filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
			x = numpy.load(filename)
			logger.info("loading orbitweights %s" % filename)
		else:
			x = x*0 + 1.
			x = random.random(len(x))
			x /= sum(x)
			u = log(x)
			
		
			
		#print rho2dmatrix.shape
		rho3dmatrix = rho3dmatrix * 1
		#rhoerror = maximum(self.rho3d_target*rho2d_error, self.rho3d_target.max() * 0.001)
		rhoerror = self.rho3d_target*rho2d_error
		s = self.rho3d_target/self.rho3d_target.max()
		rhoerror = self.rho3d_target.max() * 0.05 * maximum(0.1, s)  + self.rho3d_target * 0
		rhoerror = self.rho3d_target.max() * 0.07 * maximum(0.1/7, s)  + self.rho3d_target * 0
		rhoerror = maximum(self.rho3d_target.max() * 0.01*1.5, self.rho3d_target * 0.01*2)  + self.rho3d_target * 0
		#rhoerror = maximum(rhoerror*1e-4, rho2d_error)
		#self.opt = gdfast_schw.OptimizationProblemSchw(pmatrix, rho3dmatrix, x, self.rho3d_target, rhoerror, error_x, True, False, True)
		fit_mass_3d = False
		#fit_mass_3d = True
		if fit_mass_3d:
			mass_matrix = rho3dmatrix
			mass_target = self.rho3d_target
			mass_error = rhoerror
		else:
			rhoerror = maximum(self.rho2d_target.max() * 0.01*0.05, self.rho2d_target * 0.001)#  + self.rho3d_target * 0
			mass_matrix = rho2dmatrix
			mass_target = self.rho2d_target
			mass_error = rhoerror
		entropy_scale = 1e-20
		self.opt = gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, entropy_scale, True, True, True)
		#print "true L?", self.opt.likelihood(log(xtrue))
		self.opt_kin = gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, entropy_scale, True, False, False)
		self.opts = [
			gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, 0, True, False, False),
			gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, 0, False, True, False),
			gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, 0, False, False, True),
			gdfast_schw.OptimizationProblemSchw(pmatrix, mass_matrix, x, mass_target, mass_error, error_x, entropy_scale, False, False, False),
			]
		debug = False
		#debug = True
		if 1:
			#x = numpy.load("xlast.npy")
			#print dir(self.light_model)
			N = 250000
			light_profile = self.light_model.light_profile
			rs = light_profile.sample_r(N=N, rmax=100.)
			costheta = numpy.random.random(N) * 2 - 1
			phi = numpy.random.random(N) * 2 * pi
			eta = numpy.random.random(N) * 2 * pi
			theta = numpy.arccos(costheta)
			#sintheta = numpy.sqrt(1-costheta**2)
			sintheta = numpy.sin(theta)
			#print r.shape, sintheta.shape, phi.shape, len(dt) 
			xp = x
			x = rs * sintheta * numpy.cos(phi)
			y = rs * sintheta * numpy.sin(phi)
			Rs = sqrt(x**2+y**2)
			x = xp
			#ps = self.
			Rs = Rs[Rs<1.5]
			rs = rs[rs<1.5]
			#rs = rs[rs>0.1]
			#normal = scipy.integrate.quad(lambda R: light_profile.densityR(R,M=1.)*2*pi*R, 0, 1.5)[0]
			normal = scipy.integrate.quad(lambda r: light_profile.densityr(r,M=1.)*4*pi*r**2, 0, 1.5)[0]
			if debug:
				print "normal", normal
			#normal = 1.
			if fit_mass_3d:
				ps = [log(light_profile.densityr(r,M=1.)*4*pi*r**2/normal) for r in rs]
			else:
				ps = [log(light_profile.densityR(R,M=1.)*2*pi*R/normal) for R in Rs]
			N = len(ps)
			if debug:
				print N
				print "tot p", sum(ps)
				print "mean p", mean(ps)
				print rho3dmatrix.shape
			if fit_mass_3d:
				mass_indices = [int(r/1.5*100) for r in rs]
				mass_matrix = rho3dmatrix[:,mass_indices] * 1. / delta_r
				mass_matrixN = rho3dmatrix * 1./delta_r
				totalmass_matrix = sum(rho3dmatrix, axis=1)
				ptotalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)
				counts, bins = numpy.histogram(rs, 100, [0, 1.5], new=True)
			else:
				mass_indices = [int(R/1.5*30) for R in Rs]
				mass_matrix = self.storage_2d_losvd.masses[:,mass_indices] * 1. / delta_R
				mass_matrixN = self.storage_2d_losvd.masses * 1. / delta_R
				totalmass_matrix = ptotalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)
				counts, bins = numpy.histogram(Rs, 30, [0, 1.5])
				counts /= sum(counts)
				counts = counts / 2000
				if debug:
					print "2d, delta_R", delta_R
			#mass = dot(self.storage_2d_losvd.masses.T, xtrue)
			if debug:
				print "total 3d", sum(dot(rho3dmatrix.T, xtrue))
				print "total 2d", sum(dot(self.storage_2d_losvd.masses.T, xtrue))
				print "normal check", dot(xtrue, totalmass_matrix)
			opt_matrix_mass = gdfast_schw.OptimizationMatrix(mass_matrix, totalmass_matrix)
			counts = array(counts).astype(float64)# * 1.
			#print counts, sum(counts)
			rho3dmatrix = rho3dmatrix * 1.
			#print "-->", rho3dmatrix.shape, counts.shape
			#print rho3dmatrix.dtype, counts.dtype
			counts = self.true_rho2d
			counts /= sum(counts)
			counts *= 200000
			opt_matrix_massN = gdfast_schw.OptimizationMatrixN(mass_matrixN, counts, totalmass_matrix)
			if debug:
				print "logp", opt_matrix_mass.logp(xtrue), opt_matrix_mass.logp(x)
				print "logp", opt_matrix_massN.logp(xtrue), opt_matrix_massN.logp(x)
				print opt_matrix_mass.logp(xtrue)/N
				print sum(self.rho3d_target)
			#print (self.rho3d_target-mass)/self.rho3d_target
			#print "x =", x, sum(x)
			if 0:
				box()
				#mask = (rho3dmatrix/delta_r) > 1
				#rho3dmatrix[mask] = 1000
				I = rho3dmatrix * 1.
				I /= I.max()
				I = log10(I)
				I[I<-6] = -6
				indexedimage(I)
				draw()
			diff = log(dot(x, mass_matrix))-log(dot(xtrue, mass_matrix))
			indices = argsort(diff)

			#import pdb
			#pdb.set_trace()
			#sysdsa
			
		if 1:
			#x = xtrue
			#print "sum(x)", sum(x)
			opt_matrix_kin = gdfast_schw.OptimizationMatrix(pmatrix, ptotalmass_matrix)
			counts_kin, binsx, biny = numpy.histogram2d([star.v_index for star in stars], [star.aperture_index for star in stars], bins=[30,30], range=[(0,30),(0, 30)])
			if 0:
				counts_kin2 = numpy.zeros((30, 30))
				for star in stars:
					counts_kin2[star.v_index, star.aperture_index] += 1
				mozaic(2,2,box)
				indexedimage(counts_kin)
				select(0,1)
				indexedimage(counts_kin2)
				draw()
			
			mask = counts_kin > 0
			counts_kin = counts_kin[mask]
			counts_kin = counts_kin * 1.
			pmatrixN = losvds[:,mask]
			#debug = True
			if 1:
				pmatrixN = losvds * 1.0
				pmatrixN = pmatrixN.reshape((pmatrixN.shape[0], -1)) * 1.
				counts_kin = self.true_losvd * 1.
				counts_kin = counts_kin.reshape(-1) * 1.
				counts_kin /= sum(counts_kin)
				counts_kin *= 2000
			#@import pdb
			#pdb.set_trace()
			if debug:
				print "%d versus %d speedup: %f" % (pmatrixN.shape[1], len(stars), len(stars)*1./pmatrixN.shape[1]) 
			#pmatrixN = array(list((losvds/(self.rho2d_target/delta_R))[:,v_indices, aperture_indices]))
			#pmatrix = array(list((losvds)[:,v_indices, aperture_indices]))
			pmatrixN = pmatrixN * 1.0
			#print sum(counts_kin == 0)
			opt_matrix_kinN = gdfast_schw.OptimizationMatrixN(pmatrixN, counts_kin, ptotalmass_matrix)
			for k in [pmatrixN, counts_kin, ptotalmass_matrix]:
				print k.min(), k.max(), k.sum(), k.std()
			dsa
			opt_norm = gdfast_schw.OptimizationNormalize(1.-.000001, 0.001)
			opt_entropy = gdfast_schw.OptimizationEntropy(1.e-2)
			opt = self.opt_kin
			u = log(x)
			if debug:
				for i in [0, 1]:
					x1 = x * 1.
					x2 = x * 1.
					dx = 1e-8
					x2[i] += dx
					grad = x * 0
					for opt in [opt_matrix_kin, opt_matrix_kinN, opt_entropy]:#, opt_matrix_mass, opt_matrix_massN, opt_norm]:
						grad = x * 0
						print u.shape, grad.shape
						opt.dlogpdx(x, grad)
						#sys.exit(0)
						w1 = opt.logp(x1)
						w2 = opt.logp(x2)
						print "w", opt.logp(x)
						print "grad", grad[i]
						print "grad", (w2-w1)/dx#/x[i]
						print
					print
						#opt_matrix_kin.dlogpdx(x, grad)
						#grad *= x
						#print "grad3", grad[i]
				#print "logp", opt.likelihood(u)
				#print "logp", opt_matrix_kin.logp(x)
				print
				sys.exit(0)
		#x = x * 0 + 1
		#x /= sum(x)
		global calls
		calls = 0
		debug =  False
		debug=True
		opts = [opt_matrix_kinN, opt_matrix_massN, opt_norm]
		def f_and_g(x):
			global calls
			#if calls > 10:
			#	numpy.save("xlast.npy", x)
			#	dsa
			#calls += 1
			if 1:
				grad = x * 0
				logp = opt_matrix_kinN.logp(x)*1 +\
					opt_matrix_massN.logp(x)*0 +\
					opt_norm.logp(x) * 1
					#+\
					#opt_entropy.logp(x) * 1.
				opt_matrix_kinN.dlogpdx(x, grad)
				#opt_matrix_massN.dlogpdx(x, grad)
				opt_norm.dlogpdx(x, grad)
				#opt_entropy.dlogpdx(x, grad)
				if debug:
					print "%10f %10f %10f %10f %10f" % (logp,  sum(x), dot(totalmass_matrix, x/sum(x)), dot(ptotalmass_matrix, x/sum(x)), dot(losvds.T, x).sum() * delta_R * delta_v / dot(ptotalmass_matrix, x/sum(x)))
					#print 
				if 0:
					print ".", sum(x), dot(totalmass_matrix, x/sum(x))
					print logp
					for i in [0, 10, 100]:
						x1 = x * 1.#/sum(x)
						x2 = x * 1.#/sum(x)
						dx = 1e-7
						x2[i] += dx
						#for opt in [opt_matrix_kin, opt_matrix_massN, opt_norm]:
						for opt in [opt_matrix_massN]:
							grad = x * 0
							opt.dlogpdx(x, grad)
							w1 = opt.logp(x1)
							w2 = opt.logp(x2)
							print "grad", grad[i]
							print "grad man", (w2-w1)/dx#/x[i]
							print
				return -logp, -grad
			
			
			u = log(x)
			#print u
			w = -self.opt.likelihood(u)
			grad = u * 0
			self.opt.dfdx(u, grad)
			#print w
			if 0:
				print w
				for i in [0, 1]:
					u1 = u * 1.
					u2 = u * 1.
					du = 1e-4
					u2[i] += du
					for opt in self.opts:
						grad = u * 0
						opt.dfdx(u, grad)
						w1 = -opt.likelihood(u1)
						w2 = -opt.likelihood(u2)
						print "grad", grad[i]
						print "grad man", (w2-w1)/du#/x[i]
						print
			if 0:
				mass = numpy.tensordot(rho3dmatrix, x, axes=[(0,),(0,)])
				mass_true = numpy.tensordot(rho3dmatrix, xtrue, axes=[(0,),(0,)])
				print sum(((mass-self.rho3d_target)**2/rhoerror**2)), sum(((mass_true-self.rho3d_target)**2/rhoerror**2))
			return w, grad/x
		if 1:
			#f_and_g(x)
			#dsa
			#u = scipy.optimize.fmin_l_bfgs_b(f_and_g, log(x+1e-9), None, iprint=1,factr=1e3,maxfun=1000)[0]
			bounds = [(1e-10, 1) for i in range(len(x))]
			#x = xtrue
			#bounds = None
			approx_grad = False
			if 0: #scope["options"].mcmctest:
				filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
				x = numpy.load(filename)
				if 0:
					cdf = cumsum(dot(mass_matrixN.T, x))
					cdf /= cdf.max()
					values = cdf[mass_indices]
					borders = array([0] + list(cdf))
					counts, _ = numpy.histogram(values, borders)
					dcdf = borders[1:] - borders[:-1]
					density = counts/dcdf
					box()
					#histogram(values, binwidth=0.01)
					histogramline(borders, density, binwidth=0.01)
					draw()
					import pdb;
					pdb.set_trace()
				else:
					
					domcmc(x, opts)
					
					hessian = zeros((N, N))
					for i in range(N):
						grad1 = zeros(N)
						grad2 = zeros(N)
						dx = 1e-8
						x1 = x * 1.
						x2 = x * 1.
						x2[i] += dx
						for opt in opts:
							opt.dlogpdx(x1, grad1)
							opt.dlogpdx(x2, grad2)
						if i == 0:
							print grad1
							print grad2
							print grad2-grad1
						hessian[i] = (grad2-grad1)/dx
					print hessian[:4,:4]
					
					u, s, vd = numpy.linalg.svd(hessian)
					self.u = u
					self.vd = vd
					self.s = s
					print s
					S = zeros((len(s), len(s)))
					for i in range(50):
						S[-1-i,-1-i] = 1/s[-1-i]
						#S[i,i] = 1/s[i]
					#print S
					#self.hessian2 = dot(u, dot(S, vd))
					#self.cov2 = numpy.matrix(self.hessian2).I
					self.cov2 = numpy.matrix(vd).T * S * numpy.matrix(u).T
					mozaic(2,2,box)
					select(0, 0)
					indexedimage(x.reshape(8,20))
					select(1, 0)
					indexedimage(self.cov2)
					draw()
					import pdb
					pdb.set_trace()
			else:
				x = ones(len(x))
				x /= sum(x)
				x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=approx_grad, iprint=-1,factr=1e-2,maxfun=200000)[0]
			#x /= sum(x)
			u = log(x)
			if debug:
				print "res:", [opt.logp(x) for opt in [opt_matrix_kinN, opt_matrix_massN, opt_norm]]
			
		elif 1:
			#def myfunc(x, grad):
			def myfunc(x, grad):
				#print "%%%"
				#print x
				#print sum(x)
				u = log(x)
				#print u
				#x = exp(u)
				if grad.size > 0:
					self.opt.dfdx(u, grad)
					#grad[:] = grad#/x
					grad /= x
				w = -self.opt.likelihood(u)
				
				if 0:
					print w
					print u
					for i in [0, 1]:
						u1 = u * 1.
						u2 = u * 1.
						du = 1e-4
						u2[i] += du
						w1 = -self.opt.likelihood(u1)
						w2 = -self.opt.likelihood(u2)
						print "grad", grad[i]
						print "grad man", (w2-w1)/du#/x[i]
					
				#print sum(x), w
				return w

			epsilon = 0.01
			def normconstraint1(x, grad, *a):
				#print "test"
				if grad.size > 0:
					#print "error, no grad possible"
					grad[:] = 1.
				#print ">>>>>>>", x
				return -(sum(x) - 1.+epsilon)
			def normconstraint2(x, grad, *a):
				#print "test"
				if grad.size > 0:
					#print "error, no grad possible"
					grad[:] = 1.
				#print ">>>>>>>", x
				return -(sum(x) - 1.-epsilon)

			import nlopt
			dim = x.shape[0]
			opt = nlopt.opt(nlopt.LD_MMA, dim)
			opt = nlopt.opt(nlopt.LD_LBFGS, dim)
			opt = nlopt.opt(nlopt.LD_TNEWTON_PRECOND_RESTART, dim)
			
			#print dir(nlopt)
			#opt = nlopt.opt(nlopt.LN_COBYLA, dim)
			
			opt.set_lower_bounds([1e-10] * dim)
			#opt.set_lower_bounds([0, float('inf')] * dim)
			opt.set_min_objective(myfunc)
			#opt.add_inequality_constraint(normconstraint, 1e-5)
			#opt.add_inequality_constraint(lambda x,grad: normconstraint1(x,grad), 1e-2)
			#opt.add_equality_constraint(lambda x,grad: normconstraint1(x,grad), 1e-3)
			#opt.set_xtol_rel(1e-4)
			#opt.set_xtol_rel(1e-8)
			#opt.set_ftol_abs(1e-15)
			opt.set_ftol_abs(1e-6)
			opt.set_maxtime(60*20)
			try:
				#x = opt.optimize(x)
				x = opt.optimize(x)
			except:
				logger.error("oops, opt failed, result code = %d " % opt.last_optimize_result())
				#raise
			#x = opt.last_optimize_result()
			#import pdb; pdb.set_trace()
			#print x
			logger.debug("sum orbit weight: %f " % sum(x))
			u = log(x)
			#x = exp(u)
			minf = opt.last_optimum_value()
			#print "optimum at ", x
			#print "minimum value = ", minf
		
		else:
			#x = x * 0 + 1
			#x /= sum(x)
			u = log(x)
			#print u
			self.opt.optimize(10000, 10, u)
			x = exp(u)
			x /= sum(x)
			u = log(x)
		#numpy.save(filename, x)
		#print u
		logL = self.opt.likelihood(u)
		logL = -f_and_g(x)[0]
		#logL = self.opt_kin.likelihood(u)
		refLogL = -29244.4240491
		#refLogL = -7392.92252529
		refLogL = -7392.92252529
		refLogL = -7406.08983537-22200.0538537 #-927.128433428
		refLogL = -35247.3864887
		refLogL = -35161.7892221
		refLogL = -7495.98083661
		refLogL= -6.06900E+03
		refLogL = -5552.10049844+24210.9986633+80
		
		if isinf(logL):
			logL = -100000
		if isnan(logL):
			logL = -100000
		if debug:
			print "log L", logL, logL-refLogL, "sum(x)", sum(x), [opt.likelihood(u) for opt in self.opts]
		
		logL -= refLogL
		if debug:
			if fit_mass_3d:
				mass = numpy.tensordot(rho3dmatrix, x, axes=[(0,),(0,)])
				mass_true = numpy.tensordot(rho3dmatrix, xtrue, axes=[(0,),(0,)])
				#print sum(((mass-self.rho3d_target)**2/rhoerror**2)), sum(((mass_true-self.rho3d_target)**2/rhoerror**2))
				if debug:
					print "mass", sum(mass), sum(mass_true)
					diff = mass_true - mass
					print diff
					print argsort(diff)
					diff /= rhoerror
					print diff
					print argsort(diff)
			else:
				mass = numpy.tensordot(rho2dmatrix, x, axes=[(0,),(0,)])
				mass_true = numpy.tensordot(rho2dmatrix, xtrue, axes=[(0,),(0,)])
				#print sum(((mass-self.rho3d_target)**2/rhoerror**2)), sum(((mass_true-self.rho3d_target)**2/rhoerror**2))
				if debug:
					print "mass", sum(mass), sum(mass_true), sum(self.rho2d_target)
					diff = mass_true - mass
					print diff
					print argsort(diff)
					diff /= rhoerror
					print diff
					print argsort(diff)
				
		orbitweights = x
		orbitweights /= sum(orbitweights)
		mozaic(2,2,box)
		x = x.reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		#s = x.std(axis=0).reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		select(0,0)
		indexedimage(x)
		
		
		#draw()
		
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		numpy.save(filename, orbitweights)
		logger.info("saving orbitweights %s" % filename)
		
		
		moments_m0 = tensordot(orbitweights, self.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
		if debug:
			print "sanity check mass", sum(moments_m0[0])
		#import pdb;
		#pdb.set_trace()
		moments_m2 = tensordot(orbitweights, self.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
		moments_m4 = tensordot(orbitweights, self.storage_2d_m4.projectedmoments, axes=[(0,), (0,)])
		moments3d = tensordot(orbitweights, self.storage_3d.moments3d, axes=[(0,), (0,)])
		moments_m0[1:] /= moments_m0[0]
		moments_m2[1:] /= self.binned_data_m2.moments[0]
		moments_m4[1:] /= self.binned_data_m4.moments[0]
		moments3d[1:] /= moments3d[0]
		self.moments2d_solution_m2 = moments_m2
		self.moments2d_solution_m4 = moments_m4
		self.moments3d_solution = moments3d
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m0.npy")
		numpy.save(filename, moments_m0)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		numpy.save(filename, moments_m2)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		numpy.save(filename, moments_m4)
		logger.debug("saving %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		numpy.save(filename, moments3d)
		logger.debug("saving %s" % filename)
		
		
		chisq = -logL*2
		chisqs = [(name, chisq)  for name in "m2 m4 reg kin".split()]
		
		f = file(os.path.join(self.dirname, "results/probability" +self.postfix +".txt"), "w")
		print >>f, "%e" % exp(-0.5*chisq)
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (-0.5*chisq)
		
		logps = [(name, -0.5*chisq) for name, chisq in chisqs] 
		f = file(os.path.join(self.dirname, "results/logprobability_seperate" +self.postfix +".txt"), "w")
		print >>f, repr(logps)
		
		
		
	def load(self):
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		self.orbitweights = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m0.npy")
		self.moments_solution_m0 = numpy.load(filename)
		logger.debug("loading %s"% filename)

		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		self.moments_solution_m2 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		self.moments_solution_m4 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		self.moments3d_solution = numpy.load(filename)
		logger.debug("loading %s" % filename)

		#return self.findsolution()
	def calculate_solution_moments(self, moments):
		return tensordot(self.orbitweights, moments, axes=[(0,), (0,)])
		
def lnprob_(x):
	#u = log(x)
	u = x
	#print u
	logL = opt.likelihood(u)
	print logL
	return logL
	#return 

	return random.random()
	
class DiscretePhotometry(Discrete):
	def __init__(self, photometry, **kwargs):
		self.photometry = photometry
		super(DiscretePhotometry, self).__init__(**kwargs)
		
	def init(self):
		super(DiscretePhotometry, self).init()
		self.photometry.load()
		
	def solve(self):
		#mass_matrix = self.storage_2d_losvd.masses[:,mass_indices] * 1. / delta_R
		Rborders = arange(self.storage_2d_losvd.NR+1) / (0.0+self.storage_2d_losvd.NR) * (self.storage_2d_losvd.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		
		
		mass_matrixN = self.storage_2d_losvd.masses * 1. / delta_R
		totalmass_matrix = ptotalmass_matrix = sum(self.storage_2d_losvd.masses, axis=2)
		counts = self.photometry.grid
		print mass_matrixN.shape, counts.shape, totalmass_matrix.shape
		shape = mass_matrixN.shape
		newshape = (shape[0] * shape[1], ) + shape[2:]
		print "shape:", shape, newshape
		mass_matrixN = mass_matrixN.reshape(newshape) * 1.
		
		shape = totalmass_matrix.shape
		newshape = (shape[0] * shape[1], ) + shape[2:]
		print "shape:", shape, newshape
		totalmass_matrix = totalmass_matrix.reshape(newshape) * 1.
		
		mask = counts == 0
		print sum(mask)
		#dsa
		mass_matrixN = mass_matrixN[:,~mask] * 1
		counts = counts[~mask] * 1.
		
		print mass_matrixN.shape, counts.shape, totalmass_matrix.shape
		
		opt_matrix_massN = gdfast_schw.OptimizationMatrixN(mass_matrixN, counts, totalmass_matrix)
		
		x = zeros(newshape[0]) + 1.
		x /= sum(x)
		
		def f_and_g(x):
			grad = x * 0
			logp = opt_matrix_massN.logp(x)*1
			opt_matrix_massN.dlogpdx(x, grad)
			if debug:
				print "%10f %10f %10f" % (logp,  sum(x), dot(totalmass_matrix, x/sum(x)))
			return -logp, -grad
			u = log(x)
			#print u
			w = -self.opt.likelihood(u)
			grad = u * 0
			self.opt.dfdx(u, grad)
			#print w
			if 0:
				print w
				for i in [0, 1]:
					u1 = u * 1.
					u2 = u * 1.
					du = 1e-4
					u2[i] += du
					for opt in self.opts:
						grad = u * 0
						opt.dfdx(u, grad)
						w1 = -opt.likelihood(u1)
						w2 = -opt.likelihood(u2)
						print "grad", grad[i]
						print "grad man", (w2-w1)/du#/x[i]
						print
		bounds = [(1e-10, 1) for i in range(len(x))]
		approx_grad = False
		x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=approx_grad, iprint=-1,factr=1e-2,maxfun=200000)[0]
		print x
		u = log(x)
		#if debug:
		#	print "res:", [opt.logp(x) for opt in [opt_matrix_kinN, opt_matrix_massN, opt_norm]]		
			
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		numpy.save(filename, x)
		logger.debug("saving %s" % filename)
		
	def load(self):
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		self.orbitweights = numpy.load(filename)
		logger.debug("loading %s" % filename)



class DiscreteMCMC(Discrete):
	#def __init__(self, solution):
	#	self.solution = solution
		
		
	def _run(self, args, opts, scope):
		self.solution.init()
		self.solution.load()

		ndim = self.storage_2d_losvd.losvds.shape[0]
		nwalkers = 3*ndim
		#ivar = 1./np.random.rand(ndim)
		def x0():
			x = array([random.random() for i in range(ndim)])
			x = self.x0 + x
			x /= sum(x)
			return x
		p0 = [log(x0()) for i in xrange(nwalkers)]

		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[], threads=50)
		result = sampler.run_mcmc(p0, 100)
		print sampler.flatchain.shape
		x = sampler.flatchain
		mozaic(2,2,box)
		m = x.mean(axis=0).reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		s = x.std(axis=0).reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		select(0,0)
		indexedimage(m.T)
		select(1,0)
		indexedimage(s.T)
		draw()
		#import pdb; pdb.set_trace()


