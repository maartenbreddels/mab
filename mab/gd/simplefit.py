# -*- coding: utf-8 -*-
import scipy.optimize
from numpy import *
import mab.gd.logging as logging
logger = logging.getLogger("gd.simplefit")
import emcee
from kaplot import *

#fitter = None

def lnprob(x, fitter):
	print x, fitter
	return fitter.logL(x)

class MCMCExplore(object):
	def __init__(self, simplefit):
		self.simplefit = simplefit
		
	def run(self, *args):
		self.simplefit.load()
		
		stars = self.simplefit.observation.load()
		
		self.re = stars.re
		self.vlos = stars.vlos_helio
		self.e_vlos = stars.e_vlos
		self.N = len(self.re)
		self.cats = stars.catalogue_mask.astype(int)
		#print re.max()
		#dsa
		x0 = self.simplefit.x0
		#print self.logL(*x0)
		
		print len(x0)
		ndim, nwalkers = len(x0), len(x0)*2
		sigmas = [0.05, 0.05, 0.1, 5, 1, 5, 1] 
		def gen(x0):
			return [x+random.normal(0, sigma) for x,sigma in zip(x0, sigmas)]
		p0 = [gen(x0) for k in range(nwalkers)]
		#print "p0", p0
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[self])
		#sampler.run_mcmc(p0, 1000)
		sampler.run_mcmc(p0, 500)
		
		#hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
		bs = 10**sampler.flatchain[:,2]
		#dens = 10**sampler.flatchain[:,0]
		print len(bs)
		
		page(fontsize="16pt")
		mozaic(1,1,box)
		histogram(bs, bincount=50, datamin=0., datamax=1.)
		print bs
		#, dens
		labels("b (kpc)", "N")
		vline(0.3, color="red")
		#vline(0.16, color="red")
		#vline(0.38775660972067766, color="red") # sxt
		#vline(0.7868, color="red") # fnx
		if 0:
			select(1, 0)
			histogram(dens, bincount=50)
			labels("mu0", "N")
			select(0, 1)
			scatter(bs, dens)
			labels("b", "mu0")
		draw()
		

	def logL(self, x):
		#print "X", x
		#x = x[0]
		#print x
		self.simplefit.set_parameters(x)
		logL = 0
		for i in range(self.N):
			#p = exp(self.galaxy_velocity_model.logL(vlos[i], e_vlos[i]))*w1 \
			#	+ w2 * exp(self.foreground_velocity_model.logL(vlos[i], e_vlos[i]))
			#print p
			#if p == 0: import pdb; pdb.set_trace()
			#logL += log(p )
			logL1 = -inf
			for j in range(2):
				if (1<<j) & (self.cats[i]):
					logL1 = max(logL1, self.simplefit.logL(self.re[i], self.vlos[i], self.e_vlos[i], j))
			#logL += self.simplefit.logL(self.re[i], self.vlos[i], self.e_vlos[i])
			logL += logL1
		#print "    ", logL, x
		return logL
		

class SimpleFit(object):
	def __init__(self, filename, light_profile, light_model, observation, galaxy_velocity_model, foreground_velocity_model, const_surface_density=1e-4, fit_light=False):
		self.light_profile = light_profile
		self.light_model = light_model
		self.observation = observation
		self.galaxy_velocity_model = galaxy_velocity_model
		self.foreground_velocity_model = foreground_velocity_model
		self.filename = filename
		self.const_surface_density = const_surface_density
		self.fit_light = fit_light
		#self.const_surface_density = self.light_profile.densityR(2., M=1.)
		
	def load(self):
		self.x0 = parameters = load(self.filename)
		logger.info("parameters: %s" % parameters)
		self.set_parameters(parameters)
		
		stars = self.observation.load()
		re = stars.re
		print re.max()
		R = r = re.max()
		if 0:
			#print self.const_surface_density * self.w2 * r**2
			R = self.light_model.arcsec_to_kpc(r)
			I, err = scipy.integrate.quad(lambda r: 2*pi*r*self.light_profile.densityR(r, M=1.), 0, R)
			print I, err
			#self.const_surface_density = 0.35
			re = stars.re
			vlos = stars.vlos_helio
			e_vlos = stars.e_vlos
			N = len(re)
			print re.max()
			logL = 0
			for i in range(N):
				#p = exp(self.galaxy_velocity_model.logL(vlos[i], e_vlos[i]))*w1 \
				#	+ w2 * exp(self.foreground_velocity_model.logL(vlos[i], e_vlos[i]))
				#print p
				#if p == 0: import pdb; pdb.set_trace()
				#logL += log(p )
				logL += self.logL(re[i], vlos[i], e_vlos[i])
			print "    ", logL
			
		#dsa
		
		
		
	def set_parameters(self, x):
		#R = x[0]
		self.const_surface_density = 10**(x[0])
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				self.light_profile.b = 10**(x[1])
			else:
				self.light_profile.scale = 10**(x[1])
			offset = 2
		else:
			offset = 1
		n_params_galaxy = len(self.galaxy_velocity_model.parameters)
		n_params_foreground = len(self.foreground_velocity_model.parameters)
		for i in range(n_params_galaxy):
			value = x[offset+i]
			self.galaxy_velocity_model.parameters[i].set(value)
			#print "set galaxy", self.galaxy_velocity_model.parameters[i].name, "to", value
		for i in range(n_params_foreground):
			value = x[offset+n_params_galaxy+i]
			self.foreground_velocity_model.parameters[i].set(value)
			#print "set foreground", self.foreground_velocity_model.parameters[i].name, "to", value
		#r = abs(R)
		#self.w1 = 1/(1+r)
		#self.w2 = r/(1+r)
		#assert abs((self.w2/self.w1)- r) < 1e-6
		
	def save(self):
		parameters = [log10(self.const_surface_density)]
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				parameters.append(log10(self.light_profile.b))
			else:
				parameters.append(log10(self.light_profile.scale))
		parameters += [k.get() for k in self.galaxy_velocity_model.parameters] + [k.get() for k in self.foreground_velocity_model.parameters]
		parameters = array(parameters)
		save(self.filename, parameters)
		
	def logL(self, Re, vlos, e_vlos):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_density
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		#print w1, w2, ratio
		p =  exp(self.galaxy_velocity_model.logL(vlos, e_vlos))*w1 \
			 +exp(self.foreground_velocity_model.logL(vlos, e_vlos)) * w2
		#p = exp(self.galaxy_velocity_model.logL(vlos, e_vlos))*self.w1 \
		#	+ self.w2 * exp(self.foreground_velocity_model.logL(vlos, e_vlos))
		return log(p)
		
	def memberRe(self, Re, weight=1.):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_density * weight
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		return (w1)
	
	def non_memberRe(self, Re, weight=1.):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_density * weight
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		return (w2)
		#return self.const_surface_density * self.w2

		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		
		re = stars.re
		vlos = stars.vlos_helio
		e_vlos = stars.e_vlos
		N = len(re)
		print re.max()
		#dsa
		def f_and_g(x):
			self.set_parameters(x)
			logL = 0
			for i in range(N):
				#p = exp(self.galaxy_velocity_model.logL(vlos[i], e_vlos[i]))*w1 \
				#	+ w2 * exp(self.foreground_velocity_model.logL(vlos[i], e_vlos[i]))
				#print p
				#if p == 0: import pdb; pdb.set_trace()
				#logL += log(p )
				logL += self.logL(re[i], vlos[i], e_vlos[i])
			print "    ", logL, x
			return -logL
		x0 = [log10(self.const_surface_density)]
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				x0.append(log10(self.light_profile.b))
			else:
				x0.append(log10(self.light_profile.scale))
		
		x0 += [k.get() for k in self.galaxy_velocity_model.parameters] + [k.get() for k in self.foreground_velocity_model.parameters]
		
		bounds = [(0, 1)] + [(None,None)] * (len(x0)-1)
		bounds = None
		#x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
		x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
		print x
		self.set_parameters(x)
		self.save()
		

class SimpleMultiFit(object):
	def __init__(self, filename, light_profile, light_model, observation, galaxy_velocity_model, foreground_velocity_model, const_surface_density=1e-4, fit_light=False):
		self.light_profile = light_profile
		self.light_model = light_model
		self.observation = observation
		self.galaxy_velocity_model = galaxy_velocity_model
		self.foreground_velocity_model = foreground_velocity_model
		self.filename = filename
		self.const_surface_density = const_surface_density
		self.fit_light = fit_light
		self.catalogues = 2
		#self.const_surface_density = self.light_profile.densityR(2., M=1.)
		
	def load(self):
		self.x0 = parameters = load(self.filename)
		logger.info("parameters: %s" % parameters)
		self.set_parameters(parameters)
		
	def set_parameters(self, x):
		#self.const_surface_density = 10**(x[0])
		self.const_surface_densities = [10**k for k in x[0:self.catalogues]]
		x = x[self.catalogues:]
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				self.light_profile.b = 10**(x[0])
			else:
				self.light_profile.scale = 10**(x[0])
			x = x[1:]
		else:
			pass
		n_params_galaxy = len(self.galaxy_velocity_model.parameters)
		n_params_foreground = len(self.foreground_velocity_model.parameters)
		for i in range(n_params_galaxy):
			value = x[i]
			self.galaxy_velocity_model.parameters[i].set(value)
			#print "set galaxy", self.galaxy_velocity_model.parameters[i].name, "to", value
		for i in range(n_params_foreground):
			value = x[n_params_galaxy+i]
			self.foreground_velocity_model.parameters[i].set(value)
			#print "set foreground", self.foreground_velocity_model.parameters[i].name, "to", value
		#r = abs(R)
		#self.w1 = 1/(1+r)
		#self.w2 = r/(1+r)
		#assert abs((self.w2/self.w1)- r) < 1e-6
		
	def save(self):
		parameters = [log10(k) for k in self.const_surface_densities]
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				parameters.append(log10(self.light_profile.b))
			else:
				parameters.append(log10(self.light_profile.scale))
		parameters += [k.get() for k in self.galaxy_velocity_model.parameters] + [k.get() for k in self.foreground_velocity_model.parameters]
		parameters = array(parameters)
		save(self.filename, parameters)
		
	def logL(self, Re, vlos, e_vlos, catalogue):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_densities[catalogue]
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		#print w1, w2, ratio
		p =  exp(self.galaxy_velocity_model.logL(vlos, e_vlos))*w1 \
			 +exp(self.foreground_velocity_model.logL(vlos, e_vlos)) * w2
		#p = exp(self.galaxy_velocity_model.logL(vlos, e_vlos))*self.w1 \
		#	+ self.w2 * exp(self.foreground_velocity_model.logL(vlos, e_vlos))
		return log(p)
		
	def memberRetot(self, Re, mask):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = mean([self.light_profile.densityR(Re, M=1.) for i in range(2) if ((1<<i)&int(mask))]) / mean([self.const_surface_densities[i] for i in range(2) if ((1<<i)&int(mask))])
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		return (w1)
		
	def memberRe(self, Re, catalogue=0, weight=1.):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_densities[catalogue] * weight
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		return (w1)
	
	def non_memberRe(self, Re, catalogue=0, weight=1.):
		Re = self.light_model.arcsec_to_kpc(Re)
		ratio = self.light_profile.densityR(Re, M=1.)/self.const_surface_densities[catalogue] * weight
		w1 = ratio/(1+ratio)
		w2 = 1/(1+ratio)
		return (w2)
		#return self.const_surface_density * self.w2

		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		star_groups = [stars.filter(lambda star: (int(star.catalogue_mask) & (1<<i)) != 0) for i in range(self.catalogues)]
		for i in range(self.catalogues):
			print "length", len(star_groups[i])
		print self.catalogues
		
		re = stars.re
		vlos = stars.vlos_helio
		e_vlos = stars.e_vlos
		N = len(re)
		cats = stars.catalogue_mask.astype(int)
		print re.max()
		#dsa
		def f_and_g(x):
			self.set_parameters(x)
			logL = 0
			if 1:
				#for i in range(self.catalogues):
				#	for star in star_groups[i]:
				#		logL += self.logL(star.re, star.vlos_helio, star.e_vlos, i)
				#count1 = len(re)
				#count2 = 0
				for i in range(len(re)):
					#logL1 = -inf
					for j in range(2):
						if (1<<j) & (cats[i]):
							#logL1 = max(logL1, self.logL(re[i], vlos[i], e_vlos[i], j))
							logL += self.logL(re[i], vlos[i], e_vlos[i], j)
							#count2 += 1
							#print j, cats[i]
							#assert j == 0
					#logL += logL1
				#print count1, count2
			else:
				#for i in range(N):
				#	logL += self.logL(re[i], vlos[i], e_vlos[i], 0)
				pass
			print "    ", logL, x
			return -logL
		x0 = [log10(self.const_surface_density) for i in range(self.catalogues)]
		if self.fit_light:
			if hasattr(self.light_profile, "b"):
				x0.append(log10(self.light_profile.b))
			else:
				x0.append(log10(self.light_profile.scale))
		
		x0 += [k.get() for k in self.galaxy_velocity_model.parameters] + [k.get() for k in self.foreground_velocity_model.parameters]
		
		#bounds = [(0, 1)] + [(None,None)] * (len(x0)-1)
		bounds = None
		#x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
		x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
		print x
		self.set_parameters(x)
		self.save()
		

		
