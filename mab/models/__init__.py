# -*- coding: utf-8 -*-
from mab.functions import gaussian
import emcee
import mab.gd.logging as logging
logger = logging.getLogger("mab.models")
from numpy import *
import sys
from kaplot import *
import scipy


class TwoGaussian(object):
	def __init__(self, mu1, sigma1, mu2, sigma2, ratio, values):
		self.mu1 = mu1
		self.sigma1 = sigma1
		self.mu2 = mu2
		self.sigma2 = sigma2
		self.ratio = ratio
		self.values = numpy.array(values)
		self.x0 = [mu1, sigma1, mu2, sigma2, ratio]
		logger.info("values length: %d" % len(values)) 
		
	def __len__(self):
		return 5
		
	def __call__(self, x):
		#print ">>", p
		mu1, sigma1, mu2, sigma2, ratio = x
		ratio = 10**(ratio)
		#ratio is f1/f2 
		f1 = ratio/(1.+ratio)
		f2 = 1./(1.+ratio)
		print ".",
		sys.stdout.flush()
		ps = gaussian(self.values, mu1, sigma1)*f1 + gaussian(self.values, mu2, sigma2)*f2
		return -numpy.sum(log(ps))
		#return sum([log(gaussian(xi, mu1, sigma1)*f1 + gaussian(xi, mu2, sigma2)*f2) for xi in self.values])
	
class ModelBestFit(object):
	def __init__(self, model, ):
		self.model = model
		
	def run(self, args, opts, scope):
		res = scipy.optimize.fmin(self.model, self.model.x0, ftol=1e-7)
		print res
		#res = [  3.15544562e+01 , -1.17340052e-03 ,  1.12697458e+02 ,  1.08677279e+02,-2.40627963e-01]
		mu1, sigma1, mu2, sigma2, ratio = res
		x = mgrid[-200:200:300j]
		ratio = 10**(ratio)
		#ratio is f1/f2 
		f1 = ratio/(1.+ratio)
		f2 = 1./(1.+ratio)
		xi = x
		
		
		#box()
		mozaic(2,2,box)
		select(0,0)
		y = gaussian(xi, mu1, sigma1)*f1 + gaussian(xi, mu2, sigma2)*f2
		graph(x, y)
		y = gaussian(xi, mu1, sigma1)*f1
		graph(x, y, color="green")
		y = gaussian(xi, mu2, sigma2)*f2
		graph(x, y, color="blue")
		b = scope["besancon_data"]
		graph(b.p_helio_bins, b.p_helio, color="red")
		
		select(1,0)
		y = gaussian(xi, mu1, sigma1)*f1 + gaussian(xi, mu2, sigma2)*f2
		graph(x, cumsum(y)/sum(y))
		b = scope["besancon_data"]
		graph(b.p_helio_bins, cumsum(b.p_helio)/sum(b.p_helio), color="red")
		
		draw()
		

class ModelMCMC(object):
	def __init__(self, model, nwalkers=10, N=1000):
		self.model = model
		self.nwalkers = nwalkers
		self.N = N
		
		
	def run(self, args, opts, scope):
		dim = len(self.model)
		logger.info("dimension: %d" % dim)
		x0 = [[k+random.random() for k in self.model.x0] for i in range(self.nwalkers)]
		print x0
		#p0 = [(random.normal(0, 50),random.normal(0,50)) for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(self.nwalkers, dim, self.model, args=[], threads=2)
		pos, lnp, state = sampler.run_mcmc(x0, self.N)
		print pos
		print lnp
		print state
		#print res
		print sampler.acceptance_fraction
		print "N", sampler.flatchain.shape

		values = sampler.flatchain.T
		print values
		mozaic(dim, dim, box)
		for i in range(dim):
			for j in range(dim):
				select(i,j)
				if i == j:
					#histogram(values[i])
					pass
				else:
					scatter(values[i], values[j])
		if 0:
			select(1, 0)
			histogram(x, datamin=-10, datamax=10, bincount=150., normalize=True)
			x = mgrid[0:1:100j]
			#graph(x, gaussian(x, 0.5, 0.1), color="red")
			select(1, 1)
			histogram(y, datamin=-10, datamax=10, bincount=150., normalize=True)
			#graph(x, gaussian(x, 0.5, 0.1), color="red")
			
		draw()	
