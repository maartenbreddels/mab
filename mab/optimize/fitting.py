# -*- coding: utf-8 -*-
import numpy
import mab.gd.logging as logging
logger = logging.getLogger("gd.optimize.fitting")
from kaplot import *

class Fitting(object):
	def __init__(self, fitter, data, xattr, sigma_xattr, filename):
		self.fitter = fitter
		self.data = data
		self.xattr = xattr
		self.sigma_xattr  = sigma_xattr
		self.filename = filename
		
		
	def run(self, args, opts, scope):
		obj = self.data.load()
		xtest = x = getattr(obj, self.xattr)
		sigma_x = getattr(obj, self.sigma_xattr)
		print x
		print sigma_x
		#dsa
		model = self.fitter.get_model()
		if 1:
			self.fitter.fit(x, sigma_x)
			self.x = [param.get() for param in model.parameters]
			for param in model.parameters:
				print param.name, param.get()
			if 0:
				for x, e in zip(x, sigma_x):
					y = model.logL(x, e)
					print y, x, e
					if numpy.isnan(y):
						import pdb; pdb.set_trace()
			print "x =", self.x
			self.save()
		else:
			self.load()
		xmin = -200
		xmax = 300
		mozaic(2,1,box)
		histogram(x, binwidth=2, datamin=xmin, datamax=xmax, normalize=True)
		xs = arange(xmin, xmax, 0.5)
		logy = [model.logL(x, 0) for x in xs]
		graph(xs, exp(logy), color="red")
		
		
		select(1,0)
		N = len(sigma_x)
		M = 10000
		logLs = []
		for j in range(M):
			#print j, M
			#samples_new = array([model.sample() for i in range(len(samples))])
			logL = sum([model.logL(model.sample(), sigma_x[i]) for i in range(N)])
			logLs.append(logL)
		histogram(logLs, bincount=100, normalize=True)#, datamin=datamin, datamax=datamax, )
		logL = sum([model.logL(xtest[i], sigma_x[i]) for i in range(N)])
		vline(logL, color="green")
		
		print logL
		draw()
		
		
		
	def load(self):
		logger.info("loading fit parameters from: " + self.filename)
		values = numpy.load(self.filename)
		for param, value in zip(self.fitter.get_model().parameters, values):
			#print param, param.name, value
			param.set(value)
		#self.fitter.get_model().parameters[-1].set(1)
		#self.fitter.get_model().parameters[-1].set(5)
		#param.set(0.4)
	
	def save(self):
		logger.info("saving fit parameters to: " + self.filename)
		values = numpy.save(self.filename, self.x)
		
		
		