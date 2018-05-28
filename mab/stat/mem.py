# -*- coding: utf-8 -*-
from numpy import *
import scipy.optimize
import scipy.integrate

# Maximum entropy moment based pdf
class MEM(object):
	def __init__(self, sigma=1., gamma=3., xmax=3., n_moments=4):
		self.sigma = sigma
		self.gamma = gamma
		self.xmax = xmax
		self.n_moments = n_moments
		self.lambdas = zeros(n_moments)
		#self.lambdas[0] = 0
		#self.lambdas[1] = 1.
		self.l2 = 1.
		self.l4 = 0
		self.normalization = 1.
		self.solve()
		self.label = "mem g2=%.2f" % self.gamma
		#print ">>>>>>>>>>>>>>", self.moment(0)
		
	def solve(self):
		def f(args):
			#print "args", args
			self.normalization, self.l2, self.l4 = args
			self.normalization = abs(self.normalization)
			m0 = self.moment(0)
			m2 = self.moment(2)
			m4 = self.moment(4)
			k = m4/m2**2
			#print ">>", self.normalization, m2, m4, k
			#sprint ">>   ", self.normalization, m2, m4, k-self.
			error = (m2/self.sigma**2 - 1)**2 + (k/self.gamma - 1) **2 + (m0 - 1)**2
			#print error
			return error
		res = scipy.optimize.fmin(f, [1, 1, 0])
		print "result for mem fitting", res
		self.normalization, self.l2, self.l4 = res
		print "moment0", self.moment(0)
		print "moment2", self.moment(2)
		m4 = self.moment(4)
		print "moment4", m4
		print "gamma  ", m4/self.moment(2)**2
	def __call__(self, x):
		#return exp([x**k * self.lambdas[k] for k in range(self.n_moments)])
		return exp(x**2 * self.l2 + x**4*self.l4) * self.normalization * (abs(x) <= self.xmax)
	
	def moment(self, k=1):
		#I1, error = scipy.integrate.quad(lambda x, k: 2* self(x) * x**k, 0., 1, args=(k,)) #self.xmax)
		#I2, error = scipy.integrate.quad(lambda x, k: 2* self(x) * x**k, 1., self.xmax, args=(k,))
		#return I1+I2
		I, error = scipy.integrate.quad(lambda x: 2* self(x) * x**k, 0., self.xmax, limit=2000)
		return I
		 