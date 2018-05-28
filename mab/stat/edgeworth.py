# -*- coding: utf-8 -*-
from numpy import *
import scipy

def h3(x):
	return x**3-3*x
def h4(x):
	return x**4-6*x**2+3
def h6(x):
	return x**6-15*x**4+45*x**2-15

fac3 = 1*2*3
fac4 = fac3*4
fac6 = fac4*5*6

	
class Edgeworth(object):
	def __init__(self, sigma, gamma):
		self.sigma = sigma
		self.gamma = gamma
		self.s3 = 0
		self.kappa2 = self.sigma**2
		self.label = "edw. g2=%.2f" % self.gamma
		#kappa4 = gamma * kappa2**2
		
		#self.kappa2 = sigma**2
		self.kappa4 = gamma * self.kappa2**2
		if gamma < 3:
			pass
			#self.kappa4 = self.kappa4 * (1 - (3-gamma)**2*0.3 + (3-gamma)**3*0.07)
			#self.sigma = self.sigma * (1 - (3-gamma)**2*0.02 + (3-gamma)**3*0.03)
		self.s4 = (self.kappa4 - 3 * sigma**4)/self.sigma**6
		#print "gamma/w", gamma, self.w
		#self.fa
		
	def logpositivity(self, x):
		return log(x)
	
	def __call__(self, x):
		a = 1/self.sigma * 1/sqrt(2*pi) * exp(-0.5 * (x/self.sigma)**2)
		b1 = self.sigma*(1./fac3 * self.s3 * h3(x/self.sigma))
		b2 = self.sigma**2 * (1./fac4 * self.s4 * h4(x/self.sigma) + 10./fac6 * self.s3**2*h6(x/self.sigma))
		#return maximum(1e-10, a * (1 + b1 + b2))
		return a * (1 + b1 + b2)
			
			
	def logprob(self, x):
		a = log(1/self.sigma * 1/sqrt(2*pi) * exp(-0.5 * (x/self.sigma)**2))
		b1 = self.sigma*(1./fac3 * self.s3 * h3(x/self.sigma))
		b2 = self.sigma**2 * (1./fac4 * self.s4 * h4(x/self.sigma) + 10./fac6 * self.s3**2*h6(x/self.sigma))
		P = self.logpositivity(b1+b2)
		print b1
		print b2
		#import pdb; pdb.set_trace()
		#if any(P<0):
		#	import pdb; pdb.set_trace()
		#assert all(P>=0), P
		#print exp(a + P)
		if any(isinf(exp(a + P))):
			import pdb; pdb.set_trace()
		if any(isnan(exp(a + P))):
			import pdb; pdb.set_trace()
		return (a + P)
	 
	def cumulants(self, x, p, dx):
		mm1 = sum(x*p*dx)
		mm2 = sum(x**2*p*dx)
		mm3 = sum(x**3*p*dx)
		mm4 = sum(x**4*p*dx)
		mm5 = sum(x**5*p*dx)
		mm6 = sum(x**6*p*dx)
		c1 = mm1
		c2 = mm2
		c3 = mm3
		c4 = mm4 - 3*c2**2
		c5 = mm5 - 10*c3*c2
		c6 = mm6 - 15*c4*c2 - 10*c3**2-15*c2**3
		
		return array([c1, c2, c3, c4, c5, c5])
		
	def moment(self, k=1):
		#I1, error = scipy.integrate.quad(lambda x, k: 2* self(x) * x**k, 0., 1, args=(k,)) #self.xmax)
		#I2, error = scipy.integrate.quad(lambda x, k: 2* self(x) * x**k, 1., self.xmax, args=(k,))
		#return I1+I2
		I, error = scipy.integrate.quad(lambda x: 2* self(x) * x**k, 0., scipy.integrate.inf, limit=2000)
		n, error = scipy.integrate.quad(lambda x: 2* self(x), 0., scipy.integrate.inf, limit=2000)
		#print n
		return I/n
	
		
class EdgeworthMod(Edgeworth):
	def __init__(self, sigma, gamma):
		Edgeworth.__init__(self, sigma, gamma)
		self.w = self.sigmoid()
		
	def sigmoid(self):
		return 1./(1+exp(-(self.gamma-3)*10000000))
	
	def logpositivity(self, u):
		s = 1.0
		#return log( (exp(s*u)/s -1/s + 1)*(1-self.w)  + self.w*(1+u) )
		p = u * 1.
		mask = u > 0
		p[mask] = log(1+u[mask])
		return p
		if u < 0:
			return u
		else:
			return log(1+u)
		if (1-self.w) < 1e-10:
			return log(1+u)
		if self.w < 1e-10:
			return u
		else:
			return log( (exp(u))*(1-self.w)  + self.w*(1+u) )
		return log( (exp(u))*(1-self.w)  + self.w*(1+u) )
		return (1+ tanh(u))*(1-self.w)  + self.w*(1+u)
		
class EdgeworthCorrected(EdgeworthMod):
	def __init__(self, sigma, gamma):
		a1, a2, a3, b1, b2,b3 = [-0.10971413, -0.03265063,  0.02697869,  -0.52812814, -0.05336269, -0.0592689] # unnormalized
		if (gamma < 3.0):
			#print "%e" % gamma
			sigma_mod = sigma * (1 + (3-gamma)**2*a1 + (3-gamma)**3*a2 + (3-gamma)**4*a3)
			kappa2 = sigma**2
			kappa4 = gamma * kappa2**2
			kappa4_mod = kappa4 * (1 + (3-gamma)**2*b1 + (3-gamma)**3*b2 + (3-gamma)**4*b3)
			gamma_mod = kappa4_mod/kappa2**2
		else:
			sigma_mod = sigma
			gamma_mod = gamma
		EdgeworthMod.__init__(self, sigma_mod, gamma_mod)