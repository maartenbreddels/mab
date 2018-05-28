from interpol import *
from numpy import *

class InterpolCache(object):
	def __init__(self, f, xmin, xmax, dx, interpolator=interpolator1d_le):
		self.f = f
		self.x = arange(xmin, xmax+dx/2, dx)
		self.xmin = xmin
		self.xmax = xmax
		self.dx = dx
		self.y = f(self.x)
		self.finterpol = interpolator(self.x, self.y)
		self.call = lambda x: self.finterpol(x)
		self.callv = vectorize(self.call)
		
	def __call__(self, x):
		return self.callv(x) #self.finterpol(x)
		#if any(x >
		
		
	#	return  
	
		