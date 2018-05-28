from numpy import *
import numpy
import numpy.linalg


def gaussian(x, mu, sigma):
	y = 1/(numpy.sqrt(math.pi*2)*sigma)*numpy.exp(-0.5*((x-mu)/sigma)**2)
	return y
	
def gaussian2d(x, mu_x, sigma_x, y, mu_y, sigma_y, rho=1.):
	y = 1/(numpy.sqrt(math.pi*2)*sigma_x*sigma_y*numpy.sqrt(1-rho**2))*numpy.exp(-1/(2*(1-rho**2)) * ( (x-mu_x)**2/sigma_x**2 + (y-mu_y)**2/sigma_y**2 - 2 * rho * (x-mu_x) * (y-mu_y)/(sigma_x*sigma_y)) )
	return y

class Gaussian2d(object):
	def __init__(self, mux, muy, sigmax, sigmay, rho):
		self.mux = mux
		self.muy = muy
		self.cov = numpy.matrix([[sigmax**2, sigmax*sigmay*rho], [sigmax*sigmay*rho, sigmay**2]])
		self.covI = self.cov.I
		self.det = numpy.linalg.det(self.cov)
		
	def __call__(self, x, y=None):
		M = array(self.covI)
		if y is None:
			x, y = x
		x = x - self.mux
		y = y - self.muy
		#u = array([x, y])
		#print x
		a = (x * M[0,0] + y * M[1,0]) 
		b = (x * M[0,1] + y * M[1,1])
		chisq = a * x + b * y
		#chisq = -(a + b)
		#print u.shape, M.shape
		#a = u.T * M[0] + u.T * M[1]
		#print "a",a.shape
		
		#b = a * u.T
		#chisq = sum(b, axis=-1)
		#print chisq.shape  
		#b = dot(u, a)
		#b = tensordot(u, a, [(0,), (-1,)])
		#print b.shape
		#chisq = cum(
		N = 1/(2*pi * self.det**0.5)
		return N * exp(-0.5*chisq) 