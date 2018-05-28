# -*- coding: utf-8 -*-
from scipy.optimize import fmin
from numpy import sqrt, array

def findLmax(E, potentialr):
	def f(r, E=E):
		r = abs(r)
		return -r**2*(2*(E-potentialr(r)))
	rLmax, Lmaxsq = fmin(f, 0.1,full_output=True, disp=False)[:2]
	#Lmax = extra[0]
	#print Lmax
	return abs(rLmax[0]), sqrt(abs(Lmaxsq))


class CommandPrintVdisp0(object):
	def __init__(self, galaxy):
		self.galaxy = galaxy
		
		
	def run(self, args, opts, scope):
		jeans = self.galaxy.jeans()
		print jeans
		print jeans.sigmar(0)
		print jeans.sigma_los(array([0, 1e-7, 0.01, 0.1, 0.]))