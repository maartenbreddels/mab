# -*- coding: utf-8 -*-
#from constants import *
from numpy import *
import mab.constants
import mab.astrounits
#from numpy import *
#import numpy
#import math
#from math import pi
#import mab.cosmology
from mab.gd import gdfast
#import sys
#from scipy.optimize import fsolve, fmin, brentq
#from scipy.integrate import quad, inf
import scipy.optimize

G = mab.constants.G.asNumber(mab.astrounits.KM**2/mab.astrounits.S**2 *mab.astrounits.KPC/mab.astrounits.MSOL)

class PotentialAxi(object):
	def potentialrhoz_eff(self, rho, z, Lz):
		return self.potentialrhoz(rho, z) + Lz**2/(2*R**2)

class MiyamotoNagai(PotentialAxi):
	def __init__(self, M, a, b, G=G):
		self.M = M
		self.a = a
		self.b = b
		self.G = G
		
	def densityrhoz(self, rho, z):
		x = sqrt(z**2+self.b**2)
		f1 =self.a * rho**2 + (self.a + 3*x) *(self.a + x)**2
		f2 = (rho**2 + (self.a + x)**2)**(5./2) * (z**2+self.b**2)**(3./2)
		return (self.b**2*self.M)/4*pi * f1/f2
	
	def potentialrhoz(self, rho, z):
		return - self.G * self.M / sqrt(rho**2 + (self.a + sqrt(z**2+self.b**2))**2)
	
	
class Logarithmic(PotentialAxi):
	def __init__(self, v0=1., q=1., xc=0, G=G):
		self.v0 = v0
		self.q = q
		self.xc = xc
		self.G = G
		self.fast_ = gdfast.LogarithmicAxi(self.v0, self.q, self.xc, self.G)
		
	def fast(self):
		return self.fast_
	
	def densityxz(self, x, z):
		def wrap(x, z):
			return self.fast_.densityxz(x, z)
		return vectorize(wrap)(x, z)
		
	def potentialxz(self, x, z):
		def wrap(x, z):
			return self.fast_.potentialxz(x, z)
		return vectorize(wrap)(x, z)
		#return 0.5*self.v0**2 * log(rho**2 + z**2/self.q**2)
	
	def potentialxz_eff(self, x, z, Lz):
		def wrap(x, z, Lz):
			return self.fast_.potentialxz_eff(x, z, Lz)
		return vectorize(wrap)(x, z, Lz)
		#return 0.5*self.v0**2 * log(rho**2 + z**2/self.q**2)
	


		
		
