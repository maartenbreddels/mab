# -*- coding: utf-8 -*-
from scipy.optimize import fsolve, fmin, brentq
from scipy.integrate import quad
#from numy inf
from numpy import *
from mab.gd import gdfast 
from mab.gd.jeans.jeans import Jeans, JeansAnisotropicConstant
import inspect

class Param(object):
	def __init__(self, name, argname, info, default):
		self.name = name
		self.argname = argname
		self.info = info
		self.default = default
		
	def __str__(self):
		if self.default is not None:
			return "%s: %s"
		else:
			return "%s: %s [default: %r]" % (self.name, self.info, self.default) 
class GalaxyFactory1C_constant_anisotropy(object):
	def __init__(self, stellar_profile_class, name):
		self.stellar_profile_class = stellar_profile_class
		#argspec = inspect.getargspec(stellar_profile_class.__init__)
		#self.parameters = 0
		self.parameters = []
		self.parametersinfo = []
		for arg, description, defaultvalue in stellar_profile_class.arguments():
			print arg, description 
			self.parameters.append(name +"_" +arg)
			self.parametersinfo.append("%s: [default: %s]" % (description, defaultvalue))
		self.parameters.append("distance")
		self.parametersinfo.append(description)
	
	def create(self, parameters):
		pass
		

class Galaxy1C_constant_anisotropy(object):
	def __init__(self, stellar_profile, distance, beta):
		self.stellar_profile = stellar_profile
		self.distance = distance
		self.beta = beta
		self.stellar_profilef = self.stellar_profile.fast() # keep reference
		#stellar_profilehalf = 
		#self.stellar_profilef = self.stellar_profilehalf.fast() # keep reference
		self.null_profilef = gdfast.NullProfile()
		self._fast = self.fast()
		
		
	def __reduce__(self):
		return (Galaxy1C_constant_anisotropy, (self.stellar_profile, self.distance, self.beta))
		
	def vmax(self):
		return sqrt(-2*self.potentialr(1e-9))
		
	def jeansfast(self):
		return gdfast.JeansAnisotropicConstant(self.stellar_profilef, self.null_profilef, self.beta, 10.)
	
	def fast(self):
		return gdfast.Galaxy1C_constant_anisotropy(self.stellar_profilef, self.distance, self.beta)
	
	def jeans(self):
		return JeansAnisotropicConstant(self.stellar_profile, self.getprofiles(), self.beta)
		
	def kpc_to_arcsec(self, d): 
		return d * 1.0/(self.distance) / (pi/180/60/60)
	
	def pc_to_arcsec(self, d): 
		return d * 1e3/(self.distance) / (pi/180/60/60)
	
	def arcsec_to_kpc(self, R): 
		return R / (1.0/(self.distance) / (pi/180/60/60))
	
	def arcsec_to_pc(self, d): 
		return R / (1e3/(self.distance) / (pi/180/60/60))
	
	
	def getprofiles(self):
		return [self.stellar_profile]
	
	def densityr(self, r):
		return self.stellar_profile.densityr(r)
	
	def potentialr(self, r):
		return self.stellar_profile.potentialr(r)
	
	def dphidr(self, r):
		return self.stellar_profile.dphidr(r)
	
	def vcirc(self, r):
		return sqrt(r * self.dphidr(r))
	
	def L_and_r_max_at_E(self, E, rtry=1):
		def f(r, E=E):
			r = abs(r)
			return -r**2*(2*(E-self.potentialr(r)))
		rLmax, Lmaxsq = fmin(f, rtry, full_output=True, disp=False)[:2]
		return sqrt(abs(Lmaxsq)), abs(rLmax[0])
	 
	def Lmax_and_E_at_r(self, r):
		def f(L, r=r):
			return L**2/(2*r**2) + self.potentialr(r)
		Lmax, E = fmin(f, 1. ,full_output=True, disp=False)[:2]
		return Lmax, E #sqrt(abs(Lmaxsq)), abs(rLmax[0])
	 
	def Lmax_at_E(self, E):
		return self.L_and_r_max_at_E(E)[0]
	
	def fE(self, E):
		# distribution function for isotropic cases
		assert self.beta == 0, "fE only available for isotropic models, beta = %f" % self.beta
		return self.stellar_profile.fE(E)
		
	def fEL(self, E, L):
		return self.fE(E)
		

	def findR_at_EL(self, E, L, rtry=1., r1=None, r2=None):
		def f(r, E=E, L=L):
			r = abs(r)
			return L**2/(2*r**2) + (self.potentialr(r) - E)
		if r1 and r2:
			r = brentq(f, r1, r2)
		else:
			r = fsolve(f, rtry) #,full_output=True)[:2]
		#Lmax = extra[0]
		#print r
		return abs(r)
	
	def get_apo_peri(self, E, L, rtry=1e-5):
		def ekinrad(r): # kinetic energy in radial direction
			return (-L**2/(2*r**2) - (self.potentialr(r) - E))
		
		# just find apo or pericenter
		rstart = self.findR_at_EL(E, L, rtry)
		s = (1+1e-5) # scale factor for testing apo/peri
		r = rstart
		#print rstart, ekinrad(rstart/s), ekinrad(rstart*s), -L**2/(2*r**2) - (self.potentialr(r) - E), -L**2/(2*r**2) , (self.potentialr(r) - E)
		if (ekinrad(rstart/s) < 0) and (ekinrad(rstart*s) > 0): # we found pericenter
			rp = rstart
			ra = brentq(ekinrad, rstart*s, 1e9)
		else: # we found apocenter
			ra = rstart
			rp = brentq(ekinrad, 1e-9, rstart/s)
		
		# sanity checks
		assert ekinrad(ra*s) < 0, "available energy ar r > r_apo should be negative"
		assert ekinrad(rp/s) < 0, "available energy ar r < r_peri should be negative"
		assert ekinrad(ra/s) > 0, "available energy ar r < r_apo should be positive"
		assert ekinrad(rp*s) > 0, "available energy ar r > r_peri should be positive"
		assert ra > rp, "apocenter should be larger than pericenter" 
		return ra, rp
		
		
	def gEL(self, E, L):
		# density of states 
		return 8*pi**2*L * self.Tr(E, L) # see BT problem 4.8
		
	def gEL2(self, E, L):
		# density of states 
		return 8*pi**2*L * self.Tr2(E, L) # see BT problem 4.8
		
	def gE(self, E):
		# density of states 
		def dgE(r):
			return sqrt(2*(E-potentialr(r)))*r*r
		rmax = self.findR_at_EL(E, L)
		assert not isnan(rmax), "sanity check"
		i, err = quad(dgE, 0, rmax)
		return i
			
			
	def Tr(self, E, L, rtry=1e-5):
		# radial period
		def dg(r):
			return 2/sqrt(-L**2/r**2-2*self.potentialr(r)+2*E)
		
		#print "1"
		ra, rp = self.get_apo_peri(E, L, rtry=rtry)
		#print "2"
		# initial orbit, can be apo or pericenter
		
		
		
		i, err = quad(dg, rp, ra)
		#epsilon = 1e-2
		#print ra, rp, dg(ra*(1-epsilon)), dg(rp*(1+epsilon))
		#t = self._fast.Tr(E, L, rp*(1+epsilon), ra*(1-epsilon))
		#print t
		#return t
		if abs(err/i) > 1e-4:
			import pdb; pdb.set_trace()  
		#print "3"
		if isnan(i):
			import pdb; pdb.set_trace()  
		return i
		
	def Tphi(self, E, L, Tr, rtry=1e-5):
		# radial period
		def dphi(r):
			return 2 * L /(r*r*sqrt(-L**2/r**2-2*self.potentialr(r)+2*E))
		
		#print "1"
		ra, rp = self.get_apo_peri(E, L, rtry=rtry)
		
		deltaphi, err = quad(dphi, rp, ra)
		Tphi = 2 * pi / abs(deltaphi) * Tr
		return Tphi
		
		
	def Tr2(self, E, L, rtry=1e-5):
		# radial period
		def dg(r):
			return 2/sqrt(-L**2/r**2-2*self.potentialr(r)+2*E)
		
		#print "1"
		ra, rp = self.get_apo_peri(E, L, rtry=rtry)
		#print "2"
		# initial orbit, can be apo or pericenter
		
		
		
			
		i, err = quad(dg, rp, ra)
		print i, err
		
		epsilon = 1e-9
		print ra, rp, dg(ra*(1-epsilon)), dg(rp*(1+epsilon))
		try:
			t = self._fast.Tr(E, L, ra*(1-epsilon), rp*(1+epsilon))
		except:
			raise
		print t
		return t
		if abs(err/i) > 1e-4:
			import pdb; pdb.set_trace()  
		#print "3"
		if isnan(i):
			import pdb; pdb.set_trace()  
		return i
		
	def Jr(self, E, L, rtry=1e-5):
		# action in radial dir
		def dJr(r):
			return sqrt(2*E-2*self.potentialr(r)-L**2/r**2)
		
		#print "1"
		ra, rp = self.get_apo_peri(E, L, rtry=rtry)
		
		Jrpi, err = quad(dJr, rp, ra)
		Jr = Jrpi/pi
		
		if 0:
			print "rp", 2*E-2*self.potentialr(rp)- L**2/rp**2
			rp *= 0.999
			print "rp", 2*E-2*self.potentialr(rp)- L**2/rp**2
			rp /= 0.999
			rp /= 0.999
			print "rp", 2*E-2*self.potentialr(rp)- L**2/rp**2
			print "ra", 2*E-2*self.potentialr(ra)- L**2/ra**2
			ra *= 1.001
			print "ra", 2*E-2*self.potentialr(ra)- L**2/ra**2
			ra /= 1.001
			ra /= 1.001
			print "ra", 2*E-2*self.potentialr(ra)- L**2/ra**2
		return Jr
		
			
			
class Galaxy2C_constant_anisotropy(Galaxy1C_constant_anisotropy):
	def __init__(self, stellar_profile, dm_profile, distance, beta):
		self.dm_profile = dm_profile
		self.dm_profilef = self.dm_profile.fast() # keep reference
		Galaxy1C_constant_anisotropy.__init__(self, stellar_profile, distance, beta)
		
	def __getinitargs__(self):
		return (self.stellar_profile, self.dm_profile, self.distance, self.beta)
	
	def fast(self):
		return gdfast.Galaxy2C_constant_anisotropy(self.stellar_profilef, self.dm_profilef, self.distance, self.beta)
		
	def jeansfast(self):
		return gdfast.JeansAnisotropicConstant(self.stellar_profilef, self.dm_profilef, self.beta, 10.)
		
	def getprofiles(self):
		return [self.stellar_profile, self.dm_profile]
	
	def potentialr(self, r):
		return self.stellar_profile.potentialr(r) + self.dm_profile.potentialr(r)
	
	def dphidr(self, r):
		return self.stellar_profile.dphidr(r) + self.dm_profile.dphidr(r)
	
	def fE(self, E):
		assert False, "not implemented yet"
			