# -*- coding: utf-8 -*-
from scipy.optimize import fsolve, fmin, brentq
from scipy.integrate import quad
from numpy import *
from mab.gd import gdfast 
from mab.gd.jeans.jeans import Jeans, JeansAnisotropicConstant
import mab.gd.logging as logging
import scipy
logger = logging.getLogger("gd.schw.galaxy")


class ProfileModel1C(object):
	def __init__(self, light_profile): 
		self.light_profile = light_profile
		self.light_profile_fast = self.light_profile.fast()
		self.fast_ = gdfast.ProfileModel1C(self.light_profile_fast)
		
	def vmax(self):
		return sqrt(-2*self.potentialr(1e-5))
	
	def fast(self):
		return self.fast_
	
	def jeans(self):
		return 
	
	def total_mass(self):
		return self.light_profile.total_mass()
		
	#def fast(self):
	#	return None #gdfast.Galaxy1C_constant_anisotropy(self.stellar_profilef, self.distance, self.beta)
	
	def getprofiles(self):
		return [self.light_profile]
	
	def densityr(self, r, M=None):
		return self.light_profile.densityr(r, M)
	
	def potentialr(self, r):
		return self.light_profile.potentialr(r)
	
	def dphidr(self, r):
		return self.light_profile.dphidr(r)
	
	def vcirc(self, r):
		return sqrt(r * self.dphidr(r))
	
	def L_and_r_max_at_E(self, E, rtry=1):
		def f(r, E=E):
			r = abs(r)
			return -r**2*(2*(E-self.potentialr(r)))
		rLmax, Lmaxsq = scipy.optimize.fmin(f, rtry, full_output=True, disp=False)[:2]
		return sqrt(abs(Lmaxsq)), abs(rLmax[0])
	 
	def Lmax_and_E_at_r(self, r):
		def f(L, r=r):
			return L**2/(2*r**2) + self.potentialr(r)
		Lmax, E = fmin(f, 1. ,full_output=True, disp=False)[:2]
		return Lmax, E #sqrt(abs(Lmaxsq)), abs(rLmax[0])
	 
	def Lmax_at_E(self, E):
		return self.L_and_r_max_at_E(E)[0]
	
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
		try:
			rstart = rstart[0]
		except:
			pass
		s = (1+1e-5) # scale factor for testing apo/peri
		r = rstart
		#print rstart, ekinrad(rstart/s), ekinrad(rstart*s), -L**2/(2*r**2) - (self.potentialr(r) - E), -L**2/(2*r**2) , (self.potentialr(r) - E)
		if (ekinrad(rstart/s) < 0) and (ekinrad(rstart*s) > 0): # we found pericenter
			rp = rstart
			ra = brentq(ekinrad, rstart*s, 1e9)
		else: # we found apocenter
			ra = rstart
			rp = brentq(ekinrad, 1e-9, rstart/s)
		
		#while ekinrad(ra*s) < 0:
		#	ra 
		# sanity checks
		if 0:
			assert ekinrad(ra*s) < 0, "available energy ar r > r_apo should be negative"
			assert ekinrad(rp/s) < 0, "available energy ar r < r_peri should be negative"
			assert ekinrad(ra/s) > 0, "available energy ar r < r_apo should be positive"
			assert ekinrad(rp*s) > 0, "available energy ar r > r_peri should be positive"
			assert ra > rp, "apocenter should be larger than pericenter" 
		return ra/s, rp*s
		
	def gEL(self, E, L):
		# density of states 
		return 8*pi**2*L * self.Tr(E, L) # see BT problem 4.8
	
	def Tr(self, E, L, rtry=1e-5, ra=None, rp=None):
		def ekinrad(r): # kinetic energy in radial direction
			return (-L**2/(2*r**2) - (self.potentialr(r) - E))
		# radial period
		def dg(r):
			x = 2/sqrt(-L**2/r**2-2*self.potentialr(r)+2*E)
			#print r, x, -L**2/r**2-2*self.potentialr(r)+2*E
			if isinf(x):
				
				import pdb; pdb.set_trace()
			return x
			#return 2/sqrt(-L**2/(2*r**2)-(self.potentialr(r)-E))/sqrt(2)
		
		#print "1"
		if ra is None and rp is None:
			ra, rp = self.get_apo_peri(E, L, rtry=rtry)
		scale = 1+1e-5
		while isnan(dg(ra)):
			ra /= scale
		while isnan(dg(rp)):
			rp *= scale
		#print "2"
		# initial orbit, can be apo or pericenter
		
		
		#print E, L, ra, rp
		#return 1.
		#print "new rp/ra", rp, ra, "A", dg(rp), dg(ra)
		ret = quad(dg, rp, ra, epsabs=0, epsrel=1e-7, full_output=1)
		ravg = 0.5 * (ra+rp)
		i = ret[0]
		#import pdb; pdb.set_trace()
		logger.debug("ra=%r (%r,%r), rp=%r (%r,%r), ravg=%r (%r) %f" % (ra, dg(ra), ekinrad(ra), rp, dg(rp), ekinrad(rp), ravg, dg(ravg),i))
		#print "l", len(ret)
		#print "%g" % i
		#return 0
		err = ret[1]
		#print ret[2]
		#return 1.
		#epsilon = 1e-2
		#print ra, rp, dg(ra*(1-epsilon)), dg(rp*(1+epsilon))
		#t = self._fast.Tr(E, L, rp*(1+epsilon), ra*(1-epsilon))
		#print t
		#return t
		#if abs(err/i) > 1e-4:
		#	import pdb; pdb.set_trace()  
		#print "3"
		if isnan(i):
			logger.error("ra=%r (%r,%r), rp=%r (%r,%r), ravg=%r (%r)" % (ra, dg(ra), ekinrad(ra), rp, dg(rp), ekinrad(rp), ravg, dg(ravg)))
			import pdb; pdb.set_trace()
			raise Exception("integral is nan for E=%r L=%r" % (E, L))  
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
	
class ProfileModel2C(ProfileModel1C):
	def __init__(self, light_profile, dm_profile):
		ProfileModel1C.__init__(self, light_profile) 
		self.dm_profile = dm_profile
		self.dm_profile_fast = self.dm_profile.fast()
		self.fast_ = gdfast.ProfileModel2C(self.light_profile_fast, self.dm_profile_fast)
		 
	def vmax(self):
		return sqrt(-2*self.potentialr(1e-6))
		
	def fast(self):
		return self.fast_
		
	def enclosed_mass(self, r1, r2):
		return self.dm_profile.cumdensityr(r1, r2) +  self.light_profile.cumdensityr(r1, r2)
	
	def getprofiles(self):
		return [self.light_profile, self.dm_profile]
	
	def densityr(self, r):
		return self.light_profile.densityr(r) +  self.dm_profile.densityr(r) 
	
	def potentialr(self, r):
		return self.light_profile.potentialr(r) + self.dm_profile.potentialr(r)
	
	def dphidr(self, r):
		return self.light_profile.dphidr(r) + self.dm_profile.dphidr(r)
	
	def vcirc(self, r):
		return sqrt(r * self.dphidr(r))
	
class LightModel1C(ProfileModel1C):
	def __init__(self, light_profile, distance):
		ProfileModel1C.__init__(self, light_profile) 
		self.distance = distance
		
	def sample_r(self, N=None):
		return self.light_profile.sample_r(N=N)
		
	def cumdensityR(self, R1, R2, M=None):
		return self.light_profile.cumdensityR(R1, R2, M=M)
		
	def cumdensityr(self, r1, r2, M=None):
		return self.light_profile.cumdensityr(r1, r2, M)
		
	def densityR(self, R, M=None):
		return self.light_profile.densityR(R, M=M)
	
	def kpc_to_arcsec(self, d): 
		return d * 1.0/(self.distance) / (pi/180/60/60)
	
	def pc_to_arcsec(self, d): 
		return d * 1e3/(self.distance) / (pi/180/60/60)
	
	def arcsec_to_kpc(self, R): 
		return R / (1.0/(self.distance) / (pi/180/60/60))
	
	def arcsec_to_pc(self, d): 
		return R / (1e3/(self.distance) / (pi/180/60/60))
	
"""class Galaxy2C_constant_anisotropy(Galaxy1C_constant_anisotropy):
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
		"""
			
class Galaxy(object):
	def __init__(self, light_model, profile_model, distance):
		self.light_model = light_model
		self.profile_model = profile_model
		self.distance = distance
	#def vmax(self):
	#	return sqrt(-2*self.potentialr(1e-9))
		
	#def jeansfast(self):
	#	return gdfast.JeansAnisotropicConstant(self.stellar_profilef, self.null_profilef, self.beta, 10.)
	
	#def fast(self):
	#	return gdfast.Galaxy1C_constant_anisotropy(self.stellar_profilef, self.distance, self.beta)
	
	#def jeans(self):
	#	return JeansAnisotropicConstant(self.stellar_profile, self.getprofiles(), self.beta)
		
	def kpc_to_arcsec(self, d): 
		return d * 1.0/(self.distance) / (pi/180/60/60)
	
	def pc_to_arcsec(self, d): 
		return d * 1e3/(self.distance) / (pi/180/60/60)
	
	def arcsec_to_kpc(self, R): 
		return R / (1.0/(self.distance) / (pi/180/60/60))
	
	def arcsec_to_pc(self, d): 
		return R / (1e3/(self.distance) / (pi/180/60/60))
	
	#def getprofiles(self):
	#	return self.profile_model.getprofiles()
	
	#def densityr(self, r):
	#	return self.profile_model.densityr(r)
	
	def potentialr(self, r):
		return self.profile_model.potentialr(r)
	
	def dphidr(self, r):
		return self.profile_model.dphidr(r)
	
	def vcirc(self, r):
		return self.profile_model.vcirc(r)
		

class Galaxy_constant_anisotropy(Galaxy):
	def __init__(self, light_model, profile_model, distance, beta):
		Galaxy.__init__(self, light_model, profile_model, distance)
		self.beta = beta
		#self.stellar_profilef = self.stellar_profile.fast() # keep reference
		#self.null_profilef = gdfast.NullProfile()
		#self._fast = self.fast()
		
	def fL(self, L):
		return L**(-2*self.beta)
		
	def jeans(self):
		return JeansAnisotropicConstant(self.light_model, self.profile_model, self.beta)
	
	#def m2(self, r):
		
	
	
class Galaxy_double_anisotropy(Galaxy):
	def __init__(self, light_model, profile_model, distance, beta0, betainf, beta_rs):
		Galaxy.__init__(self, light_model, profile_model, distance)
		self.beta0 = beta0
		self.betainf = betainf
		self.beta_rs = beta_rs
		vcirc = self.vcirc(beta_rs)
		self.L0 = vcirc * beta_rs 
		
	def fL(self, L):
		power = 4.
		return (1+L**(power)/(2*self.L0**(power)))**(2/power*(-self.betainf+self.beta0))*(L**-(2*self.beta0))
	
	
class OldStuff:
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
		
			
			
"""class Galaxy2C_constant_anisotropy(Galaxy1C_constant_anisotropy):
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
"""
			