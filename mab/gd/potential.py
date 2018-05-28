# -*- coding: utf-8 -*-

"""Contain some common profiles
"""



#from constants import *
from numpy import *
import mab.constants
import mab.astrounits
import numpy
import math
from math import pi
import mab.cosmology
from mab.gd import gdfast
import sys
from scipy.optimize import fsolve, fmin, brentq
from scipy.integrate import quad
import scipy.optimize

G = mab.constants.G.asNumber(mab.astrounits.KM**2/mab.astrounits.S**2 *mab.astrounits.KPC/mab.astrounits.MSOL)



class CommandMassAtLogr(object):
	def __init__(self, profile):
		self.profile = profile
		
	def run(self, args, opts, scope):
		#print args
		for arg in args[1:]:
			logr = float(arg)
			print "logr = %f " % logr, "r = %f" % 10**logr
			mass = self.profile.enclosed_mass(10**logr)
			print "mass = %f, %e" %(mass,mass), "log mass = %f" % log10(mass)
		return len(args)
		
		
class Potential(object):
	"""Base class for potentials"""
	def findLmax(E):
		def f(r, E=E):
			r = abs(r)
			return -r**2*(2*(E-potentialr(r)))
		rLmax, Lmaxsq = fmin(f, 0.1,full_output=True, disp=False)[:2]
		#Lmax = extra[0]
		#print Lmax
		return abs(rLmax[0]), sqrt(abs(Lmaxsq))
	
	def calc_r3(self):
		def f(r):
			r = abs(r)
			s = self.logslope(r)[0]
			return (s+3.0)**2
		r3 = abs(scipy.optimize.fmin(f, self.scale*10, disp=False)[0])
		#print "r3", r3
		return r3
		
	def calc_rhalf(self):
		def f(r):
			r = abs(r)
			s = self.enclosed_mass(r)[0]/self.M
			#print r, s
			return (s-0.5)**2
		rhalf = abs(scipy.optimize.fmin(f, self.scale*10, disp=False)[0])
		#print "r3", r3
		return rhalf
			
		
	def vcirc(self, r):
		return sqrt(r * self.dphidr(r))

	def enclosed_mass_num(self, r):
		return vectorize(self._enclosed_mass_num)(r)
		
	def enclosed_binding_energy(self, r):
		return vectorize(self._enclosed_binding_energy)(r)
	def _enclosed_binding_energy(self, r):
		def dW(x):
			return -4 * pi * self.G * self.densityr(x) * self.enclosed_mass(x) * x
			#return 4 * pi * self.densityr(x) * self.potentialr(x) * x**2
		W, err = quad(dW, 0, r)
		return W
	def denclosed_binding_energy(self, r):
		return vectorize(self._denclosed_binding_energy)(r)
	def _denclosed_binding_energy(self, r):
		return -4 * pi * self.G * self.densityr(r) * self.enclosed_mass(r) * r
		
	
	def ___enclosedmass_num(self, r):
		def dM(x):
			rp = 1/x 
			return 4 * pi * rp**2 * self.densityr(rp) / x**2
		x = 1/r
		M, err = quad(dM, x, inf)
		return M 
	
	def _enclosed_mass_num(self, r):
		def dM(rp):
			return 4 * pi * rp**2 * self.densityr(rp)
		M, err = quad(dM, 0, r)
		return M 
	
	def enclosed_mass(self, r):
		return self.enclosed_mass_num(r)
	 
	def total_mass(self):
		return self.enclosed_mass_num(inf)
	
	def cumdensityr(self, r1, r2, M=None):
		if M is None:
			M = self.M
		density, error = quad(lambda r: 4*pi*self.densityr(r)*r**2 * M / self.M, r1, r2)
		return density
	
	def cumdensityR(self, R1, R2, M=None):
		if M is None:
			M = self.M
		density, error = quad(lambda R: 2*pi*self.densityR(R)*M/self.M*R, R1, R2)
		return density
	
			
	def sample_bla(self, N=1, rmax=None):
		ok = False
		rs = []
		vs = []
		#Emax = -0.0001*(self.G*self.M/self.a)
		#fmax = self.f(-(self.G*self.M/self.a)*0.99) * 150
		E0 = self.G*self.M/self.scale
		dE = E0/1000
		Es = arange(-E0*0.9, -dE/2, dE)
		fmax = max(self.f(Es))
		y = self.f(Es)
		print y, argmax(y)
		print self.f(-E0)
		from scipy.integrate import quad
		print quad(self.f, -E0, 0)
		#print Nmax
		while not ok:
			#print "dsa"
			r = self.sampler(N=1, rmax=rmax)[0]
			#print "r =", r
			Epot = self.potentialr(r)
			#print Epot
			#Emin = self.potentialr(0)
			vmax = sqrt(-2*Epot)
			#print zip(r, vmax)
			vtest = arange(0, vmax, 0.01)
			fmax = max(self.f(0.5*vtest**2+Epot)*abs(vtest)**2) * 2
			accepted_v = False
			while not accepted_v:
				# first draw a random v
				v = numpy.random.random() * vmax
				E = 0.5 * v**2 + Epot
				# what is the prop of this v
				fsample = self.f(E) * abs(v)**2
				# should we accept this v?
				frandom = numpy.random.random() * fmax
				#print frandom, fsample
				
				accepted_v = frandom < fsample
			print ".",
			vs.append(v)
			rs.append(r)
				#print sum(okindices)
				#print frandom
				#print fsample
				#rs = numpy.append(rs, r[okindices])
				#vs = numpy.append(vs, v[okindices])
			#if len(vs) > 0:
			#print len(vs)
			#lalala
			if len(rs) >= N:
				ok = True
		return numpy.array(rs[:N]), numpy.array(vs[:N])

			
	def sample_rv(self, N=1, beta=0, rmax=None):
		ok = False
		rs = []
		vrs = []
		vphis = []
		vthetas = []
		#Emax = -0.0001*(self.G*self.M/self.a)
		#fmax = self.f(-(self.G*self.M/self.a)*0.99) * 150
		E0 = self.G*self.M/self.scale
		dE = E0/1000
		Es = arange(-E0*0.9, -dE/2, dE)
		fmax = max(self.fE(Es))
		y = self.fE(Es)
		#print y, argmax(y)
		#print self.f(-E0)
		from scipy.integrate import quad
		print quad(self.fE, -E0, 0)
		#print Nmax
		while not ok:
			#print "dsa"
			r = self.sampler(N=1, rmax=rmax)[0]#  * 0 + 1.
			#print "r =", r
			Epot = self.potentialr(r)
			#print Epot
			#Emin = self.potentialr(0)
			Lmax = r*sqrt(-2*Epot)
			#vtmax =
			vmax = sqrt(-2*Epot)
			#print zip(r, vmax)
			vtest = arange(0.0, vmax+0.01/2, 0.01)
			fEmax = max(self.fE(0.5*vtest**2+Epot)) * 2
			Ltest = vtest * r
			fEmax = max(self.fE(0.5*vtest**2+Epot)*Ltest**(-2*beta)) * 2
			#print fEmax
			#fEmax = max(self.f(0.5*vtest**2+Epot)*abs(vtest)**1*((vtest/r)**(-2*beta))) * 2
			#if beta < 0:
			#	fEmax *= Lmax**(-2*beta)
			accepted_v = False
			#print Lmax, vmax, r
			while not accepted_v:
				# first draw a random v
				#L = numpy.random.random() * Lmax * 0.1 + Lmax * 0.9
				vphi = (numpy.random.random() * 2-1)*vmax
				#L = vt * r
				#L = numpy.random.random() * Lmax
				#vrmax = sqrt(-2*Epot-L**2/r**2) 
				vrmax = sqrt(-2*Epot-vphi**2)
				vr = (numpy.random.random()*2-1) * vrmax
				
				vthetamax = sqrt(-2*Epot-vphi**2-vr**2)
				vtheta = (numpy.random.random()*2-1) * vthetamax
				
				vt = sqrt(vphi**2+vtheta**2)
				L = vt*r
				#vt = L/r
				v = sqrt(vr**2+vt**2)
				assert v <= vmax
				#print L, Lmax, vmax, vt, vr, r, vmax*r, Lmax
				assert L <= Lmax
				E = 0.5 * v**2 + Epot
				assert E <= 0
				# what is the prop of this v
				fsample = self.fE(E) * L**(-2*beta)
				if fsample > fEmax:
					print fsample, fEmax
				assert fsample <= fEmax
				# should we accept this v?
				frandom = numpy.random.random() * fEmax
				#print frandom, fsample, fEmax,  L**(-2*beta), self.f(E) * abs(v)**2, Lmax**(-2*beta)
				#print " *)", fEmax, frandom, fsample
				
				accepted_v = frandom < fsample
				#print "*",
			print ".",
			sys.stdout.flush()
			#print "=" * 70
			vrs.append(vr)
			vphis.append(vphi)
			vthetas.append(vtheta)
			rs.append(r)
				#print sum(okindices)
				#print frandom
				#print fsample
				#rs = numpy.append(rs, r[okindices])
				#vs = numpy.append(vs, v[okindices])
			#if len(vs) > 0:
			#print len(vs)
			#lalala
			if len(rs) >= N:
				ok = True
		return numpy.array(rs[:N]), numpy.array(vrs[:N]), numpy.array(vphis[:N]), numpy.array(vthetas[:N])
	
	def f1a(self, E):
		e = -E
		ec = e*self.scale/(self.G*self.M)
		return 3 * ec**2/(4*pi**3*self.G*self.M*self.scale)
	
	def f1betaminhalf(self, E):
		e = -E
		ec = e*self.scale/(self.G*self.M)
		#return 3 * ec**2/(4*pi**3*self.G*self.M*self.scale)
		#print ec
		return 1./(4*pi**3*(self.G*self.M*self.scale)**2) * ( 20*ec**3/(1-ec)**2 + 20*ec**4/(1-ec)**3 + 6*ec**5/(1-ec)**4 )
	
	def fELbetaminhalf(self, E, L):
		return self.f1betaminhalf(E) * L
	
	def f1(self, E, beta):
		if beta == -0.5:
			return self.f1betaminhalf(E)
		if beta == 0.0:
			return self.fE(E)
	
	def samplev(self, rs, beta, fLz=lambda Lz,Lmax: 1):
		done = False
		
		phis = []
		thetas = []
		
		vrs = []
		vphis = []
		vthetas = []
		
		i = 0
		if beta == 0.0:
			f1 = self.fE
		elif beta == -0.5:
			f1 = self.f1betaminhalf
		else:
			raise Exception, "unsupported beta %f" % beta 
		for r in rs:
			Epot = self.potentialr(r)
			vmax = sqrt(-2*Epot)
			Lmax = vmax * r
			dv = 0.01
			vtest = arange(0.01, vmax+dv/2, dv)
			Ltest = vtest * r
			
			phi = random.random() * 2 * pi
			costheta = random.random() * 2 - 1
			theta = arccos(costheta)
			
			cosphi = cos(phi)
			sinphi = sin(phi)
			#sintheta = sqrt(1-costheta**2)
			sintheta = sin(theta)
			
			
			x = r * sintheta * cosphi
			y = r * sintheta * sinphi
			z = r * costheta
			rho = sqrt(x**2+y**2)
			
			Lzmax = rho * vmax
			vtest = arange(0.01, vmax+dv/2, dv)
			Lztest = vtest * rho
			
			
			assert all(fLz(Lztest, Lzmax) > 0)
			#print fLz(-Lztest, Lzmax)
			assert all(fLz(-Lztest, Lzmax) > 0)
			fmax1 = max(f1(0.5*vtest**2+Epot)*Ltest**(-2*beta)*fLz(Lztest, Lzmax)) * 2
			fmax2 = max(f1(0.5*vtest**2+Epot)*Ltest**(-2*beta)*fLz(-Lztest, Lzmax)) * 2
			fmax = max(fmax1, fmax2)
			accepted = False
			while not accepted:
				# sample uniformly (important!)
				vr     = (numpy.random.random()*2-1) * vmax
				vphi   = (numpy.random.random()*2-1) * vmax
				vtheta = (numpy.random.random()*2-1) * vmax
				vsquare = vr**2+vphi**2+vtheta**2
				if vsquare > vmax**2: # retry
					continue
				
				vt = sqrt(vphi**2+vtheta**2)
				v = sqrt(vr**2+vt**2)
				L = vt*r
				Lz = vphi * rho
				E = 0.5 * v**2 + Epot
				
				assert v <= vmax
				assert E <= 0
				
				# rejection sampling
				fsample = f1(E) * L**(-2*beta) * fLz(Lz, Lzmax)
				frandom = numpy.random.random() * fmax
				accepted = frandom < fsample
				#print accepted, fLz(Lz, Lmax), Lz
				
				# sanity check
				if fsample > fmax:
					print fsample, fmax
				assert fsample <= fmax
			#assert Lz > 0
			sys.stdout.write(".")
			sys.stdout.flush()
			if len(vrs) > 0 and len(vrs) % 1000 == 0:
				p = len(vrs) * 100 / len(rs)
				sys.stdout.write("\n%d samples, %.2f%%\n" % (len(vrs), p))
			vrs.append(vr)
			vphis.append(vphi)
			vthetas.append(vtheta)
			phis.append(phi)
			thetas.append(theta)
			
		return numpy.array(phis), numpy.array(thetas), numpy.array(vrs), numpy.array(vphis), numpy.array(vthetas)
				
				
				
				
	
	def sample_rv2(self, N=1, beta=0, rmax=None):
		ok = False
		rs = []
		vrs = []
		vphis = []
		vthetas = []
		#Emax = -0.0001*(self.G*self.M/self.a)
		#fmax = self.f(-(self.G*self.M/self.a)*0.99) * 150
		#E0 = self.G*self.M/self.scale * 0.9999
		E0 = self.potentialr(0) * 0.999
		dE = E0/1000
		#Es = arange(-E0*0.9, -dE/2, dE)
		#fmax = max(self.f(Es))
		#y = self.f(Es)
		#print y, argmax(y)
		#print self.f(-E0)
		#from scipy.integrate import quad
		#print quad(self.f1, -E0, -1e-5)
		#print Nmax
		while not ok:
			#print "dsa"
			r = self.sampler(N=1, rmax=rmax)[0]#  * 0 + 1.
			#print r
			#print "r =", r
			Epot = self.potentialr(r)
			#print Epot
			#Emin = self.potentialr(0)
			Lmax = r*sqrt(-2*Epot)
			#vtmax =
			vmax = sqrt(-2*Epot)
			#print zip(r, vmax)
			vtest = arange(0.01, vmax+0.01/2, 0.01)
			#fEmax = max(self.f1(0.5*vtest**2+Epot)) * 2
			Ltest = vtest * r
			fEmax = max(self.f1(0.5*vtest**2+Epot)*Ltest**(-2*beta)) * 2
			#print fEmax
			#fEmax = max(self.f(0.5*vtest**2+Epot)*abs(vtest)**1*((vtest/r)**(-2*beta))) * 2
			#if beta < 0:
			#	fEmax *= Lmax**(-2*beta)
			accepted_v = False
			#print Lmax, vmax, r
			#print fEmax
			while not accepted_v:
				# first draw a random v
				#L = numpy.random.random() * Lmax * 0.1 + Lmax * 0.9
				#vrmax = sqrt(-2*Epot-vphi**2)
				vrmax = vmax
				vphimax = vmax
				vthetamax = vmax
				vr = (numpy.random.random()*2-1) * vrmax
				vphi = (numpy.random.random() * 2-1)*vphimax
				vtheta = (numpy.random.random()*2-1) * vthetamax
				vsquare = vr**2+vphi**2+vtheta**2
				if vsquare > vmax**2:
					continue
				# 
				#vphimax = sqrt(-2*Epot-vr**2)
				#L = vt * r
				#L = numpy.random.random() * Lmax
				#vrmax = sqrt(-2*Epot-L**2/r**2) 
				
				#vthetamax = sqrt(-2*Epot-vphi**2-vr**2)
				
				vt = sqrt(vphi**2+vtheta**2)
				L = vt*r
				#vt = L/r
				v = sqrt(vr**2+vt**2)
				assert v <= vmax
				#print L, Lmax, vmax, vt, vr, r, vmax*r, Lmax
				assert L <= Lmax
				E = 0.5 * v**2 + Epot
				assert E <= 0
				# what is the prop of this v
				fsample = self.f1(E) * L**(-2*beta)
				if fsample > fEmax:
					print fsample, fEmax
				assert fsample <= fEmax
				# should we accept this v?
				frandom = numpy.random.random() * fEmax
				#print frandom, fsample, fEmax,  L**(-2*beta), self.f(E) * abs(v)**2, Lmax**(-2*beta)
				#print " *)", fEmax, frandom, fsample
				
				accepted_v = frandom < fsample
				#print "*",
			print ".",
			sys.stdout.flush()
			#print "=" * 70
			vrs.append(vr)
			vphis.append(vphi)
			vthetas.append(vtheta)
			rs.append(r)
				#print sum(okindices)
				#print frandom
				#print fsample
				#rs = numpy.append(rs, r[okindices])
				#vs = numpy.append(vs, v[okindices])
			#if len(vs) > 0:
			#print len(vs)
			#lalala
			if len(rs) >= N:
				ok = True
		return numpy.array(rs[:N]), numpy.array(vrs[:N]), numpy.array(vphis[:N]), numpy.array(vthetas[:N])
		
			
class ProjectedExponential(Potential):
	def __init__(self, M, scale, G=G):
		self.M = M
		self.scale = scale
		self.rho0 = 1.
		#sprint self.scale
		Mcurrent = self.enclosed_mass(inf)
		self.rho0 = M/Mcurrent
		self._fast = mab.gd.gdfast.ProjectedExponential(M, scale, G)
		
	def fast(self):
		return self._fast
		
		#2 \[Pi] Rs (Rs - E^(-(r/Rs)) (r + Rs))
		
	def densityR(self, r, M=None):
		if M is None:
			scaled_mass = 1
		else:
			scaled_mass = M/self.M
			
		return exp(-r/self.scale) * self.rho0 * scaled_mass
	
	def densityr(self, r, M=None):
		if M is None:
			scaled_mass = 1
		else:
			scaled_mass = M/self.M
		# kn is the (modified) bessel of second kind (integer order = 0)
		return self.rho0 * scipy.special.kn(0., r/self.scale) / (self.scale * pi) * scaled_mass
		
	def logslope(self, r):
		return -r/self.scale * scipy.special.kn(1., r/self.scale) / scipy.special.kn(0., r/self.scale)
	
	def potentialr(self, r):
		return 0
	
	def dphidr(self, r):
		return 0
		

class Plummer(Potential):
	"""Implements the Plummer profile
	 
Constructor::
	
	M - mass parameter
	b - scale parameter
		"""
	def __init__(self, M, b, G=G):
		self.M = M
		self.b = b
		self.G = G
		self.scale = b
		
	@classmethod
	def arguments(cls):
		return [("M", "total mass - solar masses", 1e6), ("b", "scale parameter - kpc", 1)]
		
	def densityR(self, R, M=None):
		"""Density at projected radius R"""
		if M is None:
			M = self.M
		return M * self.b**2 / (pi * (self.b**2+R**2)**2)
	def potentialr(self, r):
		return - self.G * self.M / numpy.sqrt(r**2 + self.b**2)
	
	def potentialxy(self, x, y):
		r = sqrt(x**2+y**2)
		return - self.G * self.M / numpy.sqrt(r**2 + self.b**2)
		
	def logslope(self, r):
		return -5*r**2/(self.b**2+r**2)
		
	
	def densityr(self, r, M=None):
		if M == None:
			M = self.M
		return 3 * self.M * self.b ** 2 / (4*pi) / (r**2 + self.b**2) ** (5./2) * M / self.M
	
	def vcircr(self, r):
		return sqrt(self.G * self.M * r ** 2 / (r**2 + self.b**2) ** (3./2))
	
	def dphidr(self, r):
		return self.G * self.M * r / (r**2 + self.b**2) ** (3./2)
	
	def c_code_dphidr(self, r="r"):
		variables = {}
		variables["G"] = self.G
		variables["M"] = self.M
		variables["b_square"] = self.b**2
		variables["r"] = r
		return "%(G)r * %(M)r * %(r)s / pow(pow(%(r)s, 2) + %(b_square)r, 3./2)" % variables
	
	def c_code_potential(self, r="r"):
		variables = {}
		variables["G"] = self.G
		variables["M"] = self.M
		variables["b_square"] = self.b**2
		variables["r"] = r
		return "-%(G)r * %(M)r / sqrt(pow(%(r)s, 2) + %(b_square)r)" % variables
	
	def sample_r(self, rmax=1e10, N=1):
		ok = False
		r = []
		if N is None:
			single_value = True
			N = 1
		else:
			single_value = False
		M = N/2
		if N == 1:
			M = 1
		while not ok:
			u = numpy.random.random(M)
			alpha = 2./3
			r = numpy.append(r, [self.b * numpy.sqrt(u**alpha / (1-u**alpha))])
			if rmax:
				r = r[r<rmax]
			#print M, u, r, rmax
			if len(r) >= N:
				ok = True
		if single_value:
			return r[0]
		else:
			return r[:N]
		
	def samplerReject(self, rmax, N=1):
		# rejection sampling
		ok = False
		samples = numpy.zeros(N) * 1.0
		i = 0
		# max of rho(r)*r**2 is at this r
		r_max_density = math.sqrt(2./3) * self.b
		M = self.densityr(r_max_density) * r_max_density**2
		#print "M =", M, r_max_density
		#r = r_max_density + 1e-2
		#print self.densityr(r) * r**2
		#r = r_max_density - 1e-2
		#print self.densityr(r) * r**2
		#M = 10
		while i < N:
			#x = random.random() * r
			x = numpy.random.random() * rmax
			u = numpy.random.random()
			fx = self.densityr(x)*x**2
			if u < (fx/M):
				samples[i] = x
				i += 1
		if N == 1:
			return samples[0]
		else:
			return samples
		
	def fE(self, E):
		phimax = self.G * self.M/self.scale
		return (-E)**(7./2) * 9/2 / (phimax**(9./2))
	
	def fast(self):
		return gdfast.Plummer(self.M, self.scale, self.G)

class Isochrone(Potential):
	def __init__(self, M, b, G=G):
		self.M = M
		self.b = b
		self.G = G
		self.scale = b
		
	def densityR(self, R, I0=1.):
		raise "unknown"
	
	def potentialr(self, r):
		return - self.G * self.M / (self.b + numpy.sqrt(r**2 + self.b**2))
	
	def potentialxy(self, x, y):
		r = sqrt(x**2+y**2)
		return self.potentialr(r)
	
	def densityr(self, r):
		a = sqrt(self.b**2+r**2)
		return self.M * ( 3*(self.b+a)*a**2 - r**2*(self.b+3*a) ) / ( 4 * pi * (self.b + a)**3*a**3)
	
	def vcircr(self, r):
		a = sqrt(self.b**2+r**2)
		return sqrt( self.G * self.M * r ** 2 / ((self.b+a)**2*a) )
	
	def dphidr(self, r):
		raise "unknown"
	
	def c_code_dphidr(self, r="r"):
		raise "unknown"
	
	def c_code_potential(self, r="r"):
		raise "unknown"
	
	def sampler(self, rmax=1e10, N=1):
		raise "unknown"
		
	def samplerReject(self, rmax, N=1):
		# rejection sampling
		ok = False
		samples = numpy.zeros(N) * 1.0
		i = 0
		# max of rho(r)*r**2 is at this r
		r_max_density = math.sqrt(2./3) * self.b
		M = self.densityr(r_max_density) * r_max_density**2
		#print "M =", M, r_max_density
		#r = r_max_density + 1e-2
		#print self.densityr(r) * r**2
		#r = r_max_density - 1e-2
		#print self.densityr(r) * r**2
		#M = 10
		while i < N:
			#x = random.random() * r
			x = numpy.random.random() * rmax
			u = numpy.random.random()
			fx = self.densityr(x)*x**2
			if u < (fx/M):
				samples[i] = x
				i += 1
		if N == 1:
			return samples[0]
		else:
			return samples
		
	def fE(self, E):
		raise "unknown"
	
	def fast(self):
		return gdfast.Isochrone(self.M, self.scale, self.G)
		
class LogarithmicProfile(Potential):
	def __init__(self, vcirc, G=G):
		self.vcirc = vcirc
		self.G = G
		
	def densityR(self, R, I0=1.):
		raise "unknown"
	
	def potentialr(self, r):
		return self.vcirc**2 * log(r)
	
	def potentialxy(self, x, y):
		r = sqrt(x**2+y**2)
		return self.potentialr(r)
	
	def densityr(self, r):
		raise "unknown"
	
	def vcircr(self, r):
		raise "unknown"
	
	def dphidr(self, r):
		raise "unknown"
	
	def c_code_dphidr(self, r="r"):
		raise "unknown"
	
	def c_code_potential(self, r="r"):
		raise "unknown"
	
	def sampler(self, rmax=1e10, N=1):
		raise "unknown"
		
	def samplerReject(self, rmax, N=1):
		raise "unknown"
		
	def fE(self, E):
		raise "unknown"
	
	def fast(self):
		return gdfast.LogarithmicProfile(self.vcirc, self.G)
		
class Hernquist(Potential):
	def __init__(self, M, a, G=G):
		self.M = M
		self.a = a
		self.rho0 = self.M / (2*pi*self.a**3)
		self.G = G
		self.scale = a
		
	def potentialr(self, r):
		return - 4 * pi * self.G * self.rho0 * self.a**2 /(2*(1+r/self.a))
	
	def densityr(self, r, M=None):
		if M is None:
			M = self.M
		m = r/self.a
		return self.rho0 / (m * (1+m)**3) * M / self.M
	
	def vcircr(self, r):
		return r * self.dphidr(r)
	
	def dphidr(self, r):
		#return self.G * self.M * r / (r**2 + self.b**2) ** (3./2)
		return 4 * pi * self.G * self.rho0 * self.a / (2*(1+r/self.a)**2)
	
	def c_code_dphidr(self, r="r"):
		variables = {}
		variables["pi"] = pi
		variables["G"] = self.G
		variables["rho0"] = self.rho0
		variables["a"] = self.a
		variables["r"] = r
		return "4 * %(pi)r * %(G)r * %(rho0)r * %(a)r / (2*pow((1+%(r)s/%(a)r), 2))" % variables
	
	def sample_r(self, rmax=1e10, N=1):
		return self.sampler(rmax, N=N)
	def sampler(self, rmax=1e10, N=1):
		if N is None:
			N = 1
		ok = False
		r = []
		M = N/2
		if N == 1:
			M = 1
		while not ok:
			u = numpy.random.random(M)
			alpha = 2./3
			r = numpy.append(r, (-self.a * numpy.sqrt(u) - self.a * u) / (u-1))
			#r = numpy.append(r, self.b * numpy.sqrt(u**alpha / (1-u**alpha)))
			if rmax:
				r = r[r<rmax]
			if len(r) >= N:
				ok = True
		if N == 1:
			return r[0]
		else:
			return r[:N]
	
	def fE(self, E):
		ec = -E * self.a / (self.G*self.M)
		return 1 / (sqrt(2)*(2*pi)**3 * (self.G * self.M * self.a) ** (3./2)) * numpy.sqrt(ec)/(1-ec)**2 * \
			((1-2*ec)*(8*ec**2-8*ec-3) + 3*numpy.arcsin(numpy.sqrt(ec)) / numpy.sqrt(ec*(1-ec)))
			
	def g(self, E):
		A = self.G * self.M / (self.a*abs(E))
		return (4*pi)**2*self.a**3 * sqrt(2*abs(E)) * \
			(numpy.sqrt(A-1)*(1./8*A**2-5./12*A-1./3) + 1./8 * A * (A**2-4*A+8)*numpy.arccos(A**-0.5))
			
	def N(self, E):
		return self.fE(E) * self.g(E)
			
	def sampleE(self, N=1):
		ok = False
		Es = []
		Emax = -0.0001*(self.G*self.M/self.a)
		Nmax = self.N(Emax)
		#print Nmax
		while not ok:
			# draw random energy
			E = -numpy.random.random(N/2) * (self.G*self.M/self.a)
			Nrandom = numpy.random.random(N/2) * Nmax
			Nsample = self.N(E)
			Edrawn = E[Nrandom<Nsample]
			Es = numpy.append(Es, Edrawn)
			#print len(Edrawn)
			#print Nrandom, Nsample
			#r = r[r<rmax]
			if len(Es) >= N:
				ok = True
		return Es[:N]
		
	def sample_rvL(self, N=1):
		ok = False
		rs = []
		vs = []
		Emax = -0.0001*(self.G*self.M/self.a)
		Nmax = self.N(Emax)
		#print Nmax
		while not ok:
			# draw random energy
			E = -numpy.random.random(N/2) * (self.G*self.M/self.a)
			Nrandom = numpy.random.random(N/2) * Nmax
			Nsample = self.N(E)
			Edrawn = E[Nrandom<Nsample]
			Es = numpy.append(Es, Edrawn)
			#print len(Edrawn)
			#print Nrandom, Nsample
			#r = r[r<rmax]
			if len(Es) >= N:
				ok = True
		return Es[:N]
	 
	def fast(self):
		return gdfast.Hernquist(self.M, self.scale, self.G)
	
	 
		
class NFWTriApprox(object):
	def __init__(self, rhoc, rs, a, b, c, G=G):
		self.M = M
		self.rs = rs
		self.a = a
		self.b = b
		self.c = c
		self.rhoc = rhoc
		self.G = G
		self.scale = rs
		
	#def potential(self, x, y, z):
	#	return - 4 * pi * self.G * self.rho0 * self.a**2 /(2*(1+r/self.a))
	
	#def densityr(self, r):
	#	m = r/self.a
	#	return self.rho0 / (m * (1+m)**3)
	
	#def vcircr(self, r):
	#	return r * self.dphidr(r)
	
	def dphi(self, x, y, z):
		r = sqrt(x*x+y*y+z*z)
		ra = self.rs
		re = sqrt(x**2/self.a**2 + y**2/self.b**2 + z**2/self.c**2)
		rtilde = (ra+r)*re / (ra+re)

		t1 = -4.0*pi*self.G*self.rhoc*self.rs**3/rtilde**2
		t2 = log(1.0 + rtilde/dm_profile_parameter)
		t3 = (rtilde/self.rs)/(1.0+rtilde/self.rs)
		
		dphidr = -t1*(t2 - t3)

		drtildedx = 1 / (ra+re)**2 * (1/re * x / self.a**2 * (ra+r) * ra + 1/r * x * (ra+re)*re)
		drtildedy = 1 / (ra+re)**2 * (1/re * y / self.b**2 * (ra+r) * ra + 1/r * y * (ra+re)*re)
		drtildedz = 1 / (ra+re)**2 * (1/re * z / self.c**2 * (ra+r) * ra + 1/r * z * (ra+re)*re)
		return	dphidr * drtildedx,\
				dphidr * drtildedy,\
				dphidr * drtildedz
		
		#return self.G * self.M * r / (r**2 + self.b**2) ** (3./2)
		#return 4 * pi * self.G * self.rho0 * self.a / (2*(1+r/self.a)**2)
	
	
def concentration2rs(M200, concentration, cosmology=mab.cosmology.wmap3y):
	rhoc = (cosmology.rho_crit * cosmology.h**2).asNumber(mab.astrounits.MSOL/mab.astrounits.KPC**3)
	r200 = (M200/200*3/(4*pi)/rhoc)**(1./3)
	#print "r200", self.r200
	return r200/concentration
	




class NFW(Potential):
	def __init__(self, M200, rs, G=G, cosmology=mab.cosmology.wmap3y):
		Potential.__init__(self)
		self.M200 = M200
		self.M = self.M200
		self.rs = rs
		self.scale = rs
		self.cosmology = cosmology 
		rhoc = (cosmology.rho_crit * cosmology.h**2).asNumber(mab.astrounits.MSOL/mab.astrounits.KM**3)
		#print rhoc
		rhoc = (cosmology.rho_crit * cosmology.h**2).asNumber(mab.astrounits.MSOL/mab.astrounits.KPC**3)
		#print rhoc
		self.rhoc = rhoc
		self.r200 = (self.M200/200*3/(4*pi)/rhoc)**(1./3)
		#print "r200", self.r200
		self.c = self.r200/self.rs
		#print "c", self.c
		x = self.r200 / self.rs
		#self.rho0 = rhoc * x * (1+x)**2
		self.rho0 = self.M200 / (4*pi*rs**3*(log(1+x)-x/(1+x)))
		
		#print "1)", 4 * pi * self.rho0 * rs ** 3 * (log(1+x) - x / (1+x))
		#print "2)", 200 * rhoc * 4/3 * pi * self.r200**3
		self.G = G
		self.fast_ = gdfast.NFW(self.M200, self.scale, self.G, self.rhoc)
		
	def run(self, *args, **kwargs):
		print "concentation", self.c
		
	@staticmethod
	def from_r3(M3, r3, slope):
		d = NFW(1., 1.)
		def f(rs):
			d.rs = abs(rs)
			d.scale = abs(rs)
			#print rs, d.logslope(r3), slope
			#import pdb
			#pdb.set_trace()
			return (d.logslope(r3)[0]-slope)**2
			
		rs = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		d.scale = rs
		d.rs = rs
		#d.normalize(r3, M3)
		Menc = d.enclosed_mass(r3)
		d.rho0 *= M3/Menc
		d.fast_ = gdfast.NFW(d.M200, d.scale, d.G, d.rhoc)
		
		#self.rho0 = M1kpc/(4*pi*rs**3*(log(1+x)-x/(1+x)))
		#self.scale = self.rs = rs
		#self.G = G
		#print self.enclosedmass(1.), M1kpc
		#assert abs(self.enclosedmass(1.) - M1kpc) < 0.01
		
		#r200 = exp(fmin(mass_diff, 10, args=(rhoc, self.rho0, rs), disp=1)[0])
		#print rhoc, self.rho0
		logr200 = (fsolve(mass_diff, 2, args=(d.rhoc, d.rho0, rs)))[0]
		
		r200 = exp(logr200)
		d.r200 = r200
		d.M200 = 200*d.rhoc*4*pi/3*r200**3
		d.M = d.M200
		d.c = r200/rs
		#print d.M200, d.scale, d.G, d.rhoc
		d.fast_ = gdfast.NFW(d.M200, d.scale, d.G, d.rhoc)
		#print d.logslope(r3)
		#print d.enclosed_mass(r3)
		
		
		
		
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		#d.fast_ = mab.gd.gdfast.TwoSlopeDensity(d.rho0, d.alpha, d.beta, d.rs, d.gamma)
		
		#d.update()
		return d

	@staticmethod
	def fromEnclosed(rs, r, M):
		nfw = NFW(1e8, rs)
		Mcurr = nfw.enclosed_mass(r)
		nfw.rho0 *= (M/Mcurr)#**-1
		nfw.M200 = nfw.enclosed_mass(200)
		nfw.r200 = (nfw.M200/200*3/(4*pi)/nfw.rhoc)**(1./3)
		nfw.c = nfw.r200/nfw.rs
		return nfw
		
	def potentialr(self, r):
		x = r/self.rs
		return - 4 * pi * self.G * self.rho0 * self.rs**2 * log(1+x)/x
	
	def c_code_potential(self, r="r"):
		variables = {}
		variables["G"] = self.G
		variables["rho0"] = self.rho0
		variables["rs"] = self.rs
		variables["r"] = r
		return "-4 * M_PI * %(G)r * %(rho0)r * %(rs)r * %(rs)r * log(1+r/%(rs)r)/%(rs)r" % variables
	
	def densityr(self, r, M=None):
		if M == None:
			M = self.M
		x = r/self.rs
		return self.rho0 / (x * (1+x)**2) * M / self.M
	
	def logslope(self, r):
		return -2 * r / (r + self.rs) -1
	
	def r_at_slope(self, slope):
		return -self.rs*(slope+1)/(3+slope)
	
	def enclosed_mass(self, r):
		return self.enclosedmass(r)
	 
	def enclosedmass(self, r):
		return vectorize(self._enclosedmass)(r)
	def _enclosedmass(self, r):
		#try:
		#	len(r)
		#	indices = argsort(r)
		#	sorted_r = 
		#except:
		x = r/self.scale
		return 4 * pi * self.rho0 * self.scale**3 * (log(1+x)-x/(1+x))
	
	#def vcircr(self, r):
	#	return r * self.dphidr(r)
	
	def dphidr(self, r):
		#r = sqrt(x*x+y*y+z*z)
		#ra = self.rs
		#re = sqrt(x**2/self.a**2 + y**2/self.b**2 + z**2/self.c**2)
		#rtilde = (ra+r)*re / (ra+re)
		x = r/self.rs
		t1 = -4.0*pi*self.G*self.rho0*self.rs**3/r**2
		t2 = log(1.0 + x)
		t3 = x/(1.0+x)
		
		dphidr = -t1*(t2 - t3)

		return	dphidr
		
		#return self.G * self.M * r / (r**2 + self.b**2) ** (3./2)
		#return 4 * pi * self.G * self.rho0 * self.a / (2*(1+r/self.a)**2)
	
	 
	def c_code_dphidr(self, r="r"):
		variables = {}
		variables["pi"] = pi
		variables["G"] = self.G
		variables["rho0"] = self.rho0
		variables["rs"] = self.rs
		variables["rs_cubed"] = self.rs**3
		variables["r"] = r
		return "4 * %(pi)r * %(G)r * %(rho0)r * %(rs_cubed)r/pow(%(r)s, 2)* (log(1.0+%(r)s/%(rs)r) - (%(r)s/%(rs)r)/(1+%(r)s/%(rs)r)) " % variables
	

	def fast(self):
		return self.fast_
	
	def sample_r(self, rmax=10**3, N=1):
		# rejection sampling
		ok = False
		samples = numpy.zeros(N) * 1.0
		i = 0
		# max of rho(r)*r**2 is at this r
		r_max_density = self.rs
		M = self.densityr(r_max_density) * r_max_density**2
		#print "M =", M, r_max_density
		#r = r_max_density + 1e-2
		#print self.densityr(r) * r**2
		#r = r_max_density - 1e-2
		#print self.densityr(r) * r**2
		#M = 10
		while i < N:
			#x = random.random() * r
			x = numpy.random.random() * rmax
			u = numpy.random.random()
			fx = self.densityr(x)*x**2
			if u < (fx/M):
				samples[i] = x
				i += 1
		if N == 1:
			return samples[0]
		else:
			return samples
	

def mass_diff(logr200, rhoc, rho0, rs):
	r200 = exp(logr200)
	#r200 = abs(r200)
	x = r200/rs
	m1 = 200*rhoc*4*pi/3*r200**3
	m2 = 4 * pi * rho0 * rs**3 * (log(1+x)-x/(1+x))
	return (m1 - m2)/(m1+m2)

class NFW1kpc(NFW):
	def __init__(self, M1kpc, rs, G=G, cosmology=mab.cosmology.wmap3y):
		self.cosmology = cosmology
		self.rhoc = rhoc = (cosmology.rho_crit * cosmology.h**2).asNumber(mab.astrounits.MSOL/mab.astrounits.KPC**3)
		x = 1./rs
		self.rho0 = M1kpc/(4*pi*rs**3*(log(1+x)-x/(1+x)))
		self.scale = self.rs = rs
		self.G = G
		#print self.enclosedmass(1.), M1kpc
		#assert abs(self.enclosedmass(1.) - M1kpc) < 0.01
		
		#r200 = exp(fmin(mass_diff, 10, args=(rhoc, self.rho0, rs), disp=1)[0])
		#print rhoc, self.rho0
		logr200 = (fsolve(mass_diff, 2, args=(rhoc, self.rho0, rs)))#[0]
		try:
			len(logr200)
			logr200 = logr200[0]
		except:
			 pass
		   
		#print logr200, exp(logr200)
		r200 = exp(logr200)
		self.r200 = r200
		self.M200 = 200*rhoc*4*pi/3*r200**3
		x = r200/rs
		#print r200
		#print 200*rhoc*4*pi/3*r200**3, 4 * pi * self.rho0 * rs**3 * (log(1+x)-x/(1+x))
		#print self.M200, self.enclosedmass(self.r200), self.M200/self.enclosedmass(self.r200)
		assert abs((self.enclosedmass(self.r200) - self.M200)/self.M200) < 0.01
		#print self.r200, self.M200
		#self.rhoc = rhoc
		#self.r200 = (self.M200/200*3/(4*pi)/rhoc)**(1./3)
		#M200 = rho0
		#NFW1kpc.__init__(M200, rs=rs, G=G, cosmology=cosmology)
		self.M1kpc = M1kpc
		self.M = M1kpc
		self.fast_ = gdfast.NFW(self.M200, self.scale, self.G, self.rhoc)
		self.c = self.r200/self.rs
	
	

	
	
class Density(object):
	def _enclosed_mass(self, r):
		M, err = quad(lambda r: 4*pi*self.densityr(r)*r**2, 0, r)
		return M
	
	def enclosed_mass(self, r):
		isvector = False
		try:
			len(r)
			isvector = True
		except:
			pass
		if isvector:
			n = len(r)
			masses = zeros(n)
			r1 = 0
			cummass = 0
			for i in range(n):
				r2 = r[i]
				dmass, err = quad(lambda r: 4*pi*self.densityr(r)*r**2, r1, r2)
				cummass += dmass
				masses[i] = cummass
				r1 = r2
			return masses
		else:
			return numpy.vectorize(self._enclosed_mass)(r)
		
	def normalize(self, r, Mtarget):
		def ftest(log10rho, density, r, Mtarget):
			self.rho0 = 10**log10rho
			Mr = self.enclosed_mass(r)
			#diff = (Mr - Mtarget) / (Mr + Mtarget)
			diff = (Mtarget-Mr) / (Mr + Mtarget)
			#print "(", `log10rho`, r, Mtarget, ")"
			#print Mtarget, Mr, density.rho0, r, diff
			return diff
			
		#print "start", `self.rho0`, self.alpha, self.beta
		self.rho0 = 10**scipy.optimize.brentq(ftest, -20, 50, args=(self, r, Mtarget))
		#print scipy.optimize.__file__
		#print "res", res
		

class NFWCut(Density):
	def __init__(self, rho0, rs, rte):
		self.rho0 = rho0
		self.rs = rs
		self.rte = rte
		self.fast_ = mab.gd.gdfast.NFWCut(self.rho0, self.rs, self.rte)
		
	@staticmethod
	def from_r3(M3, r3, slope, rte):
		d = NFWCut(1., 1., rte)
		def f(logrs):
			d.rs = 10**logrs
			return (d.logslope(r3)[0]-slope)**2
			
		rs = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		d.rs = rs
		M3now = d.enclosed_mass(r3)
		d.rho0 = d.rho0 / M3now * M3
		#d.normalize(r3, M3)
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		#d.fast_ = mab.gd.gdfast.TwoSlopeDensity(d.rho0, d.alpha, d.beta, d.rs, d.gamma)
		
		d.update()
		return d
		
		
	def update(self):
		self.fast_.rho0 = self.rho0
		self.fast_.rs = self.rs
		self.fast_.rte = self.rte
		
	def densityr(self, r, M=None):
		x = r/self.rs;
		rte = self.rte
		return self.rho0 / (x * (1+x)**2) / (1 + (r/rte)**3)
	
	def fast(self):
		return self.fast_
		
	def logslope(self, r):
		rs = self.rs
		rte = self.rte
		slope = -(2 * r**3 * (3*r+2*rs) + (3*r + rs)*rte**3) / ((r+rs)*(r**3+rte**3))
		return slope
		
	def setslope_at(self, r, slope):
		rte = self.rte
		self.rs = (-6*r**4-3*r*rte**3-r**4*slope-r*rte**3*slope) / ( 4*r**3 + rte**3 + r**3*slope + rte**3*slope)
		#self.fast_.rs = self.rs
		#return 
		

class BrokenPowerLawDensitySoft3(Density):
	def __init__(self, s1, s2, s3, gamma1, gamma2, rs1, rs2, rho0=1., M1kpc=None):
		self.rho0 = rho0
		self.s1 = s1
		self.s2 = s2
		self.s3 = s3
		self.gamma1 = gamma1
		self.gamma2 = gamma2
		self.rs1 = rs1
		self.rs2 = rs2
		self.M1kpc = M1kpc
		if M1kpc:
			self.normalize(1., M1kpc)
			Mr = self.enclosed_mass(1)
			#print "enclosed mass at %f (kpc) is %f (should be %f)" % (1., Mr, M1kpc)
			if abs((Mr-M1kpc)/M1kpc) > 0.01:
				raise "error, enclosed mass differs by more than 1%"
		self.fast_ = mab.gd.gdfast.BrokenPowerLawDensitySoft3(self.rho0, self.s1, self.s2, self.s3, self.gamma1, self.gamma2, self.rs1, self.rs2)
		
	def update(self):
		self.fast_.rho0 = self.rho0
		self.fast_.s1 = self.s1
		self.fast_.s2 = self.s2
		self.fast_.s3 = self.s3
		self.fast_.gamma1 = self.gamma1
		self.fast_.gamma2 = self.gamma2
		self.fast_.rs1 = self.rs1
		self.fast_.rs2 = self.rs2
		if self.M1kpc:
			self.normalize(1., self.M1kpc)
			Mr = self.enclosed_mass(1)
			print "enclosed mass at %f (kpc) is %f (should be %f)" % (1., Mr, self.M1kpc)
			if abs((Mr-self.M1kpc)/self.M1kpc) > 0.01:
				raise "error, enclosed mass differs by more than 1%"
			self.fast_.rho0 = self.rho0
		
	def densityr(self, r):
		x1 = r/self.rs1;
		x2 = r/self.rs2;
		return self.rho0 * x1**self.s1 * (1+x1**self.gamma1)**((self.s2-self.s1)/self.gamma1) * (1+x2**self.gamma2)**((self.s3-self.s2)/self.gamma2);
	
	def fast(self):
		return self.fast_
		
	def logslope(self, r):
		#return -(self.rs*self.alpha + r*self.beta)/(r + self.rs)
		x1 = r/self.rs1
		x2 = r/self.rs2
		return self.s3 + (self.s1-self.s2)/(1+x1**self.gamma1) + (self.s2 - self.s3) / (1+x2**self.gamma2)

class TwoSlopeDensityCut(Density):
	def __init__(self, alpha, beta, rs, gamma=1., rho0=1., rte=1):
		self.rho0 = rho0
		self.alpha = alpha
		self.beta = beta
		self.gamma = float(gamma)
		self.rs = rs
		self.rte = rte
		self.fast_ = mab.gd.gdfast.TwoSlopeDensityCut(self.rho0, self.alpha, self.beta, self.rs, self.gamma, self.rte)
		
	def update(self):
		self.fast_.alpha = self.alpha
		self.fast_.beta = self.beta
		self.fast_.rs = self.rs
		self.fast_.rho0 = self.rho0
		self.fast_.gamma = self.gamma
		self.fast_.rte = self.rte
		

	def densityr(self, r, M=None):
		x = r/self.rs;
		scale = 1
		rte = self.rte
		#if M:
		#	scale = self.M
		return self.rho0 / (x**-self.alpha * (1+x**self.gamma)**((-self.beta+self.alpha)/self.gamma)) * scale  / (1 + (r/rte)**3);
	
	def fast(self):
		return self.fast_
		
	def logslope(self, r):
		#return -(self.rs*-self.alpha + r*-self.beta)/(r + self.rs)
		return -3*r**3 / (r**3+self.rte**3)  -(-self.beta + (-self.alpha+self.beta)/(1+(r/self.rs)**self.gamma))
		
		#return -(-self.beta + (self.beta-self.alpha)/(1+(r/self.rs)**self.gamma))
	
	def setslope_at(self, r, slope):
		#rte = self.rte
		#self.rs = (-6*r**4-3*r*rte**3-r**4*slope-r*rte**3*slope) / ( 4*r**3 + rte**3 + r**3*slope + rte**3*slope)
		#self.rs = -r * (self.rte**3*(3+self.gamma) + r**3*(6+self.gamma)) / (self.rte**3*(1+self.gamma) + r**3 * (4 + self.gamma
		self.rs = r * (- (rte**3*(self.alpha+slope) + r**3*(3 + self.alpha+slope)) / (rte**3*(self.beta+slope) + r**3*(3+self.beta+slope)) ) ** (-1/self.gamma)


import scipy.special
class TwoSlopeDensityWrapper(Density):
	def __init__(self, alpha, beta, rs, gamma=1., rho0=1., M1kpc=None, G=G):
		self.rho0 = rho0
		self.alpha = alpha
		self.beta = beta
		self.gamma = float(gamma)
		#print alpha, beta
		self.rs = rs
		self.M1kpc = M1kpc
		self.G = G
		if M1kpc:
			self.normalize(1., M1kpc)
			Mr = self.enclosed_mass(1)
			#print "enclosed mass at %f (kpc) is %f (should be %f)" % (1., Mr, M1kpc)
			if abs((Mr-M1kpc)/M1kpc) > 0.01:
				raise "error, enclosed mass differs by more than 1%"
		self.fast_ = mab.gd.gdfast.TwoSlopeDensity(self.rho0, self.alpha, self.beta, self.rs, self.gamma)
		
	def enclosed_mass_(self, r):
		a = self.alpha
		b = self.beta
		g = self.gamma
		g = scipy.special.hyp2f1((3+a)/g,(-b+a)/g,(g+a+3)/g,-r**g)
		return -4 * pi * self.rho0 * r**(3+self.alpha) * g / (-3-a)
		
	def enclosed_binding_energy(self, r):
		return vectorize(self._enclosed_binding_energy)(r)
	def _enclosed_binding_energy(self, r):
		def dW(x):
			return -4 * pi * self.G * self.densityr(x) * self.enclosed_mass(x) * x
			#return 4 * pi * self.densityr(x) * self.potentialr(x) * x**2
		W, err = quad(dW, 0, r)
		return W
	def denclosed_binding_energy(self, r):
		return vectorize(self._denclosed_binding_energy)(r)
	def _denclosed_binding_energy(self, r):
		return -4 * pi * self.G * self.densityr(r) * self.enclosed_mass(r) * r
		

	@staticmethod
	def from_r3(M3, r3, slope, alpha, beta, gamma):
		d = TwoSlopeDensityWrapper(alpha, beta, 1., gamma)
		def f(rs):
			d.rs = abs(rs)
			#print rs, d.logslope(r3), slope
			#import pdb
			#pdb.set_trace()
			return (d.logslope(r3)[0]-slope)**2
			
		rs = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		d.rs = rs
		d.normalize(r3, M3)
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		#d.fast_ = mab.gd.gdfast.TwoSlopeDensity(d.rho0, d.alpha, d.beta, d.rs, d.gamma)
		
		d.update()
		return d
		
		
		
	def update(self):
		self.fast_.alpha = self.alpha
		self.fast_.beta = self.beta
		self.fast_.rs = self.rs
		self.fast_.rho0 = self.rho0
		self.fast_.gamma = self.gamma
		if self.M1kpc:
			self.normalize(1., self.M1kpc)
			Mr = self.enclosed_mass(1)
			print "enclosed mass at %f (kpc) is %f (should be %f)" % (1., Mr, self.M1kpc)
			if abs((Mr-self.M1kpc)/self.M1kpc) > 0.01:
				raise "error, enclosed mass differs by more than 1%"
			self.fast_.rho0 = self.rho0
			
	def run(self, args, opts, scope):
		print self.concentration()
		

	def concentration(self):
		cosmology=mab.cosmology.wmap3y
		rhoc = (cosmology.rho_crit * cosmology.h**2).asNumber(mab.astrounits.MSOL/mab.astrounits.KPC**3)
		def ftest(r200, self=self):
			print r200
			M200 = self.enclosed_mass(r200)
			volume = 4*pi/3 * r200**3

			rhoavg = M200/volume
			diff = log10(rhoavg) - log10(rhoc*200)
			#diff = (Mr - Mtarget) / (Mr + Mtarget)
			#diff = (Mtarget-Mr) / (Mr + Mtarget)
			#print "(", `log10rho`, r, Mtarget, ")"
			#print Mtarget, Mr, density.rho0, r, diff
			print diff
			return diff
			
		#print "start", `self.rho0`, self.alpha, self.beta
		r200 = scipy.optimize.brentq(ftest, 1e-4, 1e8, args=(self))
		print r200
		M200 = self.enclosed_mass(r200)
		print "m200 %g" % M200
		return r200/self.rs
		#self.r200 = (self.M200/200*3/(4*pi)/rhoc)**(1./3)
			
		
	def densityr(self, r, M=None):
		x = r/self.rs;
		scale = 1
		#if M:
		#	scale = self.M
		return self.rho0 / (x**-self.alpha * (1+x**self.gamma)**((-self.beta+self.alpha)/self.gamma)) * scale;
	
	def fast(self):
		return self.fast_
		
	def logslope(self, r):
		#return -(self.rs*-self.alpha + r*-self.beta)/(r + self.rs)
		return -(-self.beta + (self.beta-self.alpha)/(1+(r/self.rs)**self.gamma))
	
	def r_at_slope(self, slope):
		return ((-self.alpha-slope)/(slope+self.beta))**(1/self.gamma)
	
class ProfileNumerical1dWrapper(Potential):
	def __init__(self, n, density, G=G):
		self.G = G
		self.density = density
		self.fastdensity = density.fast()
		self.transformation = mab.gd.gdfast.TransformationSphericalLog_in_3d(10.)
		self.fast_ = mab.gd.gdfast.ProfileNumerical1dL3(n, self.fastdensity, self.transformation, G, -6, 6)
		self.n = n;
		self.rs = density.rs
		
	def update(self):
		self.density.update()
		self.fast_ = mab.gd.gdfast.ProfileNumerical1dL3(self.n, self.fastdensity, self.transformation, self.G, -6, 6)
		
	def fast(self):
		return self.fast_
		
	def potentialr(self, r):
		#return self.density.potentialr(r)
		if isinstance(r, float):
			return self.fast_.potentialr(r)
		else:
			return array([self.fast_.potentialr(k) for k in r])
	
	def densityr(self, r, M=None):
		return self.density.densityr(r, M=M)
	
	def dphidr(self, r):
		if isinstance(r, float):
			return self.fast_.dphidr(r)
		else:
			return array([self.fast_.dphidr(k) for k in r])
		#return self.density.dphidr(r)
		#return self.fast_.dphidr(r)
	def enclosed_mass(self, r):
		return self.density.enclosed_mass(r)
	
	def logslope(self, r):
		return self.density.logslope(r)
	
	def r_at_slope(self, slope):
		return self.density.r_at_slope(r)
		
	def sample_r(self, N=1, rmax=100):
		Mmax = self.enclosed_mass(rmax)
		rs = []
		while len(rs) < N:
			def ftest(r, u):
				fraction = self.enclosed_mass(r)/Mmax
				return fraction-u
			u = numpy.random.random()
			r = scipy.optimize.brentq(ftest, 1e-6, rmax, args=(u))
			rs.append(r)
		if N == 1:
			return rs[0]
		else:
			return array(rs)
		
	
# rho0_bre          rs_bre       alpha_bre          rt_bre          rd_bre
class Breddels(Density):
	def __init__(self, rho0, rs, alpha, rtidal, rdecay, G=G):
		self.rho0  = rho0
		self.rs = rs
		self.alpha = alpha
		self.rtidal = rtidal
		self.rdecay = rdecay
		self.G = G
	
	def _densityr(self, r):
		x = r / self.rs
		x1 = self.rtidal / self.rs;
		gamma = self.rdecay / self.rs;
		eps = -(1 + 3 * x1) / (1 + x1) + x1 / gamma;

		if x < x1:
			f = exp (-(2 / self.alpha) * (pow (x, self.alpha) - 1))
		else:
			f = exp (-(2 / self.alpha) * (pow (x1, self.alpha) - 1))
			f *= pow(x / x1, eps) * exp (-(x - x1) / gamma);
		return f
	
	def densityr(self, r):
		return vectorize(self._densityr)(r)
			
	def _logslope(self, r):
		x = r / self.rs
		x1 = self.rtidal / self.rs;
		gamma = self.rdecay / self.rs;
		eps = -2 * x1 ** self.alpha + x1 / gamma

		
		if x < x1:
			return -2 * x**self.alpha
		else:
			return eps - x / gamma
		
		#i = where (x lt x1, ni, complement=j, ncomplement=nj)
		#y = dblarr (n_elements (x))

		#if ni gt 0 then y[i] =  -2 * x[i]^beta
		#if nj gt 0 then begin
		#	y[j] = eps - x[j] / alpha
		#endif
		
		

	def logslope(self, r):
		return vectorize(self._logslope)(r)
  
		dr = 1e-7
		rho1 = self.densityr(r)
		return  (self.densityr(r+dr) - rho1) / dr * r/rho1
		
		
	
class Moore(Density):
	def __init__(self, rho0, rs, rtidal, rdecay, G=G):
		self.rho0  = rho0
		self.rs = rs
		self.rtidal = rtidal
		self.rdecay = rdecay
		self.G = G
	
	def _densityr(self, r):
		x = r/self.rs
		x1 = self.rtidal / self.rs
		alpha = self.rdecay / self.rs
		epsilon = - (1+3*x1)/(1+x1) + x1/alpha
		if x < x1:
			return 1/(x*(1+x)**2)
		else:
			return 1/(x1*(1+x1)**2) * (x/x1)**epsilon * exp(-(x-x1)/alpha)
	
	def densityr(self, r):
		return vectorize(self._densityr)(r)
			
	def logslope(self, r):
		dr = 1e-6
		rho1 = self.densityr(r)
		return  (self.densityr(r+dr) - rho1) / dr * r/rho1
	
class Einasto(Potential):
	def __init__(self, rho_2, r_2, alpha, G=G):
		self.rho_2 = rho_2
		self.r_2 = r_2
		self.rs = self.r_2
		self.alpha = alpha
		self.scale = r_2
		self.G = G
		self.fast_ = gdfast.Einasto(rho_2, r_2, alpha, self.G)
		
	@staticmethod
	def from_r3(M3, r3, slope, alpha):
		d = Einasto(1., 1., alpha)
		def f(r_2):
			d.r_2 = abs(r_2)
			#print rs, d.logslope(r3), slope
			#import pdb
			#pdb.set_trace()
			return (d.logslope(r3)[0]-slope)**2
			
		r_2 = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		d.scale = r_2
		d.rs = r_2
		d.r_2 = r_2
		#d.normalize(r3, M3)
		Menc = d.enclosed_mass(r3)
		d.rho_2 *= M3/Menc
		d.fast_ = gdfast.Einasto(d.rho_2, d.r_2, d.alpha, d.G)
		
		
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		#d.fast_ = mab.gd.gdfast.TwoSlopeDensity(d.rho0, d.alpha, d.beta, d.rs, d.gamma)
		
		#d.update()
		return d

	def densityr(self, r, M=None):
		#if M == None:
		#	M = 1#@self.M
		x = r/self.r_2
		return self.rho_2 * exp(-2./self.alpha * (x**self.alpha-1))
	
	def potentialr(self, r):
		return 0
	
	
	def logslope(self, r):
		x = r/self.r_2
		return -2 * x**self.alpha
	
	#def r_at_slope(self, slope):
	#	return -self.rs*(slope+1)/(3+slope)
	
	def enclosed_mass(self, r):
		return self.enclosedmass(r)
	
	def enclosed_mass_num(self, r):
		return self.enclosedmass_num(r)
	 
	def enclosedmass(self, r):
		return vectorize(self._enclosedmass)(r)
	def enclosedmass_num(self, r):
		return vectorize(self._enclosedmass_num)(r)
	def _enclosedmass(self, r):
		x = r/self.scale
		a = self.alpha
		R = scipy.special.gammainc(3./a, (2*x**a)/a) * scipy.special.gamma(3./a)
		return 2**(2.-3./a) * exp(2./a) * pi * ((1/self.r_2)**a/a)**(-3./a) * self.rho_2 * R / a
		
		
		def g(x):
			#R = scipy.special.gamma(3./self.alpha) - (1-scipy.special.gammainc(3./self.alpha, (2*x**self.alpha)/self.alpha))
			R = scipy.special.gammainc(3./self.alpha, (2*x**self.alpha)/self.alpha)
			return 8.**(1./self.alpha) * exp(2./self.alpha) * self.alpha **(-1.+3./self.alpha) * R
		return 4 * pi * self.rho_2 * self.scale**3 * g(x)
	
	def _enclosedmass_num(self, r):
		def dM(rp):
			return 4 * pi * rp**2 * self.densityr(rp)
		M, err = quad(dM, 0, r)
		return M 
	
	def dphidr(self, r):
		return  self.G * self._enclosedmass(r) / r**2
	def fast(self):
		return self.fast_
	
	@staticmethod
	def fromM1kpc(M1kpc, rs, alpha):
		p = Einasto(1, rs, alpha)
		Menc = p.enclosed_mass(1.)
		p.rho_2 *= M1kpc/Menc
		p.fast_ = gdfast.Einasto(p.rho_2, rs, p.alpha, p.G)
		return p

class Burkert(Potential):
	def __init__(self, rho0, rs, G=G):
		self.rho0 = rho0
		self.rs = rs
		self.G = G
		#self.fast_ = gdfast.Einasto(rho_2, r_2, alpha, self.G)
		self.fast_ = mab.gd.gdfast.Burkert(self.rho0, self.rs, self.G)
		
	@staticmethod
	def from_r3(M3, r3, slope):
		d = Burkert(1., 1.)
		def f(rs):
			d.rs = abs(rs)
			#print rs, d.logslope(r3), slope
			#import pdb
			#pdb.set_trace()
			return (d.logslope(r3)[0]-slope)**2
			
		rs = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		#d.scale = r_2
		d.rs = rs
		Menc = d.enclosed_mass(r3)
		d.rho0 *= M3/Menc
		
		#Menc = d.enclosed_mass(r3)
		#d.rho_2 *= M3/Menc
		#d.fast_ = gdfast.Einasto(d.rho_2, d.r_2, d.alpha, d.G)
		
		
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		d.fast_ = mab.gd.gdfast.Burkert(d.rho0, d.rs, G)
		
		#d.update()
		return d

	def densityr(self, r, M=None):
		#if M == None:
		#	M = 1#@self.M
		x = r/self.rs
		return self.rho0 * ((1+x) * (1+x**2))**-1
	
	def potentialr(self, r):
		return 0
	
	
	def logslope(self, r):
		x = r/self.rs
		slope = -3 + 1/(1 + x) + 2/(1 + x**2)
		return slope
	
	#def r_at_slope(self, slope):
	#	return -self.rs*(slope+1)/(3+slope)
	
	def enclosed_mass(self, r):
		return self.enclosedmass(r)
	
	def enclosed_mass_num(self, r):
		return self.enclosedmass_num(r)
	 
	def enclosedmass(self, r):
		return self.enclosedmass_num(r)
	def enclosedmass_num(self, r):
		return vectorize(self._enclosedmass_num)(r)
	def __enclosedmass(self, r):
		pass
	
	def _enclosedmass_num(self, r):
		def dM(rp):
			return 4 * pi * rp**2 * self.densityr(rp)
		M, err = quad(dM, 0, r)
		return M 
	
	def dphidr(self, r):
		return  self.G * self._enclosedmass(r) / r**2
	def fast(self):
		return self.fast_
	
	def sample_r(self, rmax=10**3, N=1):
		# rejection sampling
		ok = False
		samples = numpy.zeros(N) * 1.0
		i = 0
		# max of rho(r)*r**2 is at this r
		r_max_density = self.rs * 1.52
		M = self.densityr(r_max_density) * r_max_density**2 * 1.1
		#print "rmax",r_max_density , self.densityr(r_max_density), self.densityr(r_max_density-0.1), self.densityr(r_max_density+0.1)
		#print "M =", M, r_max_density
		#r = r_max_density + 1e-2
		#print self.densityr(r) * r**2
		#r = r_max_density - 1e-2
		#print self.densityr(r) * r**2
		#M = 10
		while i < N:
			#x = random.random() * r
			x = numpy.random.random() * rmax
			u = numpy.random.random()
			fx = self.densityr(x)*x**2
			if u < (fx/M):
				samples[i] = x
				i += 1
		if N == 1:
			return samples[0]
		else:
			return samples
			
class Burker(Potential):
	def __init__(self, rho_2, r_2, alpha, G=G):
		self.rho_2 = rho_2
		self.r_2 = r_2
		self.rs = self.r_2
		self.alpha = alpha
		self.scale = r_2
		self.G = G
		self.fast_ = gdfast.Einasto(rho_2, r_2, alpha, self.G)
		
	@staticmethod
	def from_r3(M3, r3, slope, alpha):
		d = Einasto(1., 1., alpha)
		def f(r_2):
			d.r_2 = abs(r_2)
			#print rs, d.logslope(r3), slope
			#import pdb
			#pdb.set_trace()
			return (d.logslope(r3)[0]-slope)**2
			
		r_2 = abs(scipy.optimize.fmin(f, 1., disp=False))[0]
		d.scale = r_2
		d.rs = r_2
		d.r_2 = r_2
		#d.normalize(r3, M3)
		Menc = d.enclosed_mass(r3)
		d.rho_2 *= M3/Menc
		d.fast_ = gdfast.Einasto(d.rho_2, d.r_2, d.alpha, d.G)
		
		
		#print >>sys.stderr, ">>>>>>>", r3, M3, d.rs
		#d.fast_ = mab.gd.gdfast.TwoSlopeDensity(d.rho0, d.alpha, d.beta, d.rs, d.gamma)
		
		#d.update()
		return d

	def densityr(self, r, M=None):
		#if M == None:
		#	M = 1#@self.M
		x = r/self.r_2
		return self.rho_2 * exp(-2./self.alpha * (x**self.alpha-1))
	
	def potentialr(self, r):
		return 0
	
	
	def logslope(self, r):
		x = r/self.r_2
		return -2 * x**self.alpha
	
	#def r_at_slope(self, slope):
	#	return -self.rs*(slope+1)/(3+slope)
	
	def enclosed_mass(self, r):
		return self.enclosedmass(r)
	
	def enclosed_mass_num(self, r):
		return self.enclosedmass_num(r)
	 
	def enclosedmass(self, r):
		return vectorize(self._enclosedmass)(r)
	def enclosedmass_num(self, r):
		return vectorize(self._enclosedmass_num)(r)
	def _enclosedmass(self, r):
		x = r/self.scale
		a = self.alpha
		R = scipy.special.gammainc(3./a, (2*x**a)/a) * scipy.special.gamma(3./a)
		return 2**(2.-3./a) * exp(2./a) * pi * ((1/self.r_2)**a/a)**(-3./a) * self.rho_2 * R / a
		
		
		def g(x):
			#R = scipy.special.gamma(3./self.alpha) - (1-scipy.special.gammainc(3./self.alpha, (2*x**self.alpha)/self.alpha))
			R = scipy.special.gammainc(3./self.alpha, (2*x**self.alpha)/self.alpha)
			return 8.**(1./self.alpha) * exp(2./self.alpha) * self.alpha **(-1.+3./self.alpha) * R
		return 4 * pi * self.rho_2 * self.scale**3 * g(x)
	
	def _enclosedmass_num(self, r):
		def dM(rp):
			return 4 * pi * rp**2 * self.densityr(rp)
		M, err = quad(dM, 0, r)
		return M 
	
	def dphidr(self, r):
		return  self.G * self._enclosedmass(r) / r**2
	def fast(self):
		return self.fast_
	
	@staticmethod
	def fromM1kpc(M1kpc, rs, alpha):
		p = Einasto(1, rs, alpha)
		Menc = p.enclosed_mass(1.)
		p.rho_2 *= M1kpc/Menc
		p.fast_ = gdfast.Einasto(p.rho_2, rs, p.alpha, p.G)
		return p


class PotentialLogarithmic2d(object):
	def __init__(self, v0, Rc, q, G=G):
		self.v0 = v0
		self.Rc = Rc
		self.q = q
		self.G = G
		self.fast_ = gdfast.Logarithmic2d(self.v0, self.q, self.Rc, G)
		
	def fast(self):
		return self.fast_
		
	def potentialxy(self, x, y):
		return 0.5 * self.v0**2 * log(self.Rc**2+x**2+y**2/self.q**2)
	
	def dphidx(self, x, y):
		return 0.5 * self.v0**2 /(self.Rc**2+x**2+y**2/q**2) * 2 * x	
	
	def dphidy(self, x, y):
		return 0.5 * self.v0**2 /(self.Rc**2+x**2+y**2/q**2) * 2 * y/q**2
	
	
	
