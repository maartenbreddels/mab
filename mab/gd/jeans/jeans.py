# -*- coding: utf-8 -*-
from scipy.integrate import quad
from numpy import *
import mab.gd.gdfast
import scipy.integrate


class Jeans(object):
	def __init__(self, stellar_potential, potentials):
		self.potentials = potentials
		self.stellar_potential = stellar_potential 
		
	
	def sigmar(self, radii, smart=True):
		potential = self.stellar_potential 
		sumint = 0
		##I, err = quad(lambda r: p1.densityr(r) * sum([p.dphidr(r) for p in potentials]), r, inf)
		res = zeros(len(radii)) * 1.0
		if smart:
			indices = argsort(radii)[::-1]
			radii = radii[indices]
			r2 = inf
			cumsumI = 0
			for i in range(len(radii)):
				r1 = radii[i] 
				#I, err = quad(lambda r: potential.densityr(r) * sum([p.dphidr(r) for p in self.potentials]), rcurr, inf)
				for p in self.potentials:
					dI, err = quad(lambda r: potential.densityr(r) * p.dphidr(r), r1, r2)
					cumsumI += dI
				#I, err = quad(lambda r: potential.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 1* p3.dphidr(r)), r, inf)
				#I, err = quad(lambda r: p1.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 0* p3.dphidr(r)), r, inf)
				res[i] = sqrt(1/potential.densityr(r1) * cumsumI)
				r2 = r1
			res2 = zeros(len(radii)) * 1.0
			res2[indices] = res
			res = res2
		else:
			for i in range(len(radii)):
				rcurr = radii[i] 
				#I, err = quad(lambda r: potential.densityr(r) * sum([p.dphidr(r) for p in self.potentials]), rcurr, inf)
				sumI = 0
				for p in self.potentials:
					I, err = quad(lambda r: potential.densityr(r) * p.dphidr(r), rcurr, inf)
					sumI += I
				#I, err = quad(lambda r: potential.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 1* p3.dphidr(r)), r, inf)
				#I, err = quad(lambda r: p1.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 0* p3.dphidr(r)), r, inf)
				res[i] = sqrt(1/potential.densityr(rcurr) * sumI)
		return array(res)

class JeansBase(object):
	def __init__(self, light_model, profile_model):
		self.light_model = light_model
		self.profile_model = profile_model
		
	def m2_sorted(self, rs):
		m2 = rs*0
		r1s = rs
		r2s = numpy.concatenate([rs[1:], [inf]])
		#m2 = array([])
		N = len(rs)
		for i in range(N):
			m2[i] = self._m2(r1s[i], r2s[i])
		return cumsum(m2[::-1])[::-1]
	def m2(self, r):
		return vectorize(self._m2)(r)
	
	def sigmar(self, r):
		return sqrt(self.m2(r))
	
	def sigma_los(self, R):
		return vectorize(self._sigma_los)(R)
	
	def _sigma_los(self, R):
		surface_density = self.light_model.densityR(R);
		def f(r): 
			sigmar = self.sigmar(r)
			sigma = sigmar * sqrt(1-self.beta(r)*R*R/(r*r))
			return 2 / (surface_density) * sigma*sigma * r * self.light_model.densityr(r) / sqrt(r*r-R*R)

		I, err = quad(f, R, inf, epsabs=1e-3, epsrel=1e-4)
		sigmas_los = sqrt(I)
		return sigmas_los
		
	def _kinetic_energy(self, s):
		def f(x):
			return 2*pi* x**2*self.light_model.densityr(x)*(3.-2*self.beta0)*self.sigmar(x)**2
		I, err = scipy.integrate.quad(f, 0, s)
		return I
		
	def kinetic_energy(self, s):
		return vectorize(self._kinetic_energy)(s)
		
			
	

class JeansAnisotropicConstant(JeansBase):
	def __init__(self, light_model, profile_model, beta):
		JeansBase.__init__(self, light_model, profile_model)
		#@self.light_model = light_model
		#self.profile_model = profile_model
		#self.allpotentials = allpotentials
		self.beta0 = beta
		self.has_cache = False
		#self.fast_ = mab.gd.gdfast.JeansAnisotropicConstant(self.light_model.fast(), self.profile_model.fast(), self.beta)
		
	def beta(self, r):
		return self.beta0 + r * 0
	def fast(self):
		return mab.gd.gdfast.JeansAnisotropicConstant(self.light_model.fast(), self.profile_model.fast(), self.beta, 1e9)
	
	def _m2(self, r):
		#import pdb;
		#pdb.set_trace()
		if self.light_model.densityr(r) < 1e-100: #== 0:
			return 0
		N = r ** (-2*self.beta0) / self.light_model.densityr(r)
		def f(rp):
			return N * rp ** (2 * self.beta0) * self.light_model.densityr(rp) * self.profile_model.dphidr(rp)
		result = quad(f, r, inf, epsabs=1e-3, epsrel=1e-6, full_output=1)
		I = result[0]
		return I
	
	def m4(self, r):
		def f(r):
			return r ** (2 * self.beta0) * self.light_model.densityr(r) * self.m2(r) * self.profile_model.dphidr(r)
		result = quad(f, r, inf, epsabs=1e-3, epsrel=1e-2, full_output=1)
		I = result[0]
		return 3 * r ** (-2*self.beta0) / self.light_model.densityr(r) * I
		
	def m4_los(self, Radii):
		surface_densities = self.light_model.densityR(Radii);
		def f(r, R, surface_density):
			m4 = self.m4(r)
			g = 1 - 2 * self.beta0 * R**2/r**2 + self.beta0*(1+self.beta0)/2*R**4/r**4
			#sigmar = self.sigmar(array([r]))[0]
			#sigma = sigmar * sqrt(1-self.beta*R*R/(r*r))
			#print sigma
			return 2 / (surface_density) * m4* r * self.light_model.densityr(r) / sqrt(r*r-R*R) * g
		#IntegratorGSL<> iGSL(integrand);
		#// integrator from R to inf
		#//cout << "integrating..." << endl;
		#double sigmasq = iGSL.integrate_to_inf(R);
		m4_los = Radii * 0
		for i in range(len(m4_los)):
			I, err = quad(f, Radii[i], inf, args=(Radii[i],surface_densities[i]), epsabs=1e-3, epsrel=1e-2)
			m4_los[i] = I
		return m4_los
		
		
	def sigma_los(self, Radii):
		surface_densities = self.light_model.densityR(Radii);
		def f(r, R, surface_density): 
			sigmar = self.sigmar(array([r]))[0]
			sigma = sigmar * sqrt(1-self.beta0*R*R/(r*r))
			#print sigma
			return 2 / (surface_density) * sigma*sigma * r * self.light_model.densityr(r) / sqrt(r*r-R*R)
		#IntegratorGSL<> iGSL(integrand);
		#// integrator from R to inf
		#//cout << "integrating..." << endl;
		#double sigmasq = iGSL.integrate_to_inf(R);
		sigmas_los = Radii * 0
		for i in range(len(sigmas_los)):
			I, err = quad(f, Radii[i], inf, args=(Radii[i],surface_densities[i]), epsabs=1e-3, epsrel=1e-2)
			sigmas_los[i] = sqrt(I)
		return sigmas_los
			
	
	def make_sigmar_cache(self, rmax, dr):
		self.cache_r = arange(0, rmax, dr)
		self.cache_integral = self.cache_r * 0.0
		r2 = inf
		N = len(self.cache_r)
		I = 0
		for i in range(N):
			r1 = self.cache_r[N-i-1]
			dI = self.integral(r1, r2)
			I += dI
			r2 = r1
			self.cache_integral[N-i-1] = I
		self.has_cache = True
		
	def sigmar_fast(self, radius):
		if not self.has_cache:
			return self.sigmar(array([radius]))[0]
		else:
			indices = ravel(where(self.cache_r > radius))
			#print len(indices), indices, len(indices) == 0
			if len(indices) == 0:
				return self.sigmar(array([radius]))[0]
			else:
				inext = indices[0]
				r1 = radius
				r2 = self.cache_r[inext]
				#print self.cache_r, inext
				#print "integrate..."
				I = self.integral(r1, r2)
				I += self.cache_integral[inext]
				#print "fast", r1, r2, I, self.cache_integral[inext]
				#self.dsigmar(array([radius]))[0]
				#res[i] = sqrt(1/(self.stellar_potential.densityr(r1) * r1**twobeta) * cumsumI)
				twobeta = self.beta0 * 2
				return sqrt(1/(self.stellar_potential.densityr(radius) * radius**twobeta) * I)
			
		
	def integral(self, r1, r2):
		I = 0
		twobeta = self.beta0 * 2
		#for p in self.allpotentials:
		#	if p.M > 0:
		#		dI, err = quad(lambda r: r**twobeta * self.stellar_potential.densityr(r) * p.dphidr(r), r1, r2)
		#		#print "--", r1, r2, dI
		#		I += dI
		I, err = quad(lambda r: r**twobeta * self.light_model.densityr(r) * self.profile_model.dphidr(r), r1, r2)
		return I
		
		
	def _sigmar(self, radii, smart=True):
		sumint = 0
		##I, err = quad(lambda r: p1.densityr(r) * sum([p.dphidr(r) for p in potentials]), r, inf)
		res = zeros(len(radii)) * 1.0
		twobeta = self.beta0 * 2
		if smart:
			indices = argsort(radii)[::-1]
			radii = radii[indices]
			r2 = inf
			cumsumI = 0
			#print radii
			for i in range(len(radii)):
				r1 = radii[i] 
				cumsumI += self.integral(r1, r2)
				#print cumsumI, r1, r2 
				res[i] = sqrt(1/(self.light_model.densityr(r1) * r1**twobeta) * cumsumI)
				r2 = r1
			res2 = zeros(len(radii)) * 1.0
			res2[indices] = res
			res = res2
		else:
			for i in range(len(radii)):
				rcurr = radii[i] 
				#I, err = quad(lambda r: potential.densityr(r) * sum([p.dphidr(r) for p in self.potentials]), rcurr, inf)
				sumI = 0
				#for p2 in self.allpotentials:
				for p2 in [self.stellar_potential]:
					for p in self.allpotentials:
						if p.M > 0:
							I, err = quad(lambda r: r**twobeta * p2.densityr(r) * p.dphidr(r), rcurr, inf)
							#print I, rcurr, inf
							sumI += I
				#I, err = quad(lambda r: potential.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 1* p3.dphidr(r)), r, inf)
				#I, err = quad(lambda r: p1.densityr(r) * (p1.dphidr(r) + p2.dphidr(r) + 0* p3.dphidr(r)), r, inf)
				#c = sum([p2.densityr(rcurr) for p2 in self.allpotentials])
				c = sum([p2.densityr(rcurr) for p2 in [self.stellar_potential]])
				res[i] = sqrt(1/(c * rcurr**twobeta) * sumI)
				#crashme
			

		return array(res)

	
		
class JeansOM(JeansBase):
	def __init__(self, light_model, profile_model, beta0, betainf, r_beta):
		JeansBase.__init__(self, light_model, profile_model)
		#self.light_model = light_model
		#self.profile_model = profile_model
		#self.allpotentials = allpotentials
		self.beta0 = beta0
		self.betainf = betainf
		self.r_beta = r_beta
		self.has_cache = False
		#self.fast_ = mab.gd.gdfast.JeansAnisotropicConstant(self.light_model.fast(), self.profile_model.fast(), self.beta)
		
	def ___fast(self):
		return mab.gd.gdfast.JeansAnisotropicConstant(self.light_model.fast(), self.profile_model.fast(), self.beta, 1e9)
	
	def beta(self, r):
		u = r**2/(r**2+self.r_beta**2)
		return self.beta0 + (self.betainf-self.beta0) * u
	
	def _m2(self, r):
		N = 1/ self.light_model.densityr(r)
		def f(rp):
			l1 = ((r**2+self.r_beta**2)*rp**2)/(self.r_beta**2+rp**2)
			l2 = (self.r_beta**2+rp**2)/(r**2+self.r_beta**2)
			u = (-2*self.beta0*log(r) + self.beta0*log(l1) + self.betainf * log(l2))
			return N * exp(u) * self.light_model.densityr(rp) * self.profile_model.dphidr(rp)
		result = quad(f, r, inf, epsabs=1e-3, epsrel=1e-6, full_output=1)
		I = result[0]
		return I
	

			
	
				