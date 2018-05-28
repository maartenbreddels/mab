# -*- coding: utf-8 -*-
import scipy.optimize
from numpy import *
from mab.constants import *
from mab.astrounits  import *

class ParticleOrbit(object):
	def __init__(self, host, apocenter, orbit_integrator, pericenter=None, eccentricity=None, Norbits=None, T=1., Tmin=None, calculate=False):
		self.host = host
		self.ra = apocenter
		print "apocenter", self.ra
		if pericenter:
			self.rp = pericenter
		else:
			e = eccentricity
			self.rp = self.ra * (1-e)/(1+e)
			print "pericenter", self.rp
		self.rmid = (self.ra + self.rp)/2
		self.orbit_integrator = orbit_integrator
		self.T = T
		self.Norbits = Norbits
		self.Tmin = Tmin
		if calculate:
			self.find_orbit()
		
	def find_orbit_EL(self, E, L):
		profile_model = self.host.profile_model
		self.E, self.L = E, L
		ra, rp = profile_model.get_apo_peri(self.E, self.L)
		#print "found ra/rp", ra, rp
		self.Torbit = profile_model.Tr(self.E, self.L)
		self.vt0 = self.L/self.ra# / 10
		#print "vt0", self.vt0
		self.x0 = ra
		self.q0 = (self.x0, 0., 0.)
		self.p0 = (0.0, self.vt0, 0.) 
		#self.q0 = (79, 0., 0.)
		#self.p0 = (160, 240., 0.) 
		if self.Norbits:
			self.Tintegrate = self.Torbit * self.Norbits
		else:
			self.Tintegrate = self.T
		if self.Tmin:
			self.Tintegrate = max(self.Tmin, self.Tintegrate)
		N = self.orbit_integrator.orbital_points
		self.ts = self.Tintegrate * arange(N) / (1+N)
		
	def run(self, args, opts, scope):
		self.find_orbit()
		print "x", self.x0
		print "vt", self.vt0
	
	def find_orbit(self):
		profile_model = self.host.profile_model
		#E = profile_model.potentialr(self.ra)
		def f(x):
			E, L = x
			#Epot = profile_model.potentialr(self.ra)
			try:
				ra, rp = profile_model.get_apo_peri(E, L)
			except:
				#print "error!"
				return 1e10
			#print ra, rp
			error = sqrt((self.ra - ra)**2 + (self.rp - rp)**2)
			#print error
			return float(error)
		E = profile_model.potentialr(self.rmid)
		L = profile_model.Lmax_at_E(E) * 0.5
		result = scipy.optimize.fmin(f, [E, L], disp=False, full_output=True)
		self.E, self.L = result[0]
		ra, rp = profile_model.get_apo_peri(self.E, self.L)
		#print "found ra/rp", ra, rp
		self.Torbit = profile_model.Tr(self.E, self.L)
		print "orbital timescale", self.Torbit
		self.vt0 = self.L/self.ra# / 10
		#print "vt0", self.vt0
		self.x0 = ra
		self.q0 = (self.x0, 0., 0.)
		self.p0 = (0.0, self.vt0, 0.) 
		#self.q0 = (79, 0., 0.)
		#self.p0 = (160, 240., 0.) 
		if self.Norbits:
			self.Tintegrate = self.Torbit * self.Norbits
		else:
			self.Tintegrate = self.T
		if self.Tmin:
			self.Tintegrate = max(self.Tmin, self.Tintegrate)
		N = self.orbit_integrator.orbital_points
		self.ts = self.Tintegrate * arange(N) / (1+N)
		
	def integrate(self, remember=False):
		
		s_to_gyr = (S/GYR).asNumber()
		N = self.orbit_integrator.orbital_points
			
		dt = self.Tintegrate/N / s_to_gyr
		#print "dt", dt
		dt = array([dt])
		#print "1", self.q0
		q0 = array([self.q0[:2]])
		p0 = array([self.p0[:2]])
		#print "2", q0
		#print "dt", dt
		#import pyublas
		#pyublas.why_not(q0)
		#print dir(pyublas)
		q, p = self.orbit_integrator.integrate(dt, q0, p0)
		if 0:
			print self.host.profile_model.potentialr(1.0)
			print self.host.profile_model.potentialr(2.0)
			print self.orbit_integrator.profile_model.potentialr(1.0)
			print self.orbit_integrator.profile_model.potentialr(2.0)
			print self.orbit_integrator.profile_model_fast.potentialr(1.0)
			print self.orbit_integrator.profile_model_fast.potentialr(2.0)
			print p.shape
		if remember:
			self.q0 = q[:,0,-1]
			self.q0 = (self.q0[0], self.q0[1], 0.)
			self.p0 = p[:,0,-1]
			self.p0 = (self.p0[0], self.p0[1], 0.)
			#print "3", self.q0
		return q[:,0,:], p[:,0,:]
		
class Orbit2d(object):
	def __init__(self, host, orbit_integrator, x0, y0, vx0, vy0, T):
		self.host = host
		self.orbit_integrator = orbit_integrator
		self.T = T
		self.q0 = (x0, y0, 0)
		self.p0 = (vx0, vy0, 0)
		
	def find_orbit(self):
		pass
	def integrate(self):
		s_to_gyr = (S/GYR).asNumber()
		N = self.orbit_integrator.orbital_points
		dt = self.T/N/s_to_gyr
		#print "dt", dt
		dt = array([dt])
		#print "1", self.q0
		q0 = array([self.q0[:2]])
		p0 = array([self.p0[:2]])
		#print "2", q0
		#print "dt", dt
		#import pyublas
		#pyublas.why_not(q0)
		#print dir(pyublas)
		q, p = self.orbit_integrator.integrate(dt, q0, p0)
		if 0:
			print self.host.profile_model.potentialr(1.0)
			print self.host.profile_model.potentialr(2.0)
			print self.orbit_integrator.profile_model.potentialr(1.0)
			print self.orbit_integrator.profile_model.potentialr(2.0)
			print self.orbit_integrator.profile_model_fast.potentialr(1.0)
			print self.orbit_integrator.profile_model_fast.potentialr(2.0)
			print p.shape
		self.q0 = q[:,0,-1]
		self.q0 = (self.q0[0], self.q0[1], 0.)
		self.p0 = p[:,0,-1]
		self.p0 = (self.p0[0], self.p0[1], 0.)
		#print "3", self.q0
		return q[:,0,:], p[:,0,:]
		
