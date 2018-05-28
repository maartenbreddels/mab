# -*- coding: utf-8 -*-
import pyublas
import numpy
from numpy import *
import mab.gd.gdfast
import sys

class OrbitIntegratorGSL2d(object):
	def __init__(self, profile_model, orbital_points, error_rel, error_abs):
		#self.galaxy = galaxy
		self.profile_model = profile_model
		self.profile_model_fast = profile_model.fast()
		self.orbital_points = orbital_points
		self.integrators = [mab.gd.gdfast.OrbitIntegratorGSL(self.profile_model_fast, error_abs, error_rel, integrator_type) for integrator_type in [mab.gd.gdfast.gsl_odeiv_step_rk8pd, mab.gd.gdfast.gsl_odeiv_step_rkck]]
		
		self.integrator_index = 0
		
	def integrate(self, dt, q0, p0, q=None, p=None):
		if q is None:
			xout, yout = q = zeros((2, len(dt), self.orbital_points))
		else:
			xout, yout = q
		if p is None:
			vxout, vyout = p = zeros((2, len(dt), self.orbital_points))
		else:
			vxout, vyout = p
		oi = self.integrators[self.integrator_index]
		#pyublas.why_not(q0)
		#pyublas.why_not(p0)
		#pyublas.why_not(xout)
		#print dt.shape, q0.shape, p0.shape
		#print xout.shape, yout.shape, vxout.shape, vyout.shape
		#print [k.dtype for k in [dt, q0, p0, xout, yout, vxout, vyout]]
		r = sqrt(q0[:,0]**2 + q0[:,1]**2)
		v = sqrt(p0[:,0]**2 + p0[:,1]**2)
		E0s = 0.5 * v**2 + self.profile_model.potentialr(r)
		oi.integrate(self.orbital_points, dt, q0, p0, xout, yout, vxout, vyout)
		
		r = sqrt(xout[:,-1]**2 + yout[:,-1]**2)
		v = sqrt(vxout[:,-1]**2 + vyout[:,-1]**2)
		Es = 0.5 * v**2 + self.profile_model.potentialr(r)
		#print E0s.shape, Es.shape
		Ediff = E0s - Es
		Ereldiff = abs(Ediff/E0s)
		Erellimit = 1e-3
		if any(Ereldiff > Erellimit):
			print sys.stderr, "energy errors in orbit integration: %s" % ",".join(["%e" % k for k in Ereldiff[Ereldiff > Erellimit]])
		return q, p
		
class OrbitIntegratorAxiGSL(object):
	def __init__(self, profile_model, orbital_points, error_rel, error_abs):
		#self.galaxy = galaxy
		self.profile_model = profile_model
		self.profile_model_fast = profile_model.fast()
		self.orbital_points = orbital_points
		self.integrators = [mab.gd.gdfast.OrbitIntegratorAxiGSL(self.profile_model_fast, error_abs, error_rel, integrator_type) for integrator_type in [mab.gd.gdfast.gsl_odeiv_step_rk8pd, mab.gd.gdfast.gsl_odeiv_step_rkck]]
		
		self.integrator_index = 0
		
	def integrate(self, dt, Lz, q0, p0, q=None, p=None):
		if q is None:
			xout, yout = q = zeros((2, len(dt), self.orbital_points))
		else:
			xout, yout = q
		if p is None:
			vxout, vyout = p = zeros((2, len(dt), self.orbital_points))
		else:
			vxout, vyout = p
		oi = self.integrators[self.integrator_index]
		#pyublas.why_not(q0)
		#pyublas.why_not(p0)
		#pyublas.why_not(xout)
		#print dt.shape, q0.shape, p0.shape
		#print xout.shape, yout.shape, vxout.shape, vyout.shape
		#print [k.dtype for k in [dt, q0, p0, xout, yout, vxout, vyout]]
		#r = sqrt(q0[:,0]**2 + q0[:,1]**2)
		#v = sqrt(p0[:,0]**2 + p0[:,1]**2)
		#E0s = 0.5 * v**2 + self.profile_model.potentialr(r)
		oi.integrate(self.orbital_points, dt, Lz, q0, p0, xout, yout, vxout, vyout)
		
		#r = sqrt(xout[:,-1]**2 + yout[:,-1]**2)
		#v = sqrt(vxout[:,-1]**2 + vyout[:,-1]**2)
		#Es = 0.5 * v**2 + self.profile_model.potentialr(r)
		#print E0s.shape, Es.shape
		#Ediff = E0s - Es
		#Ereldiff = abs(Ediff/E0s)
		#Erellimit = 1e-3
		#if any(Ereldiff > Erellimit):
		#	print sys.stderr, "energy errors in orbit integration: %s" % ",".join(["%e" % k for k in Ereldiff[Ereldiff > Erellimit]])
		return q, p
		
		
class OrbitIntegrator2dGSL_true2d(object):
	def __init__(self, profile_model, orbital_points, error_rel, error_abs):
		#self.galaxy = galaxy
		self.profile_model = profile_model
		self.profile_model_fast = profile_model.fast()
		self.orbital_points = orbital_points
		self.integrators = [mab.gd.gdfast.OrbitIntegrator2dGSL(self.profile_model_fast, error_abs, error_rel, integrator_type) for integrator_type in [mab.gd.gdfast.gsl_odeiv_step_rk8pd, mab.gd.gdfast.gsl_odeiv_step_rkck]]
		
		self.integrator_index = 0
		
	def integrate(self, dt, q0, p0, q=None, p=None, angular_velocity=0.):
		if q is None:
			xout, yout = q = zeros((2, len(dt), self.orbital_points))
		else:
			xout, yout = q
		if p is None:
			vxout, vyout = p = zeros((2, len(dt), self.orbital_points))
		else:
			vxout, vyout = p
		oi = self.integrators[self.integrator_index]
		#pyublas.why_not(q0)
		#pyublas.why_not(p0)
		#pyublas.why_not(xout)
		#print dt.shape, q0.shape, p0.shape
		#print xout.shape, yout.shape, vxout.shape, vyout.shape
		#print [k.dtype for k in [dt, q0, p0, xout, yout, vxout, vyout]]
		#r = sqrt(q0[:,0]**2 + q0[:,1]**2)
		#v = sqrt(p0[:,0]**2 + p0[:,1]**2)
		#E0s = 0.5 * v**2 + self.profile_model.potentialr(r)
		l = [dt, q0, p0, xout, yout, vxout, vyout]
		#print [k.dtype for k in l]
		#print [k.shape for k in l]
		
		#if 0:
		r = sqrt(q0[:,0]**2 + q0[:,1]**2)
		v = sqrt(p0[:,0]**2 + p0[:,1]**2)
		E0s = 0.5 * v**2 +self.profile_model.light_profile.potentialxy(q0[:,0], q0[:,1])
		
		oi.integrate(self.orbital_points, angular_velocity, dt, q0, p0, xout, yout, vxout, vyout)
		
		#r = sqrt(xout[:,-1]**2 + yout[:,-1]**2)
		#v = sqrt(vxout[:,-1]**2 + vyout[:,-1]**2)
		#Es = 0.5 * v**2 + self.profile_model.potentialr(r)
		#print E0s.shape, Es.shape
		#Ediff = E0s - Es
		#Ereldiff = abs(Ediff/E0s)
		#Erellimit = 1e-3
		#if any(Ereldiff > Erellimit):
		#	print sys.stderr, "energy errors in orbit integration: %s" % ",".join(["%e" % k for k in Ereldiff[Ereldiff > Erellimit]])
		r = sqrt(xout[:,-1]**2 + yout[:,-1]**2)
		v = sqrt(vxout[:,-1]**2 + vyout[:,-1]**2)
		Es = 0.5 * v**2 + self.profile_model.light_profile.potentialxy(xout[:,-1], yout[:,-1])
		#print E0s.shape, Es.shape
		Ediff = E0s - Es
		Ereldiff = abs(Ediff/E0s)
		Erellimit = 1e-3
		if any(Ereldiff > Erellimit):
			print sys.stderr, "energy errors in orbit integration: %s" % ",".join(["%e" % k for k in Ereldiff[Ereldiff > Erellimit]])
		return q, p
		
				