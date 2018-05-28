# -*- coding: utf-8 -*-
from numpy import *
from mab.binningtools import bingrid, binrange
import mab.gd.logging as logging

logger = logging.getLogger("gd.nbody.gadget")

		
class Component(object):
	def __init__(self, name, q, p, mass, potential=None):
		self.name = name
		self.q = q
		self.p = p
		self.mass = mass
		self.masses = None
		self.potential = potential

	def clone(self):
		return Component(self.name, self.q * 1., self.p * 1., self.mass, self.potential) # TODO check ifwe need to copy mass and potential

	def rescale(self, scale_mass, scale_size):
		self.q *= 1#0.1#scale_size
		self.p *= 1#0.1**-0.5#10**(1./3)#(scale_mass/scale_size)**0.5
		self.mass *= 1#scale_mass

	def rotate_xz(self, angle):
		x, y, z = self.q
		vx, vy, vz = self.p
		xn = x * cos(angle) + z * sin(angle)
		yn = y
		zn = x * sin(angle) - z * cos(angle)
		vxn = vx * cos(angle) + vz * sin(angle)
		vyn = vy
		vzn = vx * sin(angle) - vz * cos(angle)
		self.q = array([xn, yn, zn])
		self.p = array([vxn, vyn, vzn])

	def translate(self, q, p):
		self.q = self.q + array(q).reshape(3,1)
		self.p = self.p + array(p).reshape(3,1)
	
	def moments_polar(self, N=250):
		xp, yp, zp = self.q
		vx, vy, vz = self.p
		Nperbin = N
		r3d = sqrt(xp**2+yp**2+zp**2)
		#print "r3d min", r3d.min()
		#scale = float(i%5)/4.
		#color = Color(scale, 0, 0)
		#color = 
		#r = sqrt(xp*xp+yp*yp+zp*zp).astype(float64)
		rhosq = (xp*xp+yp*yp);
		rho = sqrt(rhosq);
		vR = (vx*xp + vy*yp)/rho;
		vphi = (vy*xp - vx*yp)/rho;
		#vtheta = (vx*zp*xp/rho + vy*zp*yp/rho - vz*rho)/r;
		#vr[isnan(vr)] = 0
		#vphi[isnan(vphi)] = 0
		#vtheta[isnan(vtheta)] = 0
		#vr[isinf(vr)] = 0
		#vphi[isinf(vphi)] = 0
		#vtheta[isinf(vtheta)] = 0
		logrho, vR, vRvar = self.moment1(log10(rho), vR, N)
		logrho, vphi, vphivar = self.moment1(log10(rho), vphi, N)
		logrho, vz, vzvar = self.moment1(log10(rho), vz, N)
		return logrho, vR, vRvar, vphi, vphivar, vz, vzvar
	
	def moments_spherical(self, N=250, moments=2):
		xp, yp, zp = self.q
		vx, vy, vz = self.p
		Nperbin = N
		r3d = sqrt(xp**2+yp**2+zp**2)
		#print "r3d min", r3d.min()
		#scale = float(i%5)/4.
		#color = Color(scale, 0, 0)
		#color = 
		r = sqrt(xp*xp+yp*yp+zp*zp).astype(float64) + 1e-10
		rhosq = (xp*xp+yp*yp);
		rho = sqrt(rhosq);
		vr = (vx*xp + vy*yp + vz*zp)/r;
		vphi = (vy*xp - vx*yp)/rho;
		vtheta = (vx*zp*xp/rho + vy*zp*yp/rho - vz*rho)/r;
		vr[isnan(vr)] = 0
		vphi[isnan(vphi)] = 0
		vtheta[isnan(vtheta)] = 0
		vr[isinf(vr)] = 0
		vphi[isinf(vphi)] = 0
		vtheta[isinf(vtheta)] = 0
		logrs, mr = self.moment1(log10(r), vr, N, moments)
		logrs, mtheta = self.moment1(log10(r), vtheta, N, moments)
		logrs, mphi = self.moment1(log10(r), vphi, N, moments)
		return logrs, mr, mphi, mtheta
	
	def moment1(self, x, v, Nperbin, moments):
		xbins = []
		momentbins = []
		def moment(x, n):
			if n == 0:
				return len(x)
			if n == 1:
				return mean(x)
			xc = x - mean(x)
			return sum(xc**n)/len(xc)
		for n, xbin, vbin in binrange(Nperbin, x, v):
			xbins.append(mean(xbin))
			momentbins.append([moment(vbin, i) for i in range(moments+1)])
		return array(xbins), array(momentbins).T
class Snapshot(object):
	def __init__(self):
		self.components = []
		self.componentmap = {}
	
	def add_component(self, component):
		self.components.append(component)
		self.componentmap[component.name] = component
		
	def center(self, name="halo", rc=0.3):
		if name not in self.componentmap:
			name = "disk"
		component = self.componentmap[name]
		if component.potential is None:
			logger.warning("component %s has no potential information, cannot center" % name)
			return array([0., 0., 0.]), array([0., 0., 0.])
		center_index = argmin(component.potential)
		center_indices = argsort(component.potential)[:10]
		q0 = mean(component.q[:,center_indices], axis=1)
		#q0 = component.q[:,center_index] * 1. # times 1 to avoid reference
		logger.debug("center: %r" % q0)
		for c in self.components:
			c.translate(-q0, [0, 0, 0])
		x, y, z = component.q
		vx, vy, vz = component.p
		r = (x**2+y**2+z**2)**0.5
		mask = r < rc
		vx0 = mean(vx[mask])
		vy0 = mean(vy[mask])
		vz0 = mean(vz[mask])
		p0 = array([vx0, vy0, vz0])
		for c in self.components:
			c.translate([0, 0, 0], -p0)
		return q0, p0
		
class SnapshotTransformed(Snapshot):
	def __init__(self, snapshot, rotate_xz_angle=0, q0=[0,0,0], p0=[0,0,0], mass_scale=1., size_scale=1.):
		super(SnapshotTransformed, self).__init__()
		self.rotate_xz_angle = rotate_xz_angle
		self.snapshot = snapshot
		self.p0 = p0
		self.q0 = q0
		self.mass_scale = mass_scale
		self.size_scale = size_scale
	
	def load(self):
		self.snapshot.load()
		for component in self.snapshot.components:
			component = component.clone()
			self.add_component(component)
			logger.debug("stds  before: %r %r "% ([std(component.p[k]) for k in range(3)], [std(component.q[k]) for k in range(3)]))
			logger.debug("means before: %r %r "% ([mean(component.p[k]) for k in range(3)], [mean(component.q[k]) for k in range(3)]))
			component.rescale(self.mass_scale, self.size_scale)
			component.rotate_xz(self.rotate_xz_angle)
			component.translate(self.q0, self.p0)
			logger.debug("stds  before: %r %r "% ([std(component.p[k]) for k in range(3)], [std(component.q[k]) for k in range(3)]))
			logger.debug("means before: %r %r "% ([mean(component.p[k]) for k in range(3)], [mean(component.q[k]) for k in range(3)]))

		
class SnapshotFiltered(Snapshot):
	def __init__(self, snapshot, component_filter):
		super(SnapshotFiltered, self).__init__()
		self.snapshot = snapshot
		self.component_filter = component_filter
	
	def load(self):
		self.snapshot.load()
		for component in self.snapshot.components:
			if self.component_filter(component):
				logger.debug("adding component: %s" % component.name)
				component = component.clone()
				self.add_component(component)
			else:
				logger.debug("skipping component: %s" % component.name)

class SnapshotCenter(Snapshot):
	def __init__(self, snapshot):
		super(SnapshotCenter, self).__init__()
		self.snapshot = snapshot
	
	def run(self, args, opts, scope):
		self.load()
		
	def load(self):
		self.snapshot.load()
		for component in self.snapshot.components:
			component = component.clone()
			x, y, z = component.q
			r = sqrt(x**2 + y**2 +z**2)
			mask = r < 1
			print "mean q         ", component.name, mean(component.q, axis=1)
			print "mean q (r<1kpc)", component.name, mean(component.q[:,mask], axis=1)
			print "mean p         ", component.name, mean(component.p, axis=1)
			print "mean p (r<1kpc)", component.name, mean(component.p[:,mask], axis=1)
			component.p -= mean(component.p[:,mask], axis=1).reshape(3,1)
			print "mean q         ", component.name, mean(component.q, axis=1)
			print "mean q (r<1kpc)", component.name, mean(component.q[:,mask], axis=1)
			print "mean p         ", component.name, mean(component.p, axis=1)
			print "mean p (r<1kpc)", component.name, mean(component.p[:,mask], axis=1)
			print 
			self.add_component(component)

class SnapshotsMerge(Snapshot):
	def __init__(self, snapshots, output_snapshot):
		super(SnapshotsMerge, self).__init__()
		self.snapshots = snapshots
		self.output_snapshot = output_snapshot
	
	def load(self):
		for snapshot in self.snapshots:
			snapshot.load()
			
		component_names = {}
		for snapshot in self.snapshots:
			for component in snapshot.components:
				component_names[component.name] = None
		names = component_names.keys()
		for name in names:
			q = zeros((3,0))
			p = zeros((3,0))
			mass = None
			masses = zeros((0))
			for snapshot in self.snapshots:
				if name in snapshot.componentmap:
					component = snapshot.componentmap[name]
					x = component.q[0]
					N = len(x)
					q = concatenate([q, component.q],1)
					#p = concatenate([p, component.p],1)
					p = concatenate([p, component.p],1)
					masses = concatenate([masses, ones(N) * component.mass])
					#if mass is None:
					#	mass = component.mass
					#assert mass == component.mass, "masses should be equal (there are %r %r)" % (mass, component.mass)
			c = Component(name, q, p, None)
			c.masses = masses
			self.add_component(c)
					
					
class Convert(object):
	def __init__(self, input, output):
		self.input = input
		self.output = output				
	def run(self, args, opts, scope):
		self.input.load()
		for component in self.input.components:
			logger.debug("stds  before: %r %r "% ([std(component.p[k]) for k in range(3)], [std(component.q[k]) for k in range(3)]))
			logger.debug("means before: %r %r "% ([mean(component.p[k]) for k in range(3)], [mean(component.q[k]) for k in range(3)]))
			self.output.add_component(component)
		self.output.save()
