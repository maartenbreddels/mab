# -*- coding: utf-8 -*-
from numpy import *
import os
import mab.gd.logging as logging
import mab.random
import sys
import time
import numpy

import mab.parallelize
import mab.utils.numpy
import scipy.integrate
import scipy.optimize

logger = logging.getLogger("gd.df")

from kaplot import *


N_SCALE_RADII = 20

class DFDelta(object):
	def __init__(self, profile_model, storage3d, storage3d_output, dfgrid, grid_logr):
		self.dfgrid = dfgrid
		self.profile_model = profile_model
		self.storage3d = storage3d
		self.storage3d_output = storage3d_output 
		self.rs = self.storage3d_output.finv(storage3d.x)
		self.r_borders = self.storage3d_output.finv(storage3d.xborders)
		self.N_r = len(self.rs)
		self.profile_model_fast = profile_model.fast()
		#self.df = mab.gd.gdfast.DF(self.profile_model_fast, 10**self.logr, 10**self.logrborders)
		self.df = mab.gd.gdfast.DF(self.profile_model_fast, self.rs, self.r_borders)
		
	def do(self):
		self.storage3d_output.init()
		rmin = 1e-6
		Nd = self.dfgrid.dither
		density = zeros(self.N_r)
		varvr = zeros(self.N_r)
		varvt = zeros(self.N_r)
		v22 = zeros(self.N_r)
		v40 = zeros(self.N_r)
		v04 = zeros(self.N_r)
		self.N_r
		for Eindex in range(self.dfgrid.n_I1):
			for i in range(self.dfgrid.dither):
				E = self.dfgrid.subgrid.Es[Eindex*self.dfgrid.dither+i]
				Lmax, rcirc = self.profile_model_fast.Lmax_and_rcirc_at_E(E)
				#print "Lmax, rcirc", Lmax, rcirc
				rmax = self.profile_model_fast.rmax_at_E(E, rcirc)
				#l = self.dfgrid.subgrid.ls[li*self.dfgrid.dither+j]
				for Lindex in range(self.dfgrid.n_I2):
					l1 = self.dfgrid.l_borders[Lindex]
					l2 = self.dfgrid.l_borders[Lindex+1]
					#L = l * Lmax
					L1 = l1 * Lmax
					L2 = l2 * Lmax
					#print l1, l2, L1, L2
					# get ra and rp from lowest angular momentum
					#print "finding ra,rp"
					ra, rp = self.profile_model_fast.get_apo_peri(E, L1, rmin, rcirc, rmax) 
					#print "found ra,rp", ra, rp
					#ra1, rp1 = self.profile_model_fast.get_apo_peri(E, L1) 
					#ra2, rp2 = self.profile_model_fast.get_apo_peri(E, L2) 
					#print "finding moments"
					self.df.momentsL2(density, varvr, varvt, v22, v40, v04, E, L1, L2, rp, ra, False)
					#print "found moments"
					
					grid = self.storage3d_output.moments3dgrid[Lindex, Eindex]
					N = sum(density)/Nd
					density[:] = 0
					varvr[:] = 0
					varvt[:] = 0
					v22[:] = 0
					v40[:] = 0
					v04[:] = 0
					grid[self.storage3d.v00index] += density/N
					grid[self.storage3d.v20index] += varvr/Nd
					grid[self.storage3d.v02index] += varvt/Nd
					grid[self.storage3d.v22index] += v22/Nd
					grid[self.storage3d.v40index] += v40/Nd
					grid[self.storage3d.v04index] += v04/Nd
		print "store"
		self.storage3d_output.save()
		
	def run(self, args, opts, scope):
		#self.do()
		#return 2
		self.storage3d.load()
		ei = int(args[1])
		li = int(args[2])
		E = self.dfgrid.Es[ei]
		l = self.dfgrid.ls[li]
		print "Lmax"
		Lmax = Lmax1 = self.profile_model_fast.Lmax_at_E(E)
		Lmax2 = self.profile_model.Lmax_at_E(E)
		print "Lmax=", Lmax1
		print "Lmax=", Lmax2
		#dsa
		print "finding apo peri"
		print Lmax
		Lmax, rcirc = self.profile_model_fast.Lmax_and_rcirc_at_E(E)
		print "Lmax, rcirc", Lmax, rcirc
		rmax = self.profile_model_fast.rmax_at_E(E, rcirc)
		print "rmax", rmax
		
		#Lmax = self.profile_model.Lmax_at_E(E)
		#print Lmax
		#sys.exit(0)
		#import pdb; pdb.set_trace()
		L = l * Lmax
		#assert Ltest
		print"E,l,L", E, l, L
		print self.storage3d.moments3dgrid.shape
		
		#mobox()
		mozaic(2,2,box)
		
		select(0, 0)
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v00index])
		labels("x", "rho")
		select(1, 0)
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v20index])
		labels("x", "v20")
		select(0, 1)
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v02index])
		labels("x", "v02")
		select(1, 1)
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v40index])
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v04index])#, linestyle="dot", linewidth="2px")
		graph(self.storage3d.moments3dgrid[li, ei, self.storage3d.v22index])
		labels("x", "v4")
		
		rmin = 1e-8
		if 0:
			density = self.rs * 0
			varvr = self.rs * 0
			varvt = self.rs * 0
			print "finding apo peri"
			ra, rp = self.profile_model_fast.get_apo_peri(E, L, 1e-6, rcirc, rmax)
			print "ra, rp", ra, rp
			def test(r):
				return (-L*L/(2*r*r) - (self.profile_model.potentialr(r) - E));
			print "test", test(rp), test(ra)
			self.df.moments(density, varvr, varvt, E, L, rp, ra)
			if 0:
				select(0, 0)
				graph(density, color="red", linestyle="dash")
				select(1, 0)
				graph(varvr, color="red", linestyle="dash")
				select(0, 1)
				graph(varvt, color="red", linestyle="dash")
			
			select(0, 0)
			density_total = self.rs * 0
			varvr_total = self.rs * 0
			varvt_total = self.rs * 0
			for i in range(self.dfgrid.dither):
				E = self.dfgrid.subgrid.Es[ei*self.dfgrid.dither+i]
				Lmax, rcirc = self.profile_model_fast.Lmax_and_rcirc_at_E(E)
				print "Lmax, rcirc", Lmax, rcirc
				rmax = self.profile_model_fast.rmax_at_E(E, rcirc)
				print "rmax", rmax
				for j in range(self.dfgrid.dither):
					l = self.dfgrid.subgrid.ls[li*self.dfgrid.dither+j]
					L = l * Lmax
					print "TEST 1"
					ra, rp = self.profile_model_fast.get_apo_peri(E, L, rmin, rcirc, rmax)
					print "TEST 2"
					density = self.rs * 0
					varvr = self.rs * 0
					varvt = self.rs * 0
					if 1:
						self.df.moments(density, varvr, varvt, E, L, rp, ra)
						density_total += density
						varvr_total += varvr
						varvt_total += varvt
			print "done"
			N = self.dfgrid.dither**2
			select(0, 0)
			graph(density_total/N, color="blue", linestyle="dot")
			select(1, 0)
			graph(varvr_total/N, color="blue", linestyle="dot")
			select(0, 1)
			graph(varvt_total/N*2, color="blue", linestyle="dot")
		if 1:
			density_total = self.rs * 0
			varvr_total = self.rs * 0
			varvt_total = self.rs * 0
			v22_total = self.rs * 0
			v40_total = self.rs * 0
			v04_total = self.rs * 0
			for i in range(self.dfgrid.dither):
				print "dither", i
				E = self.dfgrid.subgrid.Es[ei*self.dfgrid.dither+i]
				#Lmax = self.profile_model_fast.Lmax_at_E(E)
				Lmax, rcirc = self.profile_model_fast.Lmax_and_rcirc_at_E(E)
				print "Lmax, rcirc", Lmax, rcirc
				rmax = self.profile_model_fast.rmax_at_E(E, rcirc)
				print "rmax", rmax
				#l = self.dfgrid.subgrid.ls[li*self.dfgrid.dither+j]
				l1 = self.dfgrid.l_borders[li]
				l2 = self.dfgrid.l_borders[li+1]
				#L = l * Lmax
				L1 = l1 * Lmax
				L2 = l2 * Lmax
				print l1, l2, L1, L2
				# get ra and rp from lowest angular momentum
				print "getting apo/peri"
				#import pdb; pdb.set_trace() 
				ra, rp = self.profile_model_fast.get_apo_peri(E, L1, rmin, rcirc, rmax)
				print "found apo/peri", ra, rp
				#ra1, rp1 = self.profile_model_fast.get_apo_peri(E, L1) 
				#ra2, rp2 = self.profile_model_fast.get_apo_peri(E, L2) 
				density = self.rs * 0
				varvr = self.rs * 0
				varvt = self.rs * 0
				v22 = self.rs * 0
				v40 = self.rs * 0
				v04 = self.rs * 0
				try:
					self.df.momentsL2(density, varvr, varvt, v22, v40, v04, E, L1, L2, rp, ra, False)
				except:
					self.df.momentsL2(density, varvr, varvt, v22, v40, v04, E, L1, L2, rp, ra, True)
					raise
				density_total += density
				varvr_total += varvr
				varvt_total += varvt
				v40_total += v40
				v04_total += v04
				v22_total += v22
			N = sum(density_total)#self.dfgrid.dither
			print density_total
			select(0, 0)
			graph(density_total/N, color="orange", linestyle="dash")
			select(1, 0)
			graph(varvr_total/N, color="orange", linestyle="dash")
			select(0, 1)
			graph(varvt_total/N, color="orange", linestyle="dash")
			select(1, 1)
			graph(v40_total/N, color="orange", linestyle="dash")
			graph(v04_total/N, color="orange", linestyle="dash")
			graph(v22_total/N, color="orange", linestyle="dash")
		
		draw()
		
		
		return 2
		

import mab.gd.nbody.snapshot
class Snapshot(mab.gd.nbody.snapshot.Snapshot):
	def __init__(self, df_samples):
		super(Snapshot, self).__init__()
		self.df_samples = df_samples
		
	#def run(self, args, opts, scope):
	#	self.df_samples.load()
	
	def load(self):
		self.df_samples.load()
		N = len(self.df_samples.df_samples)
		q = array([self.df_samples.df_samples.x, self.df_samples.df_samples.y, self.df_samples.df_samples.z])
		p = array([self.df_samples.df_samples.vx, self.df_samples.df_samples.vy, self.df_samples.df_samples.vz])
		mass = self.df_samples.df.galaxy.light_model.light_profile.enclosed_mass(float('inf'))
		logger.info("total mass: %e, mass per particle: %e" % (mass, mass/N))
		c = mab.gd.nbody.snapshot.Component("disk", q, p, mass/N)
		self.add_component(c)
		
		N = len(self.df_samples.df_samples_dm)
		q = array([self.df_samples.df_samples_dm.x, self.df_samples.df_samples_dm.y, self.df_samples.df_samples_dm.z])
		p = array([self.df_samples.df_samples_dm.vx, self.df_samples.df_samples_dm.vy, self.df_samples.df_samples_dm.vz])
		rs = self.df_samples.df.galaxy.profile_model.dm_profile.rs
		mass = self.df_samples.df.galaxy.profile_model.dm_profile.enclosed_mass(100)
		logger.info("total mass: %e, mass per particle: %e" % (mass, mass/N))
		c = mab.gd.nbody.snapshot.Component("halo", q, p, mass/N)
		self.add_component(c)
		
class SnapshotPlus(Snapshot):
	def __init__(self, df_samples, orbits, mass):
		super(SnapshotPlus, self).__init__(df_samples)
		self.orbits = orbits
		self.mass = mass
	
	def load(self):
		super(SnapshotPlus, self).load()
		#xs = []
		#vys = []
		ps = []
		qs = []
		for orbit in self.orbits:
			logger.info("orbit: x0=%f vt0=%f" % (orbit.x0, orbit.vt0))
			qs.append(orbit.q0)
			ps.append(orbit.p0)
			#xs.append(orbit.x0)
			#vxs.append(orbit.vx0)
		ps = array(ps).T
		qs = array(qs).T
		logger.info("total mass for orbits: %e, mass per particle: %e" % (self.mass*len(qs.T), self.mass))
		c = mab.gd.nbody.snapshot.Component("bulge", qs, ps, self.mass)
		self.add_component(c)
		
class DFSamplingSeperable(object):
	def __init__(self, df, N, N_dm=None, sample_dm=False, postfix=""):
		self.df = df
		self.N = N
		self.N_dm = N_dm
		self.sample_dm = sample_dm
		if self.sample_dm and self.N_dm is None:
			self.N_dm = self.N
		
		self.dirname = os.path.join(df.modelpath, "data")
		self.filename = os.path.join(self.dirname, "3d" +postfix +".npy")
		self.filename_dm = os.path.join(self.dirname, "3d" +postfix +"_dm.npy")
		
		names = "x,y,z,r3d,phi,theta,vx,vy,vz,vr,vphi,vtheta,sigmar,sigmat,Ekin,Epot,E,Lx,Ly,Lz,L,Lmax".split(",")
		names = [k.strip() for k in names]
		types = ["f4"] * len(names)
		self.dtype = zip(names, types)
		
	def run(self, args, opts, scope):
		self.init()
		self.sample(scope)
		self.save()
			
	def init(self):
		self.df.load()
	
	def load(self):
		self.df.load()
		self.df_samples = load(self.filename).view(recarray)
		if self.sample_dm:
			self.df_samples_dm = load(self.filename_dm).view(recarray)
		
		
	def save(self):
		if not os.path.exists(self.dirname):
			logger.info("creating directory: %s" % self.dirname)
			os.makedirs(self.dirname)
		logger.info("saving samples as: %s" % self.filename)
		save(self.filename, self.df_samples)
		
		if self.sample_dm:
			logger.info("saving dark matter samples as: %s" % self.filename_dm)
			save(self.filename_dm, self.df_samples_dm)
	
	
		
	def sample(self, scope, seedoffset=None):
		self.df_samples = mab.utils.numpy.mmapzeros((self.N), dtype=self.dtype).view(recarray)
		
		def seed(id_number):
			if seedoffset is not None:
				s = id_number + seedoffset
			else:
				s = os.getpid() * 10000 + int(time.time() * 100000)
			print "seed for %d is %d (%d is pid)" % (id_number, s, os.getpid())
			numpy.random.seed(s)
		
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"], init=seed)
		def sample_wrapper(index):
			self.sample1(index)
		indices = range(self.N)
		sample_wrapper(indices)
		#self.sample1(0)
		
		if self.sample_dm:
			self.df_samples_dm = mab.utils.numpy.mmapzeros((self.N_dm), dtype=self.dtype).view(recarray)
			@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"], init=seed)
			def sample_wrapper(index):
				self.sample1(index, dark_matter=True)
			indices = range(self.N_dm)
			sample_wrapper(indices)
		
	def sample1(self, index, dark_matter=False):
		#self.xs = []
		#self.ys = []
		#self.zs = []
		#self.vrs = []
		#self.vphis = []
		#self.vthetas = []
		
		galaxy = self.df.galaxy
		Lmaxs = self.df.Lmaxs
		accepted = 0
		if dark_matter:
			fE = self.df.fE_dm
		else:
			fE = self.df.fE
		if dark_matter:
			N = self.N_dm
		else:
			N = self.N
		Es = -self.df.es
		
		if 1:
			#assert galaxy.beta < 0
			# if beta < 0, pmax is max for L is Lmax
			#print f1
			#print Lmaxs
			if dark_matter:
				pmax = max(fE)
				df_samples = self.df_samples_dm
			else:
				pmax = max(fE*galaxy.fL(Lmaxs))
				df_samples = self.df_samples
			tries = 0
			#dsas
			#while len(self.xs) < N:
			accepted_r = False
			if 1:
				while not accepted_r:
					if dark_matter:
						if hasattr(galaxy.profile_model.dm_profile, "rs"):
							scale = galaxy.profile_model.dm_profile.rs
						else:
							scale = 10.
						r = galaxy.profile_model.dm_profile.sample_r(rmax=100) #scale*N_SCALE_RADII)
						#r = galaxy.profile_model.dm_profile.sample_r(rmax=scale*100)
					else:
						r = galaxy.light_model.sample_r()
					x, y, z = mab.random.random_spherical(R=0, N=None, R0=r)
					Epot = galaxy.potentialr(r)
					vmax = sqrt(-2*Epot)
					#print "r", r, "vmax =", vmax
					
					imin = argmin(abs(Epot-Es))
					
					Lmax = sqrt(-2*Epot)*r
					if dark_matter:
						pmax = max(fE[imin:])
					else:
						pmax = max(fE[imin:]*galaxy.fL(Lmax)) #Lmax**(-2*galaxy.beta))
					#print pmax
					#u = random.random()
					#L = Lmax * u**(1.5)
					
					accepted_vel = False
					n_vel = 0
					while not accepted_vel:
					#if 0:
						vr, vphi, vtheta = mab.random.random_spherical(R=vmax, N=None)
						vt = sqrt(vphi**2+vtheta**2)
						L = vt * r
						#print L
						v = sqrt(vr**2+vt**2)
						
						E = 0.5 * v**2 + Epot 
						
						i = argmin(abs(E-Es))
						#print i, 
						
						#p = fE[i] * L ** (-2*beta)
						if dark_matter:
							p = fE[i]
						else:
							p = fE[i] * galaxy.fL(L)
						#assert beta == -0.5
						tries += 1
						u = random.random()
						#print p/pmax
						if u < p/pmax:
							accepted_vel = True
							#print ".", # retry
							#sys.stdout.flush()
						else:
							n_vel += 1
							if n_vel > 1000000:
								print "give up this r", r
								break
					if 0:
						vt = L/r
						angle = random.random() * 2 * pi
						vphi = sin(angle) * vt
						vtheta = cos(angle) * vt
						
						
						
						imin = argmin(abs(0.5*vt**2+Epot-Es))
						pmax = max(fE[imin:]) #*L**(-2*galaxy.beta)*L)
						vmax = sqrt(-2*Epot)
						if not accepted_vel:
							#v = random_spherical(R=1., N=1)[0] * (vmax-vt) + vt
							v = random.random() ** (1./3.) * (vmax-vt) + vt
							vr = sqrt(v**2-vt**2)
							#vrmax = sqrt(-2*Epot-vt**2)
							#vr = (random.random() * 2 - 1) * vrmax 
							#v = sqrt(vr**2+vt**2)
							E = 0.5 * v**2 + Epot 
							assert E < 0
							i = argmin(abs(E-Es))
							p = fE[i] #* L ** (-2*beta)*L
							assert beta == -0.5
							assert imin <= i
							tries += 1
							u = random.random()
							if u < p/pmax:
								accepted_vel = True
								print ".", # retry
								sys.stdout.flush()
							else: # retry
								#print "*", 
								n_vel += 1
								if n_vel > 1000000:
									print "give up this r", r
									break
					if accepted_vel:
						accepted_r = True
					else:
						pass # retry with other r
					
				df_samples.x[index] = x
				df_samples.y[index] = y
				df_samples.z[index] = z
				df_samples.vr[index] = vr
				df_samples.vphi[index] = vphi
				df_samples.vtheta[index] = vtheta
				
				r = sqrt(x**2+y**2+z**2)
				
				v = sqrt(vr**2+vtheta**2+vphi**2)
				Epot = galaxy.potentialr(r)
				E = 0.5 * v**2 + Epot
				Lmax, r_at_Lmax = galaxy.profile_model.L_and_r_max_at_E(E)
				#r_at_Lmax, Lmax = r*0, r*0
				
				
				vlength1 = v
				vlength2 = sqrt(vr**2+vtheta**2+vphi**2)
				vdiff = abs((vlength1-vlength2)/vlength1)
				#i = argmax(vdiff)
				#print vdiff[i], vlength1[i], vlength2[i]
				assert (vdiff) < 0.0001
				
				theta = arccos(z/r)
				phi = arctan2(y,x)
				
				#phi = random.random(N) * 2 * pi
				#costheta = random.random(N) * 2 - 1
				#theta = arccos(costheta)
				#sintheta = sqrt(1-costheta**2)
				#cosphi = cos(phi)
				#sinphi = sin(phi)
				
				# convert to cartesian
				#x = r * sintheta * cosphi
				#y = r * sintheta * sinphi
				#z = r * costheta
				
				# rotation 'matrix' to convert spherical to cartesian coordinates
				R1 = array([cos(phi) * sin(theta), sin(phi)*sin(theta), cos(theta)])
				R2 = array([-sin(phi), cos(phi), 0])
				R3 = array([cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta)])
				
				#R1 = array([cos(phi) * sin(theta), -sin(phi), cos(phi)*cos(theta)])
				#R2 = array([sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta)])
				#R3 = array([cos(theta), zeros(N)*1.0, -sin(theta)])
				
				# rotate velocities to cartesian
				#vx, vy, vz = vr, vphi,vtheta
				
				#vr, vphi, vtheta = vx * R1 + vy * R2 + vz * R3
				
				#print vr.shape, vphi.shape, vtheta.shape, R1.shape, R2.shape, R3.shape
				vx, vy, vz = vr * R1 + vphi * R2 + vtheta * R3
				
				# velocity sanity check, |v| should be the same in every coordinate system
				#print vx, vy, vz
				vlength1 = numpy.sqrt(vr**2+vphi**2+vtheta**2)
				vlength2 = numpy.sqrt(vx**2+vy**2+vz**2)
				vdiff = abs((vlength1-vlength2)/vlength1)
				#i = argmax(vdiff)
				#print vdiff[i], vlength1[i], vlength2[i]
				assert (vdiff) < 0.0001
					
				# the angular momentum
				Lx = y*vz - z*vy
				Ly = z*vx - x*vz
				Lz = x*vy - y*vx
				L = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
				
				vt = numpy.sqrt(vphi**2+vtheta**2)
				Lcheck2 = r * vt
				#assert max(abs((L-Lcheck1)/L)) < 0.0001
				#assert max(abs((L-Lcheck2)/L)) < 0.0001
				
				
				#jeans = JeansAnisotropicConstant(stellar_profile, profiles, galaxymodel.beta)
				#jeans = JeansAnisotropicConstant(dm_profile, [dm_profile], galaxymodel.beta)
				sigma_r = r * 0 #$jeans.sigmar(r, smart=False) #/ 3**0.5
				sigma_t = r * 0 #sqrt((1-galaxymodel.beta))*sigma_r
				#sigma_r = jeans.sigmar(r, smart=False) #/ 3**0.5
				#sigma_t = sqrt((1-galaxymodel.beta))*sigma_r
				
				Ekin = E - Epot
				#print Ekin
				# we're almost done, write to disk (and calculate Lmax)
				sys.stdout.flush()
				#print >>f, "x, y, z, r3d, phi, theta, vx, vy, vz, vr, vphi, vtheta, sigmar, sigmat, Ekin, Epot, E"
				
				df_samples.r3d[index] = r
				df_samples.phi[index] = phi
				df_samples.theta[index] = theta
				df_samples.vx[index] = vx
				df_samples.vy[index] = vy
				df_samples.vz[index] = vz
				df_samples.Ekin[index] = Ekin
				df_samples.E[index] = E
				df_samples.Epot[index] = Epot
				df_samples.Lmax[index] = Lmax
				df_samples.L[index] = L
				df_samples.Lx[index] = Lx
				df_samples.Ly[index] = Ly
				df_samples.Lz[index] = Lz
				
				#, phi, theta, vx, vy, vz, vr, vphi, vtheta, sigmar, sigmat, Ekin, Epot, E, Lx, Ly, Lz, L, Lmax"
				
				
				
			#	n = len(self.xs)
			#	if (n > 0) and (n % (N/100) == 0):
			#		logger.info("%.2f%%, acceptence ratio: %.2f%%" % (n*100./N, 100.*n/tries))
			#logger.info("acceptence ratio: %.2f" % (tries*100./N))
		
		

class DFSphericalSeperable(object):
	def __init__(self, modelpath, galaxy, postfix=""):
		self.modelpath = modelpath
		self.galaxy = galaxy
		self.dirname = os.path.join(modelpath, "df")
		self.filename_fE = os.path.join(self.dirname, "fE" +postfix +".npy")
		self.filename_fE_dm = os.path.join(self.dirname, "fE" +postfix +"_dm.npy")
		self.filename_E = os.path.join(self.dirname, "fE-E" +postfix +".npy")
		self.filename_Lmax = os.path.join(self.dirname, "Lmax" +postfix +".npy")
		
	def run(self, args, opts, scope):
		self.prepare()
		self.save()

	def load(self):
		self.fE = load(self.filename_fE)
		self.fE_dm = load(self.filename_fE_dm)
		self.es = -load(self.filename_E)
		self.Lmaxs = load(self.filename_Lmax)
		
	def save(self):
		if not os.path.exists(self.dirname):
			logger.info("creating directory: %s" % self.dirname)
			os.makedirs(self.dirname)
		logger.info("saving fE(E) as: %s" % self.filename_fE)
		save(self.filename_fE, self.fE)
		logger.info("saving dark matter fE(E) as: %s" % self.filename_fE_dm)
		save(self.filename_fE_dm, self.fE_dm)
		logger.info("corresponding energy as: %s" % self.filename_E)
		save(self.filename_E, -self.es)
		logger.info("corresponding Lmax as: %s" % self.filename_Lmax)
		save(self.filename_Lmax, self.Lmaxs)
	
	def prepare(self):
		logger.info("computing fE(E)...")
		galaxy = self.galaxy
		dlogr = 0.01
		logrs = arange(-3, 3, dlogr)
		rs = 10**logrs
		rho_star = galaxy.light_model.densityr(rs)
		rho_dm = galaxy.profile_model.dm_profile.densityr(rs)
		E0 = abs(galaxy.potentialr(1e-6))
		psis = -galaxy.potentialr(rs)
		es = psis[::] # go from low binding energy to high
		f_star = zeros(len(es))
		f_dm = zeros(len(es))
		rhogrid = zeros((len(es), len(rs)))
		rhogrid_dm = zeros((len(es), len(rs)))
		
		for i in range(len(es)):
			e = es[i]
			#print e
			# what is the density profile corresponding to this energy?
			i1 = 0
			i2 = i #len(rho_star) -i-1
			#i2 = argmin(abs(psis-e))
			#if psis[i2] < e:
			#	i2 -= 1
				
			
			#assert all(psis[i1:i2] >= e)
			#print sqrt(2*(psis[i1:i2]-e))
			r = rs[i1:i2]
			Lmax = sqrt(-(e-psis[i1:i2])*2)*r
			da = 0.01#/4
			for a in arange(da/2, 1-da/2, da):
				L = Lmax * a
				#print L
				#rhogrid[i][i1:i2] = 4 * pi * sqrt(2*(psis[i1:i2]-e)) * L **(-2*beta)
				assert all(psis[i1:i2]-e > 0)
				# da = const, Lda=dL
				#rhogrid[i][i1:i2] += 4 * pi / sqrt(2*(psis[i1:i2]-e)-L**2/(r**2)) * L / r**2 *  L * L **(-2*beta)
				rhogrid[i][i1:i2] += 4 * pi / sqrt(2*(psis[i1:i2]-e)-L**2/(r**2)) * L / r**2 *  L * galaxy.fL(L)
				rhogrid_dm[i][i1:i2] += 4 * pi / sqrt(2*(psis[i1:i2]-e)-L**2/(r**2)) * L / r**2 *  L # asssume DM is isotropc
			# if beta = 0
			#rhogrid[i][i1:i2] = 4 * pi *  sqrt(2*(psis[i1:i2]-e))
			if i > 0:
				if rhogrid[i][i-1] <= 0:
					import pdb; pdb.set_trace()
		
		
		rho_reproduced = rho_star * 0
		rho_reproduced_dm = rho_dm * 0
		#select(1,1)
		#indexedimage(log10(rhogrid/rhogrid.max()+1e-10), datamin=-10, datamax=0)
		#draw()
		
		# now solve the df 'coefficients' so that it generates the density distribution
		
		for i in range(1, len(es))[::-1]:
			#print rhogrid[i][i-1]
			if i > 0: 
				assert rhogrid[i][i-1] > 0
				assert rhogrid_dm[i][i-1] > 0
			
			if rhogrid[i][i-1] > 0:
				extra = 1. #/((rho_star[i-1] - rho_star[i])/(es[i-1] - es[i]))
				extra = 1./(es[i-1] - es[i])
				f_star[i] = (rho_star[i-1] - rho_reproduced[i-1])/rhogrid[i][i-1] * extra
				f_dm[i] = (rho_dm[i-1] - rho_reproduced_dm[i-1])/rhogrid_dm[i][i-1] * extra
				#f_star[i] = (rho_star[i-1]/rhogrid[i][i-1] - rho_reproduced[i-1]/rhogrid[i][i-1]) * extra
				#f_star[i] = min((rho_star[0:i]/rhogrid[i][0:i] - rho_reproduced[0:i]/rhogrid[i][0:i]) * extra)
				#print i, es[i-1]
				#assert f_star[i] > 0
				if f_star[i] < 0:
					f_star[i] = 0
				if f_dm[i] < 0:
					f_dm[i] = 0
				#print f_star[i], rhogrid[i][i-1]
				#print "extra", extra, f_star[i]
				#f_star[i] *= extra
				#print "extra", extra, f_star[i]
				#f_star[i] = rho_star[i-1]/rhogrid[i][i-1] - rho_reproduced[i-1]/rhogrid[i][i-1]
				rho_reproduced += f_star[i]  * rhogrid[i] / extra 
				rho_reproduced_dm += f_dm[i]  * rhogrid_dm[i] / extra 
		
		logger.info("computing Lmax...")
		self.Lmaxs = array([galaxy.profile_model.Lmax_at_E(E) for E in -es])
		#self.Lmaxs_dm = array([galaxy.profile_model.Lmax_at_E(E) for E in -es])
		logger.info("computing Lmax done")
		self.fE = f_star
		self.fE_dm = f_dm
		self.es = es
		
		return
		
		#box()
		mozaic(2,3,box)
		
		mask = rho_reproduced > 0
		x = logrs[mask]
		graph(x, log10(rho_reproduced[mask]), color="red")
		f = scipy.interpolate.interp1d(x, log10(rho_reproduced[mask]), kind='linear')
		labels("log r/kpc", "log10 rho/Msol/kpc<sup>3</sup>")
		select(0,1)
		print x
		x = arange(-3, 3-0.1, 0.001)
		print x
		y1 = 10**f(x)
		y2 = (galaxy.light_model.densityr(10**x))
		graph(x, (log10(y1)-log10(y2)), color="red")
		labels("log r / kpc", "&delta;rho/rho")
		s = 0.0001
		ylim(-s, s)
		
		select(0,2)
		x = logrs[mask]
		f = scipy.interpolate.interp1d(x, log10(rho_reproduced[mask]), kind='nearest')
		x = arange(-3, 3-0.1, 0.001)
		r1s = x[0:-1]
		r2s = x[1:]
		enclosed_mass_est = cumsum([scipy.integrate.quad(lambda r: (10**f(r))*(10**r)**2, r1, r2) for r1, r2 in zip(r1s, r2s)])
		enclosed_mass = cumsum([scipy.integrate.quad(lambda r: galaxy.light_model.densityr(10**r)*(10**r)**2, r1, r2) for r1, r2 in zip(r1s, r2s)])
		#import pdb
		#pdb.set_trace()
		enclosed_mass_est /= enclosed_mass_est.max()
		enclosed_mass /= enclosed_mass.max()
		graph(x, enclosed_mass_est-enclosed_mass, color="red")
		labels("log r / kpc", "&delta;Menc/Mtot")
		#graph(x, enclosed_mass_est, color="red")
		#graph(x, enclosed_mass, color="blue")
		
		
		select(1,0)
		mask = rho_reproduced_dm > 0
		graph(logrs[mask], log10(rho_reproduced_dm[mask]), color="black")
		labels("log r/kpc", "log10 rho/Msol/kpc<sup>3</sup>")
		
		select(1,1)
		x = logrs[mask]
		f = scipy.interpolate.interp1d(x, log10(rho_reproduced_dm[mask]), kind='linear')
		x = arange(-3, 3-0.1, 0.001)
		y1 = 10**f(x)
		y2 = (galaxy.profile_model.dm_profile.densityr(10**x))
		graph(x, (y1-y2)/y2, color="black")
		s = 0.001
		ylim(-s, s)
		labels("log r / kpc", "&delta;rho/rho")
		
		
		select(1,2)
		x = logrs[mask]
		f = scipy.interpolate.interp1d(x, log10(rho_reproduced_dm[mask]), kind='nearest')
		x = arange(-3, 3-0.1, 0.001)
		r1s = x[0:-1]
		r2s = x[1:]
		enclosed_mass_est = cumsum([scipy.integrate.quad(lambda r: (10**f(r))*(10**r)**2, r1, r2) for r1, r2 in zip(r1s, r2s)])
		enclosed_mass = cumsum([scipy.integrate.quad(lambda r: galaxy.profile_model.dm_profile.densityr(10**r)*(10**r)**2, r1, r2) for r1, r2 in zip(r1s, r2s)])
		#import pdb
		#pdb.set_trace()
		enclosed_mass_est /= enclosed_mass_est.max()
		enclosed_mass /= enclosed_mass.max()
		graph(x, enclosed_mass_est-enclosed_mass, color="black")
		labels("log r / kpc", "&delta;Menc/Mtot")
		
		
		draw()
		
epsabs=1e-12
epsrel=1e-12
class DFSphericalSeperable2(object):
	def __init__(self, modelpath, galaxy):
		self.modelpath = modelpath
		self.galaxy = galaxy
		self.postfix = "_2"
		self.dirname = os.path.join(modelpath, "df")
		self.filename_fE = os.path.join(self.dirname, "fE" +self.postfix +".npy")
		self.filename_E = os.path.join(self.dirname, "fE-E" +self.postfix +".npy")
		self.filename_Lmax = os.path.join(self.dirname, "Lmax" +self.postfix +".npy")
		self.filename_rhogrid = os.path.join(self.dirname, "rhogrid" +self.postfix +".npy")
		self.filename_m2_grid = os.path.join(self.dirname, "m2_grid" +self.postfix +".npy")
		Nr = 50
		NE = 40
		dlogr = 0.01
		self.logrs, self.logr_borders = mab.utils.numpy.range_and_borders(Nr, -2, 2)
		self.logrsE, self.logr_bordersE = mab.utils.numpy.range_and_borders(NE, -3, 3)
		#logrs, logr_borders = self.logrs, self.logr_borders
		self.rs = 10**self.logrs
		self.r_borders = 10**self.logr_borders
		self.rsE = 10**self.logrsE
		self.r_bordersE = 10**self.logr_bordersE
		print self.r_borders[1:] - self.r_borders[:-1] 
		#self.rho_targets = array([scipy.integrate.quad(lambda x: self.galaxy.light_model.densityr(x), r1, r2, epsabs=epsabs, epsrel=epsrel)[0]/(r2-r1) for r1, r2 in zip(self.r_borders[:-1], self.r_borders[1:])])
		
		self.rho_targets = galaxy.light_model.densityr(self.rs) * 4*pi*self.rs**2
		
		#E0 = abs(galaxy.potentialr(1e-6))
		
		
		psis = self.galaxy.potentialr(self.rsE)
		psi_borders = self.galaxy.potentialr(self.r_bordersE)
		self.es = psis[::] # go from low binding energy to high
		self.e_borders = psi_borders[::]
		self.dEs = self.e_borders[1:] - self.e_borders[:-1]
		#print "Energies", self.es
		

	def load(self):
		self.fE = load(self.filename_fE)
		#self.es = load(self.filename_E)
		self.Lmaxs = load(self.filename_Lmax)
		self.rhogrid = load(self.filename_rhogrid)
		self.m2_grid = load(self.filename_m2_grid)
		#self.rhogrida = load(self.filename_rhogrida)
		#self.rhogridv = load(self.filename_rhogridb)
		
	def save(self):
		#return
		if not os.path.exists(self.dirname):
			logger.info("creating directory: %s" % self.dirname)
			os.makedirs(self.dirname)
		logger.info("saving fE(E) as: %s" % self.filename_fE)
		save(self.filename_fE, self.fE)
		#logger.info("corresponding energy as: %s" % self.filename_E)
		#save(self.filename_E, -self.es)
		logger.info("corresponding Lmax as: %s" % self.filename_Lmax)
		save(self.filename_Lmax, self.Lmaxs)
		logger.info("saving rhogrid as: %s" % self.filename_rhogrid)
		save(self.filename_rhogrid, self.rhogrid)
		#save(self.filename_rhogridb, self.rhogridb)
		#save(self.filename_rhogrida, self.rhogrida)
		logger.info("saving m2 grid as: %s" % self.filename_m2_grid)
		save(self.filename_m2_grid, self.m2_grid)
	
	def prepare(self):
		logger.info("computing fE(E)...")
		galaxy = self.galaxy
		#logrs, logr_borders = self.logrs, self.logr_borders
		#rs = 10**logrs
		#r_borders = 10**logr_borders
		#rho_targets = [scipy.integrate.quad(lambda x: galaxy.light_model.densityr(x), r1, r2)[0]/(r2-r1) for r1, r2 in zip(r_borders[:-1], r_borders[1:])]
		
		#rho_targets = galaxy.light_model.densityr(rs)
		
		
		
		self.rhogrid = mab.utils.numpy.mmapzeros((len(self.es), len(self.rs)))
		self.rhogrida = mab.utils.numpy.mmapzeros((len(self.es), len(self.rs)))
		self.rhogridb = mab.utils.numpy.mmapzeros((len(self.es), len(self.rs)))
		self.m2_grid = mab.utils.numpy.mmapzeros((len(self.es), len(self.rs)))

		logger.info("computing Lmax...")
		self.Lmaxs = array([galaxy.profile_model.Lmax_at_E(E) for E in self.es])
		x = array([galaxy.profile_model.L_and_r_max_at_E(E) for E in self.e_borders])
		self.Lmax_borders = x[:,0]
		self.r_at_Lmax_borders = x[:,1]
		logger.info("computing Lmax done")
		
		cores=100
		info=True
		@mab.parallelize.parallelize(cores=cores, info=info)
		def do_rho(i):
			#print i+1, "of", len(es)
			e = self.es[i]
			e1 = self.e_borders[i]
			e2 = self.e_borders[i+1]
			de = self.e_borders[i+1] - self.e_borders[i]
			#print e
			# what is the density profile corresponding to this energy?
			#i1 = 0
			#i2 = i #len(rho_star) -i-1
			#i2 = argmin(abs(psis-e))
			#if psis[i2] < e:
			#	i2 -= 1
				
			
			#assert all(psis[i1:i2] >= e)
			#print sqrt(2*(psis[i1:i2]-e))
			#r = rs[i1:i2]
			#Lmax = sqrt(-(e-psis[i1:i2])*2)*r
			Lmax = self.Lmaxs[i]
			
			if 0:
				nL = 20
				for j in range(len(self.rs))[10:]:
					r = self.rs[j]
					print "r=", r, j
					e = self.es[i]
					r1 = self.r_borders[j]
					r2 = self.r_borders[j+1]
					e1 = self.e_borders[i]
					e2 = self.e_borders[i+1]
					Lmax = self.Lmaxs[i]
					Lmax1 = self.Lmax_borders[i]
					Lmax2 = self.Lmax_borders[i+1]
					r_at_Lmax1 = self.r_at_Lmax_borders[i]
					r_at_Lmax2 = self.r_at_Lmax_borders[i+1]
					#Lmax = sqrt(2*(e-galaxy.potentialr(r)))*r
					
					if 2*(e2-galaxy.potentialr(r1)) < 0: # largest energy cannot come to smallest bin
						pass
					else:
						if 2*(e2-galaxy.potentialr(r2)) < 0: # largest energy cannot cover the whole bin
							r2 = self.r_bordersE[i+1] * 0.99999 # this is r2 such that phi(r2) == e2
						def lmax(e, r1, r2):
							def f(r, E=e):
								r = abs(r)
								return -r**2*(2*(E-galaxy.potentialr(r)))
							r = abs(scipy.optimize.fminbound(f, r1, r2)[0])
							print "(",r,")"
							return (-f(r))**0.5
							
						print ">>>>>>>", 2*(e2-galaxy.potentialr(r2))*r2**2
						print "Lmax", lmax(e, r_at_Lmax1, r_at_Lmax2), Lmax, Lmax1, Lmax2
						print " ", r_at_Lmax1, r_at_Lmax2 
						for k in range(nL):
							s1 = (k+0.)/nL
							s2 = (k+1.)/nL
							L1e1 = Lmax1 * s1
							L2e1 = Lmax1 * s2
							L1e2 = Lmax2 * s1
							L2e2 = Lmax2 * s2
							r1_ = r1
							r2_ = r2
							e1_ = e1
							e2_ = e2
							def canreach(e, l, r):
								return (2*(e2-galaxy.potentialr(r1_)) - L1e2**2/r1_**2) > 0
							#if (2*(e2-galaxy.potentialr(r1_)) - L1e2**2/r1_**2) < 0: # largest energy, smallest L cannot come to smallest bin
							if (not canreach(e2, L1e2, r1)): 
								pass
							else:
								if (2*(e2-galaxy.potentialr(r2_)) - L1e2**2/r2_**2) < 0: # largest energy, smallest angular momentum cannot cover the whole bin
									def f(r, E=e):
										r = abs(r)
										return 2*(e2-galaxy.potentialr(r)) - L1e2**2/r**2
									r2_ = abs(scipy.optimize.fsolve(f, r1_)[0]) * 0.9999
									assert (2*(e2-galaxy.potentialr(r2_)) - L1e2**2/r2_**2) > 0.
								if (2*(e1-galaxy.potentialr(r1_)) - L1e1**2/r1_**2) < 0: # smallest energy, smallest angular momentum cannot reach inner radius
									e1_ = (galaxy.potentialr(r1_) + L1e1**2/r1_**2/2) * 1.00001
									print e1_, e1
									assert e1_ > e1
								print r1, r2, r2_
								def dEL(e, L, r):
									return -4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2))
								def dE(e, s1, s2):
									Lmax_ = lmax(e, r_at_Lmax1, r_at_Lmax2)
									L1 = Lmax_ * s1
									L2 = Lmax_ * s2
									r2__ = r2_
									print L1, L2, (2*(e-galaxy.potentialr(r1)) - L1**2/r1**2)
									assert (2*(e-galaxy.potentialr(r1)) - L1**2/r1**2) > 0 # the smallest angular momentum should at least reach the inner radius 
									if (2*(e-galaxy.potentialr(r2__)) - L1**2/r2__**2) < 0: # smallest angular momentum cannot cover the whole bin
										def f(r):
											r = abs(r)
											return 2*(e-galaxy.potentialr(r)) - L1**2/r**2
										r2__ = abs(scipy.optimize.fsolve(f, r1, )) * 0.9999
										print f(r2__)
										assert (2*(e-galaxy.potentialr(r2__)) - L1**2/r2__**2) > 0.
										print "."
									print dEL(e, L1, r1), dEL(e, L1, r2__), dEL(e, L2, r1), dEL(e, L2, r2__)
									return 0
									#return scipy.integrate.quad(lambda r: dEL(e, L2, r) - dEL(e, L1, r), r1_, r2_)[0]
									#return 0
									#Lmaxi = lmax(e, r_at_Lmax1, r_at_Lmax2)
								Lavg = Lmax * (k+0.5)/nL
								self.rhogrid[i][j] += scipy.integrate.quad(lambda e: dE(e, s1, s2), e1_, e2_)[0] * galaxy.fL(Lavg) 
							#self.rhogrid[i][j] += scipy.integrate.dblquad(lambda r,e: dEL(e, lmax(e,r)*s1, r)-dEL(e, lmax(e,r)*s1, r) ,\
						#e1, e2, lambda e: r1, lambda e: r2)[0] * galaxy.fL(Lavg)
					#/(e2-e1)#/(r2-r1)
						#self.m2_grid[i][j] += scipy.integrate.quad(lambda L: 4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
					#0, sqrt(2*(e-galaxy.potentialr(r)))*r)[0]#/(r2-r1)
			else:
				for j in range(len(self.rs)):
					r1 = self.r_borders[j]
					r2 = self.r_borders[j+1]
					r = self.rs[j]
					#print "r=", r, j 
					#b1 = 2*(e-galaxy.potentialr(r1))-L**2/(r1**2) 
					#b2 = 2*(e-galaxy.potentialr(r2))-L**2/(r2**2)
					#self.rhogrid[i][j] += scipy.integrate.dblquad(lambda L,x: 4 * pi / sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L),\
					# r1, r2, lambda r: 0, lambda r: sqrt(2*(e-galaxy.potentialr(r)))*r)[0]/(r2-r1)
					e1 = self.e_borders[i]
					e2 = self.e_borders[i+1]
					b = 2*(e-galaxy.potentialr(r))*r**2
					b1 = 2*(e1-galaxy.potentialr(r))*r**2
					b2 = 2*(e2-galaxy.potentialr(r))*r**2
					#if b > 0:
						#self.rhogrid[i][j] += scipy.integrate.quad(lambda L: 4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
					#0, sqrt(2*(e-galaxy.potentialr(r)))*r)[0]#/(r2-r1)
						#self.m2_grid[i][j] += scipy.integrate.quad(lambda L: 4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
					#0, sqrt(2*(e-galaxy.potentialr(r)))*r)[0]#/(r2-r1)
					#	pass
					def maxl(e):
						b = 2*(e-galaxy.potentialr(r))
						#return sqrt(b)*r
						if b > 0:
							return sqrt(b)*r
						else:
							return 0
					if b2 < 0: # if the 'orbit' with the largest energy cannot make it..
						assert b1 < 0 # ... the orbit with even less energy cannot make it at all
					else:
						if (b1 < 0): # but if only the lowest energy orbit cannot make it.. 
							e1 = galaxy.potentialr(r) * 0.9999999 # find the lowest energy that is usefull,
						# but take it slightly larger to avoid numerical problems 
						b1 = 2*(e1-galaxy.potentialr(r))*r**2
						b2 = 2*(e2-galaxy.potentialr(r))*r**2
						#print b1, b2
						if 0:
							r_at_Lmax1 = self.r_at_Lmax_borders[i]
							r_at_Lmax2 = self.r_at_Lmax_borders[i+1]
							nL = 10
							def lmax(e, r1, r2):
								def f(rc, E=e):
									rc = abs(rc)
									return -r**2*(2*(E-galaxy.potentialr(rc)))
								#rc = abs(scipy.optimize.fminbound(f, r1, r2)[0])
								rc = abs(scipy.optimize.fmin(f, r1, disp=False)[0])
								#print "(",r,")"
								return (-f(rc))**0.5
							
							def dEL(e, L, r):
								return -4*pi* 4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * r**2
							def dE(E):
								#print E
								I = 0
								e = E
								Lmax = lmax(e, r_at_Lmax1, r_at_Lmax2)
								assert Lmax > maxl(e)
								for k in range(nL):
									s1 = (k+0.)/nL
									s2 = (k+1.)/nL
									L1 = Lmax * s1
									L2 = Lmax * s2
									if L1 == 0:
										L1 = 1e-6
									Lavg = Lmax * (k+0.5)/nL
									#L1 = 0
									#L2 = maxl(e)
									#print dEL(e, L1, r), L1, 2*(e-galaxy.potentialr(r)) #-L1**2/(r**2)
									#return (0 - dEL(e, L1, r))
									if (2*(e-galaxy.potentialr(r)) - L1**2/r**2) < 0:
										pass # even no interation with low L  
									else:
										f = galaxy.fL(Lavg)
										#f = (galaxy.fL(L2)+galaxy.fL(L1))/2
										if (2*(e-galaxy.potentialr(r)) - L2**2/r**2) < 0:
											L2 = 2*(e-galaxy.potentialr(r)) * 1.0001
											#f = (galaxy.fL(L2)+galaxy.fL(L1))/2
											f = galaxy.fL(L2*0.9999)
											I += (0 - dEL(e, L1, r)) * f
											#print "1)", dEL(e, L1, r)
										else:
											I += (dEL(e, L2, r) - dEL(e, L1, r)) * f
											#print "2)", dEL(e, L1, r), dEL(e, L2, r)
									#print I
								return I
							self.rhogrid[i][j] += scipy.integrate.quad(lambda e: dE(e), e1, e2)[0]
						else:
							self.rhogrid[i][j] += scipy.integrate.dblquad(lambda L, e: 4*pi*r**2 *4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
							e1, e2, lambda e: 0, lambda e: maxl(e))[0]#/(e2-e1)#/(r2-r1)
							#self.rhogrid[i][j] += scipy.integrate.dblquad(lambda r, L, e: 4*pi*r**2 *4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
						#e1, e2, lambda e: 0, lambda e: maxl(e))[0]#/(e2-e1)#/(r2-r1)
							#def dL(e, L):
							#	return -4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L)
							#print dL(e1, maxl(e1)), 2*(e1-galaxy.potentialr(r))*r, "-", dL(e2, maxl(e2)), 2*(e2-galaxy.potentialr(r))*r, "*", dL(e1, 0), 2*(e1-galaxy.potentialr(r))*r, "-", dL(e2, 0), 2*(e2-galaxy.potentialr(r))*r
							#self.rhogrid[i][j] += scipy.integrate.quad(lambda e: 0-dL(e, 0), e1, e2)[0]
							# self.rhogrid[i][j] =  
						#	self.rhogrida[i][j] += scipy.integrate.dblquad(lambda L, e: e * 4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
						#e1, e2, lambda e: 0, lambda e: maxl(e))[0]#/(e2-e1)#/(r2-r1)
						#	self.rhogridb[i][j] += scipy.integrate.dblquad(lambda L, e: 4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
						#e1, e2, lambda e: 0, lambda e: maxl(e))[0]#/(e2-e1)#/(r2-r1)
						#self.m2_grid[i][j] += scipy.integrate.dblquad(lambda L, e: 4*pi*r**2 *4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  galaxy.fL(L),\
						#e1, e2, lambda e: 0, lambda e: maxl(e))[0]#/(e2-e1)#/(r2-r1)
						if 0:
							#print b1, br
							if (b1 > 0) and (b2 > 0): 
								self.rhogrid[i][j] += scipy.integrate.quad(lambda x: 4 * pi / sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1) * da
								self.m2_grid[i][j] += scipy.integrate.quad(lambda x: 4 * pi * sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1)* da
								pass
							#else:
			if 0:
				da = 0.01/100
				for a in arange(da/2, 1-da/2, da):
				#for a in arange(0., 1.0, da):
					L = Lmax * a
					#print L
					#rhogrid[i][i1:i2] = 4 * pi * sqrt(2*(psis[i1:i2]-e)) * L **(-2*beta)
					#assert all(psis[i1:i2]-e > 0)
					# da = const, Lda=dL
					#rhogrid[i][i1:i2] += 4 * pi / sqrt(2*(psis[i1:i2]-e)-L**2/(r**2)) * L / r**2 *  L * L **(-2*beta)
					#b = 2*(e-galaxy.potentialr(rs))-L**2/(rs**2)
					#mask = b >= 0
					for j in range(len(self.rs)):
						r1 = self.r_borders[j]
						r2 = self.r_borders[j+1]
						b1 = 2*(e-galaxy.potentialr(r1))-L**2/(r1**2) 
						b2 = 2*(e-galaxy.potentialr(r2))-L**2/(r2**2)
						if 0:
							#print b1, br
							if (b1 > 0) and (b2 > 0): 
								self.rhogrid[i][j] += scipy.integrate.quad(lambda x: 4 * pi / sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1) * da
								self.m2_grid[i][j] += scipy.integrate.quad(lambda x: 4 * pi * sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1)* da
								pass
							else:
								
								if 0 and ((b1 > 0) or (b2 > 0)):
									if b1 < 0:
										def f(r):
											return 2*(e-galaxy.potentialr(r))-L**2/(r**2)
										oldr1 = r1
										r1 = scipy.optimize.brentq(f, r1, r2)
										b1 = 2*(e-galaxy.potentialr(r1))-L**2/(r1**2)
										if r1 < oldr1:
											dr = -1e8
										else:
											dr = 1e8
										while b1 < 0:
											r1 += 1e-10
											b1 = 2*(e-galaxy.potentialr(r1))-L**2/(r1**2)
									if b2 < 0:
										def f(r):
											return 2*(e-galaxy.potentialr(r))-L**2/(r**2)
										oldr2 = r2
										r2 = scipy.optimize.brentq(f, r1, r2)
										b2 = 2*(e-galaxy.potentialr(r2))-L**2/(r2**2)
										if r2 < oldr2:
											dr = -1e8
										else:
											dr = 1e8
										while b2 < 0:
											r2 += dr
											b2 = 2*(e-galaxy.potentialr(r2))-L**2/(r2**2)
									if (not isnan(b1)) and (not isnan(b2)):
										#print b1, b2
										assert b1 >= 0, "b1=%e" % b1 
										assert b2 >= 0, "b2=%e" % b2
										self.rhogrid[i][j] += scipy.integrate.quad(lambda x: 4 * pi / sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1) * da
										self.m2_grid[i][j] += scipy.integrate.quad(lambda x: 4 * pi * sqrt(2*(e-galaxy.potentialr(x))-L**2/(x**2)) * L / x**2 *  L * galaxy.fL(L), r1, r2)[0]/(r2-r1)* da
									
							
						r = self.rs[j]
						b = 2*(e-galaxy.potentialr(r))-L**2/(r**2) 
						if b > 0:
							self.rhogrid[i][j] += 4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  L * galaxy.fL(L) * da
							#self.rhogrid[i][j] += 4 * pi / sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  L * da# / (r2-r1)
							self.m2_grid[i][j] += 4 * pi * sqrt(2*(e-galaxy.potentialr(r))-L**2/(r**2)) * L / r**2 *  L * galaxy.fL(L) * da
							#rhogrid[i][j] += 4 * pi / sqrt(b1) * L / r1**2 *  L * galaxy.fL(L)
							pass
					#rhogrid[i][mask] += 4 * pi / sqrt(b[mask]) * L / rs[mask]**2 *  L * galaxy.fL(L)
					#rhogrid[i][mask] += -4/3*pi(2*e-2*p[mask]-l**2/r**2)
				# if beta = 0
				#rhogrid[i][i1:i2] = 4 * pi *  sqrt(2*(psis[i1:i2]-e))
				#@if i > 0:
				#@	if rhogrid[i][i-1] <= 0:
				#		import pdb; pdb.set_trace()
		indices = range(len(self.es))
		do_rho(indices)
		#do_rho(12)
		
		rho_reproduced = self.rho_targets * 0
		x = rho_reproduced * 0
		
		if 0:
			for i in range(1, len(self.es))[::-1]:
				#print rhogrid[i][i-1]
				#if i > 0: 
				#	assert rhogrid[i][i-1] > 0
				dE = self.dEs[i]
				
				#if rhogrid[i][i-1] > 0:
				if 1:
					print i, self.rho_targets.shape, rho_reproduced.shape, self.rhogrid.shape
					extra = 1. #/((rho_star[i-1] - rho_star[i])/(es[i-1] - es[i]))
					#extra = 1/dE #1./(es[i-1] - es[i])
					x[i] = (self.rho_targets[i] - rho_reproduced[i])/self.rhogrid[i][i] * extra
					#f_star[i] = (rho_star[i-1]/rhogrid[i][i-1] - rho_reproduced[i-1]/rhogrid[i][i-1]) * extra
					#f_star[i] = min((rho_star[0:i]/rhogrid[i][0:i] - rho_reproduced[0:i]/rhogrid[i][0:i]) * extra)
					#print i, es[i-1]
					#assert f_star[i] > 0
					if x[i] < 0:
						x[i] = 0
					#print f_star[i], rhogrid[i][i-1]
					#print "extra", extra, f_star[i]
					#f_star[i] *= extra
					#print "extra", extra, f_star[i]
					#f_star[i] = rho_star[i-1]/rhogrid[i][i-1] - rho_reproduced[i-1]/rhogrid[i][i-1]
					rho_reproduced += x[i] * self.rhogrid[i] / extra 
		
		x = self.rho_targets * 0 + 1
		#ret = scipy.optimize.nnls(transpose(self.rhogrid), self.rho_targets)
		#x = ret[0]
		#print ret
		if 0:
			box()
			rho_reproduced = dot(transpose(self.rhogrid), x)
			graph(self.logrs, self.rho_targets, color="red")
			graph(self.logrs, rho_reproduced, linestyle="dash")
		#draw()
		#dE = e_borders[1:] - e_borders[:-1]
		self.fE = x#/self.dEs #/dE
		#self.es = es
		