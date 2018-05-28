# -*- coding: utf-8 -*-
from numpy import *
import mab.gd.gdfast
import mab.utils.numpy
import mab.random
import os
import bz2
from mab.constants import *
from mab.astrounits  import *
import numpy
import sys
import mab.gd.logging as logging

kpc_to_km = (1*KPC).asNumber(KM)
logger = logging.getLogger("gd.schw.ic")

class OrbitLibraryCommand(object):
	def __init__(self, schwmodel, modelpath, schwsetname, schwmodelname, load_compressed):
		self.schwmodel = schwmodel
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.load_compressed = load_compressed
		self.dirname = os.path.join(modelpath, "schw", schwsetname, schwmodelname)
		
	def run(self, args, opts, scope):
		filename_check = os.path.join(self.dirname, "finished.orblib")
		if os.path.exists(filename_check):
			logger.info("skipping creating orbit library (already finished)")
			return 
		nE = len(self.schwmodel.dfgrid.Es)
		nL = len(self.schwmodel.dfgrid.ls)
		i1s = [i1 for i1 in range(nE) for i2 in range(nL)]
		i2s = [i2 for i1 in range(nE) for i2 in range(nL)]
		
		
		self.schwmodel.storage_2d.init()
		self.schwmodel.storage_3d.init()
		
		
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def integrate_orbit(i1, i2):
			self.schwmodel.integrate_orbit(i1, i2)
			
		self.schwmodel.load_initial_conditions(self.schwsetname, self.schwmodelname, self.load_compressed)
	
		if scope["cores"] == 1:
			if opts.progressbar:
				pb = mab.utils.progressbar.ProgressBar(0, nE*nL-1)
			for i1 in range(nE):
				for i2 in range(nL):
					self.schwmodel.integrate_orbit(i1, i2)
					if opts.progressbar:
						pb.update(i2+i1*nL)
		else:
			integrate_orbit(i1s, i2s)
		self.schwmodel.storage_2d.save()
		self.schwmodel.storage_3d.save()
		os.system("touch %s" % filename_check)
			
class InitialConditionsCommand(object):
	def __init__(self, schwmodel, schwsetname, schwmodelname, store_compressed):
		self.schwmodel = schwmodel
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.store_compressed = store_compressed
		
	def run(self, args, opts, scope):
		#dirname = self.schwmodel.dirname
		#logger.info("changing working dir to: %s" % dirname)
		
		dirname = os.path.join(self.schwmodel.modelpath, "schw", self.schwsetname, self.schwmodelname, "intermediate")
		if not os.path.exists(dirname):
			logger.info("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.info("directory %s already exists" % dirname)
			
		if self.store_compressed:
			filename = os.path.join(self.schwmodel.modelpath, "schw", self.schwsetname, self.schwmodelname, "intermediate", "initial_conditions.npy.bz2")
		else:
			filename = os.path.join(self.schwmodel.modelpath, "schw", self.schwsetname, self.schwmodelname, "intermediate", "initial_conditions.npy")
		if os.path.exists(filename):
			logger.info("skipping creating intial conditions (%s already exists)" % filename)
			return
		
		schwmodel = self.schwmodel
		nE = len(schwmodel.dfgrid.Es)
		nL = len(schwmodel.dfgrid.ls)
		dither = schwmodel.dfgrid.dither
		
		#progressbar = mab.utils.progressbar.ProgressBar(0, nE-1)
		#progressbar.update(0)
		
		
		#names = "x, y, z, vx, vy, vz, E, Lx, Ly, Lz, Lmax, Tr".split(",")
		#names = [k.strip() for k in names]
		#types = ["f4"] * len(names)
		#dtype_initial_conditions = zip(names, types)
		#initial_conditions = mab.utils.numpy.mmapzeros((nE, nL, dither, dither), dtype=dtype_initial_conditions)
		#initial_conditions = initial_conditions.view(recarray)
		#print schwmodel.profile_model.Tr(E=-3903.0259983306432, L=0.012664602746287194)
		#sys.exit(0)
		
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def calculate_initial_conditions_parallel(i1, i2):
			schwmodel.calculate_initial_condition(i1, i2)
			
		def calculate_initial_conditions_serial():
			if opts.progressbar:
				progressbar = mab.utils.progressbar.ProgressBar(0, nE-1)
			for i1 in range(nE): #[13:]:
				for i2 in range(nL):
					#print i1, i2
					schwmodel.calculate_initial_condition(i1, i2)
				if opts.progressbar:
					progressbar.update(i1)
		
		i1s = [i1 for i1 in range(nE) for i2 in range(nL)]
		i2s = [i2 for i1 in range(nE) for i2 in range(nL)]
		logger.debug("initial conditions radii: %r" % schwmodel.dfgrid.rs)
		#print schwmodel.dfgrid.subgrid.rs
		#sys.exit(0)
		if scope["cores"] > 1:
			calculate_initial_conditions_parallel(i1s, i2s)
		else:
			calculate_initial_conditions_serial()
		#print initial_conditions.x
		
		
		schwmodel.store_initial_conditions(self.schwsetname, self.schwmodelname, compressed=self.store_compressed)
				
class SchwModelAnalytic(object):
	def __init__(self, profile_model, storage3d_output, dfgrid, grid_logr):
		self.dfgrid = dfgrid
		self.profile_model = profile_model
		self.storage3d_output = storage3d_output 
		self.rs = self.storage3d_output.finv(storage3d_output.x)
		self.r_borders = self.storage3d_output.finv(storage3d_output.xborders)
		self.profile_model_fast = profile_model.fast()
		self.df = mab.gd.gdfast.DF(self.profile_model_fast, self.rs, self.r_borders)
		self.N_r = len(self.rs)
		
	def run(self, args, opts, scope):
		self.do(scope)
		
	def do(self, scope):
		filename = self.storage3d_output.filename +".bz2"
		if os.path.exists(filename):
			logger.debug("orbit library exists: %s" % filename)
			return
		self.storage3d_output.init()
		logger.info("calculating orbit libary (analytic)")
		rmin = 1e-10
		Nd = self.dfgrid.dither
		density = zeros(self.N_r)
		varvr = zeros(self.N_r)
		varvt = zeros(self.N_r)
		v22 = zeros(self.N_r)
		v40 = zeros(self.N_r)
		v04 = zeros(self.N_r)
		
		for Eindex in range(self.dfgrid.n_I1):
			logger.debug("Eindex=%d" % Eindex)
			for i in range(self.dfgrid.dither):
				E = self.dfgrid.subgrid.Es[Eindex*self.dfgrid.dither+i]
				Lmax, rcirc = self.profile_model_fast.Lmax_and_rcirc_at_E(E)
				#print "Lmax, rcirc", Lmax, rcirc
				logger.debug("Lmax=%f rcirc=%f" % (Lmax, rcirc))
				rmax = self.profile_model_fast.rmax_at_E(E, rcirc)
				logger.debug("rmax=%f" % (rmax))
				#l = self.dfgrid.subgrid.ls[li*self.dfgrid.dither+j]
				for Lindex in range(self.dfgrid.n_I2):
					logger.debug("Lindex=%d" % Lindex)
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
					density[:] = 0
					varvr[:] = 0
					varvt[:] = 0
					v22[:] = 0
					v40[:] = 0
					v04[:] = 0
					#print "finding moments"
					try:
						self.df.momentsL2(density, varvr, varvt, v22, v40, v04, E, L1, L2, rp, ra, False)
					except:
						if 0:
							try:
								self.df.momentsL2(density, varvr, varvt, v22, v40, v04, E, L1, L2, rp, ra, True)
							except:
								pass
						print "error for", scope["model_id"], scope["model_id_path"]
						print "Eindex", Eindex
						print "Lindex", Lindex
						print "l1,l2", l1, l2
						#raise
					#print "found moments"
					
					grid = self.storage3d_output.moments3dgrid[Lindex, Eindex]
					#N = sum(density)/Nd
					grid[self.storage3d_output.v00index] += density/Nd
					grid[self.storage3d_output.v20index] += varvr/Nd
					grid[self.storage3d_output.v02index] += varvt/Nd
					grid[self.storage3d_output.v22index] += v22/Nd
					grid[self.storage3d_output.v40index] += v40/Nd
					grid[self.storage3d_output.v04index] += v04/Nd
		self.storage3d_output.save()				

class SchwModel_SphericalNonRotating(object):
	def __init__(self, light_model, profile_model, modelpath, orbitintegrator, orbital_periods, random_rotations, dfgrid, storage_2d, storage_3d):
		self.light_model = light_model
		self.profile_model = profile_model
		self.modelpath = modelpath
		self.orbitintegrator = orbitintegrator
		self.orbital_periods = orbital_periods
		self.random_rotations = random_rotations
		self.dfgrid = dfgrid
		self.n_E = self.dfgrid.n_I1
		self.n_L = self.dfgrid.n_I2
		self.n_orbits = self.n_E * self.n_L
		self.dither = self.dfgrid.dither
		#self.aperture = aperture
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		#self.storage_2d.init(self.n_orbits)
		#self.storage_3d.init(self.n_orbits)
		#self.binner = mab.gd.gdfast.Binner(self.dfgrid.logrmin, self.dfgrid.logrmin)
		
		names = "x, vy, E, Lz, Lmax, Tr".split(",")
		names = [k.strip() for k in names]
		types = ["f4"] * len(names)
		dtype_initial_conditions = zip(names, types)
		self.initial_conditions = mab.utils.numpy.mmapzeros((self.n_E, self.n_L, self.dither, self.dither), dtype=dtype_initial_conditions)
		self.initial_conditions = self.initial_conditions.view(recarray)


	def store_initial_conditions(self, schwsetname, schwmodelname, compressed=True):
		if compressed:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy.bz2")
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy")
			f = file(filename, "wb")
		save(f, self.initial_conditions)
		
	def load_initial_conditions(self, schwsetname, schwmodelname, compressed=True):
		if compressed:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy.bz2")
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy")
			f = file(filename, "rb")
		self.initial_conditions = load(f)
		self.initial_conditions = self.initial_conditions.view(recarray)
	
	
	"""def integrate_orbit2(self, i1, i2):
		orbitnr = self.index_to_orbitnr(i1, i2)
		numpy.random.seed(orbitnr)
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		dt = (self.initial_conditions.Tr[i1, i2] * self.orbital_periods/self.orbitintegrator.orbital_points).reshape(dither1*dither2).astype(float64)
		N = len(dt)
		q0 = zeros((N, 2))
		#print len(self.initial_conditions.x), N, dt.shape
		#print self.initial_conditions.x.shape, q0[:,0].shape
		q0[:,0] = self.initial_conditions.x[i1, i2].reshape(dither1*dither2)
		p0 = zeros((N, 2))
		p0[:,1] = self.initial_conditions.vy[i1, i2].reshape(dither1*dither2)
		#, Tr, q0, p0
		q, p = self.orbitintegrator.integrate(dt, q0, p0) #, xout, yout, vxout, vyout)
		
		# 'allocate' memory once
		x, y, z = zeros((3, len(q[0,0])))
		vx, vy, vz = zeros((3, len(q[0,0])))
		#vr, vphi, vtheta = zeros((3, len(dt)))
		
		for d1 in range(dither1):
			for d2 in range(dither2):
				j = d1*dither2 + d2 
				xout, yout = q[:,j,:]
				vxout, vyout = p[:,j,:]
				for r in range(self.random_rotations):
					# rotate p and q
					R = array(mab.random.randomSO3())
					mab.gd.gdfast.matmul_2d_3d(R, xout, yout, x, y, z)
					mab.gd.gdfast.matmul_2d_3d(R, vxout, vyout, vx, vy, vz)
					kpc_to_arcsec = self.light_model.kpc_to_arcsec(1.)
					x *= kpc_to_arcsec
					y *= kpc_to_arcsec
					z *= kpc_to_arcsec
					#mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, vr, vphi, vtheta)
					#r = sqrt(x*x+y*y+z*z)
					#logr = log10(r)
					
					self.storage_2d.store(x, y, vz, orbitnr, u, v)
					self.storage_2d.store(x, y, -vz, orbitnr)
					self.storage_2d.store(y, z, vx, orbitnr)
					self.storage_2d.store(y, z, -vx, orbitnr)
					self.storage_2d.store(z, x, vy, orbitnr)
					self.storage_2d.store(z, x, -vy, orbitnr)
				self.storage_3d.store(x, y, z, vx, vy, vz, orbitnr, i1, i2, d1, d2, u, v)
					#binner3d.bin3d(orbitnr, logr, vr, vphi, vtheta)
					#binner2d.bin(orbitnr, x, y, z, vx, vy, vz)
		self.storage_2d.normalize(orbitnr)
		self.storage_3d.normalize(orbitnr)"""
			
	def integrate_orbit(self, i1, i2):
		orbitnr = self.index_to_orbitnr(i1, i2)
		numpy.random.seed(orbitnr+20*1000)
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		dt = (self.initial_conditions.Tr[i1, i2] * self.orbital_periods/self.orbitintegrator.orbital_points).reshape(dither1*dither2).astype(float64)
		N = len(dt)
		q0 = zeros((N, 2))
		#print len(self.initial_conditions.x), N, dt.shape
		#print self.initial_conditions.x.shape, q0[:,0].shape
		q0[:,0] = self.initial_conditions.x[i1, i2].reshape(dither1*dither2)
		p0 = zeros((N, 2))
		p0[:,1] = self.initial_conditions.vy[i1, i2].reshape(dither1*dither2)
		#, Tr, q0, p0
		q, p = self.orbitintegrator.integrate(dt, q0, p0) #, xout, yout, vxout, vyout)
		
		# 'allocate' memory once
		x, y, z = zeros((3, len(q[0,0])))
		vx, vy, vz = zeros((3, len(q[0,0])))
		#vr, vphi, vtheta = zeros((3, len(dt)))
		
		for d1 in range(dither1):
			for d2 in range(dither2):
				j = d1*dither2 + d2
				#I1 = self.dfgrid.subgrid.I1s[i1*dither1 + d1]
				#I2 = self.dfgrid.subgrid.I2s[i2*dither2 + d2]
				xout, yout = q[:,j,:]
				vxout, vyout = p[:,j,:]
				u = float(d1 + 0.5) / dither1
				v = float(d2 + 0.5) / dither2
				M = len(xout)
				
				self.storage_3d.store(xout, yout, yout*0, vxout, vyout, vxout*0, orbitnr, u, v)
				if self.random_rotations == 0:
					#print "test"
					kpc_to_arcsec = self.light_model.kpc_to_arcsec(1.)
					
					r = numpy.sqrt(xout**2+yout**2)
					vr = (xout*vxout + yout*vyout)/r
					vt = (xout*vyout - yout*vxout)/r
					
					costheta = numpy.random.random(M) * 2 - 1
					phi = numpy.random.random(M) * 2 * pi
					eta = numpy.random.random(M) * 2 * pi
					theta = numpy.arccos(costheta)
					#sintheta = numpy.sqrt(1-costheta**2)
					sintheta = numpy.sin(theta)
					#print r.shape, sintheta.shape, phi.shape, len(dt) 
					x = r * sintheta * numpy.cos(phi)
					y = r * sintheta * numpy.sin(phi)
					vlos = vr * costheta - vt * numpy.cos(eta) * sintheta
					
					x *= kpc_to_arcsec
					y *= kpc_to_arcsec
					#z *= kpc_to_arcsec
					#mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, vr, vphi, vtheta)
					#r = sqrt(x*x+y*y+z*z)
					#logr = log10(r)
					self.storage_2d.store_test(r*kpc_to_arcsec, vr, vt, orbitnr, u, v)
					
				
				for r in range(self.random_rotations):
					#R = array(mab.random.randomSO3())
					#mab.gd.gdfast.matmul_2d_3d(R, xout, yout, x, y, z)
					#mab.gd.gdfast.matmul_2d_3d(R, vxout, vyout, vx, vy, vz)
					kpc_to_arcsec = self.light_model.kpc_to_arcsec(1.)
					
					r = numpy.sqrt(xout**2+yout**2)
					vr = (xout*vxout + yout*vyout)/r
					vt = (xout*vyout - yout*vxout)/r
					
					costheta = numpy.random.random(M) * 2 - 1
					phi = numpy.random.random(M) * 2 * pi
					eta = numpy.random.random(M) * 2 * pi
					theta = numpy.arccos(costheta)
					#sintheta = numpy.sqrt(1-costheta**2)
					sintheta = numpy.sin(theta)
					#print r.shape, sintheta.shape, phi.shape, len(dt) 
					x = r * sintheta * numpy.cos(phi)
					y = r * sintheta * numpy.sin(phi)
					vlos = vr * costheta - vt * numpy.cos(eta) * sintheta
					
					x *= kpc_to_arcsec
					y *= kpc_to_arcsec
					#z *= kpc_to_arcsec
					#mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, vr, vphi, vtheta)
					#r = sqrt(x*x+y*y+z*z)
					#logr = log10(r)
					self.storage_2d.store(x, y, vlos, orbitnr, u, v)
					#self.storage_2d.store(x, y, -vz, orbitnr, u, v)
					#self.storage_2d.store(y, z, vx, orbitnr, u, v)
					#self.storage_2d.store(y, z, -vx, orbitnr, u, v)
					#self.storage_2d.store(z, x, vy, orbitnr, u, v)
					#self.storage_2d.store(z, x, -vy, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, vx, vy, vz, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, -vx, -vy, -vz, orbitnr, u, v)
					#self.storage_2d.store(x, y, vz, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, -vx, -vy, vz, orbitnr, u, v)
					#binner3d.bin3d(orbitnr, logr, vr, vphi, vtheta)
					#binner2d.bin(orbitnr, x, y, z, vx, vy, vz)
				#self.storage_3d.store(xout, yout, xout*0, vxout, vyout, vxout*0, orbitnr, u, v)	
		self.storage_2d.normalize(orbitnr)
		self.storage_3d.normalize(orbitnr)
		
	def index_to_orbitnr(self, i1, i2):
		#return i1*self.n_L + i2
		return self.dfgrid.index_to_orbitnr(i1, i2, 0)
		#3return i1 + self.n_E * i2
		
	def calculate_initial_condition(self, i1, i2):
		kpc_to_km = (1*KPC).asNumber(KM)
		s_to_gyr = (S/GYR).asNumber()
		
		profile_model_fast = self.profile_model.fast()
		rmin=1e-6
		
		for d1 in range(self.dfgrid.dither):
			Ei = i1*self.dfgrid.dither+d1
			E = self.dfgrid.subgrid.Es[Ei]
			#Lmax = self.profile_model.Lmax_at_E(E)
			Lmax, rcirc = profile_model_fast.Lmax_and_rcirc_at_E(E)
			logger.debug("Lmax=%f rcirc=%f" % (Lmax, rcirc))
			rmax = profile_model_fast.rmax_at_E(E, rcirc)
			logger.debug("rmax=%f" % rmax)
			rs = self.dfgrid.subgrid.rs[Ei]
			#print rs
			
			for d2 in range(self.dfgrid.dither):
				Li = i2*self.dfgrid.dither+d2
				L = self.dfgrid.subgrid.ls[Li] * Lmax
				
				ra, rp = profile_model_fast.get_apo_peri(E, L, rmin, rcirc, rmax)
				#ra, rp = self.profile_model.get_apo_peri(E, L, rtry=rs)
				logger.debug("ra=%f rp=%f" % (ra, rp))
				r = ra # start at apocenter
				rc = rcirc
				#rc = self.profile_model.findR_at_EL(E, Lmax, rtry=rs)
				logger.debug("rc=%f" % rc)
				assert rc <= r
				vt = L/r
					
				#rcirc = r
				vcirc = sqrt(rc * self.profile_model.dphidr(rc))
				
				tcirc = 2*pi*(r*kpc_to_km) / vcirc
				Tr = tcirc
				Tr = self.profile_model.Tr(E, L, ra=ra, rp=rp)*kpc_to_km
				#print Tr, tcirc
				self.initial_conditions.x[i1,i2,d1,d2] = r
				self.initial_conditions.vy[i1,i2,d1,d2] = vt
				self.initial_conditions.Lz[i1,i2,d1,d2] = L
				self.initial_conditions.E[i1,i2,d1,d2] = E
				self.initial_conditions.Tr[i1,i2,d1,d2] = Tr
				self.initial_conditions.Lmax[i1,i2,d1,d2] = Lmax
				
		
		

class SchwModel_Axi(object):
	def __init__(self, light_model, profile_model, modelpath, orbitintegrator, orbital_periods, random_rotations, dfgrid, storage_2d, storage_3d):
		self.light_model = light_model
		self.profile_model = profile_model
		self.modelpath = modelpath
		self.orbitintegrator = orbitintegrator
		self.orbital_periods = orbital_periods
		self.random_rotations = random_rotations
		self.dfgrid = dfgrid
		self.n_E = self.dfgrid.n_I1
		self.n_L = self.dfgrid.n_I2
		self.n_orbits = self.n_E * self.n_L
		self.dither = self.dfgrid.dither
		#self.aperture = aperture
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		#self.storage_2d.init(self.n_orbits)
		#self.storage_3d.init(self.n_orbits)
		#self.binner = mab.gd.gdfast.Binner(self.dfgrid.logrmin, self.dfgrid.logrmin)
		
		names = "x, vy, E, Lz, Lmax, Tr".split(",")
		names = [k.strip() for k in names]
		types = ["f4"] * len(names)
		dtype_initial_conditions = zip(names, types)
		self.initial_conditions = mab.utils.numpy.mmapzeros((self.n_E, self.n_L, self.dither, self.dither), dtype=dtype_initial_conditions)
		self.initial_conditions = self.initial_conditions.view(recarray)


	def store_initial_conditions(self, schwsetname, schwmodelname, compressed=True):
		if compressed:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy.bz2")
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy")
			f = file(filename, "wb")
		save(f, self.initial_conditions)
		
	def load_initial_conditions(self, schwsetname, schwmodelname, compressed=True):
		if compressed:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy.bz2")
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = os.path.join(self.modelpath, "schw", schwsetname, schwmodelname, "intermediate", "initial_conditions.npy")
			f = file(filename, "rb")
		self.initial_conditions = load(f)
		self.initial_conditions = self.initial_conditions.view(recarray)
	
	
	"""def integrate_orbit2(self, i1, i2):
		orbitnr = self.index_to_orbitnr(i1, i2)
		numpy.random.seed(orbitnr)
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		dt = (self.initial_conditions.Tr[i1, i2] * self.orbital_periods/self.orbitintegrator.orbital_points).reshape(dither1*dither2).astype(float64)
		N = len(dt)
		q0 = zeros((N, 2))
		#print len(self.initial_conditions.x), N, dt.shape
		#print self.initial_conditions.x.shape, q0[:,0].shape
		q0[:,0] = self.initial_conditions.x[i1, i2].reshape(dither1*dither2)
		p0 = zeros((N, 2))
		p0[:,1] = self.initial_conditions.vy[i1, i2].reshape(dither1*dither2)
		#, Tr, q0, p0
		q, p = self.orbitintegrator.integrate(dt, q0, p0) #, xout, yout, vxout, vyout)
		
		# 'allocate' memory once
		x, y, z = zeros((3, len(q[0,0])))
		vx, vy, vz = zeros((3, len(q[0,0])))
		#vr, vphi, vtheta = zeros((3, len(dt)))
		
		for d1 in range(dither1):
			for d2 in range(dither2):
				j = d1*dither2 + d2 
				xout, yout = q[:,j,:]
				vxout, vyout = p[:,j,:]
				for r in range(self.random_rotations):
					# rotate p and q
					R = array(mab.random.randomSO3())
					mab.gd.gdfast.matmul_2d_3d(R, xout, yout, x, y, z)
					mab.gd.gdfast.matmul_2d_3d(R, vxout, vyout, vx, vy, vz)
					kpc_to_arcsec = self.light_model.kpc_to_arcsec(1.)
					x *= kpc_to_arcsec
					y *= kpc_to_arcsec
					z *= kpc_to_arcsec
					#mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, vr, vphi, vtheta)
					#r = sqrt(x*x+y*y+z*z)
					#logr = log10(r)
					
					self.storage_2d.store(x, y, vz, orbitnr, u, v)
					self.storage_2d.store(x, y, -vz, orbitnr)
					self.storage_2d.store(y, z, vx, orbitnr)
					self.storage_2d.store(y, z, -vx, orbitnr)
					self.storage_2d.store(z, x, vy, orbitnr)
					self.storage_2d.store(z, x, -vy, orbitnr)
				self.storage_3d.store(x, y, z, vx, vy, vz, orbitnr, i1, i2, d1, d2, u, v)
					#binner3d.bin3d(orbitnr, logr, vr, vphi, vtheta)
					#binner2d.bin(orbitnr, x, y, z, vx, vy, vz)
		self.storage_2d.normalize(orbitnr)
		self.storage_3d.normalize(orbitnr)"""
			
	def integrate_orbit(self, i1, i2):
		orbitnr = self.index_to_orbitnr(i1, i2)
		numpy.random.seed(orbitnr+20*1000)
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		dt = (self.initial_conditions.Tr[i1, i2] * self.orbital_periods/self.orbitintegrator.orbital_points).reshape(dither1*dither2).astype(float64)
		N = len(dt)
		q0 = zeros((N, 2))
		#print len(self.initial_conditions.x), N, dt.shape
		#print self.initial_conditions.x.shape, q0[:,0].shape
		q0[:,0] = self.initial_conditions.x[i1, i2].reshape(dither1*dither2)
		p0 = zeros((N, 2))
		p0[:,1] = self.initial_conditions.vy[i1, i2].reshape(dither1*dither2)
		#, Tr, q0, p0
		q, p = self.orbitintegrator.integrate(dt, q0, p0) #, xout, yout, vxout, vyout)
		
		# 'allocate' memory once
		x, y, z = zeros((3, len(q[0,0])))
		vx, vy, vz = zeros((3, len(q[0,0])))
		#vr, vphi, vtheta = zeros((3, len(dt)))
		
		for d1 in range(dither1):
			for d2 in range(dither2):
				j = d1*dither2 + d2
				#I1 = self.dfgrid.subgrid.I1s[i1*dither1 + d1]
				#I2 = self.dfgrid.subgrid.I2s[i2*dither2 + d2]
				xout, yout = q[:,j,:]
				vxout, vyout = p[:,j,:]
				u = float(d1 + 0.5) / dither1
				v = float(d2 + 0.5) / dither2
				M = len(xout)
				
				self.storage_3d.store(xout, yout, yout*0, vxout, vyout, vxout*0, orbitnr, u, v)
				
				for r in range(self.random_rotations):
					#R = array(mab.random.randomSO3())
					#mab.gd.gdfast.matmul_2d_3d(R, xout, yout, x, y, z)
					#mab.gd.gdfast.matmul_2d_3d(R, vxout, vyout, vx, vy, vz)
					kpc_to_arcsec = self.light_model.kpc_to_arcsec(1.)
					
					r = numpy.sqrt(xout**2+yout**2)
					vr = (xout*vxout + yout*vyout)/r
					vt = (xout*vyout - yout*vxout)/r
					
					costheta = numpy.random.random(M) * 2 - 1
					phi = numpy.random.random(M) * 2 * pi
					eta = numpy.random.random(M) * 2 * pi
					theta = numpy.arccos(costheta)
					#sintheta = numpy.sqrt(1-costheta**2)
					sintheta = numpy.sin(theta)
					#print r.shape, sintheta.shape, phi.shape, len(dt) 
					x = r * sintheta * numpy.cos(phi)
					y = r * sintheta * numpy.sin(phi)
					vlos = vr * costheta - vt * numpy.cos(eta) * sintheta
					
					x *= kpc_to_arcsec
					y *= kpc_to_arcsec
					#z *= kpc_to_arcsec
					#mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, vr, vphi, vtheta)
					#r = sqrt(x*x+y*y+z*z)
					#logr = log10(r)
					self.storage_2d.store(x, y, vlos, orbitnr, u, v)
					#self.storage_2d.store(x, y, -vz, orbitnr, u, v)
					#self.storage_2d.store(y, z, vx, orbitnr, u, v)
					#self.storage_2d.store(y, z, -vx, orbitnr, u, v)
					#self.storage_2d.store(z, x, vy, orbitnr, u, v)
					#self.storage_2d.store(z, x, -vy, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, vx, vy, vz, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, -vx, -vy, -vz, orbitnr, u, v)
					#self.storage_2d.store(x, y, vz, orbitnr, u, v)
					#self.storage_3d.store(x, y, z, -vx, -vy, vz, orbitnr, u, v)
					#binner3d.bin3d(orbitnr, logr, vr, vphi, vtheta)
					#binner2d.bin(orbitnr, x, y, z, vx, vy, vz)
				#self.storage_3d.store(xout, yout, xout*0, vxout, vyout, vxout*0, orbitnr, u, v)	
		self.storage_2d.normalize(orbitnr)
		self.storage_3d.normalize(orbitnr)
		
	def index_to_orbitnr(self, i1, i2):
		#return i1*self.n_L + i2
		return self.dfgrid.index_to_orbitnr(i1, i2, 0)
		#3return i1 + self.n_E * i2
		
	def calculate_initial_condition(self, i1, i2):
		kpc_to_km = (1*KPC).asNumber(KM)
		s_to_gyr = (S/GYR).asNumber()
		
		for d1 in range(self.dfgrid.dither):
			Ei = i1*self.dfgrid.dither+d1
			E = self.dfgrid.subgrid.Es[Ei]
			Lmax = self.profile_model.Lmax_at_E(E)
			rs = self.dfgrid.subgrid.rs[Ei]
			#print rs
			
			for d2 in range(self.dfgrid.dither):
				Li = i2*self.dfgrid.dither+d2
				L = self.dfgrid.subgrid.ls[Li] * Lmax
				
				ra, rp = self.profile_model.get_apo_peri(E, L, rtry=rs)
				r = ra # start at apocenter
				rc = self.profile_model.findR_at_EL(E, Lmax, rtry=rs)
				assert rc <= r
				vt = L/r
					
				#rcirc = r
				vcirc = sqrt(rc * self.profile_model.dphidr(rc))
				
				tcirc = 2*pi*(r*kpc_to_km) / vcirc
				Tr = tcirc
				Tr = self.profile_model.Tr(E, L)*kpc_to_km
				#print Tr, tcirc
				self.initial_conditions.x[i1,i2,d1,d2] = r
				self.initial_conditions.vy[i1,i2,d1,d2] = vt
				self.initial_conditions.Lz[i1,i2,d1,d2] = L
				self.initial_conditions.E[i1,i2,d1,d2] = E
				self.initial_conditions.Tr[i1,i2,d1,d2] = Tr
				self.initial_conditions.Lmax[i1,i2,d1,d2] = Lmax
				
		