import itertools
from numpy import *
import os
import mab.gd.logging as logging
try:
	import h5py
except:
	pass
from kaplot import *

logger = logging.getLogger("gd.gadget")

class GadgetObservation(object):
	def __init__(self, galaxy, galaxy_host, orbit, working_directory, cores):
		self.galaxy = galaxy
		self.galaxy_host = galaxy_host
		self.orbit = orbit
		self.working_directory = working_directory
		self.cores = cores
		
	def run(self, args, kwargs, scope):
		oldcwd = os.getcwd()
		try:
			def do(i):
				logger.info("changing working directory to: %s" % self.working_directory)
				os.chdir(self.working_directory)
				
				filename = os.path.join("output", "snapshot_%03d.hdf5" % i)
				#if not os.path.exists(filename):
				#	return
				
				print filename
				try:
					f = h5py.File(filename)
					if self.orbit:
						w = self.orbit.ra+25
					else:
						w = 200
						#w = 3.
					#rectangle(-w, -w, w, w, solid=True, color="lightgrey")
					pos = f["/PartType1/Coordinates"]
					vel = f["/PartType1/Velocities"]
					pot_dm = f["/PartType1/Potential"]
					x_dm, y_dm, z_dm = x, y, z = transpose(pos)
					vx_dm, vy_dm, vz_dm = vx, vy, vz = transpose(vel)
						
					pos = f["/PartType2/Coordinates"]
					vel = f["/PartType2/Velocities"]
					pot = f["/PartType2/Potential"]
					x, y, z = transpose(pos)
					vx, vy, vz = transpose(vel)
				except:
					print "error for", filename
					raise
				center_index = argmin(pot)
				x0 = x[center_index]
				y0 = y[center_index]
				z0 = z[center_index]
				print "zero", x0, y0, z0
				#vx0, vy0 = 0, 0
				#x0, y0, vx0, vy0 = find_center(x, y, vx, vy)
				#x0, y0, vx0, vy0 = 0, 0, 0, 0
				xc = x-x0
				yc = y-y0
				zc = z-z0
				xc_dm = x_dm - x0
				yc_dm = y_dm - y0
				zc_dm = z_dm - z0
				
				r = (xc**2+yc**2+zc**2)**0.5
				mask = r < 30
				vx0 = mean(vx[mask])
				vy0 = mean(vy[mask])
				vz0 = mean(vz[mask])
				
				p_center = array([x0, y0, z0])
				n_center = array([x0, y0, z0]) / sqrt(x0**2+y0**2+z0**2)
				x = x[mask]
				y = y[mask]
				z = z[mask]
				vx = vx[mask]
				vy = vy[mask]
				vz = vz[mask]
				n = 10000*7
				p_center = array([x0, y0, z0])
				Rs = []
				for i in range(n):
					p = array([x[i], y[i], z[i]])
					e_los = array([x[i], y[i], z[i]]) / sqrt(x[i]**2+y[i]**2+z[i]**2)
					v = array([vx[i], vy[i], vz[i]])
					v_los = dot(e_los, v)
					d = dot(p, n_center)
					p_projected = p - d * n_center
					r = sqrt(p_projected[0]**2 + p_projected[1]**2 + p_projected[2]**2)
					#print d, r
					Rs.append(r)
					
				box()
				histogram(log10(Rs), datamin=-2, datamax=2, binwidth=0.1/3)
				logRs = arange(-2, 2, 0.1/10)
				Rs = 10**logRs
				graph(logRs, self.galaxy.light_model.densityR(Rs, M=1.)*Rs**2 * log(10) * 2 * pi * n * 0.1/3)
				draw()
					
			print "do"	
			do(0)
				
			
		finally:
			os.chdir(oldcwd)		
			
	

class GadgetCreateIC(object):
	def __init__(self, filename, galaxy, galaxy_host, orbit, df_sampling, N_star, N_dm):
		self.filename = filename
		self.galaxy = galaxy
		self.galaxy_host = galaxy_host
		self.df_sampling = df_sampling
		self.orbit = orbit
		self.N_star = N_star
		self.N_dm = N_dm
	
	def run(self, args, opts, scope):
		if os.path.exists(self.filename):
			logger.info("file %s already exists, skipping" % self.filename)
			return
		else:
			logger.info("creating Gadget intial conditions: %s" % self.filename)
		if self.orbit:
			self.orbit.find_orbit()
		self.df_sampling.load()
		if self.orbit:
			x0 = self.orbit.x0
			vt0 = self.orbit.vt0
		else:
			x0 = 0
			vt0 = 0 
		
		#"x,y,z,r3d,phi,theta,vx,vy,vz,vr,vphi,vtheta,sigmar,sigmat,Ekin,Epot,E,Lx,Ly,Lz,L,Lmax".split(",")
		samples = self.df_sampling.df_samples
		samples_dm = self.df_sampling.df_samples_dm
		samples = samples[:self.N_star]
		#N = len(samples)/10
		
		#print "filename", extra.args[0]
		f = h5py.File(self.filename, "w")
		header = f.create_group("Header")
		#header.attrs["PartType2"] = 1
		
		ar = zeros(6, dtype=int32)
		ar[2] = self.N_star
		if self.df_sampling.sample_dm:
			ar[1] = self.N_dm
		
		header.attrs["NumPart_ThisFile"] = ar
		header.attrs["NumPart_Total"] = ar
		
		ar = zeros(6, dtype=int32)
		header.attrs["NumPart_Total_HighWord"] = ar
		
		ar = zeros(6, dtype=float64) + 0.0
		ar[2] = 1.e6/self.N_star/1e10
		if self.df_sampling.sample_dm:
			x, y, z = samples_dm.x, samples_dm.y, samples_dm.z
			r = sqrt(x**2+y**2+z**2)
			rmax = max(r)
			print "rmax: %e" % rmax
			print "rmax/rs: %e" % (rmax / self.galaxy.profile_model.dm_profile.rs)
			
			rs = self.galaxy.profile_model.dm_profile.rs
			rmax = 100*rs # TODO hardcoded,  also in dm sampling 
			M_dm = self.galaxy.profile_model.dm_profile.enclosed_mass(rmax)
			print "rmax: %e" % rmax
			print "rmax/rs: %e" % (rmax / self.galaxy.profile_model.dm_profile.rs)
			print "dark matter mass: %e" % M_dm
			ar[1] = M_dm/self.N_dm/1e10#/1e10 #1e10 is gadget units
		header.attrs["MassTable"] = ar
		
		
		header.attrs["Time"] = 0.0
		header.attrs["NumFilesPerSnapshot"] = 1
		header.attrs["Flag_Entropy_ICs"] = 0
		
		#import pdb; pdb.set_trace()
		
		particles = f.create_group("PartType2")
		coordinates = transpose(array([samples.x+x0, samples.y, samples.z])).astype(float32)
		velocities = transpose(array([samples.vx, samples.vy+vt0, samples.vz])).astype(float32)
		particleIDs = arange(len(samples), dtype=uint32)
		int_energy = (velocities[:,0] * 1.0).astype(float32) 
		
		particles.create_dataset("Coordinates", data=coordinates)
		particles.create_dataset("ParticleIDs", data=particleIDs)
		particles.create_dataset("Velocities", data=velocities)
		
		
		if self.df_sampling.sample_dm:
			print "including dark matter particles"
			particles = f.create_group("PartType1") # Halo particles
			samples = self.df_sampling.df_samples_dm[:self.N_dm]
			coordinates = transpose(array([samples.x+x0, samples.y, samples.z])).astype(float32)
			velocities = transpose(array([samples.vx, samples.vy+vt0, samples.vz])).astype(float32)
			particleIDs = 10**6 + arange(len(samples), dtype=uint32)
			int_energy = (velocities[:,0] * 1.0).astype(float32) 
			
			particles.create_dataset("Coordinates", data=coordinates)
			particles.create_dataset("ParticleIDs", data=particleIDs)
			particles.create_dataset("Velocities", data=velocities)
			print "use for the dark matter NFW profile, rho0=%e" % self.galaxy_host.profile_model.dm_profile.rho0
			print "M(<100kpc) = %e" % self.galaxy_host.profile_model.dm_profile.enclosed_mass(100)
			
		
		#particles.create_dataset("InternalEnergy", data=int_energy)
		
		#obj.sample(cores=int(extra.opts.cores))
		#obj.broadcast()
		#obj.save()
						
					