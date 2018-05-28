from numpy import *
import mab.parallelize
import mab.utils.numpy
import os
from mab.constants import *
from mab.astrounits  import *
import mab.gd.logging as logging

kpc_to_km = (1*KPC).asNumber(KM)
sec_in_gyr = (1*GYR).asNumber(S)

logger = logging.getLogger("gd.schw.stackel")


class OrbitLibrary(object):
	def __init__(self, dirname, orbit_integrator, initial_conditions, xy_gridder, Norbits=200, postfix=""):
		self.dirname = dirname
		self.orbit_integrator = orbit_integrator
		self.initial_conditions = initial_conditions
		self.xy_gridder = xy_gridder
		self.Norbits = Norbits
		grid = self.initial_conditions.grid
		#self.moments = [0, 1, 2]
		# 3 moments
		self.momentgrid = mab.utils.numpy.mmapzeros((grid.gridx.N, grid.gridy.N, grid.gridz.N, 5, xy_gridder.nx, xy_gridder.ny))
		
		self.filename = os.path.join(dirname, "momentgrid" + postfix +".npy")
		
		
	def save(self):
		logger.info("saving momentgrid to: %s" % self.filename)
		save(self.filename, self.momentgrid)
	
	def load(self):
		logger.info("loading momentgrid from: %s" % self.filename)
		self.momentgrid = load(self.filename)
		
		
	def do_orbit(self, i, j, k):
		x, y, vx, vy, E, i2, i3, i2_max, i3_max, Torbit = self.ic[i,j,k]
		Lz = sqrt(2*i2)
		N = 1
		q0 = zeros((N, 2))
		p0 = zeros((N, 2))
		dt = zeros((N))
		Lzs = zeros(N)
		q0[:,0] = [x]
		q0[:,1] = [y]
		p0[:,0] = [vx]
		p0[:,1] = [vy]
		Lzs[:] = [Lz]
		dt[:] = [Torbit/self.orbit_integrator.orbital_points*self.Norbits]
		q1, p1 = self.orbit_integrator.integrate(dt, Lzs, q0, p0)
		x = ravel(q1[0,:])
		y = ravel(q1[1,:])
		vx = ravel(p1[0,:])
		vy = ravel(p1[1,:])
		self.momentgrid[i,j,k,0] += self.xy_gridder(x, abs(y))
		self.momentgrid[i,j,k,1] += self.xy_gridder(x, abs(y), weights=vx**2)
		self.momentgrid[i,j,k,2] += self.xy_gridder(x, abs(y), weights=vy**2)
		#print "."
		
	def run(self, args, opts, scope):
		self.initial_conditions.load()
		self.ic = self.initial_conditions.initial_conditions
		
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def do(ij):
			#for j in range(self.ic.shape[1]):
			i = ij % self.ic.shape[0]
			j = ij / self.ic.shape[0]
			for k in range(self.ic.shape[2]):
				self.do_orbit(i, j, k)
		ijs_ = range(self.ic.shape[0]*self.ic.shape[1])
		if 0:
			for i in is_:
				do(i)
		else:
			do(ijs_)
		self.save()
		
class OrbitLibraryAnalytic(OrbitLibrary):
	def __init__(self, dirname, orbit_integrator, initial_conditions, xy_gridder, Norbits=200):
		OrbitLibrary.__init__(self, dirname, orbit_integrator, initial_conditions, xy_gridder, Norbits)
		self.filename = os.path.join(dirname, "momentgrid_analytic.npy")
		
		
	def do_orbit(self, i, j, k):
		x, y, vx, vy, E, i2, i3, i2_max, i3_max, Torbit = self.ic[i,j,k]
		profile = self.initial_conditions.profile
		x, y = self.xy_gridder.meshgrid()
		self.momentgrid[i,j,k,0] += profile.density_orbit_xy(x, y, E, i2, i3)
		#x = ravel(q1[0,:])
		#y = ravel(q1[1,:])
		#vx = ravel(p1[0,:])
		#vy = ravel(p1[1,:])
		#self.momentgrid[i,j,k,0] += self.xy_gridder(x, abs(y))

class InitialConditions(object):
	def __init__(self, dirname, profile, grid, angle=pi/2):
		self.profile = profile
		self.grid = grid
		self.angle = angle
		
		
		self.integrals = mab.utils.numpy.mmapzeros((grid.gridx.N, grid.gridy.N, grid.gridz.N, 3))
		self.integrals_dither = mab.utils.numpy.mmapzeros((grid.subgrid.gridx.N, grid.subgrid.gridy.N, grid.subgrid.gridz.N, 3))
		
		names = "x, y, vx, vy, E, i2, i3, i2_max, i3_max, Torbit, la_max".split(",")
		names = [k.strip() for k in names]
		types = ["f8"] * len(names)
		dtype_initial_conditions = zip(names, types)
		#dithered = False
		#grud
		self.initial_conditions = mab.utils.numpy.mmapzeros((grid.subgrid.gridx.N, grid.subgrid.gridy.N, grid.subgrid.gridz.N), dtype=dtype_initial_conditions)
		self.initial_conditions = self.initial_conditions.view(recarray)
		
		self.filename = os.path.join(dirname, "ic.npy")
		
		
	def save(self):
		logger.info("saving initial conditions to: %s" % self.filename)
		save(self.filename, self.initial_conditions)
	
	def load(self):
		logger.info("loading initial conditions from: %s" % self.filename)
		self.initial_conditions = load(self.filename).view(recarray)
		
	def run(self, args, opts, scope):
		#self.find_integrals(self.grid, self.integrals, scope)
		self.find_integrals(self.grid.subgrid, self.integrals_dither, scope)
		self.save()
		
		
	def find_integrals(self, grid, integrals_grid, scope):
		r = grid.gridx.centers
		x = r * cos(self.angle)
		y = r * sin(self.angle)
		Epot = self.profile.potential_xy(x, y)
		#for i in range(grid.gridx.N):
		
		#@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def do(i):
			debug = False
			#print i, grid.gridx.centers[i], Epot[i]
			E = Epot[i]
			i3_max, tau_max_i3 = self.profile.find_i3_max(Epot[i])
			if debug:
				print "E", Epot[i]
			#print "i3_max", i3_max
			if debug:
				print "i3_max", i3_max, tau_max_i3
			for k in range(grid.gridz.N):
				scale = 1.0
				i3 = grid.gridz.centers[k] * i3_max * scale
				if debug:
					print "i3", i3_max, grid.gridz.centers[k]
				i2_max, tau_max_i2 = self.profile.find_i2_max(Epot[i], i3)
				if debug:
					print "i2_max", i2_max
				for j in range(grid.gridy.N):
					i2 = grid.gridy.centers[j] * i2_max * scale
					#print integrals_grid[i,j,k].shape
					integrals = [Epot[i], i2, i3]
					#print integrals, i2
					integrals_grid[i,j,k,:] = integrals
					if debug:
						print "tau_max", tau_max_i2, self.profile.V_eff_tau(tau_max_i2, i2, i3)
					tau = tau_max_i2
					if debug:
						print tau, self.profile.alpha, self.profile.gamma
						print i2/(tau+self.profile.alpha), i3 / (tau + self.profile.gamma), self.profile.G_tau(tau)
						print "lala", self.profile.gamma,(tau+self.profile.gamma)
						print sqrt(-self.profile.gamma/(tau+self.profile.gamma)),arctan(sqrt((tau+self.profile.gamma)/-self.profile.gamma))
					
					
					diff = self.profile.V_eff_tau(tau_max_i2, i2, i3) - Epot[i]
					assert self.profile.V_eff_tau(tau_max_i2, i2, i3) <= Epot[i], "error: %g" % diff
					if debug:
						print ">",grid.gridy.centers[j], grid.gridz.centers[k], tau_max_i2, tau_max_i3
					x, y = self.profile.conf_to_xy(-self.profile.gamma+1e-10, tau_max_i2)
					#print i,j,k,x
					#if isnan(x):
					#	x = 0
					#	print "for", i, j, "x is zero"
					assert not isnan(x)
					assert not isnan(y)
					if debug:
						print x, y
					vx, vy = self.profile.v_xy(x, y, E, i2, i3)
					assert not isnan(vx)
					assert not isnan(vy)
					self.initial_conditions.x[i,j,k] = x
					self.initial_conditions.y[i,j,k] = y
					self.initial_conditions.vx[i,j,k] = vx
					self.initial_conditions.vy[i,j,k] = vy
					self.initial_conditions.E[i,j,k] = E
					self.initial_conditions.i2[i,j,k] = i2
					self.initial_conditions.i3[i,j,k] = i3
					self.initial_conditions.i2_max[i,j,k] = i2_max
					self.initial_conditions.i3_max[i,j,k] = i3_max
					self.initial_conditions.la_max[i,j,k] = tau_max_i2
					#print x, y, vx, vy, E, i2, i3, i2_max, i3_max
					r = sqrt(x**2+y**x)
					v = sqrt(vx**2+vy**2)
					T = r/v * kpc_to_km
					#print T, " = ", T/sec_in_gyr, "gyr"
					self.initial_conditions.Torbit[i,j,k] = T
		#do(range(grid.gridx.N))
		for i in range(grid.gridx.N):
			print i
			do(i)
		 
		