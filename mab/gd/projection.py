from numpy import *
import scipy.integrate
import os
import time
import numpy
import mab.parallelize
import mab.utils.numpy
import pyublas

class ProjectOrbitLibrary(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, projection_matrix, binned_data):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.projection_matrix = projection_matrix
		self.binned_data = binned_data
		
	def init(self):
		self.binned_data.load()
		self.projection_matrix.load()
		
	def project(self):
		dfgrid = self.projection_matrix.storage_3d.dfgrid
		constraints = self.binned_data.n_constraints
		self.projectedmoments = zeros((dfgrid.dof, 5, constraints), dtype=numpy.float64)
		#projectedmoments = zeros((
		#print projection_matrix.Mr.shape, solution.storage_3d.moments3d.shape, projectedmoments.shape
				
		R = self.projection_matrix.light_model.arcsec_to_kpc(self.projection_matrix.aperture.aperture_rcenters)
		Rb = self.projection_matrix.light_model.arcsec_to_kpc(self.projection_matrix.aperture.aperture_rborders)
		dR = Rb[1:] - Rb[0:-1]
		rs = 10**self.projection_matrix.storage_3d.x
		dlogr = self.projection_matrix.storage_3d.xborders[1] - self.projection_matrix.storage_3d.xborders[0]
		volume = 4*pi*rs**3*dlogr*log(10)
		
		rho3d = self.projection_matrix.storage_3d.moments3d[:,0,:] * 1.
		m2r = self.projection_matrix.storage_3d.moments3d[:,4,:] * 1.
		m2t = (self.projection_matrix.storage_3d.moments3d[:,5,:]+self.projection_matrix.storage_3d.moments3d[:,6,:])/2
		mask = self.projection_matrix.storage_3d.moments3d[:,0,:] > 0
		#m2r[mask] /=  solution.storage_3d.moments3d[:,0,:][mask]
		#m2t[mask] /=  solution.storage_3d.moments3d[:,0,:][mask]
		rho3dok = rho3d * 1
		rho3dok[rho3dok <= 0] = 1
		rho3d /= volume
		m2r /= volume
		m2t /= volume
		print self.projectedmoments.shape, rho3d.shape, self.projection_matrix.Ms.shape, (2*pi*R*dR).shape
		assert all(~isnan(rho3d))
		assert all(~isnan(self.projectedmoments))
		self.projectedmoments[:,0,:] = tensordot(rho3d, self.projection_matrix.Ms, axes=[(1,), (1,)]) *2*pi*R*dR
		self.projectedmoments[:,2,:] = tensordot(m2r, self.projection_matrix.Mr, axes=[(1,), (1,)])
		self.projectedmoments[:,2,:] += tensordot(m2t, self.projection_matrix.Mt,axes=[(1,), (1,)])
		
	def save(self):
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		filename = os.path.join(self.dirname, "intermediate", "projectedmoments.npy")
		save(filename, self.projectedmoments)
		

class ProjectionMatrix(object):
	def __init__(self, modelpath, name, storage_3d, aperture, light_model):
		self.modelpath = modelpath
		self.name = name
		self.storage_3d = storage_3d
		self.aperture = aperture
		self.light_model = light_model
	
	def init(self):
		self.aperture.load() 
	
	def load(self):
		self.aperture.load()
		self.storage_3d.load()
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		filenames = filenamebase + "_s.npy"
		filenamer = filenamebase + "_r.npy"
		filenamet = filenamebase + "_t.npy"
		self.Ms = load(filenames)
		self.Mr = load(filenamer)
		self.Mt = load(filenamet)
	
	def save(self):
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		filenames = filenamebase + "_s.npy"
		filenamer = filenamebase + "_r.npy"
		filenamet = filenamebase + "_t.npy"
		save(filenames, self.Ms)
		save(filenamer, self.Mr)
		save(filenamet, self.Mt)
		
	def calculate(self):
		def fM(i, j, c):
			Rc = self.light_model.arcsec_to_kpc(self.aperture.aperture_rcenters[i])
			#Rc1 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders[i])
			#Rc2 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders[i])
			r1 = 10**self.storage_3d.xborders[j]
			r2 = 10**self.storage_3d.xborders[j+1]
			rc = 10**self.storage_3d.x[j]
			if Rc > r2:
				return 0
			if Rc < r1:
				rmin = r1
			if Rc >= r1:
				rmin= Rc
			rmax = r2
			def rho(r):
				assert (r >= r1)
				assert (r <= r2)
				return (r > r1) * (r < r2) * self.light_model.densityr(r)/self.light_model.light_profile.M
			if c == 0:
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2)*rho(r) , rmin, rmax)
				return 2 * I/rho(rc)
			if c == 1:
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2) * (1-Rc**2/r**2) * rho(r) , rmin, rmax)
			if c == 2: 
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2) * (Rc**2/r**2) * rho(r), rmin, rmax)
			#print I, i, j
			if j == 0:
				print i
			#print galaxy.light_model.densityR(R1)
			return 2 * I / self.light_model.densityR(Rc, M=1)#/rho(rc)
					
		self.Ms = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.Mr = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.Mt = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		for i in range(len(self.aperture.aperture_rcenters)):
			for j in range(len(self.storage_3d.x)):
				self.Ms[i,j] = fM(i,j, 0)
				self.Mr[i,j] = fM(i,j, 1)
				if isnan(self.Mr[i,j]):
					import pdb
					pdb.set_trace()
				self.Mt[i,j] = fM(i,j, 2)
		print self.Mr.shape

class ProjectionMatrix2(object):
	def __init__(self, modelpath, name, storage_3d, aperture, light_model):
		self.modelpath = modelpath
		self.name = name
		self.storage_3d = storage_3d
		self.aperture = aperture
		self.light_model = light_model
	
	def init(self):
		self.aperture.load() 
	
	def load(self):
		self.aperture.load()
		self.storage_3d.load()
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		filenames = filenamebase + "_s.npy"
		filenamer = filenamebase + "_r.npy"
		filenamet = filenamebase + "_t.npy"
		filename4r = filenamebase + "_4r.npy"
		filename4t = filenamebase + "_4t.npy"
		filename4rt = filenamebase + "_4rt.npy"
		self.Ms = load(filenames)
		self.Mr = load(filenamer)
		self.Mt = load(filenamet)
		self.M4r = load(filename4r)
		self.M4r = load(filename4t)
		self.M4rt = load(filename4rt)
	
	def save(self):
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		print self.Ms.shape
		filenames = filenamebase + "_s.npy"
		filenamer = filenamebase + "_r.npy"
		filenamet = filenamebase + "_t.npy"
		filename4r = filenamebase + "_4r.npy"
		filename4t = filenamebase + "_4t.npy"
		filename4rt = filenamebase + "_4rt.npy"
		save(filenames, self.Ms)
		save(filenamer, self.Mr)
		save(filenamet, self.Mt)
		save(filename4r, self.M4r)
		save(filename4t, self.M4t)
		save(filename4rt, self.M4rt)
		
	def calculate(self):
		def fM(i, j, c):
			#R1 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rcenters[i])
			#R2 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rcenters[i])
			R1 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders[i])
			R2 = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders[i+1])
			r1 = 10**self.storage_3d.xborders[j]
			r2 = 10**self.storage_3d.xborders[j+1]
			rc = 10**self.storage_3d.x[j]
			#rmin = r1
			rmax = r2
			Rmin = R1
			Rmax = R2
			if R1 > r2: # no overlap since R>r2
				return 0
			if R1 > r2: # no overlap since R>r2
				return 0
			#if R1 < r1: # avoid zero-interation part where R<r1
			#	Rmin = r1
			if R2 >= r2: # avoid zero-interation part, where R>r2
				Rmax = r2
			
			#print R1, R2
			#print Rmin, Rmax, rmax, r1, r2
			if c == 0:
				#mass_ij, err = scipy.integrate.dblquad(lambda r,R: 2 * pi * R * 2 * r / sqrt(r**2-R**2) , Rmin, Rmax, lambda R: max(R, r1), lambda R: rmax)
				def f(R):
					rmin = max(R, r1) 
					#rmax = rmax
					return sqrt(rmax**2-R**2) - sqrt(rmin**2-R**2) 
				mass_ij, err = scipy.integrate.quad(lambda R: 4 * pi * R * f(R) , Rmin, Rmax)
				#mass_ij, err = scipy.integrate.quad(lambda R,r: 1 , Rmin, Rmax, lambda R: R, lambda R: rmax)
				return mass_ij
			if c == 1:
				#mass_ij, err = scipy.integrate.dblquad(lambda r,R: 2 * pi * R * 2 * r / sqrt(r**2-R**2) , Rmin, Rmax, lambda R: max(R, r1), lambda R: rmax)
				def f(R):
					rmin = max(R, r1) 
					#rmax = rmax
					return -R * (arctan(R/sqrt(rmax**2-R**2))  - arctan(R/sqrt(rmin**2-R**2))) 
				mass_ij, err = scipy.integrate.quad(lambda R: 4 * pi * R * f(R) , Rmin, Rmax)
				#mass_ij, err = scipy.integrate.quad(lambda R,r: 1 , Rmin, Rmax, lambda R: R, lambda R: rmax)
				return mass_ij
			if c == 2:
				#mass_ij, err = scipy.integrate.dblquad(lambda r,R: 2 * pi * R * 2 * r / sqrt(r**2-R**2) , Rmin, Rmax, lambda R: max(R, r1), lambda R: rmax)
				def f(R):
					rmin = max(R, r1) 
					#rmax = rmax
					return -R * (arctan(R/sqrt(rmax**2-R**2))  - arctan(R/sqrt(rmin**2-R**2)))/2 +\
						sqrt(rmax**2-R**2)*R**2/(2*rmax**2) - sqrt(rmin**2-R**2)*R**2/(2*rmin**2) 
				mass_ij, err = scipy.integrate.quad(lambda R: 4 * pi * R * f(R) , Rmin, Rmax)
				#mass_ij, err = scipy.integrate.quad(lambda R,r: 1 , Rmin, Rmax, lambda R: R, lambda R: rmax)
				return mass_ij
			#2 * I/rho(rc)
				
			rmax = r2
			def rho(r):
				assert (r >= r1)
				assert (r <= r2)
				return (r > r1) * (r < r2) * self.light_model.densityr(r)/self.light_model.light_profile.M
			if c == 0:
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2)*rho(r) , rmin, rmax)
				return 2 * I/rho(rc)
			if c == 1:
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2) * (1-Rc**2/r**2) * rho(r) , rmin, rmax)
			if c == 2: 
				I, err = scipy.integrate.quad(lambda r: r / sqrt(r**2-Rc**2) * (Rc**2/r**2) * rho(r), rmin, rmax)
			#print I, i, j
			if j == 0:
				print i
			#print galaxy.light_model.densityR(R1)
			return 2 * I / self.light_model.densityR(Rc, M=1)#/rho(rc)
					
		self.Ms = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.Mr = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.Mt = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.M4r = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.M4t = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		self.M4rt = zeros((len(self.aperture.aperture_rcenters), len(self.storage_3d.x)))
		for i in range(len(self.aperture.aperture_rcenters)):
			print i+1, "of", len(self.aperture.aperture_rcenters)
			for j in range(len(self.storage_3d.x)):
				scalar = fM(i,j, 0)
				two = fM(i,j, 1)
				four = fM(i,j, 2)
				#print four
				self.Ms[i,j] = scalar
				self.Mr[i,j] = scalar - two
				self.Mt[i,j] = two
				self.M4r[i,j] = scalar - 2*two + four
				self.M4t[i,j] = four
				self.M4rt[i,j] = two - four
				
				#print self.Mr[i,j]
				#if isnan(self.Mr[i,j]):
				#	import pdb
				#	pdb.set_trace()
				#self.Mt[i,j] = fM(i,j, 2)
		#print self.Mr.shape
				
				
class ProjectionMatrix1d_to_1d_MC(object):
	def __init__(self, modelpath, name, gridR, gridr, N, light_model):
		self.modelpath = modelpath
		self.name = name
		self.gridR = gridR
		self.gridr = gridr
		self.N = N
		#self.storage_3d = storage_3d
		#self.aperture = aperture
		self.light_model = light_model
	
	def run(self, argv, opts, scope):
		self.init()
		self.calculate(cores=opts.cores)
		self.save()
		
	def init(self):
		self.gridR.load() 
	
	def load(self):
		self.gridR.load()
		#sself.storage_3d.load()
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		filename = filenamebase +".npy"
		self.matrices = load(filename)
	
	def save(self):
		#return
		filenamebase = os.path.join(self.modelpath, "data", "projection_matrix_" +self.name)
		filename = filenamebase +".npy"
		save(filename, self.matrices)
		
	def calculate(self, cores=50, progressbar=True):
		scale = 50
		Nperjob = self.N / (cores*scale)
		seedoffset = 0
		def seed(id_number):
			if seedoffset is not None:
				s = id_number + seedoffset
			else:
				s = os.getpid() * 10000 + int(time.time() * 100000)
			#print "seed for %d is %d (%d is pid)" % (id_number, s, os.getpid())
			numpy.random.seed(s)
		
		n_R = self.gridR.length()
		n_r = self.gridr.length()
		n_types = 4
		
		matrices = mab.utils.numpy.mmapzeros((cores*scale, n_types, n_R, n_r), dtype=float)
		def do(i):
			uniform = random.random(Nperjob)
			r = self.gridr.uniform_transform(uniform)
			u = random.random(Nperjob)
			u = arange(Nperjob)/(0.+Nperjob)
			R = 2 * r * sqrt(u-u**2)
			R_arcsec = self.light_model.kpc_to_arcsec(R)
			rindices = self.gridr.findindex(r).astype(int32)
			Rindices = self.gridR.gridr.findindex(R_arcsec).astype(int32)
			weights = n_r * 1.0/self.N + 0 * rindices
			M = matrices[i,0]
			mab.gd.gdfast.add_to_matrix_indexed(M, Rindices, rindices, weights)
			weights = n_r * 1.0/self.N * (R**2/r**2)
			#weights = (R**2/r**2)
			M = matrices[i,1]
			mab.gd.gdfast.add_to_matrix_indexed(M, Rindices, rindices, weights)
			weights = n_r * 1.0/self.N * (R**4/r**4)
			#weights = (R**3/r**3)
			M = matrices[i,2]
			mab.gd.gdfast.add_to_matrix_indexed(M, Rindices, rindices, weights)
			weights = n_r * 1.0/self.N * (1-(R**2/r**2))**2
			#weights = (R**4/r**4)
			M = matrices[i,3]
			mab.gd.gdfast.add_to_matrix_indexed(M, Rindices, rindices, weights)
			#print rindices
			#print Rindices
			#mab.gdfast.place_in_matrix(matrices, R, 
			#sqrt(rmax**2-R**2) - sqrt(rmin**2-R**2) * 4 * pi * R * f(R) , Rmin, Rmax)
			#f = r / sqrt(r**2-Rc**2)*rho(r) , rmin, rmax)
		indices = range(cores*scale)
		if cores == 1:
			do(0)
		else:
			dopar = mab.parallelize.parallelize(cores=cores, info=progressbar, init=seed)(do)
			dopar(indices)
		self.matrices = sum(matrices, axis=0)