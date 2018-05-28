# -*- coding: utf-8 -*-
import os
import mab.gd.gdfast
import numpy
from numpy import *
import mab.gd.schw.aperture
import bz2
import mab.gd.logging as logging

logger = logging.getLogger("gd.schw.output")
scaling_factor = 1#0.999999918453/0.840580815266 #1.1345025156398363/0.961596185691/1.01356943378 #/3.79757705745/0.250000999817/0.88366506055

class MomentsIntrinsicSpherical(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, xmin, xmax, dfgrid, nr, f=lambda x: log10(x)):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.xmin = xmin
		self.xmax = xmax
		self.dfgrid = dfgrid
		self.nr = nr
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.binner = mab.gd.gdfast.Binner(self.xmin, self.xmax)
		self.n = 0
		self.filename = os.path.join(self.dirname, "intermediate", "moments3d.npy")
		self.x = (arange(self.nr) + 0.5 )/ (self.nr) * (xmax - xmin) + xmin
		self.xborders = (arange(self.nr+1.))/ (self.nr) * (xmax - xmin) + xmin
		self.firststore = True
		self.f = f
	
	def init(self):
		self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.n_orbits, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d = numpy.memmap(self.filename, dtype='float64', mode='w+', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr))
	
	def load(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.moments3d = numpy.load(f)
		#self.initial_conditions = self.initial_conditions.view(recarray)
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
		numpy.save(f, self.moments3d)
	
	def store(self, x, y, z, vx, vy, vz, orbitnr, u, v):
		if self.firststore:
			self.vr, self.vphi, self.vtheta = zeros((3, len(x)))
			self.firststore = False
		mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, self.vr, self.vphi, self.vtheta)
		r = sqrt(x*x+y*y+z*z)
		rp = self.f(r)
		self.binner.bin3d(rp, self.vr, self.vphi, self.vtheta, self.moments3d[orbitnr])
		self.n += len(r)
	
	def normalize(self, orbitnr):
		self.moments3d[orbitnr] /= self.n
		self.n = 0
	
class MomentsApertureStorage(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, aperture, moments, dfgrid, binned_data, **kwargs):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname) 
		self.aperture = aperture
		self.moments = moments
		self.dfgrid = dfgrid
		self.binned_data = binned_data
		self.filename = os.path.join(self.dirname, "intermediate", "projectedmoments.npy")
	
	def init(self):
		self.aperture.load()
		self.binned_data.load()
		constraints = self.binned_data.n_constraints
		#self.projectedmoments = numpy.memmap(self.filename, dtype='float64', mode='w+', shape=(self.dfgrid.n_orbits, self.moments, constraints))
		self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.n_orbits, self.moments, constraints), dtype=numpy.float64)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMomentsPolar_IrrReg(self.aperture.aperture_fast(), self.aperture.constraintnrgrid, self.projectedmoments)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.binned_data.aperture_index_to_constraint, self.projectedmoments)
		self.aperture_index_to_constraint = self.binned_data.aperture_index_to_constraint.astype(int32)
		#if isinstance(self.aperture, mab.gd.schw.aperture.Polar):
		#	print "polar aperture!, use optimized version"
		#	self.aperturestore = mab.gd.gdfast.ApertureStorageMoments_Polar2d_IrrReg_s(self.aperture.aperture_fast_s(), self.aperture_index_to_constraint, self.projectedmoments)
		#else: 
		
		self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.aperture_index_to_constraint, self.projectedmoments)
		self.n = 0
	
	def store(self, x, y, vz, orbitnr, u, v):
		self.aperturestore.bin(x, y, vz, orbitnr)
		self.n += len(x)
		
	def load(self, compressed=True):
		self.binned_data.load()
		constraints = self.binned_data.n_constraints
		
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.projectedmoments = numpy.load(f)
		#self.projectedmoments = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 5, constraints)))
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
		numpy.save(f, self.projectedmoments)
			
	def normalize(self, orbitnr):
		self.projectedmoments[orbitnr] /= self.n
		self.n = 0
		

class MomentsIntrinsicSphericalShaped(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, xmin, xmax, dfgrid, nr, f=lambda x: log10(x), finv=lambda x: 10**x, postfix="", next=None):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.xmin = xmin
		self.xmax = xmax
		self.dfgrid = dfgrid
		self.nr = nr
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.binner = mab.gd.gdfast.Binner(self.xmin, self.xmax)
		self.n = 0
		self.postfix = postfix
		self.filename = os.path.join(self.dirname, "intermediate", "moments3d" +self.postfix +".npy")
		self.x = (arange(self.nr) + 0.5 )/ (self.nr) * (xmax - xmin) + xmin
		self.xborders = (arange(self.nr+1.))/ (self.nr) * (xmax - xmin) + xmin
		self.firststore = True
		self.f = f
		self.finv = finv
		self.rborders = self.finv(self.xborders)
		self.rs = self.finv(self.x)
		self.weights = [0 for i in range(self.dfgrid.dof_per_cell)]
		# moments: v^powers1[m] * v^powers2[m]
		self.powers1 = numpy.array([0, 2, 0, 2, 4, 0], dtype=int32)
		self.powers2 = numpy.array([0, 0, 2, 2, 0, 4], dtype=int32)
		self.v00index = 0
		self.v20index = 1
		self.v02index = 2
		self.v22index = 3
		self.v40index = 4
		self.v04index = 5
		
		#self.momentmap = dict(["00 20 02 22 40 04"])
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		for i in range(self.dfgrid.dof_per_cell):
			for d1 in range(dither1):
				for d2 in range(dither2):
					u = float(d1 + 0.5) / dither1
					v = float(d2 + 0.5) / dither2
					weight = self.dfgrid.basis(i, u, v)
					self.weights[i] += weight
		self.next = next
		self.momentcount = 1+3+3
		if 1:
			self.momentcount = len(self.powers1)
		
	
	def init(self):
		self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.momentcount, self.nr), dtype=numpy.float64)
		self.moments3d_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.momentcount, self.nr), dtype=numpy.float64)
		self.moments3d_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.momentcount, self.nr), dtype=numpy.float64)
		#self.moments3dgrid = self.moments3d.reshape((self.dfgrid.n_I2, self.dfgrid.n_I1, self.momentcount, self.nr))
		if self.next:
			self.next.init()
		#self.moments3d = numpy.memmap(self.filename, dtype='float64', mode='w+', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr))
	
	def load(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.moments3d = numpy.load(f)
		#print self.dfgrid.n_I2, self.dfgrid.n_I1, self.dfgrid.dof, self.moments3d.shape
		#self.moments3dgrid = self.moments3d.reshape((self.dfgrid.n_I2, self.dfgrid.n_I1, self.momentcount, self.nr))
		if self.next:
			self.next.load(compressed)
		#self.initial_conditions = self.initial_conditions.view(recarray)
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		for i1 in range(self.dfgrid.n_I1):
			for i2 in range(self.dfgrid.n_I2):
				orbitnr = self.dfgrid.index_to_orbitnr(i1, i2, 0)
				for k in range(self.dfgrid.dof_per_cell):
					dof_index = self.dfgrid.dof_index(i1, i2, k)
					self.moments3d[dof_index] += self.moments3d_basis[orbitnr, k]
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
		logger.info("saving orbit library: %s" % filename) 
		numpy.save(f, self.moments3d)
		if self.next:
			self.next.save(compressed)
		
	
	def store(self, x, y, z, vx, vy, vz, orbitnr, u, v):
		if self.firststore:
			self.vr, self.vphi, self.vtheta = zeros((3, len(x)))
			self.firststore = False
		mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, self.vr, self.vphi, self.vtheta)
		r = sqrt(x*x+y*y+z*z)
		rp = self.f(r)
		#import pdb
		#pdb.set_trace()
		self.moments3d_unshaped[orbitnr] = 0
		if 1:
			self.binner.bin2d(rp, self.vr, self.vphi, self.moments3d_unshaped[orbitnr], self.powers1, self.powers2)
			self.binner.bin2d(rp, self.vr, -self.vphi, self.moments3d_unshaped[orbitnr], self.powers1, self.powers2)
		else:
			self.binner.bin3d(rp, self.vr, self.vphi, self.vtheta, self.moments3d_unshaped[orbitnr])
			self.binner.bin3d(rp, self.vr, -self.vphi, self.vtheta, self.moments3d_unshaped[orbitnr])
		for i in range(self.dfgrid.dof_per_cell):
			weight = self.dfgrid.basis(i, u, v)
			assert u > 0
			assert v > 0
			assert u < 1
			assert v < 1
			self.moments3d_basis[orbitnr,i,:,:] += self.moments3d_unshaped[orbitnr] * weight *\
			scaling_factor #* 5 /4 * 0.95 * 1.094303514738155 / 1.0972532454693291 # / self.weights[i]
		self.n += len(r) * 2
		if self.next:
			self.next.store(x, y, z, vx, vy, vz, orbitnr, u, v)
		
	
	def normalize(self, orbitnr):
		self.moments3d_basis[orbitnr] /= self.n
		self.n = 0
		if self.next:
			self.next.normalize(orbitnr)
		

class MomentsIntrinsicSphericalShapedAperture(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, aperture, postfix="", next=None):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.aperture = aperture
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.binner = mab.gd.gdfast.BinnerAperture(self.xmin, self.xmax)
		self.n = 0
		self.postfix = postfix
		self.filename = os.path.join(self.dirname, "intermediate", "moments3d" +self.postfix +".npy")
		#self.x = (arange(self.nr) + 0.5 )/ (self.nr) * (xmax - xmin) + xmin
		#self.xborders = (arange(self.nr+1.))/ (self.nr) * (xmax - xmin) + xmin
		self.firststore = True
		#self.f = f
		self.weights = [0 for i in range(self.dfgrid.dof_per_cell)]
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		for i in range(self.dfgrid.dof_per_cell):
			for d1 in range(dither1):
				for d2 in range(dither2):
					u = float(d1 + 0.5) / dither1
					v = float(d2 + 0.5) / dither2
					weight = self.dfgrid.basis(i, u, v)
					self.weights[i] += weight
		self.next = next
	
	def init(self):
		self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.dof, 1+3+3, self.nr), dtype=numpy.float64)
		self.moments3d_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, 1+3+3, self.nr), dtype=numpy.float64)
		self.moments3d_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, 1+3+3, self.nr), dtype=numpy.float64)
		if self.next:
			self.next.init()
		#self.moments3d = numpy.memmap(self.filename, dtype='float64', mode='w+', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr))
	
	def load(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.moments3d = numpy.load(f)
		if self.next:
			self.next.load(compressed)
		#self.initial_conditions = self.initial_conditions.view(recarray)
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		for i1 in range(self.dfgrid.n_I1):
			for i2 in range(self.dfgrid.n_I2):
				orbitnr = self.dfgrid.index_to_orbitnr(i1, i2, 0)
				for k in range(self.dfgrid.dof_per_cell):
					dof_index = self.dfgrid.dof_index(i1, i2, k)
					self.moments3d[dof_index] += self.moments3d_basis[orbitnr, k]
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
		numpy.save(f, self.moments3d)
		if self.next:
			self.next.save(compressed)
		
	
	def store(self, x, y, z, vx, vy, vz, orbitnr, u, v):
		if self.firststore:
			self.vr, self.vphi, self.vtheta = zeros((3, len(x)))
			self.firststore = False
		mab.gd.gdfast.velocity_cartesian_to_spherical(x, y, z, vx, vy, vz, self.vr, self.vphi, self.vtheta)
		r = sqrt(x*x+y*y+z*z)
		rp = self.f(r)
		#import pdb
		#pdb.set_trace()
		self.moments3d_unshaped[orbitnr] = 0
		self.binner.bin3d(rp, self.vr, self.vphi, self.vtheta, self.moments3d_unshaped[orbitnr])
		self.binner.bin3d(rp, self.vr, -self.vphi, self.vtheta, self.moments3d_unshaped[orbitnr])
		for i in range(self.dfgrid.dof_per_cell):
			weight = self.dfgrid.basis(i, u, v)
			assert u > 0
			assert v > 0
			assert u < 1
			assert v < 1
			self.moments3d_basis[orbitnr,i,:,:] += self.moments3d_unshaped[orbitnr] * weight *\
			scaling_factor #* 5 /4 * 0.95 * 1.094303514738155 / 1.0972532454693291 # / self.weights[i]
		self.n += len(r) * 2
		if self.next:
			self.next.store(x, y, z, vx, vy, vz, orbitnr, u, v)
		
	
	def normalize(self, orbitnr):
		self.moments3d_basis[orbitnr] /= self.n
		self.n = 0
		if self.next:
			self.next.normalize(orbitnr)

class MomentsApertureProjectionMatrix(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, projection, aperture, storage_3d, dfgrid, **kwargs):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.projection = projection
		self.aperture = aperture
		self.storage_3d = storage_3d
		#self.binned_data = binned_data
		self.dfgrid = dfgrid
		
		
	def init(self):
		pass
	
	def store(self, x, y, vz, orbitnr, u, v):
		pass

	def load(self, compressed=True):
		self.storage_3d.load()
		#self.binned_data.load()
		#self.massesprojected = self.binned_data.moments[0]
		self.projection.load()
		#self.moments_3d.init()
		#pass
		#constraints = self.binned_data.n_constraints
		#print self.storage_3d.moments3d.shape
		#print self.projection.matrices.shape
		constraints = self.projection.matrices.shape[1]
		#dsa
		self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.dof, 5, constraints))
		M = self.projection.matrices
		Ms = M[0]
		Mr = M[0] - M[1]
		Mt = M[1]
		u0, u2, u4, uw = M[0], M[1], M[2], M[3]
		
		rs = 10**self.storage_3d.x
		#rho = solution.light_model.densityr(rs)/solution.light_model.light_profile.M
		dlogr = self.storage_3d.xborders[1] - self.storage_3d.xborders[0]
		volume = 4*pi*rs**3*dlogr*log(10)
		#print self.storage_3d, self.projection.storage_3d
		#print self.storage_3d.moments3d.shape, volume.shape
		#print (self.storage_3d.moments3d[:,0,:]/volume).shape
		#print Ms.shape
		logger.debug("Ms.shape = %r" % (Ms.shape,))
		logger.debug("self.storage_3d.moments3d[:,0,:].shape = %r" % (self.storage_3d.moments3d[:,0,:].shape,))
		#print self.projectedmoments.shape
		#self.projectedmoments[:,0,:] = tensordot(self.projection.Ms, self.storage_3d.moments3d[:,0,:]/volume, axes=[(1,), (1,)]).T
		self.projectedmoments[:,0,:] = tensordot(Ms, self.storage_3d.moments3d[:,0,:], axes=[(1,), (1,)]).T
		#self.projectedmoments[:,0,:] = dot(Ms, self.storage_3d.moments3d[:,0,:])
		#self.projectedmoments[0] = tensordot(projection_matrix.Ms, self.storage_3d.moments3d[:,0,:]/volume, axes=[(1,), (1,)])
		
		#rho = solution.light_model.densityr(rs)/solution.light_model.light_profile.M
		
		#print projection_matrix.Ms.shape, pmasses.shape
		#projection_matrix.Ms.T/pmasses
		#moments2d_varr = tensordot((self.projection.Mr.T).T, self.storage_3d.moments3d[:,4,:]/volume, axes=[(1,), (1,)])
		#moments2d_vart = tensordot((self.projection.Mt.T).T, (self.storage_3d.moments3d[:,5,:]+self.storage_3d.moments3d[:,6,:])/2/volume, axes=[(1,), (1,)])
		if 1:
			v20 = self.storage_3d.moments3d[:,self.storage_3d.v20index,:]
			v02 = self.storage_3d.moments3d[:,self.storage_3d.v02index,:]/2
			v22 = self.storage_3d.moments3d[:,self.storage_3d.v22index,:]/2
			v40 = self.storage_3d.moments3d[:,self.storage_3d.v40index,:]
			v04 = self.storage_3d.moments3d[:,self.storage_3d.v04index,:]*3./8
			#moments2d_varr = tensordot(Mr, v20, axes=[(1,), (1,)])
			#moments2d_vart = tensordot(Mt, v02, axes=[(1,), (1,)])
		else:
			moments2d_varr = tensordot(Mr, self.storage_3d.moments3d[:,4,:], axes=[(1,), (1,)])
			moments2d_vart = tensordot(Mt, (self.storage_3d.moments3d[:,5,:]+self.storage_3d.moments3d[:,6,:])/2, axes=[(1,), (1,)])
		#moments2d_4r = tensordot((self.projection.Mr.T).T, self.storage_3d.moments3d[:,4,:]/volume, axes=[(1,), (1,)])
		#moments2d_4t = tensordot((self.projection.Mt.T).T, (self.storage_3d.moments3d[:,5,:]+self.storage_3d.moments3d[:,6,:])/2/volume, axes=[(1,), (1,)])
		self.projectedmoments[:,2,:] = \
				( \
					tensordot((u0 - u2), v20, axes=[(1,), (1,)]) +\
					tensordot((u2),      v02, axes=[(1,), (1,)])\
				).T
		self.projectedmoments[:,4,:] =\
				(\
					tensordot((uw), v40, axes=[(1,), (1,)]) +\
					tensordot((u4),         v04, axes=[(1,), (1,)]) + \
					tensordot(6*(u2-u4),  v22, axes=[(1,), (1,)])\
				).T
					#tensordot(6*(u2-u4),  v22, axes=[(1,), (1,)])\
		#self.projectedmoments[:,2,:] =\
		#	tensordot((self.projection.Mr.T).T, self.storage_3d.moments3d[:,4,:]/volume, axes=[(1,), (1,)]) +\
		#	tensordot((self.projection.Mt.T).T, (self.storage_3d.moments3d[:,5,:]+self.storage_3d.moments3d[:,6,:])/2/volume, axes=[(1,), (1,)])
		
		#print c.shape, moments2d0.shape
		#var_los_r = tensordot(c, moments2d_varr, axes=[(0,), (1,)])
		#var_los_t = tensordot(c, moments2d_vart, axes=[(0,), (1,)])
		
		
		
	
	def save(self, compressed=True):
		pass

	def normalize(self, orbitnr):
		pass		
		
		

class MomentsApertureStorageShaped(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, aperture, moments, dfgrid, binned_data, postfix="", next=None, **kwargs):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname) 
		self.aperture = aperture
		self.moments = moments
		self.dfgrid = dfgrid
		self.binned_data = binned_data
		self.filename = os.path.join(self.dirname, "intermediate", "projectedmoments" +postfix +".npy")
		self.weights = [0 for i in range(self.dfgrid.dof_per_cell)]
		self.kpc_to_arcsec = self.aperture.light_model.kpc_to_arcsec(1.)
		self.next = next
		dither1 = self.dfgrid.dither
		dither2 = self.dfgrid.dither
		for i in range(self.dfgrid.dof_per_cell):
			for d1 in range(dither1):
				for d2 in range(dither2):
					u = float(d1 + 0.5) / dither1
					v = float(d2 + 0.5) / dither2
					weight = self.dfgrid.basis(i, u, v)
					self.weights[i] += weight
		#print "weights", self.weights
	
	def init(self):
		self.aperture.load()
		self.binned_data.load()
		constraints = self.binned_data.n_constraints
		self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.moments, constraints))
		self.projectedmoments_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.moments, constraints))
		self.projectedmoments_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.moments, constraints))
		#self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.dof, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, 1+3+3, self.nr), dtype=numpy.float64)
		#self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self,dfgrid.dof_per_cell, self.moments, constraints), dtype=numpy.float64)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMomentsPolar_IrrReg(self.aperture.aperture_fast(), self.aperture.constraintnrgrid, self.projectedmoments)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.binned_data.aperture_index_to_constraint, self.projectedmoments)
		self.aperture_index_to_constraint = self.binned_data.aperture_index_to_constraint.astype(int32)
		#if isinstance(self.aperture, mab.gd.schw.aperture.Polar):
		#	print "polar aperture!, use optimized version"
		#	self.aperturestore = mab.gd.gdfast.ApertureStorageMoments_Polar2d_IrrReg_s(self.aperture.aperture_fast_s(), self.aperture_index_to_constraint, self.projectedmoments)
		#else: 
		
		self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.aperture_index_to_constraint, self.projectedmoments_unshaped)
		self.n = 0
		
		if self.next:
			self.next.init()
	
	def store(self, x, y, vz, orbitnr, u, v):
		self.projectedmoments_unshaped[orbitnr] = 0
		self.aperturestore.bin(x, y, vz, orbitnr)
		for i in range(self.dfgrid.dof_per_cell):
			weight = self.dfgrid.basis(i, u, v)
			#print weight, 
			assert any(isinf(x)) == False 
			assert any(isinf(y)) == False 
			assert any(isinf(vz)) == False 
			assert any(isinf(self.projectedmoments_unshaped)) == False 
			#print u, v
			assert u > 0
			assert v > 0
			assert u < 1
			assert v < 1
			self.projectedmoments_basis[orbitnr,i,:,:] += self.projectedmoments_unshaped[orbitnr] * weight *\
			 scaling_factor #1  * 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 # / self.weights[i]
		assert len(x) > 0
		self.n += len(x)
		if self.next:
			self.next.store(x, y, vz, orbitnr, u, v)
		
	def load(self, compressed=True):
		self.aperture.load()
		self.binned_data.load()
		constraints = self.binned_data.n_constraints
		
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.projectedmoments = numpy.load(f)
		#self.projectedmoments = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 5, constraints)))
		if self.next:
			self.next.load(compressed=compressed)
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		for i1 in range(self.dfgrid.n_I1):
			for i2 in range(self.dfgrid.n_I2):
				orbitnr = self.dfgrid.index_to_orbitnr(i1, i2, 0)
				for k in range(self.dfgrid.dof_per_cell):
					dof_index = self.dfgrid.dof_index(i1, i2, k)
					self.projectedmoments[dof_index] += self.projectedmoments_basis[orbitnr, k]
					assert any(isinf(self.projectedmoments_basis[orbitnr, k])) == False 
		assert any(isinf(self.projectedmoments)) == False 
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
		numpy.save(f, self.projectedmoments)
		if self.next:
			self.next.save(compressed=compressed)
			
	def normalize(self, orbitnr):
		#assert self.n > 0
		#if self.n > 0:
		self.projectedmoments_basis[orbitnr] /= self.n
		self.n = 0
		if self.next:
			self.next.normalize(orbitnr=orbitnr)
		
 	
class StorageMultiplexer(object):
	def __init__(self, **kwargs):
		self.objects = kwargs.values()
		
	def init(self):
		for obj in self.objects:
			obj.init()

	def store(self, x, y, vz, orbitnr, u, v):
		for obj in self.objects:
			obj.store(x, y, vz, orbitnr, u, v)

	def load(self, compressed=True):
		for obj in self.objects:
			obj.load(compressed=compressed)

	def save(self, compressed=True):
		for obj in self.objects:
			obj.save(compressed=compressed)

	def normalize(self, orbitnr):
		for obj in self.objects:
			obj.normalize(orbitnr)

class LosvdStorage(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, profile_model, light_model, Rmax, vmax, NR, Nv, dfgrid, rotations=1, **kwargs):
		self.modelpath = modelpath
		self.rotations = rotations
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname) 
		self.profile_model = profile_model
		self.light_model = light_model
		self.Rmax = Rmax
		if vmax is None:
			# get vmax from potential and
			# scale up vmax to avoid problems with small integration errors
			self.vmax = self.profile_model.vmax() * 1.05
		else:
			self.vmax = vmax
		self.NR = NR
		self.Nv = Nv
		self.dfgrid = dfgrid
		self.delta_v = self.vmax * 2 / self.Nv
		self.delta_R = self.Rmax / self.NR
		self.filename = os.path.join(self.dirname, "intermediate", "losvds.npy")
		self.filename_masses = os.path.join(self.dirname, "intermediate", "masses.npy")
		self.Rborders = arange(self.NR+1) / (0.0+self.NR) * (self.Rmax)
		self.Rcenters = (arange(self.NR)+0.5) / (0.0+self.NR) * (self.Rmax)
		self.vborders = arange(self.Nv+1) / (0.0+self.Nv) * (self.vmax*2)-self.vmax
		self.vcenters = (arange(self.Nv)+0.5) / (0.0+self.Nv) * (self.vmax*2)-self.vmax
		if 1:
			self.weights = [0 for i in range(self.dfgrid.dof_per_cell)]
			dither1 = self.dfgrid.dither
			dither2 = self.dfgrid.dither
			for i in range(self.dfgrid.dof_per_cell):
				for d1 in range(dither1):
					for d2 in range(dither2):
						u = float(d1 + 0.5) / dither1
						v = float(d2 + 0.5) / dither2
						weight = self.dfgrid.basis(i, u, v)
						self.weights[i] += weight
			#print "weights", self.weights
		self.gridr = mab.gd.gdfast.Grid1dRegular(0, self.light_model.kpc_to_arcsec(self.Rmax), self.NR)
		
		self.nphi = 1
		#self.gridphi = mab.gd.gdfast.Grid1dRegular(0, 2*pi, self.nphi)
		
		self.aperture = mab.gd.gdfast.AperturePolar2dNoAngle(self.gridr)
		
		def inrangefilter(stars):
			logger.info("checking if %d stars are in range of aperture and velocity range" % len(stars))
			stars_inrange = []
			for i in range(len(stars)):
				star = stars[i]
				#print star
				vindex = self.velocity_to_index(star.vlos) 
				if self.aperture.inrange(stars.xi[i], stars.eta[i]) and (vindex >= 0) and (vindex < self.Nv):
					stars_inrange.append(star)
			use_numpy = not isinstance(stars, mab.cvsfile.CsvObject)
			if use_numpy:
				stars_inrange = array(stars_inrange, dtype=stars.dtype).view(recarray)
			else:
				stars_inrange = mab.cvsfile.CsvObject(stars_inrange)
			logger.info("%d stars are in range of aperture and velocity range" % len(stars_inrange))
			return stars_inrange
		
		self.inrangefilter = inrangefilter
			
	def invrange(self, v):
		i = self.velocity_to_index(v)
		return (i >= 0) and (i < self.Nv)
		
	def velocity_to_index(self, v):
		v_index = numpy.int32(((v+self.vmax)/(2*self.vmax)) * self.Nv);
		return v_index
	
	def stars_to_apertureindices(self, stars):
		aperture_indices = []
		for i in range(len(stars)):
			#aperture_indices.append(self.aperture.findindex(stars.xi[i], stars.eta[i]))
			aperture_indices.append(self.aperture.findindex(stars[i].xi, stars[i].eta))
		aperture_indices = array(aperture_indices)
		return stars, aperture_indices
	
	def init(self):
		#self.aperture.load()
		#self.binned_data.load()
		#constraints = self.binned_data.n_constraints
		
		self.losvds = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.Nv, self.NR))
		self.masses = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.NR))
		self.losvds_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.Nv, self.NR))
		self.masses_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.NR))
		self.losvds_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.Nv, self.NR))
		self.masses_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.NR))
		
		self.aperturestore = mab.gd.gdfast.ApertureStorageLosvd(self.aperture, self.losvds_unshaped, self.masses_unshaped, self.vmax)
		
		#self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.dof, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, 1+3+3, self.nr), dtype=numpy.float64)
		#self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self,dfgrid.dof_per_cell, self.moments, constraints), dtype=numpy.float64)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMomentsPolar_IrrReg(self.aperture.aperture_fast(), self.aperture.constraintnrgrid, self.projectedmoments)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.binned_data.aperture_index_to_constraint, self.projectedmoments)
		#self.aperture_index_to_constraint = self.binned_data.aperture_index_to_constraint.astype(int32)
		#if isinstance(self.aperture, mab.gd.schw.aperture.Polar):
		#	print "polar aperture!, use optimized version"
		#	self.aperturestore = mab.gd.gdfast.ApertureStorageMoments_Polar2d_IrrReg_s(self.aperture.aperture_fast_s(), self.aperture_index_to_constraint, self.projectedmoments)
		#else: 
		
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.aperture_index_to_constraint, self.projectedmoments_unshaped)
		self.n = 0
	
	def store_test(self, r, vr, vt, orbitnr, u, v):
		if self.rotations > 0:
			self.losvds_unshaped[orbitnr] = 0
			self.masses_unshaped[orbitnr] = 0
			#if any(abs(vz) >= self.vmax):
				#print vz
				#print self.vmax
				#print vz[abs(vz) > self.vmax] 
			self.aperturestore.bin_rotate(r, vr, vt, orbitnr, self.rotations)
			if 1:
				for i in range(self.dfgrid.dof_per_cell):
					weight = self.dfgrid.basis(i, u, v)
					#print u, v
					assert u > 0
					assert v > 0
					assert u < 1
					assert v < 1
					self.losvds_basis[orbitnr,i,:,:] += self.losvds_unshaped[orbitnr] *\
					scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / self.weights[i]
					self.masses_basis[orbitnr,i,:] += self.masses_unshaped[orbitnr] *\
					scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / 
			self.n += len(r) * self.rotations
			
	def store(self, x, y, vz, orbitnr, u, v):
		if self.rotations > 0:
			return
		self.losvds_unshaped[orbitnr] = 0
		self.masses_unshaped[orbitnr] = 0
		#if any(abs(vz) >= self.vmax):
			#print vz
			#print self.vmax
			#print vz[abs(vz) > self.vmax] 
		self.aperturestore.bin(x, y, vz, orbitnr,0)
		if 1:
			for i in range(self.dfgrid.dof_per_cell):
				weight = self.dfgrid.basis(i, u, v)
				#print u, v
				assert u > 0
				assert v > 0
				assert u < 1
				assert v < 1
				self.losvds_basis[orbitnr,i,:,:] += self.losvds_unshaped[orbitnr] *\
				scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / self.weights[i]
				self.masses_basis[orbitnr,i,:] += self.masses_unshaped[orbitnr] *\
				scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / 
		self.n += len(x)
		
	def load(self, compressed=True):
		#self.binned_data.load()
		#constraints = self.binned_data.n_constraints
		
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = bz2.BZ2File(filename_masses, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.losvds = numpy.load(f)
		self.masses = numpy.load(f_masses)
		#self.projectedmoments = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 5, constraints)))
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		for i1 in range(self.dfgrid.n_I1):
			for i2 in range(self.dfgrid.n_I2):
				orbitnr = self.dfgrid.index_to_orbitnr(i1, i2, 0)
				for k in range(self.dfgrid.dof_per_cell):
					dof_index = self.dfgrid.dof_index(i1, i2, k)
					self.losvds[dof_index] += self.losvds_basis[orbitnr, k]
					self.masses[dof_index] += self.masses_basis[orbitnr, k]
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = bz2.BZ2File(filename_masses, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = file(filename_masses, "wb")
		logger.info("saving: %s" % filename)
		numpy.save(f, self.losvds)
		logger.info("saving: %s" % filename_masses)
		numpy.save(f_masses, self.masses)
			
	def normalize(self, orbitnr):
		self.losvds_basis[orbitnr] /= self.n
		self.masses_basis[orbitnr] /= self.n
		self.n = 0






class LosvdStorage2d(object):
	def __init__(self, modelpath, schwsetname, schwmodelname, profile_model, light_model, Rmax, vmax, NR, Nv, NLz, dfgrid, rotations=1, **kwargs):
		self.modelpath = modelpath
		self.rotations = rotations
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname) 
		self.profile_model = profile_model
		self.light_model = light_model
		self.Rmax = Rmax
		self.NLz = NLz
		if vmax is None:
			# get vmax from potential and
			# scale up vmax to avoid problems with small integration errors
			self.vmax = self.profile_model.vmax() * 1.05
		else:
			self.vmax = vmax
		self.NR = NR
		self.Nv = Nv
		self.dfgrid = dfgrid
		self.delta_v = self.vmax * 2 / self.Nv
		self.filename = os.path.join(self.dirname, "intermediate", "losvds_Lz.npy")
		self.filename_masses = os.path.join(self.dirname, "intermediate", "masses_Lz.npy")
		self.Rborders = arange(self.NR+1) / (0.0+self.NR) * (self.Rmax)
		self.vs = mgrid[-self.vmax:self.vmax:Nv*1j]
		if 1:
			self.weights = [0 for i in range(self.dfgrid.dof_per_cell)]
			dither1 = self.dfgrid.dither
			dither2 = self.dfgrid.dither
			for i in range(self.dfgrid.dof_per_cell):
				for d1 in range(dither1):
					for d2 in range(dither2):
						u = float(d1 + 0.5) / dither1
						v = float(d2 + 0.5) / dither2
						weight = self.dfgrid.basis(i, u, v)
						self.weights[i] += weight
			#print "weights", self.weights
		Rmax = self.light_model.kpc_to_arcsec(self.Rmax)
		self.gridx = mab.gd.gdfast.Grid1dRegular(-Rmax, Rmax, self.NR)
		self.gridy = mab.gd.gdfast.Grid1dRegular(-Rmax, Rmax, self.NR)
		self.masses_resize = ((-self.Rmax, -self.Rmax), (self.Rmax, self.Rmax))
		
		self.nphi = 1
		#self.gridphi = mab.gd.gdfast.Grid1dRegular(0, 2*pi, self.nphi)
		
		self.aperture = mab.gd.gdfast.ApertureGrid2d(self.gridx, self.gridy)
		
		def inrangefilter(stars):
			logger.info("checking if %d stars are in range of aperture and velocity range" % len(stars))
			stars_inrange = []
			for i in range(len(stars)):
				star = stars[i]
				#print star
				vindex = self.velocity_to_index(star.vlos) 
				if self.aperture.inrange(stars.xi[i], stars.eta[i]) and (vindex >= 0) and (vindex < self.Nv):
					stars_inrange.append(star)
			use_numpy = not isinstance(stars, mab.cvsfile.CsvObject)
			if use_numpy:
				stars_inrange = array(stars_inrange, dtype=stars.dtype).view(recarray)
			else:
				stars_inrange = mab.cvsfile.CsvObject(stars_inrange)
			logger.info("%d stars are in range of aperture and velocity range" % len(stars_inrange))
			return stars_inrange
		
		self.inrangefilter = inrangefilter
			
	def velocity_to_index(self, v):
		v_index = numpy.int32(((v+self.vmax)/(2*self.vmax)) * self.Nv);
		return v_index
	
	def stars_to_apertureindices(self, stars):
		aperture_indices = []
		for i in range(len(stars)):
			aperture_indices.append(self.aperture.findindex(stars.xi[i], stars.eta[i]))
		aperture_indices = array(aperture_indices)
		return stars, aperture_indices
	
	def init(self):
		#self.aperture.load()
		#self.binned_data.load()
		#constraints = self.binned_data.n_constraints
		
		self.losvds = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.NLz, self.Nv, self.aperture.length()))
		self.masses = mab.utils.numpy.mmapzeros((self.dfgrid.dof, self.NLz, self.aperture.length()))
		self.losvds_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.NLz, self.Nv, self.aperture.length()))
		self.masses_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.NLz, self.aperture.length()))
		self.losvds_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.NLz, self.Nv, self.aperture.length()))
		self.masses_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, self.NLz, self.aperture.length()))
		
		self.aperturestore = mab.gd.gdfast.ApertureStorageLosvd2d(self.aperture, self.losvds_unshaped, self.masses_unshaped, self.vmax)
		
		#self.moments3d = mab.utils.numpy.mmapzeros((self.dfgrid.dof, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_unshaped = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, 1+3+3, self.nr), dtype=numpy.float64)
		#self.moments3d_basis = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self.dfgrid.dof_per_cell, 1+3+3, self.nr), dtype=numpy.float64)
		#self.projectedmoments = mab.utils.numpy.mmapzeros((self.dfgrid.n_cells, self,dfgrid.dof_per_cell, self.moments, constraints), dtype=numpy.float64)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMomentsPolar_IrrReg(self.aperture.aperture_fast(), self.aperture.constraintnrgrid, self.projectedmoments)
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.binned_data.aperture_index_to_constraint, self.projectedmoments)
		#self.aperture_index_to_constraint = self.binned_data.aperture_index_to_constraint.astype(int32)
		#if isinstance(self.aperture, mab.gd.schw.aperture.Polar):
		#	print "polar aperture!, use optimized version"
		#	self.aperturestore = mab.gd.gdfast.ApertureStorageMoments_Polar2d_IrrReg_s(self.aperture.aperture_fast_s(), self.aperture_index_to_constraint, self.projectedmoments)
		#else: 
		
		#self.aperturestore = mab.gd.gdfast.ApertureStorageMoments(self.aperture.aperture_fast(), self.aperture_index_to_constraint, self.projectedmoments_unshaped)
		self.n = 0
	
	def store_test(self, r, vr, vt, orbitnr, u, v):
		if self.rotations > 0:
			self.losvds_unshaped[orbitnr] = 0
			self.masses_unshaped[orbitnr] = 0
			#if any(abs(vz) >= self.vmax):
				#print vz
				#print self.vmax
				#print vz[abs(vz) > self.vmax] 
			self.aperturestore.bin_rotate(r, vr, vt, orbitnr, self.rotations)
			if 1:
				for i in range(self.dfgrid.dof_per_cell):
					weight = self.dfgrid.basis(i, u, v)
					#print u, v
					assert u > 0
					assert v > 0
					assert u < 1
					assert v < 1
					self.losvds_basis[orbitnr,i,:,:] += self.losvds_unshaped[orbitnr] *\
					scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / self.weights[i]
					self.masses_basis[orbitnr,i,:] += self.masses_unshaped[orbitnr] *\
					scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / 
			self.n += len(r) * self.rotations
			
	def store(self, x, y, vz, orbitnr, u, v):
		if self.rotations > 0:
			return
		self.losvds_unshaped[orbitnr] = 0
		self.masses_unshaped[orbitnr] = 0
		#if any(abs(vz) >= self.vmax):
			#print vz
			#print self.vmax
			#print vz[abs(vz) > self.vmax] 
		self.aperturestore.bin(x, y, vz, orbitnr,0)
		if 1:
			for i in range(self.dfgrid.dof_per_cell):
				weight = self.dfgrid.basis(i, u, v)
				#print u, v
				assert u > 0
				assert v > 0
				assert u < 1
				assert v < 1
				self.losvds_basis[orbitnr,i,:,:] += self.losvds_unshaped[orbitnr] *\
				scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / self.weights[i]
				self.masses_basis[orbitnr,i,:] += self.masses_unshaped[orbitnr] *\
				scaling_factor #weight /scaling_factor #* 5 / 4 * 0.95 * 1.094303514738155 / 1.0972532454693291 / 1.18430766017# / 
		self.n += len(x)
		
	def load(self, compressed=True):
		#self.binned_data.load()
		#constraints = self.binned_data.n_constraints
		
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "rb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = bz2.BZ2File(filename_masses, "rb")
		else:
			filename = self.filename
			f = file(filename, "rb")
		self.losvds = numpy.load(f)
		self.masses = numpy.load(f_masses)
		#self.projectedmoments = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 5, constraints)))
		
	def save(self, compressed=True):
		#self.moments3d = array(memmap(self.filename, dtype='float64', mode='readonly', shape=(self.dfgrid.n_orbits, 1+3+3, self.nr)))
		for i1 in range(self.dfgrid.n_I1):
			for i2 in range(self.dfgrid.n_I2):
				orbitnr = self.dfgrid.index_to_orbitnr(i1, i2, 0)
				for k in range(self.dfgrid.dof_per_cell):
					dof_index = self.dfgrid.dof_index(i1, i2, k)
					self.losvds[dof_index] += self.losvds_basis[orbitnr, k]
					self.masses[dof_index] += self.masses_basis[orbitnr, k]
		if compressed:
			filename = self.filename +".bz2"
			f = bz2.BZ2File(filename, "wb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = bz2.BZ2File(filename_masses, "wb")
		else:
			filename = self.filename
			f = file(filename, "wb")
			filename_masses = self.filename_masses +".bz2"
			f_masses = file(filename_masses, "wb")
		logger.info("saving: %s" % filename)
		numpy.save(f, self.losvds)
		logger.info("saving: %s" % filename_masses)
		numpy.save(f_masses, self.masses)
			
	def normalize(self, orbitnr):
		self.losvds_basis[orbitnr] /= self.n
		self.masses_basis[orbitnr] /= self.n
		self.n = 0

