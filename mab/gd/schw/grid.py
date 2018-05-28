# -*- coding: utf-8 -*-
from numpy import *
from mab.gd.gdfast import *
import mab.gd.gdfast
import pyublas
import mab.gd.logging as logging

logger = logging.getLogger("gd.schw.grid")

class Grid2I_Block_Inside(object):
	def __init__(self, light_model, profile_model, logrmin, logrmax, n_I1, n_I2, dither=1, _issubgrid=False):
		self.light_model = light_model
		self.profile_model = profile_model
		self.logrmin = log10(light_model.arcsec_to_kpc(10**logrmin))
		self.logrmax = log10(light_model.arcsec_to_kpc(10**logrmax))
		self.n_I1 = n_I1
		self.n_I2 = n_I2
		self.logrs = (arange(0, n_I1,dtype=float)+0.5) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.logr_borders = (arange(0, n_I1+1,dtype=float)) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.Es = profile_model.potentialr(10**self.logrs)
		self.E_borders = profile_model.potentialr(10**self.logr_borders)
		epsilon = 0.0001
		self.ls = (arange(0, n_I2, dtype=float)+0.5)/n_I2 * (1-2*epsilon) + epsilon
		self.l_borders = (arange(0, n_I2+1, dtype=float))/n_I2 * (1-2*epsilon) + epsilon
		self.dither = dither
		self.n_orbits = self.n_I1 * self.n_I2
		self.dof_per_cell = 1
		self.dof = self.n_orbits
		
		if not _issubgrid:
			self.subgrid = Grid2I_Block_Inside(light_model, profile_model, logrmin, logrmax, n_I1*dither, n_I2*dither, dither=1, _issubgrid=True)
		else:
			self.subgrid = None
			
	def basis(self, i, u, v):
		return 1
			
	def index_to_dof(self, xi, yi, i):
		return i1*self.n_I2 + i2
	
	def solve_coordinates(self, inner_products):
		# system is orthogonal, so inner products are the coordinates themselves
		return inner_products * 1.0
	
	def index_to_orbitnr(self, i1, i2, i3):
		assert i3 == 0
		return i1*self.n_I2 + i2
	
	def __call__(self, coordinates, logr, l):
		return self.mesh.eval(coorinates, logr, l)
	
			
class Grid2I_Shaped_Inside(object):
	def __init__(self, light_model, profile_model, logrmin, logrmax, n_I1, n_I2, dither=1, order=0, _issubgrid=False, circular=False):
		self.light_model = light_model
		self.profile_model = profile_model
		self.umin = self.logrmin = log10(light_model.arcsec_to_kpc(10**logrmin))
		self.umax = self.logrmax = log10(light_model.arcsec_to_kpc(10**logrmax))
		#import pdb
		#pdb.set_trace()
		logger.debug("rmin: %s kpc rmax: %s kpc" % (light_model.arcsec_to_kpc(10**logrmin), light_model.arcsec_to_kpc(10**logrmax)))
		self.n_I1 = n_I1
		self.n_I2 = n_I2
		self.logrs = (arange(0, n_I1,dtype=float)+0.5) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.logr_borders = (arange(0, n_I1+1,dtype=float)) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.rs = 10**self.logrs
		self.r_borders = 10**self.logr_borders
		if circular:
			self.Es = array([profile_model.vcirc(r)**2/2 + profile_model.potentialr(r)  for r in 10**self.logrs])
			self.E_borders = array([profile_model.vcirc(r)**2/2 + profile_model.potentialr(r)  for r in 10**self.logr_borders])
			#import pdb;
			#pdb.set_trace()
		else:
			self.Es = profile_model.potentialr(10**self.logrs)
			self.E_borders = profile_model.potentialr(10**self.logr_borders)
		#print "Es", self.Es
		epsilon = 0.0001
		self.ls = (arange(0, n_I2, dtype=float)+0.5)/n_I2 * (1-2*epsilon) + epsilon
		self.l_borders = (arange(0, n_I2+1, dtype=float))/n_I2 * (1-2*epsilon) + epsilon
		self.dither = dither
		self.n_orbits = self.n_I1 * self.n_I2
		self.n_cells = self.n_orbits
		self.order = order
		
		if not _issubgrid:
			#self.mesh = MeshRegularNodal2dLagrange1(self.logrmin, 0, self.logrmax, 1., self.n_I1, self.n_I2)
			cls = getattr(mab.gd.gdfast, "MeshRegularNodal2dLagrange" + str(self.order)) 
			self.mesh = cls(self.logrmin, 0, self.logrmax, 1., self.n_I1, self.n_I2)
			self.dof = self.mesh.get_dof()
			self.dof_per_cell = self.mesh.get_dof_per_cell()
		else:
			self.mesh = None
		#self.xmin
		self.I1s = self.logrs
		self.I2s = self.ls
		
		
		if not _issubgrid:
			self.subgrid = Grid2I_Shaped_Inside(light_model, profile_model, logrmin, logrmax, n_I1*dither, n_I2*dither, dither=1, _issubgrid=True, circular=circular)
		else:
			self.subgrid = None
			
			
	def get_dof(self):
		return self.mesh.get_dof()
			
	def basis(self, i, u, v):
		return self.mesh.basis_uv(i, u, v)
			
	def dof_index(self, xi, yi, i):
		return self.mesh.dof_index(xi, yi, i)
	
	def solve_coordinates(self, inner_products):
		# system is not orthogonal, solve: M * coordinates = inner_products
		# where M_i,j = <b_i, b_j> (inner product of basis vectors)
		coordinates = inner_products * 0
		self.mesh.solve_coordinates(inner_products, coordinates)
		return coordinates
	
	def __call__(self, coordinates, logr, l):
		return self.mesh.eval(coordinates, logr, l)
	
	def index_to_orbitnr(self, i1, i2, i3):
		assert i3 == 0
		return i1 + i2 *self.n_I1 # + i2

class Grid2d_Shaped_Inside(object):
	def __init__(self, xmin, xmax, ymin, ymax, nx, ny, dither=1, order=0, _issubgrid=False):
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.nx = nx
		self.ny = ny
		
		self.resize = (self.xmin, self.ymin), (self.xmax, self.ymax)
		
		self.x = (arange(0, nx,dtype=float)+0.5) / (nx) * (self.xmax-self.xmin) + self.xmin
		self.x_borders = (arange(0, nx+1,dtype=float)) / (nx) * (self.xmax-self.xmin) + self.xmin
		
		self.y = (arange(0, ny, dtype=float)+0.5) / (ny) * (self.ymax-self.ymin) + self.ymin
		self.y_borders = (arange(0, ny+1, dtype=float)) / (nx) * (self.ymax-self.ymin) + self.ymin
		
		self.dither = dither
		self.n = self.nx * self.ny
		self.n_cells = self.n
		self.order = order
		
		if not _issubgrid:
			#self.mesh = MeshRegularNodal2dLagrange1(self.logrmin, 0, self.logrmax, 1., self.n_I1, self.n_I2)
			cls = getattr(mab.gd.gdfast, "MeshRegularNodal2dLagrange" + str(self.order)) 
			self.mesh = cls(self.xmin, self.ymin, self.xmax, self.ymax, self.nx, self.ny)
			self.dof = self.mesh.get_dof()
			self.dof_per_cell = self.mesh.get_dof_per_cell()
		else:
			self.mesh = None
		
		if not _issubgrid:
			self.subgrid = Grid2d_Shaped_Inside(xmin, xmax, ymin, ymax, nx*dither, ny*dither, dither=1, _issubgrid=True)
		else:
			self.subgrid = None
			
			
	def get_dof(self):
		return self.mesh.get_dof()
			
	def basis(self, i, u, v):
		return self.mesh.basis_uv(i, u, v)
			
	def dof_index(self, xi, yi, i):
		return self.mesh.dof_index(xi, yi, i)
	
	def solve_coordinates(self, inner_products):
		# system is not orthogonal, solve: M * coordinates = inner_products
		# where M_i,j = <b_i, b_j> (inner product of basis vectors)
		coordinates = inner_products * 0
		self.mesh.solve_coordinates(inner_products, coordinates)
		return coordinates
	
	def __call__(self, coordinates, logr, l):
		return self.mesh.eval(coordinates, logr, l)
	
	#def index_to_orbitnr(self, i1, i2, i3):
	#	assert i3 == 0
	#	return i1 + i2 *self.n_I1 # + i2


class Grid2I_Shaped_InsideTest(object):
	def __init__(self, light_model, profile_model, rmin, rmax, n_I1, n_I2, dither=1, order=0, _issubgrid=False):
		self.light_model = light_model
		self.profile_model = profile_model
		#self.logrmin = log10(light_model.arcsec_to_kpc(10**logrmin))
		#self.logrmax = log10(light_model.arcsec_to_kpc(10**logrmax))
		scale = 2
		self.umin = arctan(rmin*scale)
		self.umax = arctan(rmax*scale)
		self.n_I1 = n_I1
		self.n_I2 = n_I2
		#self.logrs = (arange(0, n_I1,dtype=float)+0.5) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		#self.logr_borders = (arange(0, n_I1+1,dtype=float)) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		us = (arange(0, n_I1,dtype=float)+0.5) / (n_I1) * (self.umax-self.umin) + self.umin
		uborders = (arange(0, n_I1+1,dtype=float)) / (n_I1) * (self.umax-self.umin) + self.umin
		self.rs = rs = tan(us)/scale
		self.r_borders = r_borders = tan(uborders)/scale 
		self.Es = profile_model.potentialr(rs)
		self.E_borders = profile_model.potentialr(r_borders)
		epsilon = 0.0001
		self.ls = (arange(0, n_I2, dtype=float)+0.5)/n_I2 * (1-2*epsilon) + epsilon
		self.l_borders = (arange(0, n_I2+1, dtype=float))/n_I2 * (1-2*epsilon) + epsilon
		self.dither = dither
		self.n_orbits = self.n_I1 * self.n_I2
		self.n_cells = self.n_orbits
		self.order = order
		
		if not _issubgrid:
			#self.mesh = MeshRegularNodal2dLagrange1(self.logrmin, 0, self.logrmax, 1., self.n_I1, self.n_I2)
			cls = getattr(mab.gd.gdfast, "MeshRegularNodal2dLagrange" + str(self.order)) 
			self.mesh = cls(self.umin, 0, self.umax, 1., self.n_I1, self.n_I2)
			self.dof = self.mesh.get_dof()
			self.dof_per_cell = self.mesh.get_dof_per_cell()
		else:
			self.mesh = None
		#self.xmin
		#self.I1s = self.logrs
		self.I2s = self.ls
		
		
		if not _issubgrid:
			self.subgrid = Grid2I_Shaped_InsideTest(light_model, profile_model, rmin, rmax, n_I1*dither, n_I2*dither, dither=1, _issubgrid=True)
		else:
			self.subgrid = None
			
	def basis(self, i, u, v):
		return self.mesh.basis_uv(i, u, v)
			
	def dof_index(self, xi, yi, i):
		return self.mesh.dof_index(xi, yi, i)
	
	def solve_coordinates(self, inner_products):
		# system is not orthogonal, solve: M * coordinates = inner_products
		# where M_i,j = <b_i, b_j> (inner product of basis vectors)
		coordinates = inner_products * 0
		self.mesh.solve_coordinates(inner_products, coordinates)
		return coordinates
	
	def __call__(self, coordinates, logr, l):
		return self.mesh.eval(coordinates, logr, l)
	
	def index_to_orbitnr(self, i1, i2, i3):
		assert i3 == 0
		return i1 + i2 *self.n_I1 # + i2
	
		
		
class Grid3I_Block_Inside(object):
	def __init__(self, light_model, profile_model, logrmin, logrmax, n_I1, n_I2, n_I3, ditherI1=1, ditherI2=1, ditherI3=1, _issubgrid=False):
		self.light_model = light_model
		self.profile_model = profile_model
		self.logrmin = log10(light_model.arcsec_to_kpc(10**logrmin))
		self.logrmax = log10(light_model.arcsec_to_kpc(10**logrmax))
		self.n_I1 = n_I1
		self.n_I2 = n_I2
		self.n_I3 = n_I3
		self.logrs = (arange(0, n_I1,dtype=float)+0.5) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.logr_borders = (arange(0, n_I1+1,dtype=float)) / (n_I1) * (self.logrmax-self.logrmin) + self.logrmin
		self.Es = profile_model.potentialr(10**self.logrs)
		self.E_borders = profile_model.potentialr(10**self.logr_borders)
		epsilon = 0.0001
		self.lzs = (arange(0, n_I2, dtype=float)+0.5)/n_I2 * (1-2*epsilon) + epsilon
		self.lz_borders = (arange(0, n_I2+1, dtype=float))/n_I2 * (1-2*epsilon) + epsilon
		self.i3s = (arange(0, n_I3, dtype=float)+0.5)/n_I3 * (1-2*epsilon) + epsilon
		self.i3_borders = (arange(0, n_32+1, dtype=float))/n_I3 * (1-2*epsilon) + epsilon
		self.ditherI1 = ditherI1
		self.ditherI2 = ditherI2
		self.ditherI3 = ditherI3
		self.n_orbits = self.n_I1 * self.n_I2 * self.n_I3
		self.dof_per_cell = 1
		self.dof = self.n_orbits
		
		if not _issubgrid:
			self.subgrid = Grid3I_Block_Inside(light_model, profile_model, logrmin, logrmax, n_I1*ditherI1, n_I2*ditherI2, n_I3*ditherI3, ditherI1=1, ditherI2=1, ditherI3=1, _issubgrid=True)
		else:
			self.subgrid = None
			
	def basis(self, i, u, v, w):
		return 1
			
	def index_to_dof(self, i1, i2, i3, i):
		return (i1*self.n_I2 + i2) * self.n_I3 + i3
	
	def solve_coordinates(self, inner_products):
		# system is orthogonal, so inner products are the coordinates themselves
		return inner_products * 1.0
	
	def index_to_orbitnr(self, i1, i2, i3):
		#assert i3 == 0
		#$return i1*self.n_I2 + i2
		return (i1*self.n_I2 + i2) * self.n_I3 + i3
	
	def __call__(self, coordinates, logr, l):
		return self.mesh.eval(coorinates, logr, l)
	
