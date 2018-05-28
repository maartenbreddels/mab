# -*- coding: utf-8 -*-
import mab
import numpy
from kaplot import *
from numpy import *
import mab.gd.logging as logging

logger = logging.getLogger("gd.schw.tri")

class LightModel(object):
	def __init__(self, distance, ellipticity, PA):
		self.distance = distance
		self.ellipticity = ellipticity
		self.PA = PA
		
	def kpc_to_arcsec(self, d): 
		return d * 1.0/(self.distance) / (pi/180/60/60)
	
	#def pc_to_arcsec(self, d): 
	#	return d * 1e3/(self.distance) / (pi/180/60/60)
	
	def arcsec_to_kpc(self, R): 
		return R / (1.0/(self.distance) / (pi/180/60/60))
		
		

class ApertureNative(object):
	def __init__(self, light_model, Nx, Ny, xmax_kpc, ymax_kpc, filename_aperture_descriptor):
		self.light_model = light_model
		self.Nx = Nx
		self.Ny = Ny
		self.xmax_kpc = xmax_kpc
		self.ymax_kpc = ymax_kpc
		self.filename_aperture_descriptor = filename_aperture_descriptor
		self.load()
		
	def load(self):
		self.aperture = zeros((self.Nx, self.Ny))
		index = 0
		for x in range(self.Nx):
			for y in range(self.Ny):
				self.aperture[x,y] = index
				index += 1
		x, x_borders = mab.utils.numpy.range_and_borders(self.Nx, -self.xmax_kpc, self.xmax_kpc)
		y, y_borders = mab.utils.numpy.range_and_borders(self.Nx, -self.ymax_kpc, self.ymax_kpc)
		x, y = numpy.meshgrid(x, y)
		self.xgrid, self.ygrid = x, y
				
	def run(self, args, opts, scope):
		self.load()
		logger.info("generating aperture descriptor: %s" % self.filename_aperture_descriptor)
		path = os.path.dirname(self.filename_aperture_descriptor)
		if not os.path.exists(path):
			os.makedirs(path)
		f = open(self.filename_aperture_descriptor, "w")
		print >>f, "#counter_rotation_boxed_aperturefile_version_2"
		print >>f, -self.light_model.kpc_to_arcsec(self.xmax_kpc), -self.light_model.kpc_to_arcsec(self.ymax_kpc)
		print >>f, self.light_model.kpc_to_arcsec(self.xmax_kpc*2), self.light_model.kpc_to_arcsec(self.ymax_kpc*2)
		print >>f, 0. # angle
		print >>f, self.Nx, self.Ny
		f.close()
				
		
class AperturePolar(object):
	def __init__(self, light_model, Nx, Ny, xmax_kpc, ymax_kpc, Nphi, Nr, filename_aperture_descriptor, filename_aperture):
		self.light_model = light_model
		self.Nx = Nx
		self.Ny = Ny
		self.xmax_kpc = xmax_kpc
		self.ymax_kpc = ymax_kpc
		self.Nr = Nr
		self.Nphi = Nphi
		self.filename_aperture_descriptor = filename_aperture_descriptor
		self.filename_aperture = filename_aperture
		self.load()
	
	def load(self):
		#dx = self.xmax_kpc*2/float(self.Nx)
		#dy = self.ymax_kpc*2/float(self.Ny)
		#x, y = mgrid[-self.xmax_kpc:self.xmax_kpc:dx, -self.ymax_kpc:self.ymax_kpc:dy]
		#x += dx/2 # convert left border to center
		#y += dy/2
		
		x, x_borders = mab.utils.numpy.range_and_borders(self.Nx, -self.xmax_kpc, self.xmax_kpc)
		y, y_borders = mab.utils.numpy.range_and_borders(self.Nx, -self.ymax_kpc, self.ymax_kpc)
		x, y = numpy.meshgrid(x, y)
		self.xgrid, self.ygrid = x, y
		print x
		print y
		r = sqrt(x**2+y**2)
		phi = (arctan2(y,x) + pi*2) % (2*pi)
		self.aperture = zeros_like(x) + 1
		for Ri in range(self.Nr):
			for Phii in range(self.Nphi):
				aperture_index = 2 + Ri * self.Nphi + Phii
				r1 = (self.xmax_kpc * 1.0 / self.Nr) * Ri
				r2 = (self.xmax_kpc * 1.0 / self.Nr) * (Ri+1.)
				phi1 = 2*pi/self.Nphi * Phii
				phi2 = 2*pi/self.Nphi * (Phii+1)
				mask = (r >= r1) & (r<r2) & (phi > phi1) & (phi <= phi2)
				self.aperture[mask] = aperture_index
		
	def run(self, args, opts, scope):
		self.load()
		logger.info("generating aperture descriptor: %s" % self.filename_aperture_descriptor)
		path = os.path.dirname(self.filename_aperture_descriptor)
		if not os.path.exists(path):
			os.makedirs(path)
		f = open(self.filename_aperture_descriptor, "w")
		print >>f, "#counter_rotation_boxed_aperturefile_version_2"
		print >>f, -self.light_model.kpc_to_arcsec(self.xmax_kpc), -self.light_model.kpc_to_arcsec(self.ymax_kpc)
		print >>f, self.light_model.kpc_to_arcsec(self.xmax_kpc*2), self.light_model.kpc_to_arcsec(self.ymax_kpc*2)
		print >>f, 0. # angle
		print >>f, self.Nx, self.Ny
		f.close()
		
		logger.info("generating aperture bins: %s" % self.filename_aperture)
		path = os.path.dirname(self.filename_aperture)
		if not os.path.exists(path):
			os.makedirs(path)
		f = open(self.filename_aperture, "w")
		print >>f, "#Counterrotaton_binning_version_1"
		print >>f, self.Nx*self.Ny,
		apertureflat = self.aperture.flat
		for i in range(len(apertureflat)):
			if (i % 10) == 0:
				print >>f, "\n",
			print >>f, "% 8d" % int(apertureflat[i]),
		f.close()
		
		
		
			
		#print aperture
		box()
		indexedimage(self.aperture)
		draw()
		
			
			