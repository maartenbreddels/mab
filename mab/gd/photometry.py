# -*- coding: utf-8 -*-
from numpy import *
import numpy
import os
import mab.gd.logging as logging
logger = logging.getLogger("gd.photometry")
import mab.astrometry


class FakePhotometry(object):
	def __init__(self, light_model, axis_ratio, N, aperture, dirname, name="photometry"):
		self.light_model = light_model
		self.axis_ratio = axis_ratio
		self.N = N
		self.aperture = aperture
		self.dirname = dirname
		self.name = name
		self.filename_xy = os.path.join(self.dirname, self.name + "_xy.npy")
		self.filename_grid = os.path.join(self.dirname, self.name + "_grid.npy")
		
	def run(self, args, opts, scope):
		self.x = []
		self.y = []
		logger.info("aperture: %dx%d (size=%d)" % (self.aperture.lengthx(), self.aperture.lengthy(), self.aperture.length()))
		self.grid = zeros(self.aperture.length(), dtype=float)
		count = 0
		light_profile = self.light_model.light_profile
		while count < self.N:
			r = light_profile.sample_r(N=1)[0]
			r = self.light_model.kpc_to_arcsec(r)
			costheta = numpy.random.random() * 2 - 1
			phi = numpy.random.random() * 2 * pi
			eta = numpy.random.random() * 2 * pi
			theta = numpy.arccos(costheta)
			#sintheta = numpy.sqrt(1-costheta**2)
			sintheta = numpy.sin(theta)
			x = r * sintheta * numpy.cos(phi)
			y = r * sintheta * numpy.sin(phi)
			R = sqrt(x**2+y**2)
			
			angle = random.random() * 2 * pi
			x = R * cos(angle)
			y = R * sin(angle)
			xp = x
			yp = y * self.axis_ratio
			if self.aperture.inrange(xp, yp):
				self.x.append(xp)
				self.y.append(yp)
				#print r, self.aperture.findindex(xp, yp)
				self.grid[self.aperture.findindex(xp, yp)] += 1
				count += 1
		self.x = array(self.x)
		self.y = array(self.y)
		self.grid2d = self.grid.reshape((self.aperture.lengthx(), self.aperture.lengthy()))
		self.save()
		
		
	def save(self):
		save(self.filename_xy, array([self.x, self.y]))
		save(self.filename_grid, self.grid)
		logger.info("saving photometry grid: %s (count=%d)" % (self.filename_grid, self.grid.sum()))
		
	def load(self):
		x, y = load(self.filename_xy)
		self.grid = load(self.filename_grid)
		self.grid2d = self.grid.reshape((self.aperture.lengthx(), self.aperture.lengthy()))
		logger.info("loading photometry grid: %s" % self.filename_grid)
		
class Photometry(object):
	def __init__(self, filename_input, aperture, dirname, name="photometry"):
		self.filename_input = filename_input
		self.aperture = aperture
		self.dirname = dirname
		self.name = name
		self.filename_xy = os.path.join(self.dirname, self.name + "_xy.npy")
		self.filename_grid = os.path.join(self.dirname, self.name + "_grid.npy")
		
	def read_input(self):
		pass
		
	def run(self, args, opts, scope):
		logger.info("aperture: %dx%d (size=%d)" % (self.aperture.lengthx(), self.aperture.lengthy(), self.aperture.length()))
		self.grid = zeros(self.aperture.length(), dtype=float)
		self.read_input()
		count = 0
		for x, y, in zip(self.x, self.y):
			xp, yp = x, y
			if self.aperture.inrange(xp, yp):
				self.x.append(xp)
				self.y.append(yp)
				#print r, self.aperture.findindex(xp, yp)
				self.grid[self.aperture.findindex(xp, yp)] += 1
				count += 1
		self.x = array(self.x)
		self.y = array(self.y)
		self.grid2d = self.grid.reshape((self.aperture.lengthx(), self.aperture.lengthy()))
		self.save()
		
		
	def save(self):
		save(self.filename_xy, array([self.x, self.y]))
		save(self.filename_grid, self.grid)
		logger.info("saving photometry grid: %s (count=%d)" % (self.filename_grid, self.grid.sum()))
		
	def load(self):
		x, y = load(self.filename_xy)
		self.grid = load(self.filename_grid)
		self.grid2d = self.grid.reshape((self.aperture.lengthx(), self.aperture.lengthy()))
		logger.info("loading photometry grid: %s" % self.filename_grid)
		
class PhotometryGius(Photometry):
	def __init__(self, ra0, dec0, **kwargs):
		super(PhotometryGius, self).__init__(**kwargs)
		self.ra0 = ra0
		self.dec0 = dec0
		
	def read_input(self):
		logger.info("reading photometry: %s" % self.filename_input)
		lines = file(self.filename_input).readlines()
		print len(lines)
		count = 0
		i = 0
		self.x = []
		self.y = []
		for line in lines:
			line = line.strip()
			#print line
			values = [k.strip() for k in line.split()]
			ra = str(values[0]) +" %02i" % int(values[1]) +" %06.3f" % float(values[2])
			dec = str(values[3]) +" %02i" % int(values[4]) +" %06.3f" % float(values[5])
			ra = mab.astrometry.read_ra_hour(ra)
			dec = mab.astrometry.read_dec_deg(dec)
			#vlos_helio = float(values[7])
			#e_vlos = float(values[8])
			#x = float(values[-3])
			#y = float(values[-2])
			flagI, flagV = int(values[-3]), int(values[-8])
			V, I = float(values[-5]), float(values[-10])
			
			#print flagI, flagV, values
			if (flagI in [-1,-2,-3]) and (flagV in [-1,-2,-3]) and (V>16) and (V<21.8) and (I<21.1):
				xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
				if 0:
					xip = xi * cos(angle) - eta * sin(angle)
					etap = xi * sin(angle) + eta * cos(angle)
				#PA = -(self.PA) # hmm, not sure i fully understand
				#print xi, eta
				self.x.append(xi*60*60*1.4)
				self.y.append(eta*60*60*1.4)
			
			
				