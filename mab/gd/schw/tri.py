# -*- coding: utf-8 -*-
import mab.gd.logging as logging
from numpy import *
import sys
import numpy
import mab.gd.gdfast
from kaplot import *

logger = logging.getLogger("gd.schw.tri")


class PlotOutput(object):
	def __init__(self, output):
		self.output = output
		
	def run(self, args, opts, scope):
		self.output.load()
		
class Output(object):
	def __init__(self, filename):
		self.filename = filename
		
	def load(self):
		self.orbitreader = mab.gd.gdfast.Orblib(self.filename)
		filename = self.filename +"-densities.npy"
		if 0:
			logger.info("reading: %s" % self.filename)
			self.orbitreader.read()
			grid = zeros((self.orbitreader.noI1*self.orbitreader.noI2*self.orbitreader.noI3, self.orbitreader.noConstraints, self.orbitreader.noMaxVelocityHistograms))
			self.orbitreader.fillorbit(grid)
			densities = grid.sum(axis=2)
			print densities.shape
			numpy.save(filename, densities)
			print densities.sum()
			import pdb;
			pdb.set_trace()
		else:
			self.orbitreader.readInfoOnly()
			densities = numpy.load(filename)
			#import pdb;
			#pdb.set_trace()
			print densities.nbytes
			Nx, Ny = self.orbitreader.noI2, self.orbitreader.noI3
			print Nx, Ny
			#dsadsa
			i = 8
			mozaic(Nx, Ny, box)
			for j in range(Nx):
				for k in range(Ny):
					select(j, k)
					#orbitnr = i + Nx * j + Nx*Ny*k
					orbitnr = i*Nx*Ny + j * Ny + k
					density = densities[orbitnr]
					density = density.reshape((40, 40))
					colormap = "whiterainbow"
					indexedimage(density, colormap=colormap)
			draw()
		
		
	def load_(self):
		f = open(self.filename, "rb")
		
		print numpy.fromfile(f, numpy.dtype('u4'), count=20)
		f.seek(0)
		orbits, nE1, nI2 ,nI3, ndithers = self._read(f, '(5,)u4')
		print orbits, nE1, nI2, nI3, ndithers
		
		quadrant_light, quad_lph, quad_lth, quad_lr = self._read(f, '(4,)u4')
		
		borders_radii = self._read(f, ('(%d,)f8' % (quad_lr+1)))
		borders_theta = self._read(f, ('(%d,)f8' % (quad_lth+1)))
		borders_phi = self._read(f, ('(%d,)f8' % (quad_lph+1)))
		print "borders_radii =", borders_radii
		print "borders_theta =", borders_theta
		print "borders_phi =", borders_phi
		
		
		#(i4, i4, double)
		#h_nconstr, t1, hist_basic(1,1)/hist_basic(1,3)
		ncontraints, histsizehalf, delta_v =  self._read(f, 'i4, i4, f8')
		print "ncontraints, histsizehalf, delta_v", ncontraints, histsizehalf, delta_v
		#numpy.fromfile(f, numpy.dtype('u4'), count=1
		
		for orbit in range(orbits):
			print "orbit", orbit
			orbit, E1,I2,I3, totalnotregularizable = self._read(f, '(5,)u4')
			print self._read(f, '(%d,)i4' % (ndithers**3))
			count = 16 * quad_lr * quad_lth * quad_lph
			# t(16,quad_nph,quad_nth,quad_nr
			intrinsic = self._read(f, '(%d, %d, %d, %d)f8' % (quad_lr, quad_lth, quad_lph, 16))
			
			for contraint in range(ncontraints):
				begin, end = self._read(f, 'i4, i4')
				print begin, end
				if begin <= end:
					print "read it"
					sys.exit(0)
				else:
					#print "skip"
					pass
		
		f.seek(0)
		print numpy.fromfile(f, numpy.dtype('u4'), count=20)
		import pdb
		pdb.set_trace()
		
		
	def _read(self, f, dataformat):
		pos = f.tell()
		prefix = numpy.fromfile(f, numpy.dtype('u4'), count=1)
		#type, N = dataformat
		#dataformat = "(%d,)%s" % (N, type)
		data = numpy.fromfile(f, numpy.dtype(dataformat), count=1)[0]
		postfix = numpy.fromfile(f, numpy.dtype('u4'), count=1)
		#prefix, data, postfix = d[0]
		if prefix != postfix:
			logger.error("unexpected data in fortran file: postfix(%d) != prefix(%d), dataformat = %s" % (prefix, postfix, dataformat))
			df = numpy.dtype(dataformat)
			
			import pdb
			pdb.set_trace()
			
			N = 10
			f.seek(pos)
			print numpy.fromfile(f, numpy.dtype('u4'), count=N)
			sys.exit(-1)
		return data
		


class Aperture(object):
	def __init__(self, rmax_kpc, N, light_model, output, angle=0):
		self.rmax_kpc = rmax_kpc
		self.Nx = self.Ny = self.N = N
		self.light_model = light_model
		self.output = output
		self.angle = angle
		self.rmax_arcsec = self.light_model.kpc_to_arcsec(self.rmax_kpc)
		self.x0 = self.rmax_arcsec/2
		self.y0 = self.rmax_arcsec/2
		self.width = self.rmax_arcsec
		self.height = self.rmax_arcsec
		
		
	def run(self, args, opts, scope):
		logger.info("writing to file: %s" % self.output)
		print "rmax_kpc", self.rmax_kpc
		rmax_arcsec = self.rmax_arcsec #self.light_model.kpc_to_arcsec(self.rmax_kpc)
		print "rmax_arcsec", rmax_arcsec
		
		data = ""
		data += "#counter_rotation_boxed_aperturefile_version_2\n"
		data += "%f %f\n" % (-rmax_arcsec/2, -rmax_arcsec/2)
		data += "%f %f\n" % (rmax_arcsec, rmax_arcsec)
		data += "%f\n" % (self.angle)
		data += "%d %d\n" % (self.N, self.N)
		file(self.output, "w").write(data)
		
class Bins(object):
	def __init__(self, aperture, output):
		self.aperture = aperture
		self.output = output
		
	def run(self, args, opts, scope):
		logger.info("writing to file: %s" % self.output)
		data = ""
		data += "#Counterrotaton_binning_version_1\n"
		size = self.aperture.N**2
		data += "%d\n" % (size)
		for i in range(size):
			data += "%d " % (i+1)
		file(self.output, "w").write(data)
		
	def xy_to_index(self, x, y):
		xindex = (x - self.aperture.x0)/self.aperture.width * self.aperture.Nx
		yindex = (y - self.aperture.y0)/self.aperture.height * self.aperture.Ny
		index = xindex + yindex * self.aperture.Nx
		
		
class DF(object):
	def __init__(self, n_E, n_I2, n_I3, n_dither, logrmin, logrmax):
		self.n_E = n_E
		self.n_I2 = n_I2
		self.n_I3 = n_I3
		self.n_dither = n_dither
		self.logrmin = logrmin
		self.logrmax = logrmax

class MGE(object):
	def __init__(self, values):
		self.values = values
		
		
class TriProfileNone(object):
	def __init__(self):
		self.type = 0
		self.extra = "\n"
		
kparsec_in_km = 3.08568025e16 
class TriProfileNFW(object):
	def __init__(self, profile, a, b, c):
		self.type = 3
		#dm_profile_rhoc, dm_profile_parameter, r200, dm_rho_crit, dm_profile_a, dm_profile_b, dm_profile_c
		self.rho0 = profile.rho0
		self.rs = profile.scale
		self.r200 = profile.r200
		self.extra = "%e %e %e %e %f %f %f\n" % (self.rho0 * kparsec_in_km**-3, self.rs*kparsec_in_km, self.r200*kparsec_in_km, profile.rhoc*kparsec_in_km**-3, a, b, c)
		
		r200 = self.r200*kparsec_in_km
		
		dm_profile_rhoc = self.rho0 * kparsec_in_km**-3
		dm_profile_parameter = self.rs*kparsec_in_km
		M = 4.0 * pi * dm_profile_rhoc * dm_profile_parameter**3 * (log(1+r200/dm_profile_parameter) - r200/dm_profile_parameter/(1+r200/dm_profile_parameter))
		print "mass = %g " % M
		print "m200 = %g " % profile.M200
		#sys.exit(0)
		
		
		
class Parameters(object):
	def __init__(self, mge, distance_Mpc, viewing_angles, rlogmin3d, rlogmax3d, output, df, profile_tri, M_over_L=1.):
		self.mge = mge
		self.distance_Mpc = distance_Mpc
		self.viewing_angles = viewing_angles
		self.M_over_L = M_over_L
		self.df = df
		self.rlogmin3d = rlogmin3d
		self.rlogmax3d = rlogmax3d
		self.output = output
		self.profile_tri = profile_tri
		
	def run(self, args, opts, scope):
		logger.info("writing to file: %s" % self.output)
		data = ""
		data += "%i 1\n" % len(self.mge.values)
		for value in self.mge.values:
			data += "%f %f %f %f\n" % value
		data += "%f\n" % self.distance_Mpc
		data += "%f %f %f\n" % self.viewing_angles
		data += "%f\n" % self.M_over_L
		data += "%f\n" % 0 # BH
		data += "%f\n" % 0 # BH
		data += "%f %f\n" % (self.rlogmin3d, self.rlogmax3d)
		data += "%d %f %f\n" % (self.df.n_E, self.df.logrmin, self.df.logrmax)
		data += "%d\n" % (self.df.n_I2)
		data += "%d\n" % (self.df.n_I3)
		data += "%d\n" % (self.df.n_dither)
		data += "%d\n" % (self.profile_tri.type)
		data += "%s\n" % (self.profile_tri.extra)
		f = open(self.output, "w")
		f.write(data)
		f.close() 
		