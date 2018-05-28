# -*- coding: utf-8 -*-
import mab.gd.gdfast
import os
import numpy
from numpy import *
import os
import pyublas
import mab.gd.logging as logging

logger = logging.getLogger("gd.schw.aperture")

class AperturePrepare(object):
	def __init__(self, aperture):
		self.aperture = aperture
		
	def run(self, args, opts, scope):
		self.aperture.init()
		self.aperture.process()
		self.aperture.save()

class AperturePolarBase(object):
	def __init__(self, modelpath, nphi, postfix, light_model, ellipticity, position_angle):
		self.modelpath = modelpath
		self.nphi = nphi
		self.postfix = postfix
		self.light_model = light_model
		self.ellipticity = ellipticity
		self.position_angle = position_angle
		def inrangefilter(stars):
			logger.info("checking if %d stars are in range of aperture" % len(stars))
			stars_inrange = []
			for i in range(len(stars)):
				star = stars[i]
				#print star
				if self.inrange(stars.xi[i], stars.eta[i]):
					stars_inrange.append(star)
			use_numpy = not isinstance(stars, mab.cvsfile.CsvObject)
			if use_numpy:
				stars_inrange = array(stars_inrange, dtype=stars.dtype).view(recarray)
			else:
				stars_inrange = mab.cvsfile.CsvObject(stars_inrange)
			logger.info("%d stars are in range of aperture" % len(stars_inrange))
			return stars_inrange
		
		self.inrangefilter = inrangefilter
		
		
	def aperture_filter(self, index):
		def index_filter(star, index=index):
			return self.aperture.findindex(star.xi, star.eta) == index
		return index_filter

		
	def init(self):
		pass
	
	def process(self):
		self.process_single()
		
	def process_single(self, index=None):
		self.rborders = self.calculate_bins()
		#self.rborders = self.rborders.reshape(self.rborders, (1, len(self.rborders)))
	
	def save(self):
		self.save_single()
		
	def save_single(self, index=None):
		dirname = os.path.join(self.modelpath, "schw")
		if not os.path.exists(dirname):
			os.makedirs(dirname)
		dirname = os.path.join(dirname, "aperture")
		if not os.path.exists(dirname):
			os.makedirs(dirname)
		rborders = self.rborders
		filename = os.path.join(dirname, "aperture_rborders" +self.postfix +".npy")
		logger.info("writing aperture borders to file: %s " % filename)
		numpy.save(filename, rborders)
		logger.info("borders arcsec: %r " % rborders)
		aperture_rborders_kpc = self.light_model.arcsec_to_kpc(rborders)
		logger.info("borders kpc   : %r " % aperture_rborders_kpc)
		
		filename = os.path.join(dirname, "aperture_rborders" +self.postfix +".asc")
		logger.info("writing aperture borders to file: %s (ascii version)" % filename)
		f = open(filename, "w")
		f.write(repr(list(rborders)))
		f.close()
		
		radial_surface_densities = []
		for i in range(len(rborders)-1):
			R1 = rborders[i]
			R2 = rborders[i+1]
			r1 = self.light_model.arcsec_to_kpc(R1)
			r2 = self.light_model.arcsec_to_kpc(R2)
			radial_surface_densities.append(self.light_model.light_profile.cumdensityR(r1, r2, M=1.))
			
		radial_surface_densities = array(radial_surface_densities)
		logger.info("total fractional surface density in aperture: %.3f " % sum(radial_surface_densities)) 
		
		filename = os.path.join(dirname, "radial_surface_densities" +self.postfix +".npy")
		logger.info("writing radial surface densities to file: %s" % filename)
		numpy.save(filename, radial_surface_densities)
		
		filename = os.path.join(dirname, "radial_surface_densities" +self.postfix +".asc")
		logger.info("writing radial surface densities to file: %s (ascii version)" % filename)
		f = open(filename, "w")
		f.write(repr(list(radial_surface_densities)))
		f.close()
		
	def load(self):
		dirname = os.path.join(self.modelpath, "schw", "aperture")
		if 0:
			filename = os.path.join(dirname, "constraintnrgrid" +self.postfix +".npy")
			self.constraintnrgrid = numpy.load(filename)
			assert self.constraintnrgrid.shape[1] == self.nphi, "grid on disk is not of same shape (%r vs %r) as given by parameters, please recreate aperture" % (self.constraintnrgrid.shape[1], self.nphi)
		
		filename = os.path.join(dirname, "aperture_rborders" +self.postfix +".npy")
		self.aperture_rborders = numpy.load(filename)# * 60 * 60 # deg to arcsec
		self.aperture_rcenters = (self.aperture_rborders[0:-1] + self.aperture_rborders[1:])/2
		self.aperture_rborders_kpc = self.light_model.arcsec_to_kpc(self.aperture_rborders)
		self.aperture_rcenters_kpc = self.light_model.arcsec_to_kpc(self.aperture_rcenters)
		
		filename = os.path.join(dirname, "radial_surface_densities" +self.postfix +".npy")
		logger.debug("loading radial surface densities from file: %s" % filename)
		self.radial_surface_densities = numpy.load(filename)
		
		self.gridr = mab.gd.gdfast.Grid1dIrregular(self.aperture_rborders)
		self.gridphi = mab.gd.gdfast.Grid1dRegular(0, 2*pi, self.nphi)
		self.aperture = mab.gd.gdfast.AperturePolar2d(self.gridr, self.gridphi, self.ellipticity, self.position_angle)
		
		self.gridr_s = mab.gd.gdfast.Grid1dIrregular_s(self.aperture_rborders)
		self.gridphi_s = mab.gd.gdfast.Grid1dRegular_s(0, 2*pi, self.nphi)
		self.aperture_s = mab.gd.gdfast.AperturePolar2d_IrrReg_s(self.gridr_s, self.gridphi_s, self.ellipticity, self.position_angle)
		
	def getRmax(self):
		return self.aperture_rborders_kpc[-1]
	
	def __len__(self):
		return self.length()
	
	def length(self):
		return self.aperture.length()
	
	def aperture_fast(self):
		return self.aperture
	
	def aperture_fast_s(self):
		return self.aperture_s
	
	def n_constraints(self):
		return self.constraints
	
	def mass(self, aperture_index):
		rindex = aperture_index / self.gridphi.length()
		return self.radial_surface_densities[rindex]/self.nphi
	
	def stars_to_indices(self, stars):
		aperture_indices = []
		for i in range(len(stars)):
			aperture_indices.append(self.aperture.findindex(stars.xi[i], stars.eta[i]))
		aperture_indices = array(aperture_indices)
		return stars, aperture_indices
	
		


class AperturePolarEquidistantBins(AperturePolarBase):
	def __init__(self, modelpath, nphi, Nr, postfix, light_model, rmin, rmax, ellipticity, position_angle):
		AperturePolarBase.__init__(self, modelpath, nphi, postfix, light_model, ellipticity, position_angle)
		self.Nr = Nr
		self.rmin = rmin
		self.rmax = rmax
		
	def calculate_bins(self):
		rborders_kpc = arange(self.Nr+1.) / (self.Nr+0.) * (self.rmax - self.rmin) + self.rmin
		rborders = self.light_model.kpc_to_arcsec(rborders_kpc)
		print rborders_kpc 
		print rborders 
		return rborders
	 
	def inrange(self, xi, eta):
		#r = sqrt(xi**2+eta**2)
		r = mab.astrometry.xi_eta_to_re(xi, eta, degrees(self.position_angle), self.ellipticity)
		return r < self.light_model.kpc_to_arcsec(self.rmax)
		


class AperturePolarList(AperturePolarBase):
	def __init__(self, modelpath, borders_kpc, nphi, postfix, light_model, ellipticity, position_angle):
		AperturePolarBase.__init__(self, modelpath, nphi, postfix, light_model, ellipticity, position_angle)
		self.borders_kpc = borders_kpc
		#self.rmin = rmin
		#self.rmax = rmax
		
	def calculate_bins(self):
		#rborders_kpc = arange(self.Nr+1.) / (self.Nr+0.) * (self.rmax - self.rmin) + self.rmin
		rborders_kpc = array(self.borders_kpc)
		print "kpc", self.borders_kpc
		rborders = self.light_model.kpc_to_arcsec(rborders_kpc)
		print "arcsec", rborders
		#print rborders_kpc 
		#print rborders 
		return rborders
	 
	#def inrange(self, xi, eta):
	#	r = sqrt(xi**2+eta**2)
	#	return r < self.light_model.kpc_to_arcsec(self.rmax)
		
class AperturePolarEqualStarcountBins(AperturePolarBase):
	def __init__(self, modelpath, nphi, nrperbin, postfix, light_model, ellipticity, position_angle, observation, filter=None, filters=None, rmin_kpc=None, startcounts=[]):
		AperturePolarBase.__init__(self, modelpath, nphi, postfix, light_model, ellipticity, position_angle)
		self.nrperbin = nrperbin
		self.observation = observation
		self.filter = filter 
		self.filters = filters
		self.rmin_kpc = rmin_kpc
		#self.rbins = rbins
		self.startcounts = startcounts
		#self.p_member_minimum = p_member_minimum
		#self.rmax_arcsec = rmax_arcsec
		#if self.rmax_arcsec is not None and self.p_member_minimum is not None:
		#	raise Exception("rmax_arcsec and p_member_minimum cannot both be an option, choose one")
		
	#def should_include(self, star):
	#	if self.rmax_arcsec is not None:
	#		return star.rc <= self.rmax_arcsec
	#	if self.p_member_minimum is not None:
	#		return star.p_member >= self.p_member_minimum
		
	def init(self):
		self.observation.load()
		
	def process(self):
		#self.stars_sets = self.observation.stars_sets
		#sets = self.stars_sets
		stars = self.observation.stars
		#if len(sets) == 1:
		if 1:
			self.process_single()
			self.aperture_rcenters = (self.rborders[0:-1] + self.rborders[1:])/2
			self.rcenters = self.aperture_rcenters
		else:
			print "set 0"
			rborders = self.calculate_bins(sets[0])
			self.rborders = zeros((len(sets), len(rborders)))
			self.rborders[0] = rborders
			for i in range(1, len(sets)):
				print "set", i
				self.rborders[i] = self.calculate_bins(sets[i])
			self.aperture_rcenters = (self.rborders[:,0:-1] + self.rborders[:,1:])/2
			self.rcenters = self.aperture_rcenters
		self.calculate_rest()
	
	def calculate_rest(self):
		self.radial_surface_densities = zeros_like(self.rcenters)
		#if 
		#for j in range(len(self.stars_sets)):
		if 1:
			rborders = self.rborders
			radial_surface_densities = []
			for i in range(len(rborders)-1):
				R1 = rborders[i]
				R2 = rborders[i+1]
				r1 = self.light_model.arcsec_to_kpc(R1)
				r2 = self.light_model.arcsec_to_kpc(R2)
				radial_surface_densities.append(self.light_model.light_profile.cumdensityR(r1, r2, M=1.))
				
			radial_surface_densities = array(radial_surface_densities)
			self.radial_surface_densities = radial_surface_densities
			logger.info("total fractional surface density in aperture: %.3f " % sum(radial_surface_densities)) 
		
				
	def calculate_bins(self, stars=None):
		if stars is None:
			stars = self.observation.stars
		if self.filter:
			stars = self.filter.filter(stars)
		if self.filters:
			for filter in self.filters:
				stars = filter(stars)
		#r = stars.rc # * 60 * 60
		#r = stars.xi 
		r = mab.astrometry.xi_eta_to_re(stars.xi, stars.eta, degrees(self.position_angle), self.ellipticity)
		if 0:
			logger.info("maximum radius: %f arcsec (=%f kpc)" % (max(r), self.light_model.arcsec_to_kpc(max(r))))
			if self.rmax_arcsec is not None:
				logger.info("using maximum radius of : %f (=%f kpc)" % (self.rmax_arcsec, self.light_model.arcsec_to_kpc(self.rmax_arcsec)))
			if self.p_member_minimum is not None:
				logger.info("filtering stars using a minimum p_member value of: %f" % self.p_member_minimum)
			stars = self.observation.stars.filter(lambda star: self.should_include(star))
			r = stars.rc
			logger.info("new maximum radius: %f arcsec (=%f kpc)" % (max(r), self.light_model.arcsec_to_kpc(max(r))))
		
		indices = argsort(r)
		r = r[indices]# * 1.0001
		i1 = 0
		N = len(r)
		rborders = []
		done = False
		#if self.rmin
		rborders.append(self.light_model.kpc_to_arcsec(self.rmin_kpc))
		logger.info("Rmin %f arcsec = %f kpc" % (self.light_model.kpc_to_arcsec(self.rmin_kpc), self.rmin_kpc))
		logger.info("log10 Rmin %f arcsec = %f kpc" % (log10(self.light_model.kpc_to_arcsec(self.rmin_kpc)), log10(self.rmin_kpc)))
		#rborders.append(0)
		#
		numbers_per_bin = []
		#logger.debug("maximum radius of borders: %f" % max(r))
		print N, max(r), argmax(r), "rmax", max(r)
		bin_index = 0
		while not done:
			#print i1, N, len(stars), stars
			if len(self.startcounts) > bin_index:
				i2 = i1 + self.startcounts[bin_index]
			else:
				i2 = i1 + self.nrperbin
			print i1, i2, N
			if i2 == N:
				done = True
				print len(r), i2
				rborders.append(r[-1])
				numbers_per_bin.append(i2-i1)
			elif i2 > N-1:
				i2 = N
				if (i2 - i1) < self.nrperbin*0.60:
					logger.info("merging last bin with %d stars" % (i2-i1))
					done = True
					print len(r), i2
					rborders[-1] = r[i2-1]
					print numbers_per_bin
					numbers_per_bin[-1] += (i2-i1)
				else:
					i2 = N
					#rborders.append((r[i2]+r[i2-1])/2)
					rborders.append(r[i2-1])
					numbers_per_bin.append(i2-i1)
			else:
				rborders.append((r[i2]+r[i2-1])/2)
				numbers_per_bin.append(i2-i1)
			i1 = i2
			bin_index += 1
		logger.info("growing outer r border by 0.1% so that all stars fit into the bins")
		print rborders, len(rborders)
		logger.info("borders: kpc: %s" % self.light_model.arcsec_to_kpc(array(rborders)))
		rborders[-1] *= 1.001
		#print rborders[-1]
		#print rborders[-1]
		logger.info("number of borders/bins: %d/%d" % (len(rborders), len(rborders)-1))
		rborders = array(rborders)
		assert N == sum(numbers_per_bin), "sanity check failed, not all stars are in bins?"
		print rborders
		logger.info("# stars per bin: %s" % " ".join(map(str, numbers_per_bin)))
		return rborders
			
			
	#def init(self, projectedmoments)
	#	self.aperturemoments = gdfast.ApertureGridMoments2dPolarIrrigularR(gridr, gridphi, constraintnrgrid, projectedmoments)
	#def 
		#print "aperture size(r x theta):", len(aperture_rborders)-1, "x", Nphi
		
		