# -*- coding: utf-8 -*-
import os
import mab.gd.logging as logging
logger = logging.getLogger("gd.schw.solution_discrete2")
from numpy import *
import numpy
import scipy
from mab.gd import gdfast_schw
from kaplot import *
import pyublas
import mab

class LosvdPerBin(object):
	def __init__(self, solution, storage_losvd, aperture, binned_data_m2):
		self.solution = solution
		self.storage_losvd = storage_losvd
		self.aperture = aperture
		self.binned_data_m2 = binned_data_m2
		
	def load(self):
		#self.solution.load()
		self.storage_losvd.load()
		self.binned_data_m2.load()
		self.binned_data_m2.observation.load()
		self.binned_data_m2.aperture.load()
		
		self.stars = []
		
		dirname = os.path.join(self.binned_data_m2.modelpath, "schw", "aperture")
		cached_filename= os.path.join(dirname, "observation_group" +self.binned_data_m2.postfix + "_cache" +".npy")
		stars_group = []
		no_groups = self.binned_data_m2.moments.shape[1]
		if not os.path.exists(cached_filename):
			for index in range(no_groups):
				print index
				stars = self.binned_data_m2.get_stars_from_aperture(index)
				print len(stars)
				stars_group.append(stars)
			
			logger.info("writing cache: " + cached_filename)
			numpy.save(cached_filename, stars_group)
		else:
			logger.info("reading cached: " + cached_filename)
			stars_group = numpy.load(cached_filename).tolist()
			#print stars_group
			
		for stars in stars_group:
			self.stars.extend(stars)
			
		self.stars_group = stars_group
		self.stars = mab.cvsfile.CsvObject(self.stars)
		return self.stars
		
			
		losvds = self.storage_losvd.losvds
		losvd = dot(losvds.T, self.solution.orbitweights).T
			
		self.losvd_per_bin = zeros((losvd.shape[0], no_groups))
		self.ks_D = []
		self.ks_p = []
		#mozaic(3,3,box)
		for j in range(2):
			for i in range(4):
				index = i + j * 4
				if index < no_groups:
					#select(i, j)
					#title("bin = %d" % index)
					stars = stars_group[index]
					_, Rindices = self.storage_losvd.stars_to_apertureindices(stars)
					#print Rindices
					for Rindex in Rindices:
						#print Rindex
						losvd_at_R = losvd[:,Rindex] * 1.
						delta_v = self.storage_losvd.delta_v
						sigma_v = 2. # TODO: hardcoded
						losvd_at_R = scipy.ndimage.gaussian_filter(losvd_at_R, [sigma_v/delta_v], mode='constant')
						#y = losvd_at_R/losvd_at_R.max()
						dv = self.storage_losvd.vcenters[1] - self.storage_losvd.vcenters[0]
						losvd_at_R = losvd_at_R / sum(losvd_at_R * dv)  # property normalize
						self.losvd_per_bin[:, index] += losvd_at_R
					self.losvd_per_bin[:, index] /= len(Rindices) # and get the average
					ycum = cumsum(self.losvd_per_bin[:, index])
					ycum /= max(ycum)
					#graph(self.storage_losvd.vcenters, ycum)
					#print self.storage_losvd.vcenters.shape, ycum.shape
					cdf = scipy.interpolate.interp1d(self.storage_losvd.vcenters, ycum)
					D, p = scipy.stats.kstest(stars.vlos, cdf)
					self.ks_D.append(D)
					self.ks_p.append(p)
					
					
		return self.stars
					
					
		

class SolutionKnown(object):
	def __init__(self, modelpath, dfname):
		self.modelpath = modelpath
		self.dfname = dfname
		
	def load(self):
		filename = os.path.join(self.modelpath, "df", "orbitweights_" +self.dfname +".npy")
		logger.info("using known solution from: %s" % filename)
		#filename = os.path.join(self.modelpath, "df", "orbitweights_tang.npy")
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		self.solution = c
		return c

class PhotometryFake(object):
	def __init__(self, sub_scope_true, storage_2d_losvd, solution_known, light_profile, N=10000):
		self.sub_scope_true = sub_scope_true
		self.storage_2d_losvd = storage_2d_losvd
		self.solution_known = solution_known
		self.light_profile  = light_profile
		self.N = N
		
	def load(self):
		self.sub_scope_true.load()
		self.storage_2d_losvd_true = self.sub_scope_true.subscope["storage_2d_losvd"]
		self.storage_2d_losvd_true.load()
		self.storage_2d_losvd.load()
		solution = self.solution_known.load()
		self.true_rho2d = numpy.tensordot(self.storage_2d_losvd_true.masses, solution, axes=[(0,),(0,)])
		
		self.mass_matrixN = self.storage_2d_losvd.masses * 1. / self.storage_2d_losvd.delta_R
		self.totalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)

		self.counts = self.true_rho2d
		self.counts /= sum(self.counts)
		self.counts *= self.N
		self.optimizer = gdfast_schw.OptimizationMatrixN(self.mass_matrixN, self.counts, self.totalmass_matrix)
		
		Rmax = self.storage_2d_losvd.Rborders[-1]
		self.totalmass = self.light_profile.cumdensityR(0, Rmax, M=1.)
		#print self.totalmass
		self.mass_vector = sum(self.storage_2d_losvd.masses, axis=1)
		#print "mass vector", self.mass_vector
		self.optimizer_extra = gdfast_schw.OptimizationNormalizeMass(self.mass_vector, self.totalmass, 0.0001)
		
		
	def logp(self, x):
		return self.optimizer.logp(x)
		
	def dlogpdx(self, x, gradient):
		self.optimizer.dlogpdx(x, gradient)

class PhotometryAnalytic(object):
	def __init__(self, storage_2d_losvd, light_profile, N=10000):
		self.storage_2d_losvd = storage_2d_losvd
		self.light_profile = light_profile
		self.N = N
		
	def logp(self, x):
		return self.optimizer.logp(x)*1 #+ self.optimizer_extra.logp(x)*1

	def dlogpdx(self, x, grad):
		self.optimizer.dlogpdx(x, grad)
		#self.optimizer_extra.dlogpdx(x, grad)

	def load(self):
		#	self.storage_2d_losvd.load()
		#solution = self.solution_known.load()
		#self.true_rho2d = numpy.tensordot(self.storage_2d_losvd.masses, solution, axes=[(0,),(0,)])
		Rborders = self.storage_2d_losvd.Rborders
		Rmax = self.storage_2d_losvd.Rborders[-1]
		#print "Rmax", Rmax
		#dsa
		N = len(Rborders)-1
		self.counts = array([self.light_profile.cumdensityR(Rborders[i], Rborders[i+1], M=1.) for i in range(N)])
		#print self.counts, sum(self.counts), len(self.counts)
		#dsa
		
		self.mass_matrixN = self.storage_2d_losvd.masses * 1. / self.storage_2d_losvd.delta_R
		self.totalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)

		#self.counts = self.true_rho2d
		self.counts /= sum(self.counts)
		self.counts *= self.N
		self.optimizer = gdfast_schw.OptimizationMatrixN(self.mass_matrixN, self.counts, self.totalmass_matrix)
		
		
		self.totalmass = self.light_profile.cumdensityR(0, Rmax, M=1.)
		#print self.totalmass
		self.mass_vector = sum(self.storage_2d_losvd.masses, axis=1)
		#print "mass vector", self.mass_vector
		self.optimizer_extra = gdfast_schw.OptimizationNormalizeMass(self.mass_vector, self.totalmass, 0.0001)
		
		
		#x = ones(160) / 4.
		#testgrad(self.optimizer_extra)
		#dsa
		
def testgrad(opt):
	x = ones(160)
	#x = random.random(160)
	x /= len(x) * 1.
	
	for i in [0, 10, 100]:
		dx = 1e-8
		x1 = x * 1.
		x2 = x * 1.
		x2[i] += dx
		p1 = opt.logp(x1)
		p2 = opt.logp(x2)
		grad = x * 0
		opt.dlogpdx(x, grad)
		#print p1, p2, 
		print (p2-p1)/dx, 
		print grad[i]
		
		
	
	
class KinematicsObservation(object):
	def __init__(self, storage_2d_losvd, observation, sigma_v = 2.01, v_index_clip=None, bootstrap_seed=None):
		self.storage_2d_losvd = storage_2d_losvd
		self.observation = observation
		self.sigma_v = sigma_v
		self.v_index_clip = v_index_clip
		self.bootstrap_seed = bootstrap_seed
		
	def logp(self, x):
		return self.optimizer.logp(x) - self.optimizer_extra.logp(x)
		
	def dlogpdx(self, x, grad):
		self.optimizer.dlogpdx(x, grad)
		gradex = grad * 0
		self.optimizer_extra.dlogpdx(x, gradex)
		grad -= gradex

		
	def load(self):
		self.storage_2d_losvd.load()
		stars = self.observation.load()
		print self.bootstrap_seed
		if self.bootstrap_seed:
			logger.info("bootstrapping kinematic sample with seed %d" % self.bootstrap_seed)
			import random
			random.seed(self.bootstrap_seed)
			N = len(stars)
			stars_bs = []
			for i in range(N):
				ri = int(random.random() * N)
				stars_bs.append(stars[ri])
			stars = mab.cvsfile.CsvObject(stars_bs)
			
			#import pdb
			#pdb.set_trace()
		
		logger.info("using %d stars/observations" % len(stars))
		#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		stars_inrange = stars.filter(lambda star: self.storage_2d_losvd.aperture.inrange(star.xi, star.eta))
		logger.info("stars in aperture range     : %d" % len(stars_inrange))
		logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		vmax = self.storage_2d_losvd.vmax
		#print "vmax", vmax
		delta_v = 2*vmax/self.storage_2d_losvd.Nv
		#print "res", delta_v
		
		losvds = self.storage_2d_losvd.losvds
		sigma_v = 2.01
		
		for i in range(losvds.shape[0]):
			#for j in range(losvds.shape[2]):
			#	pass
				#losvds[i,:,j] = scipy.ndimage.gaussian_filter(losvds[i,:,j], [self.sigma_v/self.storage_2d_losvd.delta_v], mode='constant')
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [self.sigma_v/self.storage_2d_losvd.delta_v, 0.0051], mode="constant")
			losvds[i] /= (self.storage_2d_losvd.delta_v * self.storage_2d_losvd.delta_R)
		
		
		if self.v_index_clip is not None:
			index_min = self.storage_2d_losvd.Nv/2-self.v_index_clip
			index_max = self.storage_2d_losvd.Nv/2+self.v_index_clip
			losvds = losvds[:,index_min:index_max,:]
			print "clip at v", self.v_index_clip * self.storage_2d_losvd.delta_v, "delta v is", self.storage_2d_losvd.delta_v, "Nv is", self.storage_2d_losvd.Nv
			#dsa
		else:
			index_min = 0
			index_max = self.storage_2d_losvd.Nv
			
		self.losvds = losvds
		
		
			
		for star in stars:
			star.aperture_index = self.storage_2d_losvd.aperture.findindex(star.xi, star.eta)
		
	
		#if 0:
		#	stars_invrange = stars.filter(lambda star: abs(star.vlos) < vmax)
		#	logger.info("stars in velocity range     : %d" % len(stars_invrange))
		#	logger.info("stars outside velocity range: %d" % (len(stars)-len(stars_invrange)))
		#	stars = stars_invrange
		self.used_stars = stars
		
		numpy.random.seed(8)
		for star in stars:
			#star.vlos = star.vlos_true# + numpy.random.normal(0, sigma_v)
			#star.vlos = star.vlos_true + numpy.random.normal(0, sigma_v)
			#print star.vlos, vmax, self.storage_2d_losvd.Nv
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d_losvd.Nv)  - index_min;
			outlier = True
			if star.v_index < 0:
				outlier = True
			elif star.v_index >= (index_max-index_min):
				outlier = True
			else:
				for losvd in losvds:
					if losvd[star.v_index, star.aperture_index] != 0:
						outlier = False
						break
			#print star.v_index, outlier
			star.is_outlier = outlier
		stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
		
		
		logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		stars = stars_no_outlier
		self.v_indices = v_indices = [star.v_index for star in stars]
		self.aperture_indices = aperture_indices = [star.aperture_index for star in stars]
		
		counts = losvds[0,:] * 0
		for vi, ai in zip(v_indices, aperture_indices):
			counts[vi, ai] += 1
		mask = counts > 0
		counts2d = counts
		
		self.counts_kin = counts[mask].reshape(-1)
		self.pmatrixN = (losvds[:,mask]).reshape((losvds.shape[0], -1)) * 1.
		#print len(self.counts_kin), self.counts_kin.sum()
		self.ptotalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)
		#self.ptotalmass_matrix = sum(sum(losvds, axis=1), axis=1)
		#print self.ptotalmass_matrix.shape
		#import pdb;
		#pdb.set_trace()
		self.pmatrixN = ascontiguousarray(self.pmatrixN).copy()
		self.counts_kin = ascontiguousarray(self.counts_kin).copy()
		self.ptotalmass_matrix = ascontiguousarray(self.ptotalmass_matrix).copy()
		self.optimizer = gdfast_schw.OptimizationMatrixN(self.pmatrixN, self.counts_kin, self.ptotalmass_matrix)
		

		counts = losvds[0,0,:] * 0
		for ai in aperture_indices:
			counts[ai] += 1
		mask = counts > 0
		
		self.counts_R = counts[mask].reshape(-1)
		#self.pmatrixRN = (losvds[:,mask]).reshape((losvds.shape[0], -1)) * 1.
		#print len(self.counts_R), self.counts_R.sum()
		
		self.pmatrixRN = sum(losvds * 1.0 * self.storage_2d_losvd.delta_v, axis=1)[:,mask] * 1.
		self.pmatrixRN = ascontiguousarray(self.pmatrixRN).copy()

		self.optimizer_extra = gdfast_schw.OptimizationMatrixN(self.pmatrixRN, self.counts_R, self.ptotalmass_matrix)
		#print self.storage_2d_losvd.vmax, self.storage_2d_losvd.delta_v, self.storage_2d_losvd.vmax/self.storage_2d_losvd.delta_v
		#ds
		if 0:
			box()
			indexedimage(counts2d)
			draw()
		#dsa
		if 0:
			print "testing gradient normal"
			testgrad(self)
			sys.exit(-1)
		
class KinematicsObservationForeground(object):
	def __init__(self, storage_2d_losvd, observation, member_ratios, model, foreground_model, sigma_v = 2.01, mean_vsys=None, simplemodel=None):
		self.storage_2d_losvd = storage_2d_losvd
		self.observation = observation
		self.sigma_v = sigma_v
		self.member_ratios = member_ratios
		self.model = model
		self.foreground_model = foreground_model
		self.mean_vsys = mean_vsys
		self.simplemodel = simplemodel
		
	def logp(self, x):
		return self.optimizer.logp(x)

	def dlogpdx(self, x, grad):
		self.optimizer.dlogpdx(x, grad)

	def load(self):
		self.storage_2d_losvd.load()
		self.simplemodel.load()
		stars = self.observation.load()
		
		logger.info("using %d stars/observations" % len(stars))
		#stars_inrange = stars.filter(lambda star: self.storage_2d.aperture.inrange(star.xi, star.eta))
		stars_inrange = stars.filter(lambda star: self.storage_2d_losvd.aperture.inrange(star.xi, star.eta))
		logger.info("stars in aperture range     : %d" % len(stars_inrange))
		logger.info("stars outside aperture range: %d" % (len(stars)-len(stars_inrange)))
		stars_invrange = stars_inrange.filter(lambda star: self.storage_2d_losvd.invrange(star.vlos))
		logger.info("stars in v range     : %d" % len(stars_invrange))
		logger.info("stars outside v range: %d" % (len(stars_inrange)-len(stars_invrange)))
		#rdsa
		vmax = self.storage_2d_losvd.vmax
		#print "vmax", vmax
		#delta_v = 2*vmax/self.storage_2d_losvd.Nv
		delta_v = self.storage_2d_losvd.delta_v
		#print "res", delta_v
		
		losvds = self.storage_2d_losvd.losvds
		#@self.losvds = losvds
		
		for i in range(losvds.shape[0]):
			#for j in range(losvds.shape[2]):
			#	pass
				#losvds[i,:,j] = scipy.ndimage.gaussian_filter(losvds[i,:,j], [self.sigma_v/self.storage_2d_losvd.delta_v], mode='constant')
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [self.sigma_v/self.storage_2d_losvd.delta_v, 0.0051], mode="constant")
			losvds[i] /= (self.storage_2d_losvd.delta_v * self.storage_2d_losvd.delta_R)
		
		
		stars = stars_invrange
		
		#sigma_v = 2.01
		numpy.random.seed(8)
		index_min = 0
		for star in stars:
			star.aperture_index = self.storage_2d_losvd.aperture.findindex(star.xi, star.eta)
			star.v_index = int(((star.vlos+vmax)/(2*vmax)) * self.storage_2d_losvd.Nv)  - index_min;
			#print star.v_index
			if 0:
				outlier = True
				if star.v_index < 0:
					outlier = True
				elif star.v_index >= (index_max-index_min):
					outlier = True
				else:
					for losvd in losvds:
						if losvd[star.v_index, star.aperture_index] != 0:
							outlier = False
							break
			#print star.v_index, outlier
			#star.is_outlier = outlier
		#stars_no_outlier = stars.filter(lambda star: not star.is_outlier)
		#logger.info("non-outlier stars : %d" % len(stars_no_outlier))
		#logger.info("outlier stars     : %d" % (len(stars)-len(stars_no_outlier)))
		#stars = stars_no_outlier
		self.v_indices = v_indices = [star.v_index for star in stars]
		self.aperture_indices = aperture_indices = [star.aperture_index for star in stars]
		
		
		
		#KinematicsObservation.load(self)
		self.model.load()
		self.member_ratios.load()
		if 0:
			self.ratios = []
			Rborders = self.storage_2d_losvd.Rborders
			for R1, R2 in zip(Rborders[:-1], Rborders[1:]):
				f1, f2 = self.member_ratios.ratios(R1, R2)
				self.ratios.append(f1/f2)
				print R1, R2, f1/f2, self.member_ratios.ratios1(R1)
			self.ratios = array(self.ratios)
			print self.ratios
		else:
			if 0:
				self.ratios = []
				for star in stars:
					Rkpc = self.storage_2d_losvd.light_model.arcsec_to_kpc(star.rc)
					ratio = self.member_ratios.ratios1(Rkpc)
					self.ratios.append(ratio)
					#print Rkpc, ratio
				self.ratios = array(self.ratios)
			
		#self.ratios = self.counts_kin * 0
		
		self.pmatrix = []
		self.pmatrixRN = []
		pmatrixRNai = sum(losvds * 1.0 * self.storage_2d_losvd.delta_v, axis=1)
		#print pmatrixRNai.shape
		for vi, ai in zip(self.v_indices, self.aperture_indices):
			if (vi >= 0) and (vi < losvds.shape[1]):
				self.pmatrix.append(losvds[:,vi,ai])
			else:
				self.pmatrix.append(losvds[:,0,ai]*0)
				assert False, "should not happen any more"
			self.pmatrixRN.append(pmatrixRNai[:,ai])
		self.pmatrix = ascontiguousarray(array(self.pmatrix).T * 1. * 2. / 2.).copy()
		self.pmatrixRN = ascontiguousarray(array(self.pmatrixRN).T * 1.)
		#print self.pmatrix.shape


		mean_sigma = 2.
		v1 = self.mean_vsys-self.storage_2d_losvd.vmax
		v2 = self.mean_vsys+self.storage_2d_losvd.vmax
		f_m_in = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v1, v2)[0]
		f_n_in = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v1, v2)[0]
		f_m_out = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v2, inf)[0]
		f_n_out = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v2, inf)[0]
		weight = f_m_in/f_n_in # additional weight, since we are interested in the vmin,vmax region

		
		#self.optimizer = gdfast_schw.OptimizationMatrixNForeground(self.pmatrixN, self.counts_kin, self.ptotalmass_matrix, self.ratios, self.p_v_non_members)
		x = ones(160) / 160.
		#test = (dot(self.losvds.T, x).sum(axis=1)*self.storage_2d_losvd.delta_v)/dot(self.pmatrixRN.T, x)
		#print test
		##numpy.testing.assert_approx_equal(test, ones(30))
		#import pdb; pdb.set_trace()
		
		#vs = [star.vlos_helio for star in stars]
		
		self.ratios = []
		self.p_v_non_members = []
		for star in stars:
			Rkpc = self.storage_2d_losvd.light_model.arcsec_to_kpc(star.rc)
			#Rindex = Rkpc * (0.0+self.storage_losvd.NR) /self.storage_losvd.Rmax
			#Rindex = int(Rindex)
			#ratio = self.member_ratios.ratios1(Rkpc)
			#self.ratios.append(ratio)
			Rarcsec = star.rc
			ratio = 0
			for i in range(self.simplemodel.catalogues):
				if int(star.catalogue_mask) & (1<<i):
					f_member = self.simplemodel.memberRe(Rarcsec, weight=weight, catalogue=i) # i)
					f_non_member = self.simplemodel.non_memberRe(Rarcsec, weight=weight, catalogue=i)#, i)
					#graph(Rarcsec, log10(f_member), color=color)
					#graph(Rarcsec, log10(f_non_member), color=color, linestyle="dash")
					ratio = max(f_member/f_non_member, ratio)
			#losvd_at_R = losvd[:,Rindex]
			#losvd_at_R = scipy.ndimage.gaussian_filter(losvd_at_R, [sigma_v/delta_v], mode='constant')
			#losvd_at_R = losvd_at_R / sum(losvd_at_R * delta_v)
			#p_member = ratio / (1+ratio)
			self.ratios.append(ratio)
			self.p_v_non_members.append(self.simplemodel.foreground_velocity_model(self.mean_vsys+star.vlos, mean_sigma)/f_n_in)
			
		self.ratios = array(self.ratios)# * 1e9 # * 0.25 # * 1e6
		
		#vs = [star.vlos for star in stars]
		#evs = [star.e_vlos for star in stars]
		
		#p_v_non_members = [exp(self.foreground_model.logL(v, e_v)) for v, e_v in zip(vs, evs)]
		self.p_v_non_members = array(self.p_v_non_members) #self.ratios * 0
		#print self.pmatrix.shape, self.pmatrixRN.shape
		#print self.ratios.shape, self.p_v_non_members.shape
		#l = [self.pmatrix, self.pmatrixRN, self.ratios, self.p_v_non_members]
		#for k in l:
		#	print k.shape, k.dtype
		#print self.p_v_non_members
		#import pdb
		#for k in l:
		#	print hex(id(k))
		#pdb.set_trace()
		#pyublas.set_trace(True)
		#self.p_v_non_members = self.p_v_non_members * 0
		#self.ratios = self.ratios * 0.1# + 1
		#return
		self.optimizer = gdfast_schw.OptimizationMatrixForegroundConditional(self.pmatrix, self.pmatrixRN, self.ratios, self.p_v_non_members)
		#print "test", self.optimizer.logp(x)
		#print "done"
		
		
		class dummy(object):
			def logp(self, x):
				return 0
			def dlogpdx(self, x, grad):
				return 0
				
		self.optimizer_extra = dummy()
		if 0:
			print "testing gradient fg"
			testgrad(self)
			sys.exit(-1)
		
		
		
		
class KinematicsFake(object):
	def __init__(self, sub_scope_true, storage_2d_losvd, solution_known, N=2000, sigma_v = 2.01, v_index_clip=None, bias_function=None):
		self.sub_scope_true = sub_scope_true
		self.storage_2d_losvd = storage_2d_losvd
		self.solution_known = solution_known
		self.N = N
		self.sigma_v = sigma_v
		self.v_index_clip = v_index_clip
		self.bias_function = bias_function
		
	def logp(self, x):
		return self.optimizer.logp(x)
		
	def dlogpdx(self, x, gradient):
		self.optimizer.dlogpdx(x, gradient)
		
		
	def load(self):
		self.sub_scope_true.load()
		self.storage_2d_losvd_true = self.sub_scope_true.subscope["storage_2d_losvd"]
		self.storage_2d_losvd.load()
		self.storage_2d_losvd_true.load()
		solution = self.solution_known.load()
		losvds = self.storage_2d_losvd.losvds * 1. # copy
		losvds_true = self.storage_2d_losvd_true.losvds * 1. # copy
		
		
		for i in range(losvds.shape[0]):
			#for j in range(losvds.shape[2]):
			#	pass
			losvds[i] = scipy.ndimage.gaussian_filter(losvds[i], [self.sigma_v/self.storage_2d_losvd.delta_v, 0.0051], mode="constant")
			losvds[i] /= (self.storage_2d_losvd.delta_v * self.storage_2d_losvd.delta_R)
			
			losvds_true[i] = scipy.ndimage.gaussian_filter(losvds_true[i], [self.sigma_v/self.storage_2d_losvd.delta_v, 0.0051], mode="constant")
			losvds_true[i] /= (self.storage_2d_losvd.delta_v * self.storage_2d_losvd.delta_R)
		
		if self.v_index_clip is not None:
			index_min = self.storage_2d_losvd.Nv/2-self.v_index_clip
			index_max = self.storage_2d_losvd.Nv/2+self.v_index_clip
			losvds = losvds[:,index_min:index_max,:]
			losvds_true = losvds[:,index_min:index_max,:]
			print "clip at v", self.v_index_clip * self.storage_2d_losvd.delta_v, "delta v is", self.storage_2d_losvd.delta_v, "Nv is", self.storage_2d_losvd.Nv
			#dsa
		else:
			index_min = 0
			index_max = self.storage_2d_losvd.Nv
		self.index_min = index_min
		self.index_max = index_max
		

		test = solution * 0 + 1
		test /= sum(test)
		self.testx = test
		# 'observation'
		if 1:
			self.true_losvd = numpy.tensordot(losvds_true, solution, axes=[(0,),(0,)])
		else:
			self.true_losvd = numpy.tensordot(self.storage_2d_losvd.losvds, solution, axes=[(0,),(0,)])
			#filename = os.path.join(self.storage_2d_losvd.modelpath, "df/losvd_tang.npy")
			#self.true_losvd = load(filename)
			self.true_losvd = scipy.ndimage.gaussian_filter(self.true_losvd, [self.sigma_v/self.storage_2d_losvd.delta_v, 0.0051], constant="constant")
			self.true_losvd /= (self.storage_2d_losvd.delta_v * self.storage_2d_losvd.delta_R)
			self.true_losvd = self.true_losvd[index_min:index_max,:]
		
		self.losvds_test = losvds
		counts_kin = self.true_losvd * 1.
		#self.counts_kin
		#for Rindex in range(self.true_losvd.shape[1]):
		#	counts_kin[:,Rindex] /= sum(counts_kin[:,Rindex])
		if self.bias_function:
			assert False
			assert len(self.storage_2d_losvd.Rcenters) == counts_kin.shape[1]
			for Rindex in range(counts_kin.shape[1]):
				# first, normalize the pdf at each radius, then weight it
				R = self.storage_2d_losvd.Rcenters[Rindex]
				counts = sum(counts_kin[:,Rindex])
				weight = self.bias_function(R, counts, Rindex)
				counts_kin[:,Rindex] *= 1./sum(counts_kin[:,Rindex]) * weight
				#if Rindex != 0:
				#	counts_kin[:,Rindex] = 0
		# always make sure it is normalized
		counts_kin /= sum(counts_kin)
		# then rescale by the fake number of observations
		counts_kin *= self.N
		#counts_kin[:,0] *= 10
		counts_kinR = sum(counts_kin, axis=0)
		counts_kin = counts_kin.reshape(-1) * 1.
		self.counts_kin = counts_kin * 1
		
		self.counts_kinR = counts_kinR * 1
		#print self.counts_kinR.sum(), self.counts_kin.sum()
		#dsa
		
		if 0:
			#a = self.true_losvd2 - self.true_losvd
			#print abs(a).max()
			box()
			#print self.true_losvd.shape
			indexedimage(counts_kin.reshape((index_max-index_min,30)))
			draw()
			dsa
		
		# probability matrix
		pmatrixN = losvds * 1.0
		pmatrixN = pmatrixN.reshape((pmatrixN.shape[0], -1)) * 1.
		self.pmatrixN = pmatrixN * 1.0
		
		#print losvds.shape
		#adsa
		self.ptotalmass_matrix = sum(sum(losvds, axis=1), axis=1)
		#self.ptotalmass_matrix = sum(self.storage_2d_losvd.masses, axis=1)
		#print self.ptotalmass_matrix.shape
		self.optimizer = gdfast_schw.OptimizationMatrixN(self.pmatrixN, self.counts_kin, self.ptotalmass_matrix)
		
		# probability matrix
		pmatrixRN = sum(losvds * 1.0 * self.storage_2d_losvd.delta_v, axis=1)
		#pmatrixRN = sum(losvds * 1.0, axis=1)
		#pmatrixRN = pmatrixRN.reshape((pmatrixRN.shape[0], -1)) * 1.
		self.pmatrixRN = pmatrixRN* 1.0
		#print self.pmatrixRN.shape
		#dsa
		
		self.optimizer_extra = gdfast_schw.OptimizationMatrixN(self.pmatrixRN, self.counts_kinR, self.ptotalmass_matrix)
		if 0:
			for k in [self.pmatrixN, self.counts_kin, self.ptotalmass_matrix]:
				print k.min(), k.max(), k.sum(), k.std()
			#dsa
		
		#opt_norm = gdfast_schw.OptimizationNormalize(1.-.000001, 0.001)
		#opt_entropy = gdfast_schw.OptimizationEntropy(1.e-2)
		#import pdb
		#pdb.set_trace()
		#debugger()
		#print self.pmatrixRN.shape
		x = ones((self.pmatrixRN.shape[0]))
		#c1 = self.optimizer.logp(x)
		#c2 = self.optimizer_extra.logp(x)
		#print c1, c2, c1-c2
		#sys.exit(0)

class Regularization(object):
	def __init__(self, regularizationQP, delta=1e-1*2):
		self.regularizationQP = regularizationQP
		self.delta = delta
		
	def init(self):
		Greg, hreg = self.regularizationQP.assemble_system()
		self.Greg = array(Greg)
		self.hreg = array(hreg)
		self.P = tensordot(self.Greg/self.delta, self.Greg/self.delta, axes=([1],[1]))
		self.q = -tensordot(array(hreg)/self.delta, self.Greg/self.delta, axes=([0], [1]))
		#n_nodes = self.Greg.shape
		self.optimizer = gdfast_schw.OptimizationQP(-self.P, -self.q)
		#testgrad(self.optimizer)
		#dsa
		
	def logp(self, x):
		return self.optimizer.logp(x)
		
	def dlogpdx(self, x, gradient):
		self.optimizer.dlogpdx(x, gradient)
		
class QP(object):
	def __init__(self, storage_2d_m2, storage_2d_m4, binned_data_m2, binned_data_m4, regularizationQP, delta=1e-1*2, use_fourth_moment=True):
		self.storage_2d_m2 = storage_2d_m2
		self.storage_2d_m4 = storage_2d_m4
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.regularizationQP = regularizationQP
		self.delta = delta
		self.use_fourth_moment = use_fourth_moment
		
	def load(self):
		self.init()
		
	def init(self):
		self.storage_2d_m2.load()
		self.storage_2d_m4.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		
		densityprojected_target_m2 = self.binned_data_m2.moments[0]
		densityprojected_target_m4 = self.binned_data_m4.moments[0]

		
		kinematic_m2 = self.binned_data_m2.moments[2]
		schw_m2 = self.storage_2d_m2.projectedmoments[:,2,:]
		e_kinematic_m2 = self.binned_data_m2.e_moments[2] #* 10
		
		kinematic_m4 = self.binned_data_m4.moments[4]
		schw_m4 = self.storage_2d_m4.projectedmoments[:,4,:]
		e_kinematic_m4 = self.binned_data_m4.e_moments[4] #* 10

		Greg, hreg = self.regularizationQP.assemble_system()
		self.Greg = array(Greg)
		self.hreg = array(hreg)
		
		Pt_list = []
		if self.delta:
			Pt_list.append(self.Greg/self.delta)
		
		Pt_list.append((kinematic_m2- schw_m2/densityprojected_target_m2)/(e_kinematic_m2))
		if self.use_fourth_moment:
			Pt_list.append((kinematic_m4- schw_m4/densityprojected_target_m4)/(e_kinematic_m4))

		Pt = hstack(tuple(Pt_list))
		self.P = tensordot(Pt, Pt, axes=([1],[1]))

		x_list = []
		if self.delta:
			x_list.append(array(self.hreg)/self.delta)
		x_list.append(kinematic_m2/e_kinematic_m2*0)
		if self.use_fourth_moment:
			x_list.append(kinematic_m4/e_kinematic_m4*0)
		x = hstack(tuple(x_list))
		self.q = -tensordot(x, Pt, axes=([0], [1]))

		#self.P = tensordot(self.Greg/self.delta, self.Greg/self.delta, axes=([1],[1]))
		#self.q = -tensordot(array(hreg)/self.delta, self.Greg/self.delta, axes=([0], [1]))
		#n_nodes = self.Greg.shape
		self.optimizer = gdfast_schw.OptimizationQP(-self.P, -self.q)
		#testgrad(self.optimizer)
		#dsa
		
	def logp(self, x):
		return self.optimizer.logp(x)
		
	def dlogpdx(self, x, gradient):
		self.optimizer.dlogpdx(x, gradient)
		
		
		
		

class Discrete(object):
	def __init__(self, kinematics, photometry, storage_2d_m0, storage_2d_m2, storage_2d_m4, binned_data_m2, binned_data_m4, storage_3d, modelpath, schwsetname, schwmodelname, postfix="", regularization=None):
		self.kinematics = kinematics
		self.photometry = photometry
		self.storage_2d_m0 = storage_2d_m0
		self.storage_2d_m2 = storage_2d_m2
		self.storage_2d_m4 = storage_2d_m4
		self.storage_3d = storage_3d
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.postfix = postfix
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.regularization = regularization
		
		
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		#self.debug = True
		self.debug = False
		
	def run(self, args, opts, scope):
		debug = self.debug
		self.init()
		#return
		x = self.solve(scope)
		#self.logL = -1
		#x = ones((160))/160.
		orbitweights = x
		#print sum(orbitweights)
		#orbitweights /= sum(orbitweights)
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		numpy.save(filename, orbitweights)
		
		moments_m0 = tensordot(orbitweights, self.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
		if debug:
			print "sanity check mass", sum(moments_m0[0])
		#import pdb;
		#pdb.set_trace()
		moments_m2 = tensordot(orbitweights, self.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
		moments_m4 = tensordot(orbitweights, self.storage_2d_m4.projectedmoments, axes=[(0,), (0,)])
		moments3d = tensordot(orbitweights, self.storage_3d.moments3d, axes=[(0,), (0,)])
		moments_m0[1:] /= moments_m0[0]
		moments_m2[1:] /= moments_m2[0] #self.binned_data_m2.moments[0]
		moments_m4[1:] /= moments_m4[0] #self.binned_data_m4.moments[0]
		moments3d[1:] /= moments3d[0]
		self.moments2d_solution_m2 = moments_m2
		self.moments2d_solution_m4 = moments_m4
		self.moments3d_solution = moments3d
		#print "2"
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m0.npy")
		numpy.save(filename, moments_m0)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		numpy.save(filename, moments_m2)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		numpy.save(filename, moments_m4)
		logger.debug("saving %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		numpy.save(filename, moments3d)
		logger.debug("saving %s" % filename)
		#print "3"
		
		
		chisq = -self.logL*2
		chisqs = [(name, chisq)  for name in "m2 m4 reg kin".split()]
		#print chisq, self.logL
		
		f = file(os.path.join(self.dirname, "results/probability" +self.postfix +".txt"), "w")
		#print f
		print >>f, "%e" % exp(-0.5*chisq)
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (-0.5*chisq)
		f.close()
		#print "%e" % (-0.5*chisq)
		
		logps = [(name, -0.5*chisq) for name, chisq in chisqs] 
		f = file(os.path.join(self.dirname, "results/logprobability_seperate" +self.postfix +".txt"), "w")
		#print "4"
		print >>f, repr(logps)
		
		#dsas
		
		
		
	def load(self):
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		self.orbitweights = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m0.npy")
		self.moments_solution_m0 = numpy.load(filename)
		logger.debug("loading %s"% filename)

		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		self.moments_solution_m2 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		self.moments_solution_m4 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		self.moments3d_solution = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		
		
	def init(self):
		self.kinematics.load()
		#return
		self.photometry.load()
		self.storage_2d_m0.init()
		self.storage_2d_m0.load()
		self.storage_2d_m2.init()
		self.storage_2d_m2.load()
		self.storage_2d_m4.init()
		self.storage_2d_m4.load()
		self.storage_3d.init()
		self.storage_3d.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		if self.regularization:
			self.regularization.init()
		return 
		
		self.observation.load()
		#self.storage_2d_binned.init()
		#self.storage_2d_binned.load()
		self.storage_2d_losvd.load()
		#self.storage_2d.aperture.load()
		self.aperture_light.load()
		
	def solve(self, scope):
		debug = self.debug
		debug = False
		#debug = True
		opt_norm = gdfast_schw.OptimizationNormalize(1.-.000001, 0.001)
		opt_ent = gdfast_schw.OptimizationEntropy(0.1)
		#print "fitting"
		def f_and_g(u):
			x = exp(u)
			grad = x * 0
			gradn = x * 0
			logp = 0
			logp += self.kinematics.logp(x)
			logp += self.photometry.logp(x)
			logp += opt_norm.logp(x)
			if self.regularization:
				logp += self.regularization.logp(x)
			#logp += opt_ent.logp(x)
			#logp = self.kinematics.logp(x) * 1 + \
			#	0 +\
			#	self.photometry.logp(x) * 1
				#opt_norm.logp(x) * 1+ \
				#self.kinematics.optimizer_extra.logp(x) * -1 
			#logp = opt_matrix_kinN.logp(x)*1 +\
			#	opt_matrix_massN.logp(x)*1 +\
			#	opt_norm.logp(x) * 1
			#	#+\
			#	#opt_entropy.logp(x) * 1.
			self.kinematics.dlogpdx(x, grad)
			#self.kinematics.optimizer_extra.dlogpdx(x, gradn)
			#grad -= gradn
			self.photometry.dlogpdx(x, grad)
			opt_norm.dlogpdx(x, grad)
			if self.regularization:
				self.regularization.dlogpdx(x, grad)
			#opt_ent.dlogpdx(x, grad)
			if debug:
				#print "%10f %10f %10f %10f %10f" % (logp,  sum(x), dot(totalmass_matrix, x/sum(x)), dot(ptotalmass_matrix, x/sum(x)), dot(losvds.T, x).sum() * delta_R * delta_v / dot(ptotalmass_matrix, x/sum(x)))
				#print "%5.15f %10f" % (logp,  sum(x)), self.photometry.optimizer.logp(x), self.kinematics.optimizer.logp(x)
				print "%5.15f %10f" % (logp,  sum(x)),  self.regularization.optimizer.logp(x)
				pass
			grad = grad*x
			return -logp, -grad
			#return -logp
		
		Nparams = 160 #self.kinematics.pmatrixRN.shape[0]
		x = ones(Nparams) * 1.
		x /= sum(x)
		logger.info("%d free parameters" % Nparams)
		#bounds = [(1e-10, 1) for i in range(Nparams)]
		bounds = None
		u = log(x)
		
		u = scipy.optimize.fmin_l_bfgs_b(f_and_g, u, None, bounds=bounds, approx_grad=False, iprint=-1,factr=1e-2,maxfun=200000)[0]
		x = exp(u)
		#x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=False, iprint=-1,factr=1e-2,maxfun=200000)[0]
		#x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
		#print "end fitting"
		#x /= sum(x)
		#print "f,g", f_and_g(x)[0]
		#sself.logL = -f_and_g(x)
		self.logL = -f_and_g(u)[0]
		#u = log(x)
		logger.info("logL = %f" % self.logL)
		
		#x = self.solution.orbitweights
		if 0:
			
			xlow = x * 1.
			xhigh = x * 1.
			xhigh = xhigh.reshape((8, 20))
			xlow = xlow.reshape((8, 20))
			xhigh[0:3,:]  = 0
			xlow[3:,:]  = 0
			xhigh = xhigh.reshape(8*20)
			xlow = xlow.reshape(8*20)
			mlow = tensordot(xlow, self.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
			mhigh = tensordot(xhigh, self.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])

			dr = self.storage_2d_m0.projection.gridR.rmax * 1. / mlow.shape[1]
			R = 0.100
			imax = int(R/dr)
			totlow = sum(mlow[0,:imax])
			tothigh = sum(mhigh[0,:imax])
			ratio = totlow/tothigh
			print "ratio", ratio, "inside", R
			import pdb
			pdb.set_trace()
			x = xhigh



		#self.#solution.orbitweights = c
		
		#print self.logL
		#import pdb
		#pdb.set_trace()
		#if debug:
		#	print "res:", [opt.logp(x) for opt in [opt_matrix_kinN, opt_matrix_massN, opt_norm]]
		return x

