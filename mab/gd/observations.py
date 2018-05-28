# -*- coding: utf-8 -*-
import mab.cvsfile
import mab.utils.numpy
import mab.parallelize
import os
from numpy import *
import numpy
import pyublas
import mab.gd.logging as logging
import mab.astrometry
import mab.utils.progressbar
import sys
from kaplot import *
import mab.asciifile
import scipy.ndimage
#import emcee

logger = logging.getLogger("gd.obs")

class ObservationFiltered(object):
	def __init__(self, observation, filters=[]):
		self.observation = observation
		self.filters = filters
		
	def load(self):
		stars = self.observation.load()
		for i, filter in enumerate(self.filters):
			if i == 0:
				logger.info("filter: before: %d stars" % len(stars))
			stars = filter(stars)
			logger.info("filter: after: %d stars" % len(stars))
		self.stars = stars
		return stars
		
class MockObservation(object):
	def __init__(self, modelpath, input, filters=[], vlos_helio_offset=None, **kwargs):
		self.modelpath = modelpath
		#self.moments = moments
		self.input = input
		self.filename = os.path.join(self.modelpath, "data", self.input)
		self.vlos_helio_offset = vlos_helio_offset
		self.filters = filters
		
	def load(self):
		logger.info("loading stars/observations from filename: %s" % self.filename)
		if self.filename.endswith(".csv"):
			self.stars = mab.cvsfile.readcsv(self.filename)
		else:
			self.stars_sets = load(self.filename).view(recarray)
			self.stars = self.stars_sets[0]
		logger.info("loading %d stars/observations" % len(self.stars))
		if self.vlos_helio_offset:
			logger.warning("USING OFFSET %f" % vlos_helio_offset)
			self.stars.vlos_helio += self.vlos_helio_offset
		#logger.info("%d sets" % len(self.stars_sets))
		for filter in self.filters:
			self.stars = self.stars.filter(filter)
		return self.stars
	
	def save(self, stars=None):
		if stars is not None:
			self.stars = stars
		logger.info("writing stars/observations to filename: %s ..." % self.filename)
		if self.filename.endswith(".csv"):
			mab.cvsfile.writecsv(self.filename, self.stars)
		else:
			raise "fixme"
		logger.info("wrote %d stars/observations)" % len(self.stars))
		
		#self.aperture = aperture
		
Observation = MockObservation

class ObservationWalker(object):
	def __init__(self, datapath,  input):
		self.filename = os.path.join(datapath, input)
		
	def load(self):
		logger.info("reading asci file %s" % self.filename)
		self.stars = mab.asciifile.read(self.filename, 10, 26, 30)
		return self.stars
		#import pdb; pdb.set_trace()
		
class ObservationGius(object):
	def __init__(self, datapath,  input):
		self.filename = os.path.join(datapath, input)
		
	def load(self):
		logger.info("reading asci file %s" % self.filename)
		self.lines = file(self.filename).readlines()
				
class CommandRminmax(object):
	def __init__(self, observations, light_model):
		self.observations = observations
		self.light_model = light_model
		
	def run(self, args, opts, scope):
		stars = self.observations.load()
		rs = sort(stars.rc)
		#rmin = stars.Rdeg.min()
		#rmax = stars.Rdeg.max()
		#print "Rmin=", rmin, "'' =", rmin/60, "' = ", rmin/3600, " deg", self.light_model.arcsec_to_kpc(rmin), " kpc"
		print self.light_model.arcsec_to_kpc(rs).tolist()
		print self.light_model.arcsec_to_kpc(rs[-1])
		box()
		rs = self.light_model.arcsec_to_kpc(rs)
		histogram(log10(rs), bincount=100)
		draw()
		
class CommandRmax(object):
	def __init__(self, simplemodel, ratio, member_filter):
		self.simplemodel = simplemodel
		self.ratio = ratio
		self.member_filter = member_filter
		
		
	def load(self):
		self.simplemodel.load()
		
		mean_sigma = 2.
		v1 = self.member_filter.v - 3*self.member_filter.sigma
		v2 = self.member_filter.v + 3*self.member_filter.sigma
		
		
		f_m_in = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v1, v2)[0]
		f_n_in = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v1, v2)[0]
		
		weight = f_m_in/f_n_in
		Rarcsec = 3000.

		self.Rmaxs = []
		for i in range(self.simplemodel.catalogues):
			def f(Rarcsec):
				f_member = self.simplemodel.memberRe(Rarcsec, weight=weight, catalogue=i) # i)
				f_non_member = self.simplemodel.non_memberRe(Rarcsec, weight=weight, catalogue=i)#, i)
				ratio = f_member/f_non_member
				return (ratio-self.ratio)**2
			Rarcsec = scipy.optimize.fmin(f, 1000., disp=False)[0]
			#print Rarcsec
			self.Rmaxs.append(Rarcsec)
		
	def run(self, args, opts, scope):
		self.load()
		print "arcsec", self.Rmaxs
		light_model = scope["light_model"]
		print "kpc", light_model.arcsec_to_kpc(array(self.Rmaxs))
		

		
		
		
		
class HelioToRestFrameSubtractGradient(object):
	def __init__(self, observation_helio, observation_systemic, vr, vsigma, filters=[]):
		self.observation_helio = observation_helio
		self.observation_systemic = observation_systemic
		self.vr = vr
		self.vsigma = vsigma
		self.filters = filters
		
	def run(self, args, opts, scope):
		stars = self.observation_helio.load()
		allstars = stars
		for filter in self.filters:
			stars = filter(stars)
		print len(allstars)
		mozaic(2,2,box)
		
		print "TEST"
		Nangles = 50
		Ngrads = 50
		angles = arange(0, Nangles) / (Nangles+0.) * 2 * pi
		grads = 10**((arange(0, Ngrads) / (Ngrads+0.) * 3 - 4.5))
		print grads
		#grads *= 0 
		x = stars.xi
		y = stars.eta
		logpgrid = zeros((Nangles, Ngrads))
		vlos_helio = stars.vlos_helio# - self.vr
		vr = mean(vlos_helio)
		print "vr", vr, self.vr
		for i, angle in enumerate(angles):
			print i
			for j, grad in enumerate(grads):
				xp = cos(angle) * x - sin(angle) * y
				yp = sin(angle) * x + cos(angle) * y
				vadd = grad * xp
				v = (vlos_helio) + vadd
				sigmas = sqrt(stars.e_vlos**2 + self.vsigma**2)
				ps = gaussian(v, vr, sigmas)
				#print mean(v), mean(sigmas)
				#print ps
				
				#print sum(log(ps)),
				logpgrid[i,j] = sum(log(ps))
		grid = exp(logpgrid - logpgrid.max())
		print grid.max()
		probimage2d(grid, 0, 1, angles, log10(grads))
		select(0, 1)
		probgraph(grid, 0, angles)
		select(1, 0)
		probgraph(grid, 1, log10(grads))
		index = argmax(grid)
		i,j = unravel_index(index, grid.shape)
		angle, grad = angles[i], grads[j]
		
		x = allstars.xi
		y = allstars.eta
		vlos_helio = allstars.vlos_helio# - self.vr
		xp = cos(angle) * x - sin(angle) * y
		yp = sin(angle) * x + cos(angle) * y
		vadd = grad * xp
		print len(xp), len(allstars)
		v = (vlos_helio) + vadd
		
		allstars.vlos = v - vr
		print "mean vlos", stars.vlos.mean()
		print "used mean vr", vr
		print "correction for proper motion: gradient =", grad, " angle =", angle, " =",degrees(angle),"degrees"
		for star in allstars:
			star.attributes.append("vlos")
		self.observation_systemic.save(allstars)
		print vr
		print angle, grad, grad*3600 
		#draw()

class HelioToRestFrameUsingProperMotions(object):
	def __init__(self, observation_helio, observation_systemic, vr, pm_ra, pm_dec, distance):
		self.observation_helio = observation_helio
		self.observation_systemic = observation_systemic
		self.pm_ra = pm_ra
		self.pm_dec = pm_dec
		self.vr = vr
		self.distance = distance
		
		
	def run(self, args, opts, scope):
		zero_ref_frame = mab.gd.gdfast.ZeroReferenceFrame()
		
		stars = self.observation_helio.load()
		
		cartesian = mab.gd.gdfast.Cartesian()
		spherical = mab.gd.gdfast.SphericalGalactic()
		cylindrical = mab.gd.gdfast.Cylindrical()
		
		eq_reference_frame = mab.gd.gdfast.EqReferenceFrame(radians(mab.astrometry.theta0), radians(mab.astrometry.al_NGP), radians(mab.astrometry.d_NGP), zero_ref_frame)
		
		k = 4.74057
		v_sun_lsr = array([10.0, 5.25, 7.17])
		v_lsr_gsr = array([0., 220, 0])
		v_sun_gsr = v_lsr_gsr + v_sun_lsr
		def vrpm_eq_to_car(vr, pm_ra, pm_dec, r, alpha, delta):
			coordinate = mab.gd.gdfast.Coordinate(1, alpha, -delta+pi/2, spherical)
			position = mab.gd.gdfast.Position(coordinate, eq_reference_frame)
			d = delta
			#print cos(d), math.degrees(d), math.degrees(d)
			velocity_coordinate_radial = mab.gd.gdfast.VelocityCoordinate(vr, pm_ra*1e-5 * k * r, -pm_dec*1e-5 * k * r, spherical)
			velocity = mab.gd.gdfast.Velocity(velocity_coordinate_radial, eq_reference_frame, position, eq_reference_frame)
			#print velocity.to_coordinate_system(zero_ref_frame, spherical)
			v_spherical = velocity.to_coordinate_system(zero_ref_frame, spherical)
			return velocity.to_coordinate_system(zero_ref_frame, cartesian) + v_sun_gsr
			#return velocity.to_coordinate_system(zero_ref_frame, cartesian) + v_sun_gsr
		
		def vcar_to_spherical(v, r, alpha, delta):
			coordinate = mab.gd.gdfast.Coordinate(1, alpha, -delta+pi/2, spherical)
			position = mab.gd.gdfast.Position(coordinate, eq_reference_frame)
			
			#velocity_coordinate_radial = VelocityCoordinate(vr, pm_ra*1e-5 * k *d, -pm_dec*1e-5 * k *d, spherical)
			velocity_coordinate = mab.gd.gdfast.VelocityCoordinate(v[0], v[1], v[2], cartesian)
			velocity = mab.gd.gdfast.Velocity(velocity_coordinate, zero_ref_frame, position, zero_ref_frame)
			#print velocity.to_coordinate_system(zero_ref_frame, spherical)
			#v = velocity.to_coordinate_system(zero_ref_frame, spherical)
			v = velocity.to_coordinate_system(zero_ref_frame, spherical)
			return v
		
		if 0:
			alpha, delta = mab.astrometry.read_ra_hour("02 39 53.1"), mab.astrometry.read_dec_deg("-34 30 16.0")
			alpha, delta = math.radians(alpha), math.radians(delta)
			
			coordinate = mab.gd.gdfast.Coordinate(self.distance, alpha, -delta+pi/2, spherical)
			position = mab.gd.gdfast.Position(coordinate, eq_reference_frame)
			print "x,y,z", position.to_global()
			x, y, z = position.to_global()
			x, y, z = x-8, y, z
			print "gc position", x, y, z
			coordinate_gc = mab.gd.gdfast.Coordinate(x, y, z, cartesian)
			p = position.to_coordinate_system(zero_ref_frame, spherical)
			print "r,angles", p, "l,b",(math.degrees(p[1])+360), -(math.degrees(p[2])-90)
			b = math.radians(-(math.degrees(p[2])-90))
			cosb = cos(b)
			print "b,cosb",math.degrees(b), cosb
			#print position.to_coordinate_system(zero_ref_frame, cartesian) 
			alpha, delta = math.radians(stars.ra[0]), math.radians(stars.dec[0])
			#print stars.ra[0], stars.dec[0]
			#print alpha, delta
			#dsa
			#alpha, delta = 0, math.radians(90)
			v = vrpm_eq_to_car(self.vr, self.pm_ra, self.pm_dec, self.distance * 1000, alpha, delta)
			v = vrpm_eq_to_car(53.3, 47.6, -36., 138. * 1000, alpha, delta)
			#v = vrpm_eq_to_car(53.3, 54.1, -27.5, 138. * 1000, alpha, delta)
			#v = vrpm_eq_to_car(self.vr, 0, 0, self.distance * 1, math.radians(stars.ra[0]), math.radians(stars.dec[0])) - v_sun_gsr
			print "vel cartesian", v, self.vr
			
			if 0:
				vr, pm_ra, pm_dec = 53.3, 54.1, -27.5
				r = self.distance * 1000
				velocity_coordinate_radial = mab.gd.gdfast.VelocityCoordinate(vr, pm_ra*1e-5 * k * r, -pm_dec*1e-5 * k * r, spherical)
				velocity = mab.gd.gdfast.Velocity(velocity_coordinate_radial, eq_reference_frame, position, eq_reference_frame)
				#print velocity.to_coordinate_system(zero_ref_frame, spherical)
				v_spherical = velocity.to_coordinate_system(zero_ref_frame, spherical)
				return velocity.to_coordinate_system(zero_ref_frame, cartesian) + v_sun_gsr
			
			
			
			
			d = self.distance * 1000
			velocity_coordinate = mab.gd.gdfast.VelocityCoordinate(v[0], v[1], v[2], cartesian)
			#coordinate = mab.gd.gdfast.Coordinate(1, alpha, -delta+pi/2, spherical)
			#position = mab.gd.gdfast.Position(coordinate, eq_reference_frame)
			#position = mab.gd.gdfast.Position(coordinate, zero_ref_frame)
			velocity = mab.gd.gdfast.Velocity(velocity_coordinate, zero_ref_frame, position, zero_ref_frame)
			pm = velocity.to_coordinate_system(zero_ref_frame, spherical)
			print "vr", pm[0], pm[1], pm[2], sqrt(pm[1]**2 + pm[2]**2)
			print "pm", pm[1]/1e-5/k/d#/cos(b)
			print "pm", -pm[2]/1e-5/k/d
			print "d", d
			
			position = mab.gd.gdfast.Position(coordinate_gc, zero_ref_frame)
			velocity = mab.gd.gdfast.Velocity(velocity_coordinate, zero_ref_frame, position, zero_ref_frame)
			pm = velocity.to_coordinate_system(zero_ref_frame, spherical)
			print "vr,vt", pm[0], sqrt(pm[1]**2 + pm[2]**2)
		
		
		vs_gc = [vrpm_eq_to_car(vr_helio, self.pm_ra, self.pm_dec, self.distance * 1000, math.radians(a), math.radians(d)) for vr_helio, a, d in zip(stars.vlos_helio, stars.ra, stars.dec)]
		
		vs_spherical = [vcar_to_spherical(v_gc, self.distance * 1000, math.radians(a), math.radians(d)) for v_gc, a, d in zip(vs_gc, stars.ra, stars.dec)]
		vs_spherical = array(vs_spherical)
		print "vr", mean(vs_spherical[:,0])#/(1e-5 * k * r)
		print mean(vs_spherical[:,1])/(1e-5 * k * self.distance * 1000)
		print mean(-vs_spherical[:,2])/(1e-5 * k * self.distance * 1000)
		vs_gc = array(vs_gc)
		print vs_gc.shape
		print [mean(vs_gc[:,k]) for k in range(3)]
		
		
		#stars = self.observation_systemic.load()
		print vs_spherical.shape
		print len(stars), len(vs_spherical[:,0])
		stars.vlos = vs_spherical[:,0]
		for star in stars:
			star.attributes.append("vlos")
		self.observation_systemic.save(stars)
		
		
class CommandAddRe(object):
	def __init__(self, input, output, ra0, dec0, PA=0., ellipticity=1.):
		self.input = input
		self.output = output
		self.ra0 = ra0
		self.dec0 = dec0
		self.PA = PA
		self.ellipticity = ellipticity
		
	def run(self, args, opts, scope):
		stars = self.input.load()
		count = 0
		for star in stars:
			if 0:
				values = star.ra.split()
				ra = str(values[1]) +" %02i" % int(values[2]) +" %06.3f" % float(values[3])
				values = star.dec.split()
				dec = str(values[1]) +" %02i" % int(values[2]) +" %06.3f" % float(values[3])
				ra = star.RAh/24.*360 + star.RAm/60. + star.RAs/60.**2
				dec = star.DEd + star.DEm/60. + star.DEs/60.**2
			
			ra, dec = mab.astrometry.read_ra_hour(star.ra), mab.astrometry.read_dec_deg(star.dec)
			xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
			xin, etan = mab.astrometry.xi_eta_to_xin_etan(xi, eta, self.PA, self.ellipticity)
			xi *= 60.**2
			eta *= 60.**2
			re = mab.astrometry.xi_eta_to_re(xi, eta, self.PA, self.ellipticity)
			star.re = re
			star.attributes.append("re")
			star.ra = ra
			star.dec = dec
			count += 1
		self.output.save(stars)
		logger.info("wrote %d observations to filename %s" % (count, self.output.filename))		
		
		
class CommandObservationWalkerConvert(object):
	def __init__(self, input, output, ra0, dec0, PA=0., ellipticity=1., Vhb=None):
		self.input = input
		self.output = output
		self.ra0 = ra0
		self.dec0 = dec0
		self.PA = PA
		self.ellipticity = ellipticity
		self.Vhb = Vhb
		print "l,b", mab.astrometry.eq_to_gal(ra0, dec0)
		
	def load(self):
		self.input.load()
		
	def convert(self):
		print len(self.input.stars)
		count = 0
		f = open(self.output.filename, "w")
		print >>f,"id,ra, dec, vlos_helio, e_vlos, xi, eta, re, rc, Rdeg, V, I, SigMg, e_SigMg, p_member, FeH "#, SNR, SigMg, e_SigMg, V, I, FeH"
		for star in self.input.stars:
			ra = star.RAh/24.*360 + star.RAm/60. + star.RAs/60.**2
			dec = star.DEd + star.DEm/60. + star.DEs/60.**2
			
			ra = str(star.RAh) +" %02i" % int(star.RAm) +" %06.3f" % float(star.RAs)
			dec = star["DE-"] + str(star.DEd) +" %02i" % int(star.DEm) +" %06.3f" % float(star.DEs)
			ra, dec = mab.astrometry.read_ra_hour(ra), mab.astrometry.read_dec_deg(dec)
			
			
			
			#print ra, dec
			#print star.RAh, star.RAm, star.RAs
			#print star["DE-"], star.DEd, star.DEm, star.DEs
			#if star["DE-"] == "-":
			#	dec *= -1
			xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
			
			PA = (self.PA) # hmm, not sure i fully understand
			#self.ellipticity = 0.
			xin, etan = mab.astrometry.xi_eta_to_xin_etan(xi, eta, PA, self.ellipticity)
			#re = astrometry.xi_eta_to_re(xi, eta, PA, e)
			
			#print xi, eta
			xi *= 60.**2
			eta *= 60.**2
			re = mab.astrometry.xi_eta_to_re(xi, eta, PA, self.ellipticity)
			#re = 0
			rc = sqrt(xi**2 + eta**2)
			Rdeg = rc / 60**2 
			V = star.Vmag
			I = V - star["V-I"]
			SigMg_prime = star.SigMg + 0.079 * (V-self.Vhb)
			FeH = 1.76 * SigMg_prime - 2.11
			values = [star.Target, ra, dec, star.VHel, star.e_VHel, xi, eta, re, rc, Rdeg, V, I, star.SigMg, star.e_SigMg, star.PM, FeH]
			values = [str(k) for k in values] 
			#star["V-I"], star.VHel
			#f = sys.stdout
			#print FeH
			#if FeH < -1.5:
			print >>f, ",".join(values)
			count += 1
		f.close()
		logger.info("wrote %d observations to filename %s" % (count, self.output.filename))
		
		
	def run(self, *args):
		self.load()
		self.convert()

class Besancon(object):
	def __init__(self, filename, deltav, vmin, vmax, smooth_sigmav, l, b, uvw_sun):
		self.filename = filename
		self.deltav = deltav
		self.vmin = vmin
		self.vmax = vmax
		self.smooth_sigmav = smooth_sigmav
		self.Nbins = int(round((vmax - vmin)/deltav))
		self.l = l
		self.b = b
		self.uvw_sun = uvw_sun
		
	def load(self):
		self.data = mab.asciifile.readsimple(self.filename, 95, 6, False)
		self.vlos_helio = self.data.Vr
		
		# calculate the radial component of the solar velocity
		cartesian = mab.gd.gdfast.Cartesian()
		spherical = mab.gd.gdfast.SphericalGalactic()
		zero_ref_frame = mab.gd.gdfast.ZeroReferenceFrame()
		coordinate = mab.gd.gdfast.Coordinate(1, self.l, -self.b+pi/2, spherical)
		position = mab.gd.gdfast.Position(coordinate, zero_ref_frame)
		u,v,w = self.uvw_sun
		velocity_coordinate = mab.gd.gdfast.VelocityCoordinate(u,v,w, cartesian)
		velocity = mab.gd.gdfast.Velocity(velocity_coordinate, zero_ref_frame, position, zero_ref_frame)
		vr = velocity.to_coordinate_system(zero_ref_frame, spherical)[0]
		
		self.vlos = self.vlos_helio + vr
		count, bins= numpy.histogram(self.vlos_helio, bins=self.Nbins, range=(self.vmin, self.vmax))#, new=True)
		count = scipy.ndimage.gaussian_filter1d(count, self.smooth_sigmav/self.deltav)
		count = count / float(sum(count))
		self.p_helio = count / self.deltav
		self.p_helio_bins = bins
		
		count, bins= numpy.histogram(self.vlos, bins=self.Nbins, range=(self.vmin, self.vmax))#, new=True)
		count = scipy.ndimage.gaussian_filter1d(count, self.smooth_sigmav/self.deltav)
		count = count / float(sum(count))
		self.p = count / self.deltav
		self.p_bins = bins
		return self.data
		
		
	def run(self, *args):
		self.load()
		print self.data.Vr
		
	def probability(self, v):
		index = (v - self.vmin) / (self.vmax - self.vmin) * self.Nbins
		if ((index >= 0) and (index < self.Nbins)):
			return self.p[index]
		else:
			return 0
		
	def probability_helio(self, v):
		index = (v - self.vmin) / (self.vmax - self.vmin) * self.Nbins
		if ((index >= 0) and (index < self.Nbins)):
			return self.p_helio[index]
		else:
			return 0
		
		

		
class CommandObservationGiusConvert(object):
	def __init__(self, input, output, ra0, dec0, PA=0., ellipticity=1.):
		self.input = input
		self.output = output
		self.ra0 = ra0
		self.dec0 = dec0
		self.PA = PA
		self.ellipticity = ellipticity
		
	def load(self):
		self.input.load()
		
	def convert(self):
		print len(self.input.lines)
		count = 0
		i = 0
		
		f = open(self.output.filename, "w")
		#f = sys.stdout
		print >>f,"id,ra, dec, vlos_helio, e_vlos, x_gius, y_gius, re_gius, xi, eta, xin, etan, re, rc, Rdeg, SNR, SigMg, e_SigMg, V, I, FeH"
		def ok(vlos_helio, e_vlos, SNR):
			return (SNR >= 10) & (e_vlos <= 5)
			
		x1 = []
		x2 = []
		for line in self.input.lines:
			line = line.strip()
			values = [k.strip() for k in line.split()]
			ra = str(values[1]) +" %02i" % int(values[2]) +" %06.3f" % float(values[3])
			dec = str(values[4]) +" %02i" % int(values[5]) +" %06.3f" % float(values[6])
			ra = mab.astrometry.read_ra_hour(ra)
			dec = mab.astrometry.read_dec_deg(dec)
			#ra = mab.astrometry.read_ra_hour("0 59 31.68")
			#dec = mab.astrometry.read_dec_deg("-32 25  2.4")
			vlos_helio = float(values[7])
			e_vlos = float(values[8])
			x = float(values[-3])
			y = float(values[-2])
			xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
			PA = (self.PA) # hmm, not sure i fully understand
			#self.ellipticity = 0.
			xin, etan = mab.astrometry.xi_eta_to_xin_etan(xi, eta, PA, self.ellipticity)
			#re = astrometry.xi_eta_to_re(xi, eta, PA, e)
			re = mab.astrometry.xi_eta_to_re(xi, eta, PA, self.ellipticity)
			
			rc = sqrt(xi**2+eta**2)
			#print re, values[-5]
			#print self.ra0, self.dec0
			#print xi, eta, rc, re
			#print rc, re, values, len(values), ra, dec
			#print values
			#import pdb
			#pdb.set_trace()
			#x1.append(dec)
			#x2.append(rc-float(values[-3]))
			#x1.append(re/(1-self.ellipticity))
			x1.append(re)
			x2.append(float(values[-3]))
			#print re/float(values[-3])
			dec_center = self.dec0
			ra_center = self.ra0
			#print 1/(sin(radians(dec_center)) * sin(radians(dec)) + cos(radians(dec_center)) * cos(radians(ra - ra_center)) * cos(radians(dec)))
			cosc = sin(radians(dec_center)) * sin(radians(dec)) + cos(radians(dec_center)) * cos(radians(dec)) * cos(radians(ra-ra_center))
			#print arccos(cosc) * 180/pi
			
			Rdeg = rc / 60**2
			re_gius = float(values[-3])
			SNR = float(values[11])/float(values[12])
			V = float(values[15])
			I = float(values[16])
			CaTEW = float(values[13])
			Vhb = 20.13
			W_T01p = CaTEW + 0.64 * (V-Vhb) # Eq. 2.11
			FeH_low = (-2.81) + 0.44 * W_T01p   # Eq  2.16
			
			id = "gius%04d" % i
			deg_to_arcsec = 60**2
			if ((V-I)>0.5) and ():
				continue
			#v = [ra, dec, vlos_helio, e_vlos, x*deg_to_arcsec, y*deg_to_arcsec, re_gius*deg_to_arcsec, xi*deg_to_arcsec, eta*deg_to_arcsec, re*deg_to_arcsec, re_gius, rc*deg_to_arcsec, Rdeg, SNR, 99999., 99999., V, I]
			v = [ra, dec, vlos_helio, e_vlos, x*deg_to_arcsec, y*deg_to_arcsec, re_gius*deg_to_arcsec, xi*deg_to_arcsec, eta*deg_to_arcsec, xin*deg_to_arcsec, etan*deg_to_arcsec, re*deg_to_arcsec,  rc*deg_to_arcsec, Rdeg, SNR, float(values[-2]), float(values[-1]), V, I, FeH_low]
			i += 1
			if ok(vlos_helio, e_vlos, SNR):
				print >>f, id + "," + ",".join([repr(k) for k in v])
				count += 1
		if 0:
			box()
			scatter(x1, x2)
			kaplot.line(0, 0, 1.5, 1.5)
			draw()
		#f.close()
		logger.info("wrote %d observations to filename %s" % (count, self.output.filename))
		
		
	def run(self, *args):
		self.load()
		self.convert()


class PhotometryModel(object):
	def __init__(self, data, photometry_profile, parameters_shape, parameter_profile):
		self.data = data
		self.photometry_profile = photometry_profile
		self.parameters_shape = parameters_shape
		self.parameter_profile = parameter_profile
		
	def run(self, *args):
		stars = self.data.load()
		cache = "cache.npy"
		if 1: #not os.path.exists(cache):
			for star in stars:
				values = star.ra.split()
				ra = str(values[0]) +" %02i" % int(values[1]) +" %06.3f" % float(values[2])
				values = star.dec.split()
				dec = str(values[0]) +" %02i" % int(values[1]) +" %06.3f" % float(values[2])
				ra = mab.astrometry.read_ra_hour(ra)
				dec = mab.astrometry.read_dec_deg(dec)
				ras.append(ra)
				decs.append(dec)
			ras = array(ras)
			decs = array(decs)
			#save(cache, [res, xs, ys] )
		else:
			res, xs, ys = load(cache)
			
			
		def lnprob(x, ivar=None):
			for i in range(len(self.parameters_shape)):
				setattr(self, self.parameters_shape[i].name, x[i])
			offset = len(self.parameters_shape)
			for i in range(offset, offset+len(parameters_profile)):
				setattr(self.photometry_profile, self.parameters_profile[i-offset].name, x[i])
			
				
				
			global model
			fom = 0
			model = dfL * 0
			for g in range(Ng):
				w, mu, sigma = x[g*3:g*3+3]
				w = (arctan(w)*2/pi+1)*50 # between 0 and 5
				sigma = exp(sigma)
				model += w * kaplot.gaussian(Ls, mu, sigma)
				#print "%10.4f %10.4f %10.4f %10.4f %10.4f" %(mux, muy, sigmax, sigmay, rho)
				print w, mu, sigma
			print ""
				
			#for e in range(nE):
			#	for l in range(nL):
			#		fom += df[l,e] * log(mvg(e, l, mux, muy, sigmax, sigmay, rho))
			#print fom
			chisq = sum((dfL - model)**2)
			print chisq
			#if isnan(fom):
			#	sys.exit(0)
			return chisq*10000000
			#print x
			#print ivar
			#dsa
			#print x
			#return -0.5
			
		n = 1
		ndim, nwalkers = Ng*3, 4*3
		ivar = 1. / numpy.random.rand(ndim)
		p0 = [numpy.random.rand(ndim) for i in range(nwalkers)]			
	
		
		
		

class PhotometryGius(object):
	def __init__(self, observations, position_angle, ellipticity, ra0, dec0, light_profile, light_model, gridder):
		self.observations = observations
		self.position_angle = position_angle
		self.ellipticity = ellipticity
		self.ra0 = ra0
		self.dec0 = dec0
		self.light_profile = light_profile
		self.light_model = light_model
		self.gridder = gridder
		
	def run(self, *args):
		stars = self.observations.load()
		stars.VminI = stars.V - stars.I
		#stars = stars.filter(lambda star: (star.V < 20.5) & (star.V > 20) & (star.VminI > 0.5) & (star.VminI < 0.7))
		#stars = stars.filter(lambda star: (star.V < 20.8) & (star.V > 19.8) & (star.VminI > -0.2) & (star.VminI < 0.3))
		#stars = stars.filter(lambda star: (star.V < 20.) & (star.V > 16.5) & (star.VminI > 0.7) & (star.VminI < 2))
		print len(stars)
		res = []
		xs = []
		ys = []
		cache = "cache.npy"
		if 1: #not os.path.exists(cache):
			for star in stars:
				values = star.ra.split()
				ra = str(values[0]) +" %02i" % int(values[1]) +" %06.3f" % float(values[2])
				values = star.dec.split()
				dec = str(values[0]) +" %02i" % int(values[1]) +" %06.3f" % float(values[2])
				ra = mab.astrometry.read_ra_hour(ra)
				dec = mab.astrometry.read_dec_deg(dec)
				xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
				xin, etan = mab.astrometry.xi_eta_to_xin_etan(xi, eta, self.position_angle, self.ellipticity)
				re = mab.astrometry.xi_eta_to_re(xi, eta, self.position_angle, self.ellipticity)
				res.append(re)
				xs.append(xin)
				ys.append(etan)
			res = array(res)
			xs = array(xs)
			ys = array(ys)
			save(cache, [res, xs, ys] )
		else:
			res, xs, ys = load(cache)
		print res
		#box()
		mozaic(3,3,box)
		h = histogram(res, binwidth=0.1/5, datamin=0, datamax=1.)
		
		select(2,0)
		xs = self.light_model.arcsec_to_kpc(xs*60*60)
		ys = self.light_model.arcsec_to_kpc(ys*60*60)
		#density2d(xs, ys, binc
		I = self.gridder(xs, ys)
		I = scipy.ndimage.gaussian_filter(I, [1., 1.])
		I /= I.max()
		I = log10(I)
		level = -3
		I[I<level] = level
		indexedimage(I, colormap="whiteblack")
		mask = I > level
		
		if 1:
			ResMax = 2.5
			res2 = arange(0, ResMax, 0.01)
			weights = res2 * 0
			Nangles = 400
			angles = arange(0, 2*pi, 2*pi/Nangles)
			for i in range(len(res2)):
				Re = res2[i]
				xtest = Re*cos(angles)
				ytest = Re*sin(angles)
				xindices = ((xtest - self.gridder.xmin) / (self.gridder.xmax-self.gridder.xmin) * self.gridder.nx).astype(int)
				yindices = ((ytest - self.gridder.ymin) / (self.gridder.ymax-self.gridder.ymin) * self.gridder.ny).astype(int)
				indexmask = (xindices < self.gridder.nx) & (xindices >= 0) & (yindices < self.gridder.ny) & (yindices >=0)
				values = mask[xindices[indexmask], yindices[indexmask]]
				weights[i] = 1.*sum(values) / Nangles
				
				#for angle in angles
			select(2,1)
			graph(weights)

		if 0:
			sum_scl = 0
			sum_const = 0
			for i in range(len(res2)-1):
				I, _ = scipy.integrate.quad(lambda r: 2*pi*r *self.light_profile.densityR(r, M=1), res2[i], res2[i+1])
				sum_scl += I * weights[i]
				I, _ = scipy.integrate.quad(lambda r: 2*pi*r, res2[i], res2[i+1])
				sum_const += I * weights[i]
			print sum_scl, sum_const
			#dsa
				
			res_kpc = self.light_model.arcsec_to_kpc(res*60*60)
			
			def f(x):
				b, ratio = 10**x[0], 10**x[1]
				self.light_profile.b = b
				sum_scl = 0
				for i in range(len(res2)-1):
					I, _ = scipy.integrate.quad(lambda r: 2*pi*r *self.light_profile.densityR(r, M=1), res2[i], res2[i+1])
					sum_scl += I * weights[i]
				print x
				#ratio = I0
				w1 = ratio/(1+ratio)
				w2 = 1/(1+ratio)
				ps = self.light_profile.densityR(res_kpc, M=1)/sum_scl * w1 + w2* 1./sum_const
				print ratio, w1, w2, w1+w2
				resindices = (res/ResMax * len(res2)).astype(int)
				return -sum(log(ps*weights[resindices]))
			x0 = [log10(0.3), log10(100)]
			bounds = None
			#x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
			x = scipy.optimize.fmin_l_bfgs_b(f, x0, None, bounds=bounds, approx_grad=True, iprint=-1,factr=1e-2,maxfun=200000)[0]
			print x
			print "values", 10**x[0], 10**x[1]
			print self.light_model.kpc_to_arcsec(10**x[0])/60
			#dsa

		select(1,0)
		scatter(xs, ys)
		select(0, 1)
		counts = h.data
		rs = h.bins[:-1] * 60
		area = (rs[1:]**2 - rs[:-1]**2) * 2 * pi
		print counts.shape, area.shape, rs.shape
		rcenters = (rs[1:] + rs[:-1])/2
		#graph(log10(rcenters), log10(counts/area))
		graph((rcenters), (counts/area))
		#import pdb
		#pdb.set_trace()
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters*60)
		dens = self.light_profile.densityR(rcenters_kpc)
		#graph(log10(rcenters), log10(dens/1e6)+0.5+0.25, color="red", linestyle="dash")
		#graph(log10(rcenters), log10(dens/1e6)+0.5+0.15-0.25+0.05, color="red", linestyle="dash")
		graph((rcenters), (dens/1e6)*10**(+0.5+0.15-0.25+0.05), color="red", linestyle="dash")
		#graph(log10(rcenters), log10(dens/1e6)+0.45, color="red", linestyle="dash")
		labels("R/kpc" , "surf. density")
		select(0,2)
		graph(log10(rcenters), log10(counts/area))
		graph(log10(rcenters), log10(dens/1e6)+0.5+0.15-0.25+0.05, color="red", linestyle="dash")
		labels("log R/kpc" , "log surf. density")
		
		
		print log10(rcenters), log10(dens)
		print self.light_profile
		#import pdb; pdb.set_trace()
		select(1,1)
		density2d(xs, ys, contourlevels=30, bincount=200)#, drawimage=False, drawcontour=True)
		draw()
			
		
		
		
class CommandObservationElseConvert(object):
	def __init__(self, input, output, ra0, dec0, PA=0., ellipticity=1.):
		self.input = input
		self.output = output
		self.ra0 = ra0
		self.dec0 = dec0
		self.PA = PA
		self.ellipticity = ellipticity
		
	def load(self):
		self.input.load()
		
	def convert(self):
		print len(self.input.lines)
		count = 0
		i = 0
		f = open(self.output.filename, "w")
		#f = sys.stdout
		print >>f,"id,ra, dec, vlos_helio, e_vlos, xi, eta, xin, etan, re, rc, Rdeg, SNR, V, I, FeH"
		def ok(vlos_helio, e_vlos, SNR):
			return (SNR >= 10) & (e_vlos <= 5)
		for line in self.input.lines[1:]:
			line = line.strip()
			values = [k.strip() for k in line.split()]
			ra = str(values[2]) +" %02i" % int(values[3]) +" %06.3f" % float(values[4])
			dec = str(values[5]) +" %02i" % int(values[6]) +" %06.3f" % float(values[7])
			ra = mab.astrometry.read_ra_hour(ra)
			dec = mab.astrometry.read_dec_deg(dec)
			vlos_helio = float(values[8])
			e_vlos = float(values[9])
			#x = float(values[-3])
			#y = float(values[-2])
			xi, eta = mab.astrometry.ra_dec_to_xi_eta(ra, dec, self.ra0, self.dec0)
			PA = (self.PA) # hmm, not sure i fully understand
			#self.ellipticity = 0.
			xin, etan = mab.astrometry.xi_eta_to_xin_etan(xi, eta, PA, self.ellipticity)
			#re = astrometry.xi_eta_to_re(xi, eta, PA, e)
			re = mab.astrometry.xi_eta_to_re(xi, eta, PA, self.ellipticity)
			rc = sqrt(xi**2+eta**2)
			Rdeg = rc / 60**2
			#import pdb;
			#pdb.set_trace()
			#dsas
			#re_gius = float(values[-3])
			#import pdb
			#pdb.set_trace()
			SNR = float(values[12])/float(values[13])
			V = float(values[16])
			I = float(values[17])
			#$CaTEW = float(values[13])
			#Vhb = 20.13
			#W_T01p = CaTEW + 0.64 * (V-Vhb) # Eq. 2.11
			#FeH_low = (-2.81) + 0.44 * W_T01p   # Eq  2.16
			FeH_low = float(values[20])
			
			#id = "gius%04d" % i
			id = values[1]
			deg_to_arcsec = 60**2
			#if ((V-I)>0.5) and ():
			#	continue
			#v = [ra, dec, vlos_helio, e_vlos, x*deg_to_arcsec, y*deg_to_arcsec, re_gius*deg_to_arcsec, xi*deg_to_arcsec, eta*deg_to_arcsec, re*deg_to_arcsec, re_gius, rc*deg_to_arcsec, Rdeg, SNR, 99999., 99999., V, I]
			v = [ra, dec, vlos_helio, e_vlos, xi*deg_to_arcsec, eta*deg_to_arcsec, xin*deg_to_arcsec, etan*deg_to_arcsec, re*deg_to_arcsec,  rc*deg_to_arcsec, Rdeg, SNR, V, I, FeH_low]
			i += 1
			if ok(vlos_helio, e_vlos, SNR):
				print >>f, id + "," + ",".join([repr(k) for k in v])
				count += 1
		
		#f.close()
		logger.info("wrote %d observations to filename %s" % (count, self.output.filename))
		
		
	def run(self, *args):
		self.load()
		self.convert()
		
		
class MembershipTest(object):
	def __init__(self, observation, light_profile, aperture, member_estimate_filter, foreground, output):
		self.observation = observation
		self.light_profile = light_profile
		self.aperture = aperture
		self.member_estimate_filter = member_estimate_filter
		self.foreground = foreground
		self.output = output
		
	def run(self, args, opts, scope):
		all_stars = self.observation.load()
		def ok(star):
			return (star.vlos > -200) & (star.vlos < 50)
		all_stars = all_stars.filter(ok)
		all_member_stars = self.member_estimate_filter(all_stars)
		all_nomember_stars = [star for star in all_stars if star not in all_member_stars]
		logger.info("member: %d non-member: %d total: %d" % (len(all_member_stars), len(all_nomember_stars), len(all_stars)))
		
		#self.light_profile.load()
		self.aperture.load()
		
		#print scipy.integrate.quad(lambda x: self.foreground(x, 10.), -1000, 1000)


		borders = self.aperture.aperture_rborders
		borders_kpc = self.aperture.aperture_rborders_kpc
		
		box()
		outputstars = []
		#@mab.parallelize.parallelize(cores=opts["cores"], info=opts["show_progressbar"])
		def do(i, draw=False):
			r1 = borders[i]
			r2 = borders[i+1]
			stars = [star for star in all_stars if ((star.rc >= r1) and (star.rc < r2))]
			member_stars = [star for star in all_member_stars if ((star.rc >= r1) and (star.rc < r2))]
			nomember_stars = [star for star in all_nomember_stars if ((star.rc >= r1) and (star.rc < r2))]
			
			N = len(stars)
			N_members = len(member_stars)
			N_nomembers = len(nomember_stars)
			ratio = N_members*1./N_nomembers
			fraction_members =  N_members*1. / N
			print i, "fraction", fraction_members, N
			return
			
			
			#print "logps", logps
			logsigmas = arange(log(1.), log(100), 0.05)
			sigmas = 10**logsigmas
			logps = array([logp(mean_v, sigma) for sigma in sigmas])
			graph(logsigmas, exp(logps))
			ps = exp(logps-logps.max())
			mean_sigma = sum(sigmas*ps)/sum(ps)
			print "mean_sigma", mean_sigma
			
			def p_member_j(v, e_v, mean_v, sigma):
				p_no_member = (1-fraction_members) * self.foreground(v, e_v)
				p_is_member = fraction_members * gaussian(v, mean_v, sqrt(sigma**2+e_v**2))
				p_total =  p_is_member + p_no_member
				return p_is_member / (p_total)
			ps = [p_member_j(star.vlos, star.e_vlos, mean_v, mean_sigma) for star in stars]
			for star, p in zip(stars, ps):
				star.attributes.append("p_member")
				star.p_member = p
				outputstars.append(star)
			ps = array(ps)
			
			#if i == 6:
			if draw:
				vs = [star.vlos for star in stars]
				histogram(vs, datamin=-200, datamax=100, bincount=100, normalize=True)
				v = arange(-200, 100)
				#fraction_members = 0.0
				#mean_sigma = 10.
				pv = fraction_members * gaussian(v, 0, mean_sigma) +  (1-fraction_members) * self.foreground(v, 0)
				graph(v, pv)
				print max(vs), min(vs)
			#print ps
			#print ps.mean(), ps.min(), ps.max()
			#def p_member_j(v, e_v):
			#	ps = array([
				
			#member_stars = [star for star in all_member_stars if ((star.re >= r1) and (star.re < r2))]
			#stars = [star for star in member_stars if ((star.re_gius >= r1) and (star.re_gius < r2))]
			#counts.append(len(stars))
			#r1 = borders_kpc[i]
			#r2 = borders_kpc[i+1]
		
		if 1:
			for i in range(len(borders)-1):
				do(i)
			#do(12)
		self.output.save(outputstars)
		#draw()

		
		


class CommandObservationsMerge(object):
	def __init__(self, inputs, output, blacklist):
		self.inputs = inputs
		self.output = output
		self.blacklist =blacklist
		
	def load(self):
		for input in self.inputs:
			input.load()
		
	def merge(self):
		allstars = []
		for i, input in enumerate(self.inputs):
			allstars += list(input.stars)
			input.stars.catalogue_mask = ones(len(input.stars), dtype=int) << i
				
		allstars_copy = list(allstars)
		unique_stars = []
		logger.info("total # stars: %d" % len(allstars))
		
		def similar_ra_dec(star1, star2):
			distance1 = sqrt((star1.eta-star2.eta)**2)
			distance2 = sqrt((star1.xi-star2.xi)**2)
			distance = sqrt(distance1**2+distance2**2)
			maxdistance = 1
			#return (distance1 < maxdistance) and (distance2 < maxdistance)
			return distance < maxdistance
		
		def similar_vlos(star1, star2):
			sigma = sqrt(star1.e_vlos**2 + star2.e_vlos**2)
			sigmas = (star1.vlos_helio-star2.vlos_helio)/sigma
			return abs(sigmas) < 3
		
		def similar(star1, star2):
			return similar_ra_dec(star1, star2) and similar_vlos(star1, star2)
			
		
		i = 0
		counts = {}
		for i in range(5):
			counts[i] = 0
		while allstars:
			print ".",
			sys.stdout.flush()
			current_star = allstars[0]
			similar_stars = [current_star]
			for j in range(1, len(allstars)):
				if similar(current_star, allstars[j]):
					similar_stars.append(allstars[j])
			for star in similar_stars:
				allstars.remove(star)
			if len(similar_stars) > 1:
				logger.debug("duplicates found: %d" % len(similar_stars))
				counts[len(similar_stars)] += 1
			if len(similar_stars) > 1:
				problem = False
				for star in similar_stars[1:]:
					sigma = sqrt(current_star.e_vlos**2 + star.e_vlos**2)
					sigmas = (current_star.vlos_helio-star.vlos_helio)/sigma
					if abs(sigmas) > 3:
						distance1 = sqrt((current_star.eta-star.eta)**2)
						distance2 = sqrt((current_star.xi-star.xi)**2)
						distance = sqrt(distance1**2+distance2**2)
						print "*** WARNING, vlos doesn't correspond", sigmas, "avg sigma =", sigma, "distance =", distance
						problem = True
				if problem: 
					for star in similar_stars:
						print star.id, star.eta, star.xi, star.vlos_helio, star.e_vlos
					print
					
			star = current_star.clone()
			star.catalogue_mask = 0
			for similar_star in similar_stars:
				star.catalogue_mask |= similar_star.catalogue_mask
				star.attributes.append("catalogue_mask")
			star.vlos_helio = mean([k.vlos_helio for k in similar_stars])
			star.e_vlos = sqrt(mean([k.e_vlos**2 for k in similar_stars]))
			#sigmgs = [k.SigMg for k in similar_stars if k.SigMg != 99999.0]
			#if len(sigmgs) > 0:
			#	star.SigMg = mean(sigmgs)
			#	star.e_SigMg = sqrt(mean([k.e_SigMg**2 for k in similar_stars if k.SigMg != 99999.0]))
			#else:
			#	star.SigMg = 99999.
			#	star.e_SigMg = 99999.
			unique_stars.append(star)
		print "unique stars:", len(unique_stars), "out of", len(allstars_copy)
		for i in range(5):
			print counts[i]
		#blacklist = ["Scl-0280", "Scl-0307"]
		unique_stars = [k for k in unique_stars if k.id not in self.blacklist]
		print "unique stars (after blacklist):", len(unique_stars), "out of", len(allstars_copy)
		
		outputfilename = self.output.filename #os.path.join(datapath,  "vlos_helio_fnx_gw.csv")
		
		stars = mab.cvsfile.CsvObject(unique_stars)
		f = file(outputfilename, "w")
		columns = "id ra dec xi eta rc re Rdeg vlos_helio e_vlos catalogue_mask".split() # SigMg e_SigMg
		mab.cvsfile.writecolumns(f, [stars.id, stars.ra, stars.dec, stars.xi, stars.eta, stars.rc, stars.re, stars.Rdeg, stars.vlos_helio, stars.e_vlos, stars.catalogue_mask], columns) # stars.SigMg, stars.e_SigMg
		count = len(stars)
		logger.info("wrote %d observations to filename %s" % (count, outputfilename))
		
	def run(self, *args):
		self.load()
		self.merge()
		
class CommandObservationsMerge2(object):
	def __init__(self, inputs, output, blacklist):
		self.inputs = inputs
		self.output = output
		self.blacklist =blacklist
		
	def load(self):
		for input in self.inputs:
			input.load()
		
	def merge(self):
		allstars = []
		for i, input in enumerate(self.inputs):
			allstars += list(input.stars)
			input.stars.catalogue_mask = ones(len(input.stars), dtype=int) << i
				
		allstars_copy = list(allstars)
		unique_stars = []
		logger.info("total # stars: %d" % len(allstars))
		
		def similar_ra_dec(star1, star2):
			distance1 = sqrt((star1.eta-star2.eta)**2)
			distance2 = sqrt((star1.xi-star2.xi)**2)
			distance = sqrt(distance1**2+distance2**2)
			maxdistance = 1
			#return (distance1 < maxdistance) and (distance2 < maxdistance)
			return distance < maxdistance
		
		def similar_vlos(star1, star2):
			sigma = sqrt(star1.e_vlos**2 + star2.e_vlos**2)
			sigmas = (star1.vlos_helio-star2.vlos_helio)/sigma
			return abs(sigmas) < 3
		
		def similar(star1, star2):
			return similar_ra_dec(star1, star2) and similar_vlos(star1, star2)
			
		
		i = 0
		while allstars:
			print ".",
			sys.stdout.flush()
			current_star = allstars[0]
			similar_stars = [current_star]
			for j in range(1, len(allstars)):
				if similar(current_star, allstars[j]):
					similar_stars.append(allstars[j])
			for star in similar_stars:
				allstars.remove(star)
			if len(similar_stars) > 1:
				logger.debug("duplicates found: %d" % len(similar_stars))
			if len(similar_stars) > 1:
				problem = False
				for star in similar_stars[1:]:
					sigma = sqrt(current_star.e_vlos**2 + star.e_vlos**2)
					sigmas = (current_star.vlos_helio-star.vlos_helio)/sigma
					if abs(sigmas) > 3:
						distance1 = sqrt((current_star.eta-star.eta)**2)
						distance2 = sqrt((current_star.xi-star.xi)**2)
						distance = sqrt(distance1**2+distance2**2)
						print "*** WARNING, vlos doesn't correspond", sigmas, "avg sigma =", sigma, "distance =", distance
						problem = True
				if problem: 
					for star in similar_stars:
						print star.id, star.eta, star.xi, star.vlos_helio, star.e_vlos
					print
					
			star = current_star.clone()
			star.catalogue_mask = 0
			for similar_star in similar_stars:
				star.catalogue_mask |= similar_star.catalogue_mask
				star.attributes.append("catalogue_mask")
			star.vlos_helio = mean([k.vlos_helio for k in similar_stars])
			star.vlos = mean([k.vlos for k in similar_stars])
			star.e_vlos = sqrt(mean([k.e_vlos**2 for k in similar_stars]))
			star.p_member = mean([k.p_member for k in similar_stars])
			#sigmgs = [k.SigMg for k in similar_stars if k.SigMg != 99999.0]
			#if len(sigmgs) > 0:
			#	star.SigMg = mean(sigmgs)
			#	star.e_SigMg = sqrt(mean([k.e_SigMg**2 for k in similar_stars if k.SigMg != 99999.0]))
			#else:
			#	star.SigMg = 99999.
			#	star.e_SigMg = 99999.
			unique_stars.append(star)
		print "unique stars:", len(unique_stars), "out of", len(allstars_copy)
		
		#blacklist = ["Scl-0280", "Scl-0307"]
		unique_stars = [k for k in unique_stars if k.id not in self.blacklist]
		print "unique stars (after blacklist):", len(unique_stars), "out of", len(allstars_copy)
		
		outputfilename = self.output.filename #os.path.join(datapath,  "vlos_helio_fnx_gw.csv")
		
		stars = mab.cvsfile.CsvObject(unique_stars)
		f = file(outputfilename, "w")
		columns = "id ra dec xi eta rc re Rdeg vlos_helio vlos e_vlos catalogue_mask p_member".split() # SigMg e_SigMg
		mab.cvsfile.writecolumns(f, [stars.id, stars.ra, stars.dec, stars.xi, stars.eta, stars.rc, stars.re, stars.Rdeg, stars.vlos_helio, stars.vlos, stars.e_vlos, stars.catalogue_mask, stars.p_member], columns) # stars.SigMg, stars.e_SigMg
		count = len(stars)
		logger.info("wrote %d observations to filename %s" % (count, outputfilename))
		
	def run(self, *args):
		self.load()
		self.merge()

class BesanconGaussianFit(object):
	def __init__(self, besancon_data, observation, parameterset, vmin, vmax, vdelta):
		self.besancon_data = besancon_data
		self.observation = observation
		self.parameterset = parameterset
		self.vmin = vmin
		self.vmax = vmax
		self.vdelta = vdelta
		self.Nbins = int(round((vmax - vmin)/vdelta))
		self.filename = os.path.join(self.parameterset.modelpath, "data", self.parameterset.type +"_" +self.parameterset.name+".npy")
		
	def run(self, *args):
		self.besancon_data.load()
		self.observation.load()
		def filter(star):
			return (star.vlos >= self.vmin) and (star.vlos < self.vmax)
		#def filter(star):
		#	return (star.vlos >= 0) and (star.vlos < 100)
		stars = self.observation.stars.filter(filter)
		
		count_besancon, bins_becanson = numpy.histogram(self.besancon_data.vlos, bins=self.Nbins, range=(self.vmin, self.vmax), new=True)
		count_stars, bins = numpy.histogram(stars.vlos, bins=self.Nbins, range=(self.vmin, self.vmax), new=True)
		bin_centers = bins[0:-1] + self.vdelta/2.
		#print data, bins
		#print mean(self.observation.stars.e_vlos)
		
		besancon_indices = ((stars.vlos - self.vmin) * self.Nbins / (self.vmax - self.vmin)).astype(int) 
		
		count_besancon_norm = count_besancon / float(sum(count_besancon))
		count_besancon_norm = scipy.ndimage.gaussian_filter1d(count_besancon_norm, 2./self.vdelta)
		count_besancon_norm = count_besancon_norm / float(sum(count_besancon_norm))
		pdf_besancon_discrete = count_besancon_norm / self.vdelta
		
		document(size="30cm,25cm")
		page(fontsize="11pt")
		mozaic(3,3,box)
		#print len(bins), len(count_stars)
		histogramline(bins, count_besancon, color="red", )
		histogramline(bins, pdf_besancon_discrete*1.0/sum(pdf_besancon_discrete)*len(self.besancon_data.vlos), color="orange", )
		histogramline(bins, count_stars, color="black", )
		labels("v (km/s)", "count")
		
		select(0, 1)
		labels("v (km/s)", "cumulative count")

		cumcount_besancon = cumsum(count_besancon)
		cumcount_stars = cumsum(count_stars)
		histogramline(bins, cumcount_besancon/float(max(cumcount_besancon)), color="red", )
		x = pdf_besancon_discrete/float(max(pdf_besancon_discrete))
		x = cumsum(x)
		x = x / float(max(x))
		histogramline(bins, x, color="orange", )
		histogramline(bins, cumcount_stars/float(max(cumcount_stars)), color="black", )
		print len(self.parameterset.parameter_values)
		#grid = zeros((self.parameterset.n,) * 3)
		i = 0
		errors = []
		logps = []
		for values in self.parameterset.parameter_values:
			#print values
			mu, sigma, scale = values.dict["mu"], values.dict["sigma"], values.dict["scale"]
			sigmas = sqrt(stars.e_vlos**2 + sigma**2)# * 0 + sigma
			ps_gaussian = gaussian(stars.vlos, mu, sigmas)
			ps_basancon = pdf_besancon_discrete[besancon_indices]
			#scale = 1.0
			ps = ps_gaussian * scale + ps_basancon* (1-scale)
			logp = sum(log(ps))
			logps.append(logp)
			#model = (gaussian(bin_centers, mu, sigma)*scale + count_besancon_norm*(1-scale))
			#cum = cumsum(model)
			#error = max(abs(cum - cumcount_stars/float(max(cumcount_stars))))
			#errors.append(error)
			#grid[i,j,k] = error
		best_index = argmax(logps)
		best_values = self.parameterset.parameter_values[best_index]
		mu, sigma, scale = best_values.dict["mu"], best_values.dict["sigma"], best_values.dict["scale"]
		#scale = 1.
		
		sigmas = sqrt(stars.e_vlos**2 + sigma**2)
		best_model = (gaussian(bin_centers, mu, sigma)*scale + pdf_besancon_discrete*(1-scale)) * self.vdelta
		best_model_gaussian = (gaussian(bin_centers, mu, sigma)) * self.vdelta
		#ps_gaussian = gaussian(stars.vlos, mu, sigmas)
		#ps_basancon = 
		#scale = 1.0
		#ps = ps_gaussian * scale + ps_basancon* (1-scale)
		
		histogramline(bins, cumsum(best_model), color="green", )
		select(0,0)
		histogramline(bins, best_model* len(stars), color="green", )
		vline(mu)
		vline(mu-sigma*3)
		vline(mu+sigma*3)
		select(1,0)
		labels("v (km/s)", "log count")
		vline(mu)
		vline(mu-sigma*3)
		vline(mu+sigma*3)
		nz = count_stars > 0
		histogramline(bins[nz], log10(count_stars[nz]), color="black", )
		histogramline(bins, log10(best_model* len(stars)), color="green", )
		histogramline(bins, log10(best_model_gaussian* len(stars)), color="orange", )
		ylim(-2, 2.5)
		select(2,0)
		nonzero = pdf_besancon_discrete > 0
		print "nz", nonzero.sum(), len(pdf_besancon_discrete)
		#histogramline(bins[nonzero], log10(pdf_besancon_discrete[nonzero]), color="red", )
		#histogramline(bins, (pdf_besancon_discrete), color="red", )
		
		ps = exp(logps - max(logps))
		print ps.max()
		grid = self.parameterset.makegrid(ps)
		print grid.max()
		select(2, 0)
		probgraph(grid, 0, resize=self.parameterset.range1d(0))
		labels(self.parameterset.parameter_range_list[0].label, "prob. dens.")
		select(2, 1)
		probgraph(grid, 1, resize=self.parameterset.range1d(1))
		labels(self.parameterset.parameter_range_list[1].label, "prob. dens.")
		select(2, 2)
		probgraph(grid, 2, resize=self.parameterset.range1d(2))
		labels(self.parameterset.parameter_range_list[2].label, "prob. dens.")
		select(0, 2)
		probimage2d(grid, 0, 1, resize=self.parameterset.range2d(0,1), drawcontourlines=True)
		labels(self.parameterset.parameter_range_list[0].label, self.parameterset.parameter_range_list[1].label)
		select(1, 2)
		probimage2d(grid, 0, 2, resize=self.parameterset.range2d(0,2), drawcontourlines=True)
		labels(self.parameterset.parameter_range_list[0].label, self.parameterset.parameter_range_list[2].label)
		logger.info("writing grid to: %s" % self.filename)
		numpy.save(self.filename, grid)
		print mu, sigma, scale
		draw()
		 
	def load(self):
		self.grid = numpy.load(self.filename)
		return self.grid
		 
class BesanconFitMerge(object):
	def __init__(self, fit1, fit2):
		self.fit1 = fit1
		self.fit2 = fit2
		self.parameterset = self.fit1.parameterset
		
	def run(self, *args):
		g1 = self.fit1.load()
		g2 = self.fit2.load()
		print "max", g1.max()
		print "max", g2.max()
		
		document(size="25cm,25cm")
		mozaic(4,3,box)
		
		select(0, 0)
		probgraph(g1, 0, resize=self.parameterset.range1d(0))
		select(1, 0)
		probgraph(g1, 1, resize=self.parameterset.range1d(1))
		select(2, 0)
		probgraph(g1, 2, resize=self.parameterset.range1d(2))
		
		select(0, 1)
		probgraph(g2, 0, resize=self.parameterset.range1d(0))
		select(1, 1)
		probgraph(g2, 1, resize=self.parameterset.range1d(1))
		select(2, 1)
		probgraph(g2, 2, resize=self.parameterset.range1d(2))
		if 1:
			select(3, 0)
			probimage2d(g1, 0, 1, resize=self.parameterset.range2d(0,1), drawcontourlines=True)
			select(3, 1)
			probimage2d(g2, 0, 1, resize=self.parameterset.range2d(0,1), drawcontourlines=True)
			
		N = g1.shape[-1]
		merge = zeros((g1.shape + (N,)))
		for i in range(N):
			for j in range(N):
				merge[:,:,i,j] = g1[:,:,i] * g2[:,:,j]
		select(3, 2)
		p_mu_sigma = probimage2d(merge, 0, 1, resize=self.parameterset.range2d(0,1), drawcontourlines=True)
		
		select(0, 2)
		probgraph(merge, 0, resize=self.parameterset.range1d(0))
		select(1, 2)
		probgraph(merge, 1, resize=self.parameterset.range1d(1))
		select(2, 2)
		probgraph(merge, 3, resize=self.parameterset.range1d(2))
		
		best_index = argmax(merge)
		
		indices = unravel_index(best_index, merge.shape)
		uniform_values = [k / float(n-1.) for k, n in zip(indices, merge.shape)]
		parameter_range_list = self.parameterset.parameter_range_list + [self.parameterset.parameter_range_list[-1]]
		values = [parameter_range.scale_from_uniform(uniform_value) for  parameter_range,uniform_value in zip(parameter_range_list, uniform_values)]
		
		#self.parameterset.scale_from_uniform
		print values
		
		print indices 
		
		draw()
		
		

class Filter(object):
	
	def __call__(self, stars):
		return self.filter(stars)
	
	
class FilterRmaxMulti(Filter):
	def __init__(self, rmax):
		self.rmax = rmax
		self.rmax.load()
		
	def filter(self, stars):
		logger.info("cutting on radii (elliptical) %s" % (self.rmax.Rmaxs))
		n1 = len(stars)
		def f(star):
			ok = False
			for i in range(self.rmax.simplemodel.catalogues):
				if int(star.catalogue_mask) & (1<<i):
					ok = ok or (star.re < self.rmax.Rmaxs[i])
			return ok
		stars = stars.filter(f)
		n2 = len(stars)
		logger.info("(before n=%d after n=%d)" % (n1, n2))
		return stars

class FilterCatalogue(Filter):
	def __init__(self, catalogue_number):
		self.catalogue_number = catalogue_number
		
	def filter(self, stars):
		logger.info("making a cut on catalogue: number %s" % self.catalogue_number)
		n1 = len(stars)
		stars = stars.filter(lambda star: int(star.catalogue_mask) & (1<<self.catalogue_number))
		n2 = len(stars)
		logger.info("(before n=%d after n=%d)" % (n1, n2))
		return stars

class FilterList(Filter):
	def __init__(self, *filters):
		self.filters = filters
		
	def load(self):
		for filter in self.filters:
			filter.load()
			
	def __call__(self, stars):
		return self.filter(stars)
	
	def filter(self, stars):
		for i, filter in enumerate(self.filters):
			if i == 0:
				logger.info("filter: before: %d stars" % len(stars))
			stars = filter.filter(stars)
			logger.info("filter: after: %d stars" % len(stars))
		return stars

class FilterCutOnR(Filter):
	def __init__(self, light_model, rmax_arcsec):
		self.light_model = light_model
		self.rmax_arcsec = rmax_arcsec
		
	def load(self):
		pass
		
	def filter(self, stars):
		logger.info("making a cut on radius: %f (=%f kpc)" % (self.rmax_arcsec, self.light_model.arcsec_to_kpc(self.rmax_arcsec)))
		n1 = len(stars)
		stars = stars.filter(lambda star: star.rc <= self.rmax_arcsec)
		n2 = len(stars)
		logger.info("(before n=%d after n=%d)" % (n1, n2))
		return stars

class FilterLambda(Filter):
	def __init__(self, f, name="no name"):
		self.f = f
		self.name = name
		
	def load(self):
		pass
		
	def filter(self, stars):
		logger.info("performing a '%s' filter (lambda function)" % self.name)
		n1 = len(stars)
		stars = stars.filter(self.f)
		n2 = len(stars)
		logger.info("(before n=%d after n=%d)" % (n1, n2))
		return stars
	
class FilterCutThreeSigmaHelio(Filter):
	def __init__(self, v, sigma):
		self.v = v
		self.sigma = sigma
		
	def load(self):
		pass
		
	def filter(self, stars, invert=False):
		logger.info("making a three sigma cut: %f (+/ 3*%f)" % (self.v, self.sigma))
		n1 = len(stars)
		if invert:
			stars = stars.filter(lambda star: not(abs(star.vlos_helio - self.v) <= 3*self.sigma))
		else:
			stars = stars.filter(lambda star: abs(star.vlos_helio - self.v) <= 3*self.sigma)
		n2 = len(stars)
		logger.info("(before n=%d after n=%d" % (n1, n2))
		return stars
	
class FilterCutThreeSigmaSys(Filter):
	def __init__(self, v, sigma):
		self.v = v
		self.sigma = sigma
		
	def load(self):
		pass
		
	def filter(self, stars, invert=False):
		logger.info("making a three sigma cut: %f (+/ 3*%f)" % (self.v, self.sigma))
		n1 = len(stars)
		if invert:
			stars = stars.filter(lambda star: not( abs(star.vlos - self.v) <= 3*self.sigma))
		else:
			stars = stars.filter(lambda star: abs(star.vlos - self.v) <= 3*self.sigma)
		n2 = len(stars)
		logger.info("(before n=%d after n=%d" % (n1, n2))
		return stars
	

class ProcessHelioToSys(object):
	def __init__(self, correction, input_observation, output_observation):
		#self.modelpath = modelpath
		self.correction = correction
		self.input_observation = input_observation
		self.output_observation = output_observation
		#self.output = output
		#self.filename = os.path.join(modelpath, self.output)
		#self.vlsr_gsr, self.vsun_lsr, self.vsys_gsr = array(vlsr_gsr), array(vsun_lsr), array(vsys_gsr)
		#self.distance = distance
		
	#def read(self):
	#	self.stars = mab.cvsfile.readcsv(os.path.join(self.modelpath, self.input))
	#	logger.info("read %d heliocentric observations" % len(self.stars))
	 
	def run(self, args, opts, scope):
		self.load()
		self.process()
		self.save()
		
	def load(self):
		self.correction.load()
		self.stars = self.input_observation.load()
		
	def save(self):
		self.output_observation.save(self.stars)
		
	def process(self):
		stars = self.stars
		logger.info("calculating systemic velocities")
		#self.correction.read()
		#self.correction.observation_helio.load()
		#stars = self.correction.observation_helio.stars
		vlos_helio = stars.vlos_helio
		e_vlos = stars.e_vlos
		l, b = mab.astrometry.eq_to_gal(stars.ra, stars.dec)
		e_x = cos(radians(l)) * cos(radians(b))
		e_y = sin(radians(l)) * cos(radians(b))
		e_z = sin(radians(b))
		e_los = transpose(array([e_x, e_y, e_z]))
		
		#v_sun_gsr = (self.correction.observation_helio.vsun_lsr + self.correction.observation_helio.vlsr_gsr)
		v_sun_gsr = self.correction.v_sun_gsr
		
		v_sys_gsr = self.correction.get_vsys_gsr()
		logger.info("systemic velocity of Sculptor (GSR): %r" % v_sys_gsr)
		
		# TODO: double check with paper!
		vlos_gsr = vlos_helio - dot(e_los, v_sys_gsr-v_sun_gsr)
		stars = stars.clone()
		stars.vlos = vlos_gsr
		stars[0].attributes.append("vlos")
		#logger.info("writing to %s" % self.filename)
		#mab.cvsfile.writecsv(self.filename, stars)
		self.stars = stars
		#return stars
		
		#vhelio = self.stars.vlos
		#vsun_gsr = self.vsun_lsr + self.vlsr_gsr
		# calculate the line of sight vectors of the sun towards all the stars
		#l, b = mab.astrometry.eq_to_gal(self.stars.ra, self.stars.dec)
		# TODO: we don't need the distance for this right? we can set it to unit distance
		#x_star_sun, y_star_sun, z_star_sun = mab.astrometry.lbr_to_xyz(l, b, self.distance)
		#ex_los = x_star_sun / self.distance
		#ey_los = y_star_sun / self.distance
		#ez_los = z_star_sun / self.distance
		#self.e_los = transpose(array([ex_los, ey_los, ez_los]))
		#vsun_gsr_los = dot(self.e_los, vsun_gsr)
		#self.vstar_scl = vhelio + vsun_gsr_los
		
	
	#def write(self):
	#	print self.vstar_scl
	
	
class ProcessObservation(object):
	def __init__(self, modelpath, input_observation, output_observation, filter=None):
		self.input_observation = input_observation
		self.filter = filter
		self.output_observation = output_observation
		
	def load(self):
		self.filter.load()
		logger.info("loading input observations")
		self.stars = self.input_observation.load()
	
	def process(self, stars=None):
		if stars is None:
			stars = self.stars
		logger.info("# stars = %d" % len(stars))
		if self.filter:
			self.stars = self.filter.filter(self.stars)
		logger.info("# stars = %d" % len(self.stars))
		return self.stars
		
	def save(self, stars=None):
		logger.info("saving output")
		if stars is None:
			stars = self.stars
		self.output_observation.save(stars)
		

class ObservationHelio(MockObservation):
	def __init__(self, modelpath, input, vlsr_gsr, vsun_lsr):
		self.modelpath = modelpath
		self.input = input
		self.vlsr_gsr, self.vsun_lsr = array(vlsr_gsr), array(vsun_lsr)
		self.filename = os.path.join(self.modelpath, "data", self.input)
		
	#def load(self):
	#	self.stars = mab.cvsfile.readcsv(self.filename)
	#	logger.info("read %d heliocentric observations" % len(self.stars))
	#	return self.stars
		 
class LikelihoodvGSR(object):
	def __init__(self, observation_helio, output, sigma, vxrange, vyrange, vzrange, filter=None, filters=None):
		self.observation_helio = observation_helio
		self.output = output
		self.sigma = sigma
		self.filter = filter
		self.filters = filters
		self.vxrange = vxrange
		self.vyrange = vyrange
		self.vzrange = vzrange
		self.vxs = arange(*self.vxrange)
		self.vys = arange(*self.vyrange)
		self.vzs = arange(*self.vzrange)
		self.filename = os.path.join(self.observation_helio.modelpath, "data", self.output)
		self.v_sun_gsr = (self.observation_helio.vsun_lsr + self.observation_helio.vlsr_gsr)
		
	def load(self):
		self.gridL = load(self.filename)
		
	def get_vsys_gsr(self):
		i = self.gridL.argmax()
		ix, iy, iz = unravel_index(i, self.gridL.shape)
		return array([self.vxs[ix], self.vys[iy], self.vzs[iz]])
	
	def process(self, show_progressbar=True, cores=8):
		self.observation_helio.load()
		stars = self.observation_helio.stars
		filters = []
		if self.filter:
			filters.append(self.filter)
		filters += self.filters
		for filter in filters:
			stars = self.filter(stars)
		vlos_helio = stars.vlos_helio
		e_vlos = stars.e_vlos
		l, b = mab.astrometry.eq_to_gal(stars.ra, stars.dec)
		e_x = cos(radians(l)) * cos(radians(b))
		e_y = sin(radians(l)) * cos(radians(b))
		e_z = sin(radians(b))
		e_los = transpose(array([e_x, e_y, e_z]))
		
		v_sun_gsr = (self.observation_helio.vsun_lsr + self.observation_helio.vlsr_gsr)
		
		def logL(v_sys_gsr):
			vlos_gsr = vlos_helio - dot(e_los, v_sys_gsr-v_sun_gsr)
			return -0.5*sum(vlos_gsr**2/(self.sigma**2+e_vlos**2))
		
		progressbar = mab.utils.progressbar.ProgressBar(0, len(self.vxs)-1)
		v_sys_gsr = zeros(3)
		#self.gridlogL = zeros((len(self.vxs), len(self.vys), len(self.vzs)))
		self.gridlogL = mab.utils.numpy.mmapzeros((len(self.vxs), len(self.vys), len(self.vzs)))
		@mab.parallelize.parallelize(cores=cores, info=show_progressbar)
		def do(i):
			#for i, vx in enumerate(self.vxs):
			if 1:
				vx = self.vxs[i]
				#if show_progressbar:
				#	progressbar.update(i)
				for j, vy in enumerate(self.vys):
					for k, vz in enumerate(self.vzs):
						v_sys_gsr[0] = vx
						v_sys_gsr[1] = vy
						v_sys_gsr[2] = vz
						logLi = logL(v_sys_gsr)
						self.gridlogL[i,j,k] = logLi
						#if logLi > logLmax:
						#	logLmax = logLi
						#	v = array(v_scl_gsr)
						#print logL(v_scl_gsr)
			#print logLmax, v
		do(range(len(self.vxs)))
		self.gridlogL -= self.gridlogL.max()
		self.gridL = exp(self.gridlogL)
	
	def write(self):
		logger.info("saving to: %s" % self.filename)
		save(self.filename, self.gridL)


class SubMean_(object): # depr, see commands
	def __init__(self, modelpath, input_observation, output_observation, vmean=None):
		self.modelpath = modelpath
		self.input_observation = input_observation
		self.output_observation = output_observation
		self.vmean = vmean
		#self.output = output
		#self.filename = os.path.join(modelpath, "data", self.output)
		#self.vlsr_gsr, self.vsun_lsr, self.vsys_gsr = array(vlsr_gsr), array(vsun_lsr), array(vsys_gsr)
		#self.distance = distance
		
	 
	def load(self):
		#self.correction.load()
		self.stars = self.input_observation.load()
		
	def save(self):
		self.output_observation.save(self.stars)
	
	
	#def load(self):
	#	self.observation.read()
	#def read(self):
	#	self.stars = mab.cvsfile.readcsv(os.path.join(self.modelpath, self.input))
	#	logger.info("read %d heliocentric observations" % len(self.stars)) 
		
	def process(self):
		#self.observation.read()
		stars = self.input_observation.stars
		#vlos = stars.vlos
		#e_vlos = stars.e_vlos
		#stars = stars.clone()
		if self.vmean is None:
			vmean = stars.vlos_helio.mean()
		else:
			vmean = self.vmean
		print "vmean", vmean
		stars.vlos = stars.vlos_helio - vmean
		stars[0].attributes.append("vlos")
		logger.info("subtracted mean velocity")
		#logger.info("writing to %s" % self.filename)
		#mab.cvsfile.writecsv(self.filename, stars)
		

def moment(x, k):
	return sum(x**k)/len(x)

def calcmoments(x, ex):
	rawmoments = [moment(x, i) for i in range(4*2+1)]
	#print rawmoments[0]
	N = len(x)
	sigma_noise = mean(ex) # TODO: check what if the data has no heteroscedasticity
	moments = zeros(4*2+1)
	moments[0] = rawmoments[0]
	moments[1] = rawmoments[1]
	moments[2] = rawmoments[2] - sigma_noise**2
	moments[3] = rawmoments[3] - 3*moments[1]*sigma_noise**2
	moments[4] = rawmoments[4] - 3*sigma_noise**4 - 6*moments[2]*sigma_noise**2 
	moments[5] = rawmoments[5] - 10*moments[3]*sigma_noise**2 - 5*moments[1]*3*sigma_noise**4 
	moments[6] = rawmoments[6] - 15*sigma_noise**6 - 15*moments[2]*3*sigma_noise**4 - 15*moments[4]*sigma_noise**2
	moments[7] = rawmoments[7] - 35*moments[3]*3*sigma_noise**4 - 21*moments[5]*sigma_noise**2 - 7*moments[1]*15*sigma_noise**6 
	moments[8] = rawmoments[8] - 105*sigma_noise**8 - 70*moments[2]*3*sigma_noise**4 - 28*moments[6]*sigma_noise**2 - 28*moments[2]*15*sigma_noise**6
	exps = [] # expectation values
	vars = [] # variances
	for k in [1,2,3,4]:
		exps.append(moments[k])
		if k == 3:
			varm1 = (rawmoments[1*2]-rawmoments[1]**2)/N
			vars.append((rawmoments[k*2]-rawmoments[k]**2)/N + 9*sigma_noise**4*varm1)
		elif k == 4:
			varm2 = (rawmoments[2*2]-rawmoments[2]**2)/N
			vars.append((rawmoments[k*2]-rawmoments[k]**2)/N + 36*sigma_noise**4*varm2)
		else:
			#exps.append(moment(samples,k))
			vars.append((rawmoments[k*2]-rawmoments[k]**2)/N)
			
	return exps, (vars)

class IdealizeMoments(object):
	def __init__(self, binned_data, galaxy):
		self.binned_data = binned_data
		self.galaxy = galaxy
		
	def run(self, args, opts, scope):
		self.binned_data.load()
		self.binned_data.aperture.load()
		
		jeans = self.galaxy.jeans()
		r = self.binned_data.aperture.aperture_rcenters_kpc
		sigma_los = jeans.sigma_los(r)
		m4_los = jeans.m4_los(r)
		self.binned_data.moments[2] = sigma_los**2
		self.binned_data.moments[3] *= 0
		self.binned_data.moments[4] = m4_los
		self.binned_data.save()
		
		
class BinningSimple(object):
	def __init__(self, modelpath, aperture, observation, filter=None, filters=None, filter_aperture=None, mean_v=0, postfix=""):
		self.modelpath = modelpath
		self.observation = observation
		self.aperture = aperture
		self.postfix = postfix
		self.filter = filter
		self.filters = filters
		self.filter_aperture = filter_aperture
		self.logger = logger = logging.getLogger("gd.schw.obs.BinningSimple")
		self.mean_v = mean_v
		
		
	def get_stars_from_aperture(self, index):
		stars = self.observation.stars
		if self.filter:
			stars = self.filter.filter(stars)
		if self.filters:
			for filter in self.filters:
				stars = filter(stars)
		stars = stars.filter(self.aperture.aperture_filter(index))
		#aperture = self.aperture.aperture
		if self.filter_aperture is None:
			return stars
		else:
			aperture_stars = self.filter_aperture.filter_aperture(stars, index)
		return aperture_stars
		
	
	def prepare(self, cores=48, show_progressbar=True):
		self.observation.load()
		self.aperture.load()
		aperture = self.aperture.aperture
		n_aperture = aperture.length()
		
		#print self.aperture
		#dsa
		
		logger.info("maximum aperture index/length: %d/%d" % (n_aperture-1, n_aperture))
		
		self.aperture_index_to_constraint = zeros(n_aperture, dtype=int)-1
		aperture_stars = [list() for k in range(n_aperture)]
		aperture_indices = [list() for k in range(n_aperture)]
		aperture_velocities = [list() for k in range(n_aperture)]
		aperture_velocity_errors = [list() for k in range(n_aperture)] 
		
		stars = self.observation.stars
		if self.filter:
			stars = self.filter.filter(stars)
		if self.filters:
			for filter in self.filters:
				stars = filter(stars)
		x = stars.xi
		y = stars.eta
		#x = stars.re * 60 * 60
		#y = x * 0
		#angles = (arctan2(y, x) + 2 * pi) % (2*pi)
		r = sqrt(x*x+y*y)
		
		vlos = stars.vlos
		e_vlos = stars.e_vlos
		#print len(vlos)
		
		# look up in which aperture the star belongs, and add the velocity and its error to that specific list
		for i in range(len(x)):
			#print r[i], stars.rc[i]
			if not aperture.inrange(x[i], y[i]):
				#msg = "star :%d not in aperture (has p_member value of %.3f, r=%f)" % (i, stars[i].p_member, r[i])
				msg = "star :%d not in aperture (r=%f)" % (i, r[i])
				self.logger.info(msg)
				#sys.exit(-1)
				#else: # if not, it's and error
				#	self.logger.info(msg)
			aperture_index = aperture.findindex(x[i], y[i])
			#print aperture_index,
			aperture_stars[aperture_index].append(stars[i])
			aperture_indices[aperture_index].append(i)
			aperture_velocities[aperture_index].append(vlos[i])
			aperture_velocity_errors[aperture_index].append(e_vlos[i])
		constraint = 0
		# now count the stars in each aperture, and give each non empty aperture and constraint number
		for i in range(n_aperture):
			n = len(aperture_velocities[i])
			if n > 0:
				self.aperture_index_to_constraint[i] = constraint
				constraint += 1
			else:
				self.logger.warning("warning: aperture index %d is empty" % i)
		self.n_constraints = constraint 
		self.moments = zeros((5, self.n_constraints)) 
		self.e_moments = zeros((5, self.n_constraints)) 
		constraint = 0
		
		self.moments, self.e_moments = self.calculate_moments(aperture_velocities, aperture_velocity_errors, stars, aperture_indices)
		#self.moments[2] -= self.moments[1]**2
		#print std(aperture_velocities[0]), self.moments[2][0]**0.5
		self.logger.info("average velocity dispersion: %.2f" % mean(self.moments[2]**0.5))
		
		
		self.logger.info("number of bins: %d" % self.moments.shape[1])
		#print mean(self.moments[2]**0.5)
		for i in [2,4]:
			S = mean(self.moments[i])
			N = mean(self.e_moments[i])
			SNR = S/N
			self.logger.info("moment %i SNR: %f (%f/%f)" % (i, SNR, S, N))
		
		sigma = sqrt(self.moments[2])
		#print "sigma", sigma
		#error = sqrt(self.e_moments[2]**2/(2*sigma)**2)
		error = 1/(2*sigma) * self.e_moments[2]
		S = mean(sigma)
		N = mean(error)
		SNR = S/N
		self.logger.info("sigma SNR: %f (%f/%f)" % (SNR, S, N))
		
		m2 = self.moments[2]
		m4 = self.moments[4]
		e_m2 = self.e_moments[2]
		e_m4 = self.e_moments[4]
		kappa = m4/m2**2
		e_kappa = sqrt(1/m2**4 * e_m4**2 + (2 * m4/m2**4)**2*e_m2**2)
		S = mean(kappa)
		N = mean(e_kappa)
		SNR = S/N
		self.logger.info("kappa SNR: %f (%f/%f)" % (SNR, S, N))
		
	def calculate_moments(self, aperture_velocities, aperture_velocity_errors, stars, aperture_indices):
		n_aperture = len(aperture_velocities)
		moments_constraints = zeros((5, self.n_constraints))
		e_moments_constraints = zeros((5, self.n_constraints))
		constraint = 0 
		print "filter", self.filter_aperture
		for i in range(n_aperture):
			n = len(aperture_velocities[i])
			#sda
			if self.filter_aperture:
				aperture_stars = self.filter_aperture.filter_aperture(stars, i)
			else:
				#raise "old"
				aperture_stars = None
			if n > 0:
				if aperture_stars is None:
					vlos = array(aperture_velocities[i])
					e_vlos = array(aperture_velocity_errors[i])
				else:
					vlos = aperture_stars.vlos
					e_vlos = aperture_stars.e_vlos
				print n,  "stars in aperture", i
				#vlos = array(aperture_velocities[i]) - self.mean_v
				#e_vlos = array(aperture_velocity_errors[i])
				moments, errors = calcmoments(vlos, e_vlos)
				#print std(vlos), moments
				#print i
				#print aperture_velocities[i]
				#print "sigma", sqrt(moments[1])
				if i == 0:
					#print "vlosses", list(vlos)
					#print moments, errors
					pass
				#self.moments[0, constraint] = n
				#self.e_moments[0, constraint] = sqrt(n)
				#print constraint, moments_constraints.shape 
				moments_constraints[0, constraint] = self.aperture.mass(i)
				e_moments_constraints[0, constraint] = self.aperture.mass(i)
				moments_constraints[1:5, constraint] = moments
				#moments_constraints[2, constraint] = var(vlos)
				import pdb;
				#pdb.set_trace()
				moments_constraints[2, constraint] = moments_constraints[2,constraint] - moments_constraints[1,constraint]**2
				#v = vlos - mean(vlos)
				#w = 1/(0**2+e_vlos**2) * 0 + 1
				#w = 1/(0**2+e_vlos**2)* 0 + 1
				#w = v * 0 + 1.
				#self.moments[2, constraint] = sum(v**2*w)/sum(w)#std(vlos)**2
				#self.moments[2, constraint] = sigmas_ml[i]**2
				#print self.moments[2, constraint]
				#self.moments[2, constraint] = std(vlos)**2
				e_moments_constraints[1:5, constraint] = sqrt(errors) # go from var to std
				#self.e_moments[2, constraint] = sqrt(2*sqrt(self.moments[2, constraint]) * sigmas_errors[i]**2)
				#print n, moments#, vlos
				constraint += 1
		return moments_constraints, e_moments_constraints
		
	def save(self):
		self.store()			
	def store(self):
		dirname = os.path.join(self.modelpath, "schw", "aperture")
		
		filename = os.path.join(dirname, "observation_moments" +self.postfix +".npy")
		moments = array([self.moments, self.e_moments])
		self.logger.info("saving moments to file: %s" % filename)
		numpy.save(filename, moments)
		
		filename = os.path.join(dirname, "aperture_index_to_constraint" +self.postfix +".npy")
		self.logger.info("saving aperture index to constraint array to file: %s" % filename)
		numpy.save(filename, self.aperture_index_to_constraint)
		
		filename = os.path.join(dirname, "aperture_index_to_constraint" +self.postfix +".asc")
		self.logger.info("saving aperture index to constraint array to file: %s (ascii version)" % filename)
		f = open(filename, "w")
		f.write(repr(list(self.aperture_index_to_constraint)))
		f.close()
		#self.logger.info("saving observations to disk again since they now include p_member values")
		#self.observation.save()
	
	def load(self):
		dirname = os.path.join(self.modelpath, "schw", "aperture")
		filename = os.path.join(dirname, "observation_moments" +self.postfix +".npy")
		moments = numpy.load(filename)
		self.moments, self.e_moments = moments
		
		filename = os.path.join(dirname, "aperture_index_to_constraint" +self.postfix +".npy")
		self.aperture_index_to_constraint = numpy.load(filename)
		self.n_constraints = self.aperture_index_to_constraint.max() + 1
		 
		 
class BinningWithMembership(BinningSimple):
	def __init__(self, modelpath, aperture, observation, membership, postfix=""):
		BinningSimple.__init__(self, modelpath, aperture, observation, postfix=postfix)
		self.logger = logger = logging.getLogger("gd.schw.obs.BinningWithMembership")
		self.membership = membership
	
	def prepare(self, cores=48, show_progressbar=True):
		self.observation.load()
		self.aperture.load()
		aperture = self.aperture.aperture
		n_aperture = aperture.length()
		
		logger.debug("maximum aperture index/length: %d/%d" % (n_aperture-1, n_aperture))
		
		self.aperture_index_to_constraint = zeros(n_aperture, dtype=int)-1
		aperture_stars = [list() for k in range(n_aperture)]
		aperture_indices = [list() for k in range(n_aperture)]
		aperture_velocities = [list() for k in range(n_aperture)]
		aperture_velocity_errors = [list() for k in range(n_aperture)] 
		
		stars = self.observation.stars
		x = stars.xi
		y = stars.eta
		#x = stars.re * 60 * 60
		#y = x * 0
		#angles = (arctan2(y, x) + 2 * pi) % (2*pi)
		r = sqrt(x*x+y*y)
		
		vlos = self.observation.stars.vlos_helio
		e_vlos = self.observation.stars.e_vlos
		
		# look up in which aperture the star belongs, and add the velocity and its error to that specific list
		for i in range(len(x)):
			#print r[i], stars.rc[i]
			if not aperture.inrange(x[i], y[i]):
				#msg = "star :%d not in aperture (has p_member value of %.3f, r=%f)" % (i, stars[i].p_member, r[i])
				msg = "star :%d not in aperture (r=%f)" % (i, r[i])
				self.logger.info(msg)
				#sys.exit(-1)
				#else: # if not, it's and error
				#	self.logger.info(msg)
			else:
				aperture_index = aperture.findindex(x[i], y[i])
				aperture_stars[aperture_index].append(stars[i])
				aperture_indices[aperture_index].append(i)
				aperture_velocities[aperture_index].append(vlos[i])
				aperture_velocity_errors[aperture_index].append(e_vlos[i])
		constraint = 0
		# now count the stars in each aperture, and give each non empty aperture and constraint number
		for i in range(n_aperture):
			n = len(aperture_velocities[i])
			if n > 0:
				self.aperture_index_to_constraint[i] = constraint
				constraint += 1
			else:
				self.logger.warning("warning: aperture index %d is empty" % i)
		self.n_constraints = constraint 
		self.moments = zeros((5, self.n_constraints)) 
		self.e_moments = zeros((5, self.n_constraints)) 
		constraint = 0
		#print aperture_velocities[-1]
		#print aperture_velocity_errors[-1]
		if 1:
			document(size="35cm,15cm")
			sigmas = arange(1, 40,0.25)
			#box()
			vsplit(box,3)
			sigmas_ml = []
			sigmas_errors = []
			vs_ml = []
			ai = -1 #-8
			vsys = 110.6
			vsigma = 10.1
			v0 = self.membership.member_distribution.v
			def member(star):
				return (star.vlos_helio < (vsys+3*vsigma)) & (star.vlos_helio > (vsys-3*vsigma))
			def non_member(star):
				return not member(star)
			@mab.parallelize.parallelize(cores=cores, info=show_progressbar)
			def do(i):
				#vs_helio = array([star.vlos_helio for star in aperture_stars[i]])
				#selection = (vs_helio < (vsys+3*vsigma)) & (vs_helio > (vsys-3*vsigma))
				
				#logps = array([self.membership.logp_sigma(array(aperture_stars[i])[selection], sigma) for sigma in sigmas])
				maxlogp = -1e1000
				logps = None 
				vbest = v0
				Ntot = len(aperture_stars[i])
				Nmw = sum([1 for star in aperture_stars[i] if non_member(star)])
				Ngal = Ntot-Nmw
				for v in (v0+arange(-1.5, 1.5, 0.25)):
				#for v in [v0]
					newlogps = array([self.membership.logp_sigma(aperture_stars[i], sigma, v, Ngal, Nmw) for sigma in sigmas])
					if newlogps.max() > maxlogp:
						logps = newlogps
						maxlogp = logps.max()
						vbest = v
					#logpss.append(logps)
				#print vbest, len(logps)
				logps -= max(logps)
				p = exp(logps)
				p /= p.max()
				#print p
				v_sigma = sigmas[argmax(p)]
				v_mean = vbest
				#sigmas_ml.append(v_sigma)
				#vs_ml.append(vbest)
				sigma_mean = sum(sigmas * p) / sum(p)
				sigma_error = sqrt(sum((sigmas-sigma_mean)**2 * p) / sum(p))
				#print logps, p
				#print "sigma", v_sigma, "+/-", sigma_error
				sigmas_errors.append(sigma_error)
				graph(sigmas,p) # + len(aperture_stars[i]))
				p_member_list = []
				for star in aperture_stars[i]:
					#star.p_member = self.membership.p_member(star, v_mean, v_sigma)
					p_member_list.append(self.membership.p_member(star, v_mean, v_sigma, Ngal, Nmw))
					if "p_member" not in star.attributes: 
						star.attributes.append("p_member")
				return (v_mean, v_sigma, sigma_error, p_member_list,newlogps)
			results = do(range(len(aperture_stars)))
			sigmas_ml = [v_sigma for (v_mean, v_sigma, sigma_error, p_member_list, newlogps) in results]
			vs_ml = [v_mean for (v_mean, v_sigma, sigma_error, p_member_lis, newlogpst) in results]
			sigmas_errors = [sigma_error for (v_mean, v_sigma, sigma_error, p_member_list, newlogps) in results]
			p_member_lists = [p_member_list for (v_mean, v_sigma, sigma_error, p_member_list, newlogps) in results]
			newlogpss =  [newlogps for (v_mean, v_sigma, sigma_error, p_member_list, newlogps) in results]
			for i in range(len(aperture_stars)):
				for star, p_member in zip(aperture_stars[i], p_member_lists[i]):
					star.p_member = p_member
					if "p_member" not in star.attributes: 
						star.attributes.append("p_member")
			
			print "r1)", self.aperture.aperture_rborders[ai]
			print "r1)",self.aperture.aperture_rborders[ai-1]
			print [star.vlos_helio for star in aperture_stars[ai]]
			print [star.e_vlos for star in aperture_stars[ai]]
			print len(aperture_stars[ai])
			print self.aperture.aperture_rborders
			#vline(10)
			select(0)
			graph(sigmas, exp(newlogpss[ai]))
			select(1)
			if 1:
				vs = arange(-50, 350)
				r = aperture_stars[ai][0].re
				ngalaxy = self.membership.density_ratio(r)
				ngalactic = 1-ngalaxy
				self.membership.member_distribution.sigma = sigmas_ml[ai]
				self.membership.member_distribution.v = vs_ml[ai]
				
				Ntot = len(aperture_stars[ai])
				Nmw = sum([1 for star in aperture_stars[ai] if non_member(star)])
				Ngal = Ntot-Nmw
				
				N_tot = Ngal + Nmw
				f_galaxy = Ngal / float(N_tot)
				f_mw = Nmw / float(N_tot)
				
				p1 = array([self.membership.member_distribution(v, 0.) * f_galaxy for v in vs])
				p2 = array([self.membership.galatic_distribution(v, 0.) * f_mw for v in vs])
				
				vs_helio = array([star.vlos_helio for star in aperture_stars[ai]])
				histogram(vs_helio, datamin=-50, datamax=350, binwidth=4, normalize=True, differential=True)
				graph(vs, p1)
				graph(vs, p2, color="red")
				graph(vs, p1+p2, color="blue", linewidth="2pt")
				
				if 1:
					selection = (vs_helio < (vsys+3*vsigma)) & (vs_helio > (vsys-3*vsigma))
					print "corr", std(vs_helio[selection]), std(vs_helio), sigmas_ml[ai]
					print len(vs_helio[selection]), len(vs_helio)
					histogram(vs_helio[selection], datamin=-50, datamax=350, binwidth=4, normalize=True, differential=True, color="red", linewidth="2pt")
					self.membership.member_distribution.sigma = std(vs_helio[selection])
					print "sigmas", self.membership.member_distribution.sigma, sigmas_ml[ai]
					p1 = array([self.membership.member_distribution(v, 0.) for v in vs])
					graph(vs, p1, color="green", linestyle="dot")
			select(2)
			#print sigmas_ml
			print std(vs_helio)
			hline(std(vs_helio))
			hline(sigmas_ml[ai], color="red")
			graph(sigmas_ml)
			errorbars(arange(len(sigmas_ml)), sigmas_ml, yerr=sigmas_errors)
			ylim(0, 20)
			
			draw()
		
		
		for i in range(n_aperture):
			n = len(aperture_velocities[i])
			if n > 0:
				vlos = array(aperture_velocities[i])
				e_vlos = array(aperture_velocity_errors[i])
				moments, errors = calcmoments(vlos, e_vlos)
				#self.moments[0, constraint] = n
				#self.e_moments[0, constraint] = sqrt(n)
				self.moments[0, constraint] = self.aperture.mass(i)
				self.e_moments[0, constraint] = self.aperture.mass(i)
				self.moments[1:5, constraint] = moments
				v = vlos - mean(vlos)
				#w = 1/(0**2+e_vlos**2) * 0 + 1
				#w = 1/(0**2+e_vlos**2)* 0 + 1
				w = v * 0 + 1.
				self.moments[2, constraint] = sum(v**2*w)/sum(w)#std(vlos)**2
				self.moments[2, constraint] = sigmas_ml[i]**2
				self.moments[1, constraint] = vs_ml[i]
				#print self.moments[2, constraint]
				#self.moments[2, constraint] = std(vlos)**2
				self.e_moments[1:5, constraint] = sqrt(errors) # go from var to std
				#self.e_moments[2, constraint] = sqrt(2*sqrt(self.moments[2, constraint]) * sigmas_errors[i]**2)
				self.e_moments[2, constraint] = sqrt(4*self.moments[2, constraint] * sigmas_errors[i]**2)				
				#print n, moments#, vlos
				constraint += 1


class CommandBinData(object):
	def __init__(self, binned_data):
		self.binned_data = binned_data
		
	def run(self, args, opts, scope):
		self.binned_data.prepare()
		self.binned_data.store()
		
		
class BinningWithForeground(BinningSimple):
	def __init__(self, modelpath, aperture, observation, foreground, fractions, mean_v, filters, postfix=""):
		BinningSimple.__init__(self, modelpath, aperture, observation, filters=filters, postfix=postfix)
		self.foreground = foreground
		self.fractions = fractions
		self.mean_v = mean_v
		
	def calculate_moments(self, aperture_velocities, aperture_velocity_errors, stars, aperture_indices):
		n_aperture = len(aperture_velocities)
		moments_constraints = zeros((5, self.n_constraints))
		e_moments_constraints = zeros((5, self.n_constraints))
		constraint = 0 
		no_catalogues = len(self.fractions)
		borders = self.aperture.aperture_rborders
		print borders
		plot = False
		if plot:
			#box()
			mozaic(2,2,box)
		test = []
		Ntot = len(stars)
		if 0:
			for i in range(n_aperture):
				n = len(aperture_velocities[i])
				Nmem_i = Nmem * self.aperture.radial_surface_densities[i] / sum(self.aperture.radial_surface_densities)
				Nfg_i = Nfg * (borders[i+1]**2 - borders[i]**2)/borders[-1]**2
				#area_fraction = A_fg / A_gal
				test.append(Nfg_i)
			print test
			print sum(test), Ntot, Nfg, Nmem
		#if 1:
		for i in range(n_aperture):
			n = len(aperture_velocities[i])
			#print area_fraction
			if n > 0:
				#print aperture_indices
				vlos = array(stars[aperture_indices[i]].vlos)
				e_vlos = array(stars[aperture_indices[i]].e_vlos)
				catalogue_masks = stars.catalogue_mask[aperture_indices[i]] 
				sigmas = arange(0, 30, 0.25)
				logps_sigma = []
				print len(vlos), mean(array(stars[aperture_indices[i]].rc)), 
				for sigma in sigmas:
					#for v, e in zip(vlos, e_vlos):
					logps_sigma.append(0)
					for j in range(len(vlos)):
						p_gaussian = gaussian(vlos[j], self.mean_v, sqrt(e_vlos[j]**2+sigma**2))
						ps_foreground = []
						ps = []
						for k in range(no_catalogues):
							if ((1<<k) & int(catalogue_masks[j])):
								#ps_foreground.append(self.foreground.probability(v))
								fraction = self.fractions[k]
								fraction = 0.9
								Nmem = Ntot * self.fractions[k]
								Nfg = Ntot - Nmem
								Nmem_i = Nmem * self.aperture.radial_surface_densities[i] / sum(self.aperture.radial_surface_densities)
								Nfg_i = Nfg * (borders[i+1]**2 - borders[i]**2)/borders[-1]**2
								fraction = Nmem_i / (Nmem_i + Nfg_i)
								#print fraction
								ps.append(p_gaussian * fraction + (1-fraction) * self.foreground.probability(vlos[j]))
						logps_sigma[-1] += log(mean(ps))
				ps = 10**(logps_sigma-max(logps_sigma))
				ps /= sum(ps)
				mean_sigma = sum(ps*sigmas)
				std_sigma = sqrt(sum(ps*(sigmas-mean_sigma)**2))
				index = argmax(logps_sigma)
				sigma = sigmas[index]
				 
				print i, sigma, mean_sigma, std_sigma
				group = [18,23]#,70,74]
				if i in group and plot:
					
					color = nicecolors[group.index(i)]
					select(0, 0)
					#graph(sigmas, logps_sigma)
					select(1,1)
					histogram(vlos, datamin=-300, datamax=200, binwidth=4, color=color, normalize=False)
					select(1,0)
					histogram(vlos, datamin=-300, datamax=200, binwidth=4, color=color, normalize=True)
					vs = arange(-200, 100)
					p = gaussian(vs, self.mean_v, sigma)
					graph(vs, p, color=color)
				#moments, errors = calcmoments(vlos, e_vlos)
				moments_constraints[0, constraint] = self.aperture.mass(i)
				e_moments_constraints[0, constraint] = self.aperture.mass(i)
				moments_constraints[2, constraint] = sigma**2
				e_moments_constraints[2, constraint] = 2*sigma*std_sigma
				#moments_constraints[1:5, constraint] = moments
				#e_moments_constraints[1:5, constraint] = sqrt(errors) # go from var to std
				constraint += 1
		if plot:
			draw()
		return moments_constraints, e_moments_constraints
		
		
	def prepare(self):
		self.foreground.load()
		BinningSimple.prepare(self)
		
	def run(self, *args):
		self.calculate()
		