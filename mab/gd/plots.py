# -*- coding: utf-8 -*-
from kaplot import *
import scipy.ndimage
from mab.binningtools import bingrid, binrange
import mab.gd.logging as logging
logger = logging.getLogger("gd.plot")
import mab.asciifile
import numpy as np

class Sigmar(object):
	def __init__(self, galaxy):
		self.galaxy = galaxy
		
	def run(self, args, opts, scope):
		logr = arange(-4, 1, 0.1)
		r = 10**logr
		jeans = self.galaxy.jeans()
		sigmar = jeans.sigmar(r)
		graph(logr, log10(sigmar), color="red", linestyle="dash")
		labels("log r/kpc", "log vel disp (km/s)")
		draw()
		
class PlotStars2D(object):
	def __init__(self, observation, filters=[]):
		self.observation = observation
		self.filters = filters
		
	def load(self):
		self.allstars = self.observation.load()
		self.stars = self.allstars
		for filter in self.filters:
			self.stars = filter(self.stars)
			
			
	def run(self, args, opts, scope):
		self.load()
		self.plot(scope)

class Vlos(object):
	def __init__(self, observation, galatic_distribution, vmin, vmax, vmean, vsigma, fraction):
		self.observation = observation
		self.galatic_distribution = galatic_distribution
		self.vmin = vmin
		self.vmax = vmax
		self.vmean = vmean
		self.vsigma = vsigma
		self.fraction = fraction
		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		
		x = mgrid[self.vmin:self.vmax:200j]
		box()
		y = gaussian(x, self.vmean, self.vsigma) * (1.-self.fraction)
		graph(x, self.galatic_distribution(x, x*0) * self.fraction + y)
		print self.fraction
		histogram(stars.vlos, datamin=self.vmin, datamax=self.vmax, binwidth=1., normalize=True)
		draw()
		
		
		

class CMD(object):
	def __init__(self, photometry, filters):
		self.photometry = photometry
		self.filters = filters
		
	def run(self, args, opts, scope):
		stars = self.photometry.load()
		print len(stars)
		for filter in self.filters:
			stars = filter(stars)
		
		box()
		scatter(stars.V - stars.I, stars.I)
		#flipx()
		flipy()
		labels("V-I", "I")
		draw()
		
	
		
		
		

class PlotStars2DKDE(PlotStars2D):
	def __init__(self, observation, filters=[], xlim=None, ylim=None, paper=False):
		super(PlotStars2DKDE, self).__init__(observation, filters)
		self.xlim = xlim
		self.ylim = ylim
		self.paper = paper
		
	def plot2d(self, x, y, scope):
		#box()
		if self.paper:
			page(fontsize="20pt")
			mozaic(1,1, box)
		else:
			mozaic(3,2, box)
		#x, y = self.stars[p1], self.stars[p2]
		import scipy.stats.kde
		
		y = y - numpy.mean(y)
		
		data = array([x, y])
		
		kde = scipy.stats.kde.gaussian_kde(data)
		#print res
		#import pdb
		#pdb.set_trace()
		
		#xmin, xmax = min(x), max(x)
		#ymin, ymax = min(y), max(y)
		xmin, xmax = self.xlim
		ymin, ymax = self.ylim
		xg, yg = mgrid[xmin:xmax:100j, ymin:ymax:100j]
		xrange, yrange = mgrid[xmin:xmax:100j], mgrid[ymin:ymax:100j]
		dens = kde([xg.ravel(), yg.ravel()]).reshape(xg.shape)

		select(0, 0)
		#indexedimage(dens.T, resize=[(xmin, ymin), (xmax, ymax)])
		dens2 = dens# / dens.sum(axis=0)
		if 1:
			for i in range(dens2.shape[0]):
				dens2[i] /= dens2[i].max()
			
		resize = (xmin, ymin), (xmax, ymax)
		indexedimage(dens2.T, colormap="rainbow", resize=resize)
		contour(dens.T, levels=10, resize=resize, color="red")
		scatter(x, y, symbolsize="30pt", color="black")
		if self.xlim:
			xlim(*self.xlim)
		if self.ylim:
			ylim(*self.ylim)
			
		if not self.paper:
				
			dens2 = dens# / dens.sum(axis=0)
			for i in range(dens.shape[0]):
				dens2[i] /= dens[i].max()
			#ylim(-3, 0)
			select(1, 0)
			indexedimage(dens2.T, resize=[(xmin, ymin), (xmax, ymax)])
			contour(dens2.T, levels=10, resize=[(xmin, ymin), (xmax, ymax)])
			scatter(x, y)
			if self.xlim:
				xlim(*self.xlim)
			if self.ylim:
				ylim(*self.ylim)
			#ylim(-3, 0)
			select(0, 1)
			dens2 = dens
			for i in range(dens.shape[1]):
				dens2[:,i] /= dens[:,i].max()
			indexedimage(dens2.T, resize=[(xmin, ymin), (xmax, ymax)])
			scatter(x, y)
			if self.xlim:
				xlim(*self.xlim)
			if self.ylim:
				ylim(*self.ylim)
			#ylim(-3, 0)
			select(1,1)
			cumhistogram(x, normalize=True)
			if self.xlim:
				xlim(*self.xlim)
				
			select(2,0)
			
			dens2 = dens# / dens.sum(axis=0)
			sigmas = []
			means = []
			for i in range(dens.shape[0]):
				#dens2[i] /= dens[i].max()
				p = dens2[i]/sum(dens2[i])
				
				mean = numpy.sum(yrange*p)
				print mean
				var = numpy.sum((yrange-mean)**2 *p)
				sigmas.append(var**0.5)
				means.append(mean)
			graph(xrange, sigmas)
			#ylim(0, 30)
			if self.xlim:
				xlim(*self.xlim)
			select(2,1)
			graph(xrange, means)
			if self.xlim:
				xlim(*self.xlim)
			#ylim(0, 30)
			
		
		
		#labels("rc", "Fe/H")
		#draw()		
		
class PlotFeHvsR(PlotStars2DKDE):
	def __init__(self, observation, filters=[]):
		super(PlotStars2DKDE, self).__init__(observation, filters)
		
	def plot(self, scope):
		self.plot2d(self.stars.rc, self.stars.FeH, scope)
		select(0, 0)
		labels("rc", "Fe/H")
		draw()
		
class PlotXYKDE(PlotStars2DKDE):
	def __init__(self, observation, property1, property2, filters=[], xlim=None, ylim=None, paper=False, label1=None, label2=None):
		super(PlotXYKDE, self).__init__(observation, filters, xlim=xlim, ylim=ylim, paper=paper)
		self.property1, self.property2 = property1, property2
		self.label1 = self.property1 if label1 is None else label1
		self.label2 = self.property2 if label2 is None else label2
		
	def plot(self, scope):
		print getattr(self.stars, self.property1)
		self.plot2d(getattr(self.stars, self.property1), getattr(self.stars, self.property2), scope)
		select(0, 0)
		#labels(self.property1, self.property2)
		labels(self.label1, self.label2)
		draw()		

class Test(object):
	def __init__(self, observation, light_profile, aperture, filter_velocity):
		self.observation = observation
		self.light_profile = light_profile
		self.aperture = aperture
		self.filter_velocity = filter_velocity
		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		member_stars = self.filter_velocity(stars)
		nomember_stars = [star for star in stars if star not in member_stars]
		logger.info("member: %d non-member: %d total: %d" % (len(member_stars), len(nomember_stars), len(stars)))
		
		#self.light_profile.load()
		self.aperture.load()
		
		def plot(stars, **kwargs):
			r = [star.rc for star in stars]
			vlos = [star.vlos_helio for star in stars]
			scatter(vlos, r, **kwargs)
		#box()
		mozaic(3,2,box)
		select(0, 0)
		plot(member_stars, color="black")
		plot(nomember_stars, color="red")
		xlim(-100, 200)
		
		
		select(1, 0)
		borders = self.aperture.aperture_rborders
		borders_kpc = self.aperture.aperture_rborders_kpc
		scale = 1.
		counts = []
		exected_counts = []
		Rmax = borders_kpc[-1]
		normalization = self.light_profile.cumdensityR(0, Rmax, 1.)
		print "light normalization", normalization
		#for r1, r2 in zip(borders[:-1], borders[1:]):
		self.light_profile.b 
		for i in range(len(borders)-1):
			r1 = borders[i]
			r2 = borders[i+1]
			stars = [star for star in member_stars if ((star.rc >= r1) and (star.rc < r2))]
			stars = [star for star in member_stars if ((star.re >= r1) and (star.re < r2))]
			#stars = [star for star in member_stars if ((star.re_gius >= r1) and (star.re_gius < r2))]
			counts.append(len(stars))
			r1 = borders_kpc[i]
			r2 = borders_kpc[i+1]
			exected_counts.append(self.light_profile.cumdensityR(r1, r2, 1.)/normalization * len(member_stars))
		#do(member_stars, color="black")
		counts = array(counts)
		exected_counts = array(exected_counts)
		#graph(counts * scale, color="black")
		#graph(exected_counts * scale, color="purple")
		graph(counts/exected_counts * scale, color="purple")
		hline(1)
		rescale = counts/exected_counts
		rescale /= mean(rescale)
		print "check norm", mean(rescale)
		print self.light_profile.b
		#select(1, 1)
		counts = []
		for r1, r2 in zip(borders[:-1], borders[1:]):
			stars = [star for star in nomember_stars if ((star.rc >= r1) and (star.rc < r2))]
			counts.append(len(stars))
		counts = array(counts)
		#graph(counts * scale, color="red")
		#do(nomember_stars, scale=1, color="red")
		
		areas = []
		scale = len(nomember_stars)
		for r1, r2 in zip(borders_kpc[:-1], borders_kpc[1:]):
			#counts.append(self.light_profile.cumdensityR(r1, r2, 1.))
			area = (r2**2-r1**2) / Rmax**2
			areas.append(area)
		logger.info("total: %f (%f)" % (sum(counts), sum(counts)*scale))
		areas = array(areas)
		#graph(areas * scale, color="green")
		#graph(areas * scale * rescale, color="orange")
		
		print counts.shape, areas.shape
		graph(counts / (areas*scale), color="blue")
		hline(1)
		
		select(2,0)
		graph(counts *rescale/areas, color="orange")
		
		draw()
		
	

class Profile1d(object):
	def __init__(self, profile):
		self.profile = profile
		
	def run(self, args, opts, scope):
		#box()
		document(size="20cm,20cm")
		#print self.profile.rs
		mozaic(2,2,box)
		
		logrs = arange(-2, log(200), 0.1)
		select(0, 0)
		self.plot_density(logrs)
		select(1, 0)
		self.plot_enclosed_mass(logrs)
		select(0, 1)
		self.plot_logslope(logrs)
		ylim(-6, 1)
		draw()
	
	def plot_logslope(self, logrs):
		slope = self.profile.logslope(10**logrs)
		ok = ~isnan(slope)
		graph(logrs[ok], slope[ok]) 
		#graph(logrs, self.profile.fast().densityr(10**logrs)), color="blue", linestyle="dash" )
		labels("log r", "dlog &rho;(r) / dlogr")
		#xlim(0, 1.5)
	def plot_density(self, logrs):
		rho = log10(self.profile.densityr(10**logrs)) 
		ok = ~(isnan(rho) | isinf(rho))
		graph(logrs[ok], rho[ok])
		ymax = (max(rho[ok]))
		ylim(ymax-10, ymax+1)
		#graph(logrs, log10(self.profile.fast().densityr(10**logrs)), color="blue", linestyle="dash" )
		labels("log r", "log &rho;(r)")
		
	def plot_enclosed_mass(self, logrs):
		rs = 10**logrs
		graph(logrs, log10(self.profile.enclosed_mass(rs)))
		#graph(logrs, log10(self.profile.enclosed_mass_num(rs)), color="red", linestyle="dot")
		#graph(logrs, log10(self.profile.fast().enclosed_mass(rs)), color="blue", linestyle="dash")
		labels("log r", "M_enc(r)")
		
class Orbit(object):
	def __init__(self, orbit):
		self.orbit = orbit
		
	def run(self, args, opts, scope):
		document(size="45cm,30cm")
		#box()
		mozaic(4,3,box)
		q, p = self.orbit.integrate()
		xmax = q.max()
		vmax = p.max()
		graph(q[0], q[1])
		print q
		xlim(-xmax, xmax)
		ylim(-xmax, xmax)
		
		select(1, 0)
		graph(q[0], p[1])
		xlim(-xmax, xmax)
		ylim(-vmax, vmax)
		labels("x", "vy")
		
		select(1, 1)
		graph(q[1], p[1])
		xlim(-xmax, xmax)
		ylim(-vmax, vmax)
		labels("y", "vy")
		
		select(0, 1)
		graph(q[0], p[0])
		xlim(-xmax, xmax)
		ylim(-vmax, vmax)
		labels("x", "vx")
		
		x, y = q
		vx, vy = p
		r = sqrt(x**2+y**2)
		phi = arctan2(y,x)
		vr = (x * vx + y * vy) / r
		vphi = (x * vy - y * vx) / r
		
		select(2, 0)
		graph(r, vr)
		#ylim(-xmax, xmax)
		#ylim(-vmax, vmax)
		labels("r", "vr")
		
		select(2, 1)
		#graph(vr, vphi-220)
		#ylim(-xmax, xmax)
		#ylim(-vmax, vmax)
		#labels("vr", "vphi")
		graph(phi, vphi*r)
		#ylim(-xmax, xmax)
		#ylim(-vmax, vmax)
		labels("phi", "p_phi=vphi*r")
		
		select(0, 2)
		#graph(vr, vphi-220)
		#ylim(-xmax, xmax)
		#ylim(-vmax, vmax)
		#labels("vr", "vphi")
		graph(vphi)
		#ylim(-xmax, xmax)
		#ylim(-vmax, vmax)
		labels("phi", "p_phi=vphi*r")
		
		draw()
		
		

class Data1d(object):
	def __init__(self, binned_data, light_model, maxR, aperture, galaxy=None):
		self.binned_data = binned_data
		self.light_model = light_model
		self.aperture = aperture
		self.maxR = maxR
		self.galaxy = galaxy
		
	def load(self):
		self.binned_data.load()
		self.aperture.load()
		
	def run(self, args, opts, scope):
		self.load()
		box()
		
		self.plot_data()
		xlim(0, self.maxR)
		labels("projected radius - kpc", "&sigma;<sub>los</sub> - km/s")
		ylim(0, 20)
		draw()
		
	def plot_data(self, **kwargs):
		moments = self.binned_data.moments
		e_moments = self.binned_data.e_moments
		
		
		r_borders = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders)
		rcenters_arcsec = (self.aperture.aperture_rborders[0:-1]+self.aperture.aperture_rborders[1:])/2
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters_arcsec)
		
		r = rcenters_kpc * 1e3
		#r = rcenters_arcsec / (60*60)
		sigma = sqrt(moments[2])
		print "r=",rcenters_arcsec
		print "sigma=", sigma
		print len(sigma), len(rcenters_arcsec)
		#errors = sqrt(e_moments[2]**2/(2*sigma))
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		
		#scatter(r, sigma)
		print sigma
		print errors
		scatter(rcenters_kpc, sigma, symbolName="squaresolid", **kwargs)
		errorbars(rcenters_kpc, sigma, yerr=errors, **kwargs)
			
		galaxy = None
		if galaxy:
			R = arange(0, 1.6, 0.05)
			jeans = galaxy.jeans()
			sigmaR = jeans.sigma_los(R)
			graph(R, sigmaR, color="red", linestyle="dash")
		#xlim(0, max(r) * 1.2)
		#xlim(0, 1.6)
		if 0: #self.galaxy:
			jeans = self.galaxy.jeans()
			n = 100
			R = arange(n)/(n+1.) * self.maxR
			sigma_los = jeans.sigma_los(R)
			graph(R, sigma_los, color="red")
		
		if 0:
			select(1,0)
			N = len(rcenters_kpc)
			x = arange(N)/(N-1.) * N
			scatter(x, sigma/20*100, symbolName="squaresolid")

class ApertureHistograms(object):
	def __init__(self, observation, aperture, binned_data, vmin, vmax, vmean, filters=[], v_mirror=None): #=None):
		self.observation = observation
		self.aperture = aperture
		self.binned_data = binned_data
		self.vmin = vmin
		self.vmax = vmax
		self.vmean = vmean
		self.filters = filters
		self.v_mirror = v_mirror
		
	def run(self, args, opts, scope):
		self.load()
		#v_center = 55.
		N = len(self.binned_data.moments[0])
		Nx = int(sqrt(N))
		Ny = int(math.ceil(float(N)/Nx))
		logger.info("%d apertures, %dx%d" % (N, Nx, Ny))
		mozaic(Nx, Ny, box)
		for i in range(N):
			stars = self.stars.filter(self.aperture.aperture_filter(i))
			for filter in self.filters:
				stars = filter(stars)
				
			if 1:
				vlos = stars.vlos
				vlos = vlos[(vlos > self.vmin) & (vlos < self.vmax)]
				print "mean/sigma", mean(vlos), std(vlos), self.binned_data.moments[2,i]**0.5 #min(vlos), max(vlos)
			#import pdb
			#pdb.set_trace()
			logger.info("%d stars" % (len(stars)))
			select(i%Nx, i/Nx)
			#vs = stars.vlos_helio
			#vs -= vs.mean()
			#print vs.var()**0.5
			normalize = True
			if self.v_mirror is None:
				#vs = stars.vlos_helio
				vs = stars.vlos
				#histogram(vs, binwidth=2, datamin=-40+v_center, datamax=40+v_center, normalize=normalize)
				histogram(vs, binwidth=2.5, datamin=self.vmin, datamax=self.vmax, normalize=normalize)
				#print stars.attributes[0]
				#scatter(vs, stars.rc/60/60/10)
			else:
				vs = abs(stars.vlos-self.v_mirror)
				histogram(vs, binwidth=2, datamin=0, datamax=40, normalize=normalize)
			#x = arange(-40+v_center, 40+v_center, 0.5)
			x = arange(self.vmin, self.vmax, 0.5)
			graph(x, gaussian(x, self.vmean, self.binned_data.moments[2,i]**0.5), color="red")
			print self.binned_data.moments[2,i]**0.5
		draw()
			
		
	def load(self):
		self.stars = self.observation.load()
		self.aperture.load()
		self.binned_data.load()

class VelocityMoments(object):
	def __init__(self, binned_data, light_model, Rmax, modelpath, galaxy=None):
		self.binned_data = binned_data
		self.light_model = light_model
		self.Rmax = Rmax
		self.galaxy = galaxy
		
		self.dirname = os.path.join(modelpath, "plots")
		if not os.path.exists(self.dirname):
			os.makedirs(self.dirname)
		print self.dirname
		
	def load(self):
		for binned_data in self.binned_data:
			binned_data.load()
			binned_data.aperture.load()
		#self.aperture.load()
		
	def run(self, args, opts, scope):
		self.load()
		#box()
		document(size="20cm,28cm")
		mozaic(1,4,box)
		
		for binned_data, color in zip(self.binned_data, goodcolors)[:]:
			self.plot_data(binned_data, color=color)
		select(0, 0)
		xlim(0, self.Rmax)
		labels("projected radius - kpc", "v<sub>los</sub> - km/s")
		hline(0)
		ylim(-10, 10)
		
		select(0, 1)
		xlim(0, self.Rmax)
		labels("projected radius - kpc", "&sigma;<sub>los</sub> - km/s")
		ylim(0, 15)
		
		select(0, 2)
		xlim(0, self.Rmax)
		labels("projected radius - kpc", "k3")
		ylim(-3, 3.)
		hline(0)
		
		select(0, 3)
		xlim(0, self.Rmax)
		labels("projected radius - kpc", "k4")
		ylim(1, 5.)
		hline(3)
		
		draw()
		
	def plot_vdisp(self, binned_data, **kwargs):
		moments = binned_data.moments
		e_moments = binned_data.e_moments
		aperture = binned_data.aperture
		
		r_borders = self.light_model.arcsec_to_kpc(aperture.aperture_rborders)
		rcenters_arcsec = (aperture.aperture_rborders[0:-1]+aperture.aperture_rborders[1:])/2
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters_arcsec)
		
		r = rcenters_kpc * 1e3
		sigma = sqrt(moments[2])
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		
		scatter(rcenters_kpc, sigma, symbolName="squaresolid", **kwargs)
		errorbars(rcenters_kpc, sigma, yerr=errors, **kwargs)

	def plot_kurtosis(self, binned_data, **kwargs):
		moments = binned_data.moments
		e_moments = binned_data.e_moments
		aperture = binned_data.aperture
		
		r_borders = self.light_model.arcsec_to_kpc(aperture.aperture_rborders)
		rcenters_arcsec = (aperture.aperture_rborders[0:-1]+aperture.aperture_rborders[1:])/2
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters_arcsec)
		
		r = rcenters_kpc * 1e3
		sigma = sqrt(moments[2])
		k4 = moments[4]/moments[2]**2
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		#k4_errors = ( (1/moments[2] **2)**2 * e_moments[4]**2 + (-2*moments[4]/moments[2]**3)**2*e_moments[2]**2) **0.5
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		rho = 0.867819669542 # numerical result for a gaussian
		k4_errors = ( (1/moments[2] **2)**2 * e_moments[4]**2 + (-2*moments[4]/moments[2]**3)**2*e_moments[2]**2\
			-2 * (1/moments[2] **2)*2*moments[4]/moments[2]**3 * rho * e_moments[4] * e_moments[2]) **0.5
		
		scatter(rcenters_kpc, k4, symbolName="squaresolid", **kwargs)
		#graph(rcenters_kpc, k4, **kwargs)
		errorbars(rcenters_kpc, k4, yerr=k4_errors, **kwargs)

	def plot_data(self, binned_data, **kwargs):
		moments = binned_data.moments
		e_moments = binned_data.e_moments
		aperture = binned_data.aperture
		
		r_borders = self.light_model.arcsec_to_kpc(aperture.aperture_rborders)
		rcenters_arcsec = (aperture.aperture_rborders[0:-1]+aperture.aperture_rborders[1:])/2
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters_arcsec)
		
		r = rcenters_kpc * 1e3
		#r = rcenters_arcsec / (60*60)
		sigma = sqrt(moments[2])
		k4 = moments[4]/moments[2]**2
		k3 = moments[3]/sigma**3
		print "r=",rcenters_arcsec, rcenters_kpc
		print "sigma=", sigma
		print len(sigma), len(rcenters_arcsec)
		#errors = sqrt(e_moments[2]**2/(2*sigma))
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		rho = 0.867819669542 # numerical result for a gaussian
		k4_errors = ( (1/moments[2] **2)**2 * e_moments[4]**2 + (-2*moments[4]/moments[2]**3)**2*e_moments[2]**2\
			-2 * (1/moments[2] **2)*2*moments[4]/moments[2]**3 * rho * e_moments[4] * e_moments[2]) **0.5
		print "correction", -2 * (1/moments[2] **2)*2*moments[4]/moments[2]**3 * rho * e_moments[4] * e_moments[2]
		k3_errors = ( (1/sigma**3)**2 * e_moments[3]**2 + (-3*moments[3]/sigma**4)**2*e_moments[2]**2) **0.5
		
		#scatter(r, sigma)
		print sigma
		print errors
		print binned_data.postfix, kwargs
		select(0, 1)
		scatter(rcenters_kpc, sigma, symbolName="squaresolid", **kwargs)
		errorbars(rcenters_kpc, sigma, yerr=errors, **kwargs)
			
		#xlim(0, max(r) * 1.2)
		#xlim(0, 1.6)
		print self.galaxy, hasattr(self.galaxy, "jeans")
		if self.galaxy and hasattr(self.galaxy, "jeans"):
			filename = os.path.join(self.dirname, "sigma_los.npy")
			print "cache:", filename
			n = 100
			R = arange(n)/(n-1.) * self.Rmax
			if os.path.exists(filename):
				sigma_los = load(filename) 
			else:
				jeans = self.galaxy.jeans()
				sigma_los = jeans.sigma_los(R)
				save(filename, sigma_los)
			graph(R, sigma_los, color="blue")
		
		select(0,3)
		if self.galaxy and hasattr(self.galaxy, "jeans"):
			filename = os.path.join(self.dirname, "m4_los.npy")
			print "cache:", filename
			n = 100
			R = arange(n)/(n-1.) * self.Rmax
			if os.path.exists(filename):
				m4_los = load(filename) 
			else:
				jeans = self.galaxy.jeans()
				m4_los = jeans.m4_los(R)
				save(filename, m4_los)
			graph(R, m4_los/sigma_los**4, color="blue")
		scatter(rcenters_kpc, k4, symbolName="squaresolid", **kwargs)
		graph(rcenters_kpc, k4, **kwargs)
		errorbars(rcenters_kpc, k4, yerr=k4_errors, **kwargs)
		select(0, 2)
		scatter(rcenters_kpc, k3, symbolName="squaresolid", **kwargs)
		graph(rcenters_kpc, k3, **kwargs)
		errorbars(rcenters_kpc, k3, yerr=k3_errors, **kwargs)
		#errorbars(rcenters_kpc, k4, yerr=k4_errors, **kwargs)
		
		select(0, 0)
		scatter(rcenters_kpc, moments[1], symbolName="squaresolid", **kwargs)
		graph(rcenters_kpc, moments[1], **kwargs)
		#errorbars(rcenters_kpc, k3, yerr=k3_errors, **kwargs)
		
		if 0:
			select(1,0)
			N = len(rcenters_kpc)
			x = arange(N)/(N-1.) * N
			scatter(x, sigma/20*100, symbolName="squaresolid")

	def plot_grid(self, grid):
		#colours = [Color.blue.blend(Color.white, a) for a in [0.2, 0.35, 0.5]][::-1]
		alphas = [0.45, 0.3, 0.15]
		#indexedimage(grid.grid.T, resize=grid.resize, colormap="whiterainbow")
		for i, level in list(enumerate(grid.levels))[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			ymin = grid.contours[:,i,0]
			ymax = grid.contours[:,i,1]
			#print len(grid.profile.u_centers), len(ymin)
			print len(grid.profile.u_centers), len(ymin), len(ymax)
			fillrange(grid.profile.u_centers, ymin, ymax, color=color)
		graph(grid.profile.u_centers, grid.median, color="green")
		labels(grid.profile.label, grid.label)
		(x1, y1), (x2, y2) = grid.resize
		xlim(x1, x2)
		ylim(y1, y2)

class Moments24(VelocityMoments):
	def __init__(self, binned_m2, binned_m4, light_model, Rmax, modelpath, vdisp_grid=None, kurtosis_grid=None, xlabel=None, ylabel1=None, ylabel2=None, title=None, scopebestfit=None, scopes=None, colors=None, legend=False):
		VelocityMoments.__init__(self, binned_data=[binned_m2, binned_m4], light_model=light_model, Rmax=Rmax, modelpath=modelpath)
		self.vdisp_grid = vdisp_grid
		self.kurtosis_grid = kurtosis_grid
		self.xlabel = xlabel
		self.ylabel1 = ylabel1
		self.ylabel2 = ylabel2
		self.title = title
		self.scopebestfit = scopebestfit
		self.scopes=scopes
		self.colors = colors
		self.legend = legend
		#self.Rmax = Rmax
		
	def run(self, args, opts, scope_):
		if self.scopes:
			for scope, color in zip(self.scopes, self.colors):
				scope.load()
				print scope.subscope["modelpath"]
			for scope, color in zip(self.scopes, self.colors):
				scope.load()
				print scope.subscope["modelpath"]
		#sys.exit(0)
		self.load()
		document(size="20cm,20cm")
		page(fontsize="20pt")
		#container(viewport=[(0, 0), (0.8, 0.8)]) 
		mozaic(1,2,box)
		if self.scopebestfit:
			self.scopebestfit.load()

		select(0,1)
		left = 0.15
		top = 0.05
		right = 0.05
		current.container.drawOutsideLeft = True
		current.container.viewport = (left, 0.5-top/2), (1.-right, 1.0-top)
		color = "black"
		if self.vdisp_grid and not self.scopes:
			self.vdisp_grid.load()
			self.plot_grid(self.vdisp_grid)
		else:
			labels(self.xlabel, self.ylabel1)
			#labels(self.vdisp_grid.profile.label, self.vdisp_grid.label)
		self.plot_vdisp(self.binned_data[0], color=color)
		if self.scopebestfit:
			storage_2d = self.scopebestfit.subscope["storage_2d_m4"]
			solution = self.scopebestfit.subscope["solution"]
			storage_2d.load()
			solution.load()
			allmoments = storage_2d.projectedmoments
			moments = solution.calculate_solution_moments(allmoments)
			moments[1:] /= moments[0]
			m4 = moments[4]
			m2 = moments[2]
			kappa = m4/m2**2
			assert len(self.vdisp_grid.profile.u_centers) == len(m2)
			graph(self.vdisp_grid.profile.u_centers, m2**0.5, color="orange")
			
		xlim(0, self.Rmax)
		#labels("projected radius - kpc", "&sigma;<sub>los</sub> - km/s")
		#hline(0)
		#ylim(-10, 10)
		ylim(0, 20)
		
		if self.title:
			select(0, 1)
			title(self.title)
		select(0, 0)
		current.container.drawOutsideLeft = True
		current.container.viewport = (left, 0), (1.-right, 0.5-top/2)
		if self.kurtosis_grid and not self.scopes:
			self.kurtosis_grid.load()
			self.plot_grid(self.kurtosis_grid)
		else:
			labels(self.xlabel, self.ylabel2)
			#labels(self.kurtosis_grid.profile.label, self.kurtosis_grid.label)

		#hline(3)
		self.plot_kurtosis(self.binned_data[1], color=color)
		if self.scopebestfit:
			graph(self.kurtosis_grid.profile.u_centers, kappa, color="orange")
		xlim(0, self.Rmax)
		#labels("projected radius - kpc", "kurtosis")
		ylim(1, 5.)
		clearautolegend()
		if self.scopes:
			for scope, color in zip(self.scopes, self.colors):
				print color
				#scope.load()
				storage_2d = scope.subscope["storage_2d_m4"]
				solution = scope.subscope["solution"]
				storage_2d.load()
				solution.load()
				allmoments = storage_2d.projectedmoments
				moments = solution.calculate_solution_moments(allmoments)
				moments[1:] /= moments[0]
				m4 = moments[4]
				m2 = moments[2]
				kappa = m4/m2**2
				assert len(self.vdisp_grid.profile.u_centers) == len(m2)
				select(0,1)
				graph(self.vdisp_grid.profile.u_centers, m2**0.5, color=color)
				select(0,0)
				graph(self.kurtosis_grid.profile.u_centers, kappa, color=color, addlegend=False)
			if self.legend:
				names = [k.subscope["dm_name"] for k in self.scopes]
				autolegend(*names)
			
		draw()

class Losvds(object):
	def __init__(self, observation, aperture, binned_data):
		self.observation = observation
		self.aperture = aperture
		self.binned_data = binned_data
		
	def run(self, args, opts, scope):
		self.aperture.load()
		self.binned_data.load()
		stars = self.observation.load()
		
		if self.binned_data.filters:
			for filter in self.binned_data.filters:
				stars = filter(stars)
		
		N = len(self.aperture.aperture_rcenters_kpc)
		print "N", N
		nx = 4
		ny = N/nx+1
		aperture_stars = [list() for k in range(N)]
		x = stars.xi
		y = stars.eta
		for i in range(len(x)):
			if self.aperture.aperture.inrange(x[i], y[i]):
				aperture_index = self.aperture.aperture.findindex(x[i], y[i])
				aperture_stars[aperture_index].append(stars[i])
		
		index = 0
		mozaic(nx,ny,box)
		#import pdb;
		#pdb.set_trace()
		#self.binned_data.moments[2,0] = 10.**2
		for j in range(ny):
			for i in range(nx):
				select(i, j)
				if index < N:
					stars = aperture_stars[index]
					print len(stars)
					vlos = [star.vlos for star in stars]
					sigma = self.binned_data.moments[2,index]**0.5
					mean = self.binned_data.moments[1,index]
					print sigma
					#mean = 0
					if 0:
						histogram(vlos, datamin=-30, datamax=30, bincount=30, normalize=True)
						x = arange(-30, 30, 0.5)
						y = gaussian(x, mean, sigma)
						#y = cumsum(y)
						#y /= y.max()
					else:
						cumhistogram(vlos, datamin=-30, datamax=30, bincount=60, normalize=True)
						x = arange(-30, 30, 0.5)
						y = gaussian(x, mean, sigma)
						y = cumsum(y)
						y /= y.max()
					graph(x, y, color="red")
					index += 1
		draw()
		

class Data1dRvsVlos(object):
	def __init__(self, observation, binned_data, light_model, maxR, aperture, filters=[], filter_r=None, filter_three_sigma=None, vmin=None, vmax=None, paper=False, title=None, filter_rmax_multi=None):
		self.observation = observation
		self.binned_data = binned_data
		self.light_model = light_model
		self.aperture = aperture
		self.maxR = maxR
		self.filter_three_sigma = filter_three_sigma
		self.filters = filters
		self.filter_r = filter_r
		self.vmin = vmin
		self.vmax = vmax
		self.paper = paper
		self.title = title
		self.filter_rmax_multi = filter_rmax_multi
		
	def run(self, args, opts, scope):
		if self.paper:
			document(size="7cm,15cm")
			page(fontsize="13pt")
		else:
			document(size="10cm,15cm")
		box()
		if self.paper:
			print dir(current.container.xaxis)
			current.container.xaxis.interval = 200
		self.observation.load()
		self.binned_data.load()
		self.aperture.load()
		if self.title:
			title(self.title)
		
		stars = self.observation.stars
		for filter in self.filters:
			stars = filter(stars)
		values = None
		if hasattr(stars[0], "p_member"):
			values = stars.p_member
		
		for r in self.aperture.aperture_rborders:
			if self.paper:
				hline(self.light_model.arcsec_to_kpc(r), color="grey")
			else:
				hline(r, color="grey")
		if self.paper:
			rs = self.light_model.arcsec_to_kpc(stars.rc) #/60**2
			scatter(stars.vlos_helio, self.light_model.arcsec_to_kpc(stars.rc), colors=values, colormap="rainbow", symbolsize="10pt")
		else:
			rs = stars.rc
			scatter(stars.vlos_helio, stars.rc, colors=values, colormap="rainbow", symbolsize="20pt")
		#scatter(stars.vlos, stars.rc)
		color = "grey"
		if self.filter_three_sigma:
			sigma = self.filter_three_sigma.sigma
			vline(self.filter_three_sigma.v, color=color)
			vline(self.filter_three_sigma.v+3*sigma, color=color)
			vline(self.filter_three_sigma.v-3*sigma, color=color)
		if self.filter_r:
			if self.paper:
				hline(self.light_model.arcsec_to_kpc(self.filter_r.rmax_arcsec), color="red")
			else:
				hline(self.filter_r.rmax_arcsec, color="red")
		if self.filter_rmax_multi:
			for i, color in zip(range(self.filter_rmax_multi.rmax.simplemodel.catalogues), nicecolors[1:]):
				hline(self.filter_rmax_multi.rmax.Rmaxs[i]/60**2, color=color)
				
			
		if self.paper:
			labels("v<sub>los,helio</sub> (km/s)", "r (kpc)")
		else:
			labels("v<sub>los,helio</sub>", "r - arcsec")
		
		if self.vmin is not None and self.vmax is not None:
			xlim(self.vmin, self.vmax)
			
		ylim(0, max(rs))
		
		
		draw()
		return
		
		moments = self.binned_data.moments
		e_moments = self.binned_data.e_moments
		
		
		r_borders = self.light_model.arcsec_to_kpc(self.aperture.aperture_rborders)
		rcenters_arcsec = (self.aperture.aperture_rborders[0:-1]+self.aperture.aperture_rborders[1:])/2
		rcenters_kpc = self.light_model.arcsec_to_kpc(rcenters_arcsec)
		
		r = rcenters_kpc * 1e3
		#r = rcenters_arcsec / (60*60)
		sigma = sqrt(moments[2])
		print "r=",rcenters_arcsec
		print "sigma=", sigma
		print len(sigma), len(rcenters_arcsec)
		#errors = sqrt(e_moments[2]**2/(2*sigma))
		errors = sqrt(e_moments[2]**2/(2*sigma)**2)
		
		#scatter(r, sigma)
		print sigma
		print errors
		scatter(rcenters_kpc, sigma, symbolName="squaresolid")
		errorbars(rcenters_kpc, sigma, yerr=errors)
			
		galaxy = None
		if galaxy:
			R = arange(0, 1.6, 0.05)
			jeans = galaxy.jeans()
			sigmaR = jeans.sigma_los(R)
			graph(R, sigmaR, color="red", linestyle="dash")
		ylim(0, 20)
		#xlim(0, max(r) * 1.2)
		xlim(0, 1.6)
		
		labels("projected radius - kpc", "&sigma;<sub>los</sub> - km/s")
		if 0:
			select(1,0)
			N = len(rcenters_kpc)
			x = arange(N)/(N-1.) * N
			scatter(x, sigma/20*100, symbolName="squaresolid")
		draw()


class DataSpatial(object):
	def __init__(self, observation, maxR, aperture=None, velocityfilter=None, filters=[]):
		self.observation = observation
		#self.light_model = light_model
		self.maxR = maxR
		self.aperture = aperture
		self.filters = filters
		self.velocityfilter = velocityfilter
		
	def run(self, args, opts, scope):
		document(size="24cm,10cm")
		#box()
		mozaic(2,1,box)
		self.observation.load()
		self.aperture.load()
		
		stars = self.observation.stars
		for filter in self.filters:
			stars = filter(stars)
		
		if self.velocityfilter:
			vlos = stars.vlos#_helio
			#scale = self.velocityfilter.scale
			#self.velocityfilter.scale = 2500*2#/2.
			#vsub = self.velocityfilter.filter(stars.xi, stars.eta, vlos)
			#self.velocityfilter.scale = scale
			import random
			if 0:
				for i in range(100):
					indices = range(len(stars))
					random.shuffle(indices)
					vlos = vlos[indices]
					v = self.velocityfilter.filter(stars.xi, stars.eta, vlos)
					print v.min(), v.max(), v.mean(), v.std()
				
			randomize = True
			randomize = False
			if randomize:
				indices = range(len(stars))
				random.shuffle(indices)
				vlos = vlos[indices]
			mean_vsys=-31.673469387755102
			#vlos = vlos - mean_vsys
			#vlos = abs(vlos)
			v = vlos- mean_vsys #
			#v = self.velocityfilter.filter(stars.xi, stars.eta, v)# - vsub
			v = abs(v)
			print v.min(), v.max(), v.mean(), v.std()
			mid, sigma = median(v), v.std()
			#colormap = ColorMap.redwhiteblue
			#colormap = ColorMap(["red", "black", "blue"])#.redblackblue
			#colormap = ColorMap(["red", "blue"])#.redblackblue
			#colormap = ColorMap(["white", "blue"])#.redblackblue
			colormap = ColorMap(["red", "blue"])#.redblackblue
			#scatter(stars.xi/3600, stars.eta/3600, colors=v, datamin=mid-sigma*1.5, datamax=mid+sigma*1.5, colormap=colormap, symbolsize="25pt")
			#scatter(stars.xi/3600, stars.eta/3600, colors=v, datamin=0, datamax=20, colormap=colormap, symbolsize="25pt")
			scatter(stars.xi/3600, stars.eta/3600, colors=v, datamin=0, datamax=15., colormap=colormap, symbolsize="20pt", alpha=1.0)
		else:
			scatter(stars.xi/3600, stars.eta/3600)
		if self.maxR:
			maxR = self.maxR
		else:
			maxR = max(stars.xi.max(), stars.eta.max())
		vline(0)
		hline(0)
		squarelim(-maxR, maxR)
		labels("xi", "eta")
		select(1,0)
		#histogram(vlos, datamin=0., datamax=10., binwidth=1.)
		if self.velocityfilter:
			histogram(v, datamin=-30., datamax=30., binwidth=1.)
		draw()
		
		

class DataVlosHistogram(object):
	def __init__(self, observation, foregrounddata=None, binwidth=1., helio=False, vmin=-100, vmax=100, filters=[]):
		self.observation = observation
		self.foregrounddata = foregrounddata
		self.filters = filters
		self.binwidth = binwidth
		self.helio = helio
		self.vmin = vmin
		self.vmax = vmax
		
	def run(self, *args):
		document(size="35cm,20cm")
		mozaic(2,1,box)
		stars = self.observation.load()
		
		for filter in self.filters:
			print filter
			stars = filter(stars)
		if self.foregrounddata:
			self.foregrounddata.load()
		#box()
		#histogram(stars.vlos_helio, binwidth=2.)
		
		#mean_vsys=-31.673469387755102
		if self.helio:
			vlos = stars.vlos_helio
		else:
			vlos = stars.vlos
		#vlos = stars.vlos - mean_vsys
		#vlos = abs(vlos)
		select(0,0)
		histogram(vlos, binwidth=self.binwidth)
		if self.foregrounddata:
			histogram(self.foregrounddata.data.Vr, binwidth=2., color="red")
		xlim(self.vmin, self.vmax)
		grow(top=1.1)
		select(1,0)
		histogram(vlos, binwidth=self.binwidth, function=lambda x: log10(maximum(1, x)))
		grow(top=1.1)
		#histogram(stars.vlos_helio-86, binwidth=1., function=lambda x: log10(maximum(1, x)), color="red")
		#xlim(-100, 100)
		xlim(self.vmin, self.vmax)
		if self.foregrounddata:
			histogram(self.foregrounddata.data.Vr, binwidth=2., color="red", function=lambda x: log10(maximum(1, x)))
		#histogram(stars.vlos, binwidth=1.)
		draw()
		
		
class ParameterSet(object):
	def __init__(self, parameterset, smooth_param=0.5):
		self.parameterset = parameterset
		self.smooth_param = smooth_param
		
	def load(self):
		self.parameterset.load()
		
		if 1:
			keys = self.parameterset.valuemap.keys()
			values = self.parameterset.valuemap.values()
			values = array(values)
			indices = argsort(values)[::-1]
			for i in range(10):
				index = indices[i]
				value = values[index]
				key = keys[index]
				print value, key, ["%g" % v for v in key], [log10(v) % v for v in key]
			#print values[indices]
			import pdb;
			#pdb.set_trace()
		
		
		#self.parameterset.grid()
		grid = self.parameterset.probability_grid
		#print grid.max(), grid.min()
		self.p = (grid)
		#self.p += self.p.min()
		scale = 6000/200.
		scale = 1.
		self.p = exp(self.p*scale)
		self.p[isnan(grid)] = 0
#<<<<<<< Updated upstream
#		self.smooth = self.p
#=======
		self.smooth = self.p #scipy.ndimage.gaussian_filter(self.p, 0.005)
		#self.smooth = scipy.ndimage.gaussian_filter(self.p, 10.5)
		
#>>>>>>> Stashed changes
		index = argmax(self.smooth.flat)
		index = unravel_index(index, self.smooth.shape)
		#scales = parametersetcollect.probability_range
		#print scales
		#print index
		parameter_range_list = self.parameterset.parameter_range_list
		
		self.max_L = [parameter_range_list[k].min + (parameter_range_list[k].max - parameter_range_list[k].min) * (index[k]+0.5)/(self.smooth.shape[k]-1.) for k in range(len(index))]
		print "max Likelihood", self.max_L
		#self.smooth = self.p #
		if 1:
			self.smooth = self.p * 1.
			#self.smooth = scipy.ndimage.gaussian_filter(self.p, 1.5)
			#self.smooth = scipy.ndimage.gaussian_filter(self.p, 0.5)
			print "self.smooth_param =", self.smooth_param
			self.smooth = scipy.ndimage.gaussian_filter(self.p, self.smooth_param)
			#self.smooth = self.p
			print "..."
			#dsa
			ranges = self.makeresize(0,1)
			import pickle
			p = self.smooth
			fn ="scl.pickle"
			with open(fn, "w") as f:
				pickle.dump((ranges, p), f)
			dsa
		
		
	def run(self, argv, opts, scope):
		self.load()
		dim = self.parameterset.dimension
		w = 7 * dim
		document("%fcm,%fcm" % (w,w))
		mozaic(dim, dim, box)
		
		print "minmax", self.smooth.max(), self.smooth.min()
		#indexedimage(self.smooth.T)
		#draw()
		
		for i in range(dim):
			for j in range(dim):
				select(i, dim-1-j)
				self.plot_parameters(i, j, scope)
		draw()
				
	def makeresize(self, i, j):
		return (self.parameterset.parameter_range_list[i].min, self.parameterset.parameter_range_list[j].min), (self.parameterset.parameter_range_list[i].max, self.parameterset.parameter_range_list[j].max)

	def plot_parameters(self, i, j, scope, prior_grid=None, addlabels=True, fill=True, color="black", linewidth="0.2pt", sigmas=3, **kwargs):
		parameterset = self.parameterset
		def __makeresize(i, j):
			#return (parameterset.ranges_org[i][0], parameterset.ranges_org[i][1]), (parameterset.ranges_org[j][0], parameterset.ranges_org[j][1])
			return (parameterset.ranges_org[i][0], parameterset.ranges_org[j][0]), (parameterset.ranges_org[i][1], parameterset.ranges_org[j][1])
		p = self.smooth
		if prior_grid is not None:
			p *= prior_grid
		if i == j:
			#probgraph(p, i, resize=parameterset.parameter_ranges[i])
			print "min/max", parameterset.parameter_range_list[i].min, parameterset.parameter_range_list[i].max
			print p, i
			probgraph(p, i, resize=(parameterset.parameter_range_list[i].min, parameterset.parameter_range_list[i].max))
			#labels(parameterset.parameter_labels[i], "p")
			if addlabels:
				labels(parameterset.parameter_range_list[i].label, "prob. density")
			delta = parameterset.parameter_range_list[i].max-parameterset.parameter_range_list[i].min
			d = delta / 2. / self.p.shape[i]
			xlim(parameterset.parameter_range_list[i].min-d, parameterset.parameter_range_list[i].max+d)
			pr = parameterset.parameter_range_list[i]
			name = pr.name
			obj = vline(self.max_L[i], color="blue")
			if name in scope:
			  value = scope[name]
			  print name, scope[name], pr.finv(value)
			  vline(pr.finv(value), color="red", linestyle="dash")
			#if name in scope:
			#xs = [pv.values_org[i] for pv in parameterset.parameter_values]
			#for x in xs:
			#	vline(x, alpha=0.2)
		else:
			#print self.p
			#print self.p.shape, i, j, self.parameterset.dimension
			#self.p
			#robimage2d(p, i, j, resize=self.makeresize(i,j), drawcontourlines=True, fill=False, colormap=None, premultiply=True, color=color, **kwargs) #(0, 1), (0,1))
			if fill:
				#probimage2d(p, i, j, resize=self.makeresize(i,j), drawcontourlines=False, colormap=None, premultiply=True, color=color, **kwargs) #(0, 1), (0,1))
				if 0:
					print p.T
					print p.min()
					mask = p.T == 0
					I = log10(p.T)
					I -= I.max()
					I[mask] = -3
					#mask = I < -7
					#I[mask] = -7
					
					indexedimage(I, resize=self.makeresize(i,j), )
					contour(I, levels=20, resize=self.makeresize(i,j), )
					contour(I, levels=[-3,-2,-1], resize=self.makeresize(i,j), )
					contour(I, levels=[-1,-0.5,-0.25], resize=self.makeresize(i,j), )
					#contour(I, levels=[-3,-6,-9], resize=self.makeresize(i,j), color="white")
				obj = probimage2d(p, i, j, resize=self.makeresize(i,j), drawcontourlines=True, fill=False, colormap=None, premultiply=True, color=color, sigmas=sigmas, linewidth=linewidth, **kwargs) #(0, 1), (0,1))
				#probimage2d(p, i, j, resize=self.makeresize(i,j), drawcontourlines=False, fill=False, colormap="whiteblack", premultiply=True, color=color, **kwargs) #(0, 1), (0,1))
			else:
				obj = probimage2d(p, i, j, resize=self.makeresize(i,j), drawcontourlines=True, colormap=None, fill=False, color=color, linewidth=linewidth, sigmas=sigmas, **kwargs) #(0, 1), (0,1))
			options = scope["options"]
			if hasattr(options, "show_grid") and options.show_grid:
				valuegrid = parameterset.valuegrid
				#x = valuegrid[i]
				#y = valuegrid[j]
				x = [pv.values_org[i] for pv in parameterset.parameter_values]
				y = [pv.values_org[j] for pv in parameterset.parameter_values]
				scatter(x, y)
				#print x
			if addlabels:
				if 0:
					o = labels(parameterset.parameter_range_list[i].label, parameterset.parameter_range_list[j].label)
					print o
					x, y = o
					y.context.fontname="CM Typewriter Greek"
				else:
					xlabel(parameterset.parameter_range_list[i].label)
					obj = ylabel(parameterset.parameter_range_list[j].label)
					if "gamma" in parameterset.parameter_range_list[j].label:
						obj.context.fontname="CM Typewriter Greek"
			pr_i = parameterset.parameter_range_list[i]
			pr_j = parameterset.parameter_range_list[j]
			#symbol(self.max_L[i], self.max_L[j], color="blue", symbolsize="40pt", addlegend=False)
			symbol(self.max_L[i], self.max_L[j], color=color, symbolsize="20pt")
			#if pr_i.name in scope:
			#	symbol(pr_i.finv(scope[pr_i.name]), pr_j.finv(scope[pr_j.name]), color="red", symbolsize="20pt")
			#labels(parameterset.parameter_labels[i], parameterset.parameter_labels[j])
			#labels(parameterset.paramlabels[i], parameterset.paramlabels[j])
		return obj
class ParameterSet2d(ParameterSet):
	def __init__(self, subscopes=[], colors=[], smooth_param=0.5, plotlegend=False):
		ParameterSet.__init__(self, None, smooth_param=smooth_param)
		self.subscopes = subscopes
		self.colors = colors
		self.plotlegend = plotlegend
		
	def run(self, argv, opts, scope, nodraw=False, init=True):
		#self.load()
		#dim = self.parameterset.dimension
		#w = 7 * dim
		if init:
			w = 7
			document("%fcm,%fcm" % (w,w))
			mozaic(x, y, box)
		
		#print "minmax", self.smooth.max(), self.smooth.min()
		#indexedimage(self.smooth.T)
		#draw()
		
		#for i in range(dim):
		#	for j in range(dim):
		#		select(i, dim-1-j)
		#self.plot_parameters(0, 1, scope)
		addlabels = True
		objs = []
		names = []
		n = len(self.subscopes)
		for color, subscope in zip(self.colors, self.subscopes):
			print scope
			subscope.load()
			scope = subscope.subscope
			self.parameterset = scope["parameterset"]
			self.load()
			obj = self.plot_parameters(0, 1, scope, color=color, addlabels=addlabels, sigmas=2, linewidth="0.5pt")
			#obj = self.plot_parameters(1, 1, scope, color=color, addlabels=addlabels, sigmas=2)
			objs.append(obj)
			if self.plotlegend:
				names.append(scope["dm_name"])
				print scope["dm_name"]
			addlabels = False
		if self.plotlegend:
			n = len(names)
			legend(["line"] * n, names, objs, edgespacing="3mm", location="right,bottom")
			
		rhalf = scope["rhalf"]
		
		#M = 2300 * (rhalf*1000)**1.4
		#vline(log10(M))
		#print M, log10(M), rhalf
		
		#print scope["dm_profile"].enclosed_mass(0.3)
		#print scope["dm_profile"].enclosed_massrhalf)
		if not nodraw:
			draw()

class ParameterSet2dBinding(ParameterSet2d):
	def __init__(self, subscopes=[], colors=[], smooth_param=0.5, plotlegend=False, subhaloslist=None):
		ParameterSet2d.__init__(self, subscopes=subscopes, colors=colors, smooth_param=smooth_param, plotlegend=plotlegend)
		self.subhaloslist = subhaloslist
		
	def run(self, argv, opts, scope):
		
		halos = self.subhaloslist.load()
		try:
			mainhalo = scope["mainhalo"]
		except:
			mainhalo = None
		try:
			plothalos= scope["plothalos"]
		except:
			plothalos = True
		try:
			plothalotransform = scope["plothalotransform"]
		except:
			plothalotransform = True
		try:
			masscontours= scope["masscontours"]
		except:
			masscontours = True
			
		try:
			allow_massloss = scope["allow_massloss"]
		except:
			allow_massloss = False
			
		if mainhalo:
			halos = halos[halos.halo == "Aq-%s-2" % mainhalo]
		#halos = halos[halos.rho0 != -1] # filter out halos that don't fit
		halos = halos[halos.rho2 != -1] # filter out halos that don't fit
		#halos = halos[halos.rho0_hay_p != -1]
		#halos = halos[halos.rho0_hay != -1]
		halos = halos[halos.dist < 280] # radial cut (as in Carlos' paper)
		#Mvmin = -13.3
		#halos = halos[halos.Mv > Mvmin] 
		#draw()
		import mab
		
		trackhalos = []
		
		try:
			energies = scope["energies"]
		except:
			energies = False
		try:
			energy_histogram = scope["histogram"]
		except:
			energy_histogram = False
		
		#print energy_histogram, energies
		#dsa
		
		r3 = scope["r3"]
		Mv = scope["Mv"]
		Mvmin = Mv - 1
		Mvmax = Mv + 1
		Mvmin = scope["Mvmin"]
		Mvmax = scope["Mvmax"]
		print r3

		mask_disrupted = halos.rtid < 0.002 * r3  #(halos.rt_hay/halos.rs_hay) < 1.
		#mask_badfit = halos.rt_hay/halos.rs_hay > 50
		mask_badfit = log10(halos.Mtid) < 1e-10  # not in use
		mask_ok = ~(mask_disrupted | mask_badfit)
		if any(mask_badfit):
			print halos[mask_badfit].halo
			print halos[mask_badfit].treepos
		#import pdb
		#pdb.set_trace()
		
		

		if 1:
			m3s = []
			slope3s = []
			#dsa
			#halos = mab.cvsfile.CsvObject(trackhalos)
			#halos = self.subhaloslist.halos
			for halo in halos:
				nfw = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
				halo.M3 = nfw.enclosed_mass(r3)
				halo.slope3 = nfw.logslope(r3)
				if 1:
					einasto = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
					halo.M3 = einasto.enclosed_mass(r3)
					halo.slope3 = einasto.logslope(r3)
				if 1:
					r3 = scope["r3"]
					print "r3 =", r3
					name = scope["object_name"]
					if name == "sxt": # different use of names
						name = "sex"
					print name
					r3key = "g" + name +"(r3)"# name for the slope column
					r3_nearest = r3
					print r3key
					
					print halos
					print halos.g300
					
					non_param_slopes = getattr(halo, r3key)
					halo.slope3 = non_param_slopes
					halo.M3 = getattr(halo, "m" + name + "(r3)")
					
					
			indices = argsort(halos.slope3)
			print indices
			halos = halos[indices]
			#for halo in halos:
			#	if (halo.Mv > Mv-0.5) & (halo.Mv < Mv+0.5):
			#		print halo.halo, halo.treepos, log10(halo.M3), halo.slope3
			#		trackhalos.append(halo)
			#trackhalos = mab.cvsfile.CsvObject(trackhalos)
			
		if name == "sex": # different use of names
			name = "sxt"
		
		
		short = True
		w = 7
		if short:
			if energies:
				document("%fcm,%fcm" % (w*(2+0.2),w*(1+0.2)))
				mozaic(2, 1, box)
			else:
				document("%fcm,%fcm" % (w*(1+0.2),w*(1+0.2)))
				page()
				mozaic(1, 1, box)
		else:
			if len(trackhalos) > 7:
				Nx = 4
			else:
				Nx = 3
			document("%fcm,%fcm" % (w*(Nx+0.2),w*(Nx+0.2)))
			mozaic(Nx, Nx, box)
		
		ParameterSet2d.run(self, argv, opts, scope, nodraw=True, init=False)
		#ylabel("&gamma;<sub>-3</sub>")
		
		#Mvmin = Mv - 2
		#Mvmax = min(max(halos.Mv), Mv + 2)
		#mask = (halos.Mv 
		if plothalos:
			Mvminplot, Mvmaxplot = -8, -16
			#Mvmin, Mvmax = -10, -14
			#color = halos.Mv / 
			#scatter(log10(halos.M3), halos.slope3, color=color)
			clearautolegend()
			dotsbw = scatter(log10(halos.M3), halos.slope3, color="black")
			#print color[mask]
			colormapbase = kaplot.colormaps["rainbow"]
			#colormapbase = kaplot.ColorMap([kaplot.Color(1., 0.2, 0.2), kaplot.Color.lightgreen, kaplot.Color.royalblue], interpolate=True)
			colormapbase = kaplot.ColorMap([kaplot.Color(1., 0.2, 0.2), kaplot.Color.green, kaplot.Color.cyan, kaplot.Color.royalblue], interpolate=True)
			#mapcolors = [colormapbase(f) for f in arange(0.2, 0.9, 0.1)]
			mapcolors = [colormapbase(f) for f in arange(0.0, 1.0+0.001, 1/8.)]
			kaplot.colormaps["heat2"] = kaplot.ColorMap(mapcolors, interpolate=False)
			colormap = "heat2"
			indices = argsort(halos.Mv)[::-1]
			halos = halos[indices]
			del indices
			
			color = halos.Mv #(halos.Mv - Mvmin) / (Mvmax - Mvmin)
			print color
			mask = (halos.Mv < Mvmin) & (halos.Mv > Mvmax)
			#scatter(log10(halos.M3[mask]), halos.slope3[mask], colors=color[mask], colormap=colormap, symbolsize="25pt", datamin=Mvmin, datamax=Mvmax)
			
			magnituderangemask = (halos.Mv > Mvmin) & (halos.Mv < Mvmax)# | (halos.Mv < 1000)
			
			print sum(mask_badfit & magnituderangemask)
			#dsa
			#import pdb
			#pdb.set_trace()
			if sum(mask_badfit & magnituderangemask):
				scatter(log10(halos.M3)[mask_badfit & magnituderangemask], halos.slope3[mask_badfit & magnituderangemask], color="blue", symbolName="circle")
			if sum(mask_disrupted & magnituderangemask):
				scatter(log10(halos.M3)[mask_disrupted & magnituderangemask], halos.slope3[mask_disrupted & magnituderangemask], color="red", symbolName="circle")
			trackhalos = halos[magnituderangemask]
			print "*******"
			print "halo, Mv, treepos, sunfind, vmax, distance, rs, rt/rs, log m3, slope3"
			for halo in trackhalos:
				print halo.halo, halo.Mv, halo.treepos, halo.subfind, halo.Vmax, halo.dist, halo.rs_hay, halo.rt_hay/halo.rs_hay, log10(halo.M3), halo.slope3
			print "*******"
			if 1:
				dotscolor = scatter(log10(halos.M3[magnituderangemask]), halos.slope3[magnituderangemask], colors=color[magnituderangemask], colormap=colormap, symbolsize="5pt", datamin=Mvminplot, datamax=Mvmaxplot, symbolName="trianglesolid")
				#scatter(log10(halos.M3[magnituderangemask]), halos.slope3[magnituderangemask], color="black", symbolsize="10pt", symbolName="plus", linewidth="0.2pt", datamin=Mvmin, datamax=Mvmax)
				#scatter(log10(halos.M3[magnituderangemask]), halos.slope3[magnituderangemask], color="black", symbolsize="35pt", symbolName="circlesmall", linewidth="0.2pt", datamin=Mvmin, datamax=Mvmax)
			#legend(["symbol"]*2, "all halos", "halos matching Mv bin", [dotsbw, dotscolor])
			if name == "scl":
				autolegend("all halos", "halos matching M<sub>V</sub> bin", location="right,bottom", edgespacing="4mm", spacing="2mm", fontsize="8pt")
			#for x, y in zip(log10(halos.M3[magnituderangemask]), halos.slope3[magnituderangemask]):
			#	vline(x, linewidth="0.1pt")
			#	hline(y, linewidth="0.1pt")
			print log10(halos.M3[magnituderangemask])
			print halos.slope3[magnituderangemask]
			print halos.Mv[magnituderangemask]
			print magnituderangemask
			print "%" * 10
			#dsa
			print arange(Mvmin, Mvmax, -0.5)
			print halos.Mv[mask]
			print mask
			#sys.exit(0)
			print color[mask]
			innercolorbar(label="M<sub>V</sub>", colormap=colormap, direction='right', location='right, top', edgespacing="5mm", levels=arange(Mvminplot, Mvmaxplot-0.1, -1.0), datamin=Mvminplot, datamax=Mvmaxplot, size="4cm,0.5cm", markers=[Mv,  Mvmax,   Mvmin], markercolor="orange")
			
		Name = name[0].upper() + name[1:]
		title("%s, M<sub>V</sub> = %.2f, r<sub>-3</sub>=%0.2f kpc " % (Name, Mv, r3))
		#show()
			#scatter(log10(trackhalos.M3), trackhalos.slope3, color="red", symbolsize="20pt")
			
		if plothalotransform or energies:
			#bestfit_density = scope["dm_density"]
			#print bestfit_density
			#symbol(log10(bestfit_density.enclosed_mass(r3)), bestfit_density.logslope(r3), symbolName="plus", color="blue")
			efficiencies = zeros((len(trackhalos), len(self.subscopes))) * 1.0
			for i, halo in enumerate(trackhalos):
				density_halo = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
				#r_end = 1e5
				r_end = halo.rtid
				#r_end = halo.rt_hay * 2
				halo_mass = density_halo.enclosed_mass(r_end)
				halo_mass_r3 = density_halo.enclosed_mass(r3)
				#if log10(halo_mass_r3) > 7.25:
				#	continue
				#if density_halo.logslope(r3) > -2.2:
				#	continue
				#halo_mass *= 0.5
				#halo_mass = 1.0

				pot = mab.gd.potential.ProfileNumerical1dWrapper(200, density_halo)
				to_erg = 1.9891e43
				E0 = pot.enclosed_binding_energy(r_end) * to_erg
				print E0

				if not short:
					index = 2 + i
					#print i, index%3, index/3
					select(index%Nx, index/Nx)
					logr = arange(-2, 2, 0.1)
					rs = arange(0.1/2, 2, 0.1/2)
					#mass = density_halo.enclosed_mass(10**logr)
					mass = density_halo.enclosed_mass(rs)
					G = pot.G
					#graph(logr, sqrt(G*mass/10**logr))
					graph(rs, sqrt(G*mass/rs))
					#vline(log10(r3))
					vline(r3)
					select(1, 0)
					graph(rs, sqrt(G*mass/rs))
					ylim(0, 45)
					xlim(0, 2.)
					select(0, 0)
					
					
					#continue
				

				Mstar = 10**(-0.4*(Mv-4.84))
				Esn = 1e51
				Mmean = 0.4
				fsn = 0.0037
				Nsn = Mstar/Mmean * fsn
				for j, (color, subscope) in enumerate(zip(self.colors, self.subscopes)[:]):
					subscope.load()
					#import pdb
					#pdb.set_trace()
					density = subscope.subscope["dm_density"]
					target_mass_r3 = density.enclosed_mass(r3)
					target_slope_r3 = density.logslope(r3)
					
					if isinstance(density, mab.gd.potential.NFW):
						if 0:
							density_fan = mab.gd.potential.NFWCut(density.rho0, density.rs, density.rs*i)
							rs_start = density_fan.rs
							Mtarget = density_fan.enclosed_mass(r_end)
							xs = []
							ys = []
							for k in range(10):
								density_fan.rs = rs_start * 10**(0-k/5.)
								density_fan.rho0 *= Mtarget/density_fan.enclosed_mass(r_end)
								density_fan.update()
								xs.append(log10(density_fan.enclosed_mass(r3)))
								ys.append(density_fan.logslope(r3))
								print ">>> mass", Mtarget, density_fan.enclosed_mass(r_end)
							graph(xs, ys, color=color, alpha=0.8, linestyle="dash")
							print xs, ys
						
						
					if isinstance(density, mab.gd.potential.NFW):
						print "nfw"
						density_mod = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
						#density_mod = None
					elif isinstance(density, mab.gd.potential.TwoSlopeDensityWrapper):
						print "two slope"
						#alpha, beta, rs, gamma=1., rho0=1., rte=1):
						density_mod =  mab.gd.potential.TwoSlopeDensityCut(density.alpha, density.beta, density.rs, density.gamma, density.rho0, halo.rt_hay)
						print "tidal", halo.rt_hay, halo.rs_hay
						#import pdb
						#pdb.set_trace()
					else:
						print "no supported"
						sys.exit(0)
						
					#log10m3target = log10(bestfit_density.enclosed_mass(r3))
					if density_mod:
						if not allow_massloss:
							# convert the density mod to have the same mass at r_-3 while keeping the total mass the same
							def f(logrs):
								#print "logrs", logrs
								rs = 10**logrs[0]
								density_mod.rs = rs
								density_mod.rho0 *= halo_mass/density_mod.enclosed_mass(r_end)
								#print " ", log10(target_mass_r3), log10(density_mod.enclosed_mass(r3)), rs
								#print "  ", log10(halo_mass), log10(density_mod.enclosed_mass(r_end))
								return (log10(target_mass_r3)- log10(density_mod.enclosed_mass(r3)))**2
							res = scipy.optimize.fmin(f, log10(density.rs))
							#continue
							density_mod.rs = 10**res[0]
							#print density_mod.rs, density_mod.rte, density_mod.rs/density_mod.rte, density_mod.rs/halo.rs_hay
						else:
							def f(x):
								logrs, u = x
								massfraction = arctan(u)/pi + 0.5
								rs = 10**logrs
								#massfraction = 1.
								density_mod.rs = rs
								density_mod.rho0 *= abs(halo_mass/density_mod.enclosed_mass(r_end) * massfraction)
								#print rs, massfraction, target_slope_r3, density_mod.logslope(r3)
								#print " ", log10(target_mass_r3), log10(density_mod.enclosed_mass(r3)), rs
								#print "  ", log10(halo_mass), log10(density_mod.enclosed_mass(r_end))
								return (log10(target_mass_r3)- log10(density_mod.enclosed_mass(r3)))**2 +\
										((target_slope_r3)- (density_mod.logslope(r3)))**2
							res = scipy.optimize.fmin(f, [log10(density.rs), 10])
							#continue
							density_mod.rs = 10**res[0]
							massfraction = arctan(res[1])/pi + 0.5
							density_mod.rho0 *= abs(halo_mass/density_mod.enclosed_mass(r_end) * massfraction)
							
							#import pdb;
							#pdb.set_trace()
							
						density_mod.update()
						#import pdb
						#pdb.set_trace()
						#print log10(density_mod.enclosed_mass(r3)), density_mod.logslope(r3)
						#line(log10(halo.M3), halo.slope3, log10(density_mod.enclosed_mass(r3)), density_mod.logslope(r3), color="lightgrey")
						if plothalotransform:
							if density_mod.enclosed_mass(r3) < 0:
								import pdb;
								pdb.set_trace()
							if density_mod.rs > 1e200:
								import pdb;
								pdb.set_trace()
							print "CATCH", r3, density_mod.rs, (halo.M3), halo.slope3, (density_mod.enclosed_mass(r3)), density_mod.logslope(r3)
							line(log10(halo.M3), halo.slope3, log10(density_mod.enclosed_mass(r3)), density_mod.logslope(r3), color=color, linestyle="dot", linewidth="0.1pt")
							symbol(log10(density_mod.enclosed_mass(r3)), density_mod.logslope(r3), symbolName="circle", color=color, symbolsize="3pt", linewidth="0.1pt")
						if energies:
							#print density_mod.rs
							#print density_mod.rho0
							pot = mab.gd.potential.ProfileNumerical1dWrapper(200, density_mod)
							to_erg = 1.9891e43
							E1 = pot.enclosed_binding_energy(r_end) * to_erg
							dE = (E1 - E0)/2
							efficiency = dE/(Esn*Nsn)
							select(1, 0)
							#print "CATCH", density_mod.rs, density_mod.rho0, E1, efficiency
							print ">>>>>>>>>>>>>>", energy_histogram, efficiency
							if energy_histogram:
								#efficiencies.append(efficiency)
								efficiencies[i,j] = efficiency
							else:
								if dE > 0:
									#symbol(log10(halo_mass_r3), log10(dE), color=color)
									symbol(log10(halo_mass_r3), log10(efficiency), color=color)
								else:
									#symbol(log10(halo_mass_r3), log10(-dE), color=color, symbolName="triangle", symbolsize="4pt")
									symbol(log10(halo_mass_r3), log10(-efficiency), color=color, symbolName="triangle", symbolsize="4pt")
								
							print ">", log10(halo_mass_r3)
							select(0, 0)
						if not short:
							index = 2 + i
							select(index%Nx, index/Nx)
							#mass = density_mod.enclosed_mass(10**logr)
							mass = density_mod.enclosed_mass(rs)
							#graph(logr, sqrt(G*mass/10**logr), color=color)
							graph(rs, sqrt(G*mass/rs), color=color)
							
							vline(log10(r3))
							
							#mass = density.enclosed_mass(10**logr)
							mass = density.enclosed_mass(rs)
							#graph(logr, sqrt(G*mass/10**logr), color=color, linestyle="dash")
							graph(rs, sqrt(G*mass/rs), color=color, linestyle="dash")
							ylim(0, 45)
							xlim(0, 2.)
							select(0, 0)
				if energies:
					select(1, 0)
					select(1, 0)
					#labels("log M(r<sub>-3</sub>)", "Energy (ergs)")
					#labels("log M(r<sub>-3</sub>)", "log SN efficiency")
					select(0, 0)
					
					select(0, 0)
					
							
					
				print
				if 0:
					nfw = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
					slope = nfw.logslope(r3)
					m3 = nfw.enclosed_mass(r3)
					xs = [m3]
					ys = [slope]
					#symbol(log10(m3), slpoe, symbolName="dot", color="blue")
					print "->", 
					#r_end = halo.rt_hay
					r_end = 10**3
					pot = mab.gd.potential.ProfileNumerical1dWrapper(200, nfw)
					to_erg = 1.9891e43
					E0 = pot.enclosed_binding_energy(r_end) * to_erg
					mtotal = nfw.enclosed_mass(r_end)
					
					log10m3target = log10(bestfit_density.enclosed_mass(r3))
					def f(logrs):
						#print "logrs", logrs
						rs = 10**logrs
						nfw.rs = rs
						nfw.rho0 *= mtotal/nfw.enclosed_mass(r_end)
						#print " ", log10m3target, log10(nfw.enclosed_mass(r3))
						return (log10m3target - log10(nfw.enclosed_mass(r3)))**2
					res = scipy.optimize.fmin(f, log10(bestfit_density.rs))
					nfw.rs = 10**res[0]
					print nfw.rs, nfw.rte, nfw.rs/nfw.rte, nfw.rs/halo.rs_hay
					nfw.update()
					symbol(log10(nfw.enclosed_mass(r3)), nfw.logslope(r3), symbolName="plus", color="orange")
					
					pot = mab.gd.potential.ProfileNumerical1dWrapper(200, nfw)
					E1 = pot.enclosed_binding_energy(r_end) * to_erg
					print log10(-E0), log10(-E1)
					dE = E1 - E0
					select(1, 0)
					if dE > 0:
						print "energy needed", dE, log10(dE)
						symbol(log10(m3), log10(dE))
					else:
						print "TOO MUCH energy", dE, log10(-dE)
					select(0, 0)
					if slope > -2.5:
						i = 0
						if 0:
							for s in arange(slope, -1.1, 0.1):
								#print s,
								nfw.setslope_at(r3, s)
								print s, nfw.logslope(r3)
								nfw.rho0 *= mtotal/nfw.enclosed_mass(r_end)
								print mtotal/nfw.enclosed_mass(r_end)
								nfw.update()
								slope = nfw.logslope(r3)
								m3 = nfw.enclosed_mass(r3)
								#symbol(log10(m3), slope, symbolName="dot", color="blue")
								xs.append(m3)
								ys.append(slope)
								if 0:
									try:
										pot = mab.gd.potential.ProfileNumerical1dWrapper(200, nfw)
										E1 = pot.enclosed_binding_energy(r_end) * to_erg
										print i, log10(-E0), log10(-E1)
										dE = E1 - E0
										if dE > 0:
											print "energy needed", dE, log10(dE)
										else:
											print "TOO MUCH energy", dE, log10(-dE)
									except:
										print "cannot do potential?"
								i += 1
							graph(log10(xs), ys, color="blue")
				
		if energies:
			select(1, 0)
			if energy_histogram:
				labels("log SN efficiency", "N")
				for j, (color, subscope) in enumerate(zip(self.colors, self.subscopes)[:]):
					x = efficiencies[:,j]
					mask = x > 0
					if len(x[mask]) > 0:
						print x[mask]
						print log10(x[mask])
						histogram(log10(x[mask]), datamin=-4., datamax=1., bincount=20, color=color)
					if len(x[~mask]) > 0:
						print -x[~mask]
						histogram(log10(-x[~mask]), datamin=-4., datamax=1., bincount=20, color=color, linestyle="dot")
					print x
					xlim(-3, 1)
			else:
				labels("log M(r<sub>-3</sub>)", "log SN efficiency")
			select(0, 0)

		for color, subscope in zip(self.colors, self.subscopes)[:]:
			##subscope.load()
			scope = subscope.subscope
			self.parameterset = scope["parameterset"]
			p1 = self.parameterset.parameter_range_list[0]
			p2 = self.parameterset.parameter_range_list[1]
		#xlim(p1.scale_from_uniform_to_range(-0.8), )
		xcorner = p1.scale_from_uniform_to_range(2.8) - 0.15
		#print self.subscopes[0].subscope["parameterset"].parameter_ranges[0]
		errorbars([xcorner], [-2.5], yerr=[0.2 if name == "sxt" else 0.1])
		symbol(xcorner, -2.5)
				
		if masscontours: # plot constant mass contours
			for color, subscope in zip(self.colors, self.subscopes)[:]:
				subscope.load()
				scope = subscope.subscope
				#scope.reset()
				scope.init()
				assert len(self.parameterset.parameter_range_list) == 2
				#v1 = p1.scale_from_uniform(0.5)
				#v2 = p2.scale_from_uniform(0.5)
				
				du = 0.1
				u = arange(0, 1, du) + du/2
				u1 = arange(-0.5, 2., du) + du/2
				grid = zeros((len(u), len(u1)), dtype=float64)
				
				for i, u1 in enumerate(u1):
					for j, u2 in enumerate(u):
						self.parameterset = scope["parameterset"]
						p1 = self.parameterset.parameter_range_list[0]
						p2 = self.parameterset.parameter_range_list[1]
						v1 = p1.scale_from_uniform(u1)
						v2 = p2.scale_from_uniform(u2)
						dmscope = scope.clone()
						dmscope[p1.name] = v1
						dmscope[p2.name] = v2
						#print p1.name, p2.name
						#import pdb
						#pdb.set_trace()
						#del dmscope.dict["dm_density"]
						density = dmscope["dm_density"]
						
						
						r_tidal = 10.
						#r_tidal = 0.4
						if isinstance(density, mab.gd.potential.NFW):
							#print "nfw"
							#r_tidal = density.rs * 1
							density_mod = mab.gd.potential.NFWCut(density.rho0, density.rs, r_tidal) #halo.rho0_hay, halo.rs_hay, halo.rt_hay)
							#density_mod = None
						elif isinstance(density, mab.gd.potential.TwoSlopeDensityWrapper):
							#print "two slope"
							#alpha, beta, rs, gamma=1., rho0=1., rte=1):
							#r_tidal = density.rs * 1
							density_mod =  mab.gd.potential.TwoSlopeDensityCut(density.alpha, density.beta, density.rs, density.gamma, density.rho0, r_tidal)
							#import pdb
							#pdb.set_trace()
						else:
							print "no supported"
							sys.exit(0)
						grid[j,i] = log10(density_mod.enclosed_mass(inf))
						#print grid[j,i], density.rho0, density.rs, v1
				
				resize = (p1.scale_from_uniform_to_range(-0.5), p2.min), (p1.scale_from_uniform_to_range(2.), p2.max)
				levels = arange(6, 12, log10(2.))
				contour(grid, levels, resize=resize, color=color, linewidth="0.5pt") #, linestyle="dot")
						
						
			
		if 0:
			for color, subscope in zip(self.colors, self.subscopes)[:]:
				##subscope.load()
				scope = subscope.subscope
				assert len(self.parameterset.parameter_range_list) == 2
				v1 = p1.scale_from_uniform(0.5)
				v2 = p2.scale_from_uniform(0.5)
				dmscope = scope.clone()
				dmscope[p1.name] = v1
				dmscope[p2.name] = v2
				refprofile = dmscope["dm_profile"]
				
				du = 0.1
				u = arange(0, 1, du) + du/2
				#u = arange(0, 1, du) + du/2
				u1 = arange(-0.5, 1.5, du) + du/2
				grid = zeros((len(u), len(u1)), dtype=float64)
				grid_virial_ratio = zeros((len(u), len(u)), dtype=float32)
				for i, u1 in enumerate(u1):
					for j, u2 in enumerate(u):
						#p1 = self.parameterset.parameter_range_list[0]
						#p2 = self.parameterset.parameter_range_list[1]
						v1 = p1.scale_from_uniform(u1)
						v2 = p2.scale_from_uniform(u2)
						dmscope = scope.clone()
						dmscope[p1.name] = v1
						dmscope[p2.name] = v2
						profile = dmscope["dm_profile"]
						#print ">>>>>>>", profile.M, v1, v2
						def f(r):
							return 4 * pi * r*  profile.G * (profile.densityr(r) * profile.enclosed_mass(r))
						r3 = dmscope["r3"]
						r3 = 10000
						binding_energy = scipy.integrate.quad(f, 1e-5, r3)[0]
						to_erg = 1.9891e43
						binding_energy = binding_energy * to_erg
						grid[j,i] = binding_energy # * to_erg
						#jeans = dmscope["jeans"]
						#K = jeans.kinetic_energy(r3) * to_erg
						#R = 2 * K / abs(binding_energy)
						#grid_virial_ratio[j,i] = R

						#binding_energy *= 1e13 # Joule to erg
						print v1, v2, "%e" % binding_energy #, R #, f(1e-6), f(1)
				#resize = 
				resize = (p1.min, p2.min), (p1.max, p2.max)
				resize = (p1.scale_from_uniform_to_range(-0.5), p2.min), (p1.scale_from_uniform_to_range(1.5), p2.max)
				
				print resize
				print "g", grid
				print log10(grid)
				levels = [11, 11.5, 12, 12.5, 13]
				#grid[:,3] = 1
				#levels = [11.5, 12, 12.5]
				levels = [9.5, 10, 10.5]
				
				indicator = 55 #.8
				delta=0.2
				
				levels = arange(52, indicator-delta/2, delta)
				contour(log10(grid), levels, resize=resize, color=color)
				levels = [indicator]
				contour(log10(grid), levels, resize=resize, color=color, linestyle="dash")
				levels = arange(indicator+delta, 60, delta)
				contour(log10(grid), levels, resize=resize, color=color)
				
				if 0:
					levels = [1, 3, 4]
					contour(grid_virial_ratio, levels, resize=resize, color=color, linestyle="dash")
					levels = [2]
					contour(grid_virial_ratio, levels, resize=resize, color=color, linestyle="dot")
					print grid_virial_ratio
				#contour(log10(grid), levels, resize=resize, color=color)
				#for 
				#self.load()
		#xlim(7.7, 8.4)
		#xlim(7.1, 7.5)
		select(0, 0)
		xlim(p1.scale_from_uniform_to_range(-0.8), p1.scale_from_uniform_to_range(2.8))
		ylim(-3, -0.5)
		if energies:
			select(1, 0)
			if not energy_histogram:
				xlim(p1.scale_from_uniform_to_range(-0.8), p1.scale_from_uniform_to_range(2.8))
				hline(0.0)
				ylim(-3, 1)
			grow(y=1.1)
		#select(1, 0)
		#grow(1.2)
		#flimp
		draw()
		
class PlotAquariusHalo(object):
	def __init__(self, filename=None):
		self.filename = filename
		
	def run(self, argv, opts, scope):
		if self.filename == None:
			filename = argv[1]
		else:
			filename = self.filename
			
		G = mab.gd.potential.G

		print "reading", filename
		halo = mab.asciifile.readsimple(filename)
		rs = halo.r
		rs = np.sort(rs)
		print rs
		M = (arange(len(rs))+1)*2.87e5/0.73 # mass or particle over h0 (hubble const)
		vcirc = (G * M / rs)**0.5
		#box()
		mozaic(2,1,box)
		graph(log10(rs), log10(vcirc))
		
		logrs = arange(-3, 3, 0.1)
		rs = 10**logrs
		if 1:
			pot1 = mab.gd.potential.NFW(1e10, 1.)
			pot2 = mab.gd.potential.NFW(1e12, 1.)
			pot3 = mab.gd.potential.NFW(1e10, 10.)
			pot4 = mab.gd.potential.NFW(1e12, 10.)

			pots = [pot1, pot2, pot3, pot4]
			for pot, color in zip(pots, nicecolors):
				M = pot.enclosed_mass(rs)
				vcirc = (G * M/rs)**0.5
				graph(log10(rs), log10(vcirc), color=color)

		#pot = mab.gd.potential.NFW(1.5e12, 6.)
		#M = pot.enclosed_mass(rs)
		#vcirc = (G * M/rs)**0.5
		#graph(log10(rs), log10(vcirc), color="orange")
		
		rs = halo.r
		rs = np.sort(rs)[10000::10]
		Mhalo = M = (arange(len(rs))+1)*2.87e5/0.73 # mass or particle over h0 (hubble const)
		def f(x):
			pot.rs = 10**x[0]
			pot.rho0 = 10**x[1]
			err = sum( (log10(pot.enclosed_mass(rs)) - log10(M))**2 )
			print x, err
			return err
			#return sum((log10(pot.enclosed_mass(rs)) - log10(M))**2)
			#return sum((log10(pot.densityr(rs)) - log10(M))**2)
			
		x = scipy.optimize.fmin(f, [log10(pot.rs), log10(pot.rho0)])
		pot.rs = 10**x[0]
		pot.rho0 = 10**x[1]
	
		M = pot.enclosed_mass(rs)
		vcirc = (G * M/rs)**0.5
		color="orange"
		graph(log10(rs), log10(vcirc), color=color)
		
		select(1,0)
		graph(log10(rs), log10(Mhalo), color="black")
		graph(log10(rs), log10(M), color=color)
		
		draw()
		return 1
		
		
		
		
class PlotAquariusHalos(object):
	def __init__(self, subhaloslist):
		self.subhaloslist = subhaloslist
		
		
	def gethalos(self, scope):
		halos = self.subhaloslist.load()
		halos = halos[halos.dist < 280]
		#halos = halos[halos.rho0 != -1]
		#halos = halos[halos.rho0_hay != -1]
		halos = halos[halos.rs_hay > 1e-10]
		#halos = halos[halos.rho0_hay_l != -1]
		try:
			mainhalo = scope["mainhalo"]
		except:
			mainhalo = None
		if mainhalo:
			halos = halos[halos.halo == "Aq-%s-2" % mainhalo]
		#print `halos.halo`
		#print (halos.halo == "Aq-B-2") | (halos.halo == "Aq-A-2")
		#halos = halos[(halos.halo == "Aq-B-2") | (halos.halo == "Aq-A-2") | (halos.halo == "Aq-D-2")]
		#halos = halos[(halos.halo == "Aq-A-2")]
		for halo in halos:
			halo.density_hay = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
			halo.density_nfw = mab.gd.potential.NFW(1., halo.rs)
			halo.density_nfw.rho0 = halo.rho0
			halo.density_ein = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
			halo.density_moo = mab.gd.potential.Moore(halo.rho0_moo, halo.rs_moo, halo.rt_moo, halo.rd_moo)
			halo.density_bre = mab.gd.potential.Breddels(halo.rho0_bre, halo.rs_bre, halo.alpha_bre, halo.rt_bre, halo.rd_bre)
			#densities = [density_nfw, density_hay, density_ein]
			#densities = [density_nfw, density_hay, density_ein]
		return halos
		
		
class PlotAquariusHalosEnergies(PlotAquariusHalos):
	def __init__(self, subhaloslist):
		PlotAquariusHalos.__init__(self, subhaloslist)

	def run(self, argv, opts, scope):
		halos = self.gethalos(scope)
		Mvmin = scope["Mvmin"]
		Mvmax = scope["Mvmax"]
		magnituderangemask = (halos.Mv > Mvmin) & (halos.Mv < Mvmax)
		r3 = scope["r3"]
		size = "10cm,16cm"
		document(size=size)
		mozaic(1,2,box)
		#mapping = {300: "g300", 500:"g500", 1000:"g1000"}
		#r3s = array(mapping.keys())
		#index = argsort(abs(r3s-r3*1000))[0]
		#r3_nearest = r3s[index]/1000.
		#print "nearest r3 to ", r3, "is", r3s[index]
		#r3key = mapping[r3s[index]]
		#print r3key
		#print halos
		#print halos.g300
		
		#non_param_slopes = halos[r3key]
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		vfill(Mvmin, Mvmax, color=color)
		#		E1 = density_halo .enclosed_binding_energy(halo.rtid)# * to_erg
		#		return -E1
		#	Ws = calc_W(halos)
		#	#scatter(log10(halos.Wtid*to_erg)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
		#	scatter(log10(halos.Mtid)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
		
		to_erg = 1.9891e43
		
		scatter(halos.Mv, log10(halos.Wtid * to_erg) )
		labels("Mv", "log W/erg")
		
		select(0,1)
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		vfill(Mvmin, Mvmax, color=color)
		hline(0)
		for color, j in zip(nicecolors[1:], range(3))[:]:
			@mab.parallelize.timed
			@mab.parallelize.parallelize()
			def do(i):
				halo = halos[i]
				density = [halo.density_nfw, halo.density_hay, halo.density_ein][j]#, halo.density_moo] #, halo.density_bre]
				#print density, r3_nearest
				obj = density
				# not all density objects have analytic energy expressions
				if not hasattr(density, "enclosed_binding_energy"):
					obj = mab.gd.potential.ProfileNumerical1dWrapper(200, density)
				return -obj.enclosed_binding_energy(halo.rtid)
			Ws = array(do(range(len(halos)))) * to_erg
			refWs = halos.Wtid * to_erg
			#print refWs, Ws
			if any(isnan(log10(Ws))):
				import pdb
				pdb.set_trace()
			converged = True
			symbolsize = "1.2pt"
			
			
			symbols(halos.Mv, log10(refWs)-log10(Ws), color=color, symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
			print log10(refWs)-log10(Ws)
		
		ylim(-0.15, 0.15)
		labels("Mv", "log Wref/W ")
		#title("r3 = %.2f, r3_near = %.2f" % (r3, r3_nearest))
		draw()




class PlotAquariusHalosSlopes(PlotAquariusHalos):
	def __init__(self, subhaloslist):
		PlotAquariusHalos.__init__(self, subhaloslist)
		
	def run(self, argv, opts, scope):
		halos = self.gethalos(scope)
		Mvmin = scope["Mvmin"]
		Mvmax = scope["Mvmax"]
		magnituderangemask = (halos.Mv > Mvmin) & (halos.Mv < Mvmax)
		r3 = scope["r3"]
		print "r3 =", r3
		
		size = "10cm,16cm"
		#box(size=size)
		document(size=size)
		mozaic(1,2,box)
		name = scope["object_name"]
		
		if name == "sxt": # different use of names
			name = "sxt" 
		print name
		#dsa
		#mapping = {300: "g300", 500:"g500", 1000:"g1000"}
		#mapping = {300: "g300", 500:"g500", 1000:"g1000"}
		#r3s = array(mapping.keys())
		#index = argsort(abs(r3s-r3*1000))[0]
		#r3_nearest = r3s[index]/1000.
		#print "nearest r3 to ", r3, "is", r3s[index]
		#r3key = mapping[r3s[index]]
		#r3key = {"fnx":"gfnx", "scl":"gscl", 
		r3key = "g" + name # name for the slope column
		r3_nearest = r3
		print r3key
		
		print halos
		print halos.g300
		
		non_param_slopes = getattr(halos, r3key)
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		vfill(Mvmin, Mvmax, color=color)
		scatter(halos.Mv, non_param_slopes)
		labels("Mv", "slope")
		
		select(0,1)
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		vfill(Mvmin, Mvmax, color=color)
		labels("Mv", "delta slope")
		xs = []
		ys = []
		hline(0)
		for halo in halos:
			#densities = [halo.density_nfw, halo.density_hay, halo.density_ein, halo.density_moo] #, halo.density_bre]
			densities = [halo.density_nfw, halo.density_hay, halo.density_ein]#, halo.density_moo] #, halo.density_bre]
			refslope = halo[r3key]
			symbolsize = "2.2pt"
			for color, density in zip(nicecolors[1:], densities[:]):
				converged = True
				#print density, r3_nearest
				y = refslope-density.logslope(r3_nearest)
				symbol(halo.Mv, y, color=color, symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
				if density == halo.density_ein:
					if abs(y) < 3:
						xs.append(halo.Mv)
						ys.append(y)
		N = len(xs)
		Nwindow = 15
		xsmooth = []
		ysmooth = []
		errors = []
		xs = array(xs)
		ys = array(ys)
		indices = argsort(xs)
		xs = xs[indices]
		ys = ys[indices]
		for i in range(N-Nwindow):
			meanx = mean(xs[i:i+Nwindow])
			meany = mean(ys[i:i+Nwindow])
			error = std(ys[i:i+Nwindow])
			xsmooth.append(meanx)
			ysmooth.append(meany)
			errors.append(error)
		#select(1, 0)
		ysmooth = array(ysmooth)
		xsmooth = array(xsmooth)
		error = array(error )
		graph(xsmooth, ysmooth, color="blue")
		graph(xsmooth, ysmooth+errors, color="blue", linestyle="dot")
		graph(xsmooth, ysmooth-errors, color="blue", linestyle="dot")
		ylim(-3, 3)
		title("r3 = %.2f, r3_near = %.2f" % (r3, r3_nearest))
		
		
		
		
		
		draw()
		
		
		
class PlotAquariusHalosEnergiesTransform(PlotAquariusHalos):
	def __init__(self, subhaloslist, subscopes, colors):
		PlotAquariusHalos.__init__(self, subhaloslist)
		self.subscopes = subscopes
		self.colors = colors

	def run(self, argv, opts, scope):
		halos = self.gethalos(scope)
		r3 = scope["r3"]
		M3sum = 0
		for subscope in self.subscopes:
			subscope.load()
			density = subscope.subscope["dm_density"]
			M3sum += density.enclosed_mass(r3)
		M3 = M3sum/len(self.subscopes)
		halos = halos[halos.Mtid > M3]
		
		Mvmin = scope["Mvmin"]
		Mvmax = scope["Mvmax"]
		Mv = scope["Mv"]
		#Mvmin = Mv+1
		#Mvmax = Mv-1
		#ewrwe
		try:
			plot_eff = scope["plot_eff"]
		except:
			plot_eff = False

		name = scope["object_name"]
		if name == "sxt": # different use of names
			name = "sxt" 
		haxis = name == "sxt"
		
		
		magnituderangemask = (halos.Mv > Mvmin) & (halos.Mv < Mvmax)
		
		size = "10cm,7cm" if haxis else "10cm,6cm"
		document(size=size)
		box(viewport=((0, 1./7 if haxis else 0.), (1.,1.)) )
		current.container.drawOutsideRight = True #scope["doleft"]
		current.container.drawOutsideTop = True #scope["doleft"]
		#current.container.drawOutsideLeft = True #scope["doleft"]
		current.container.drawOutsideBottom = True #not haxis #scope["dobottom"]
		#current.container.leftaxes[0].drawLabels = scope["doleft"]
		current.container.bottomaxes[0].drawLabels = haxis #scope["dobottom"]
		
		#mozaic(1,2,box)
		#select(0,1)
		
		stripped = halos.rtid < 2 * r3
		
		#mapping = {300: "g300", 500:"g500", 1000:"g1000"}
		#r3s = array(mapping.keys())
		#index = argsort(abs(r3s-r3*1000))[0]
		#r3_nearest = r3s[index]/1000.
		#print "nearest r3 to ", r3, "is", r3s[index]
		#r3key = mapping[r3s[index]]
		#print r3key
		#print halos
		#print halos.g300
		
		#non_param_slopes = halos[r3key]
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		#if not plot_eff:
		vfill(Mvmin, Mvmax, color=color)
		#if not plot_eff:
		hline(log10(1.0))
		#		E1 = density_halo .enclosed_binding_energy(halo.rtid)# * to_erg
		#		return -E1
		#	Ws = calc_W(halos)
		#	#scatter(log10(halos.Wtid*to_erg)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
		#	scatter(log10(halos.Mtid)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
		
		to_erg = 1.9891e43 # m_sol * km^2 / sec^2 in ergs
		symbolsize = "3pt"
		#scatter(halos.Mv, log10(halos.Wtid * to_erg), symbolName="circlesolid", symbolsize=symbolsize )
		if haxis:
			xlabel("M<sub>V</sub>")
		if plot_eff:
			ylabel("log &epsilon;<sub>SN</sub>")
		else:
			ylabel("log W<sub>dSph</sub>/W<sub>sub</sub>")
		
		#select(0,0)
		color = kaplot.Color.white.blend(kaplot.Color.blue, 0.8)
		#vfill(Mvmin, Mvmax, color=color)
		#hline(0)
		Mstar = 10**(-0.4*(Mv-4.84))
		Esn = 1e51
		Mmean = 0.4
		fsn = 0.0037
		Nsn = Mstar/Mmean * fsn
		#hline(log10(Nsn*Esn))
		#legend_types, labels, current.legend_objects
		if name == "scl" and not plot_eff:
			fontsize = "6pt"
			green = text("", -1000, -1000, color="green", symbolsize=fontsize)
			red = text("", -1000, -1000, color="red", symbolsize=fontsize)
			black = text("", -1000, -1000, color="black", symbolsize=fontsize)
			legend(["circlesolid", "trianglesolid", "circlesolid", "trianglesolid", "circlesolid", "circle"], ["NFW (r<sub>95</sub>)", "NFW (r<sub>50</sub>)", "core13 (r<sub>95</sub>)", "core13 (r<sub>50</sub>)"], [red, red, green, green, red, red, black, black], edgespacing="4.5mm", borderspacing="0mm", symbolsize=fontsize, fontsize=fontsize, location="right,bottom")
		#draw()
		#return
		
		for j, (color, subscope) in enumerate(zip(self.colors, self.subscopes)[:]):
			subscope.load()
			#import pdb
			#pdb.set_trace()
			#for symbolName, Wname, rname in zip(["circlesolid", "squaresolid", "trianglesolid"], ["W(r95)", "W(r80)", "W(r50)"], ["r95", "r80", "r50"]):
			for symbolName, Wname, rname in zip(["circlesolid", "trianglesolid"], ["W(r95)", "W(r50)"], ["r95", "r50"]):
				density = subscope.subscope["dm_density"]
				@mab.parallelize.timed
				@mab.parallelize.parallelize()
				def do(i):
					halo = halos[i]
					rmax = getattr(halo, rname)
					#density = [halo.density_nfw, halo.density_hay, halo.density_ein][j]#, halo.density_moo] #, halo.density_bre]
					#print density, r3_nearest
					obj = density
					# not all density objects have analytic energy expressions
					if not hasattr(density, "enclosed_binding_energy"):
						obj = mab.gd.potential.ProfileNumerical1dWrapper(200, density)
					return -obj.enclosed_binding_energy(rmax)
				Ws = array(do(range(len(halos)))) * to_erg
				refWs = getattr(halos, Wname) * to_erg
				print refWs, Ws
				if any(isnan(log10(Ws))):
					import pdb
					pdb.set_trace()
				converged = True
				print halos[0][rname], Ws[0], refWs[0]
				print halos[-1][rname], Ws[-1], refWs[-1]
				#symbolsize = "2.2pt"
				for i in range(len(Ws)):
					dE = (Ws[i] - refWs[i])/2.
					import pdb
					#pdb.set_trace()
					f = -dE / (Esn*Nsn)			
					sign = 1
					if dE > 0:
						#symbol(halos.Mv[i], log10(dE), color=color, symbolName="circlesolid", symbolsize=symbolsize)
						pass
					else:
						sign = -1
						if plot_eff:
							symbol(halos.Mv[i], log10(abs(f)), color=color, symbolName=symbolName, symbolsize=symbolsize, linewidth="0.4pt")
						#if stripped[i]:
						#if plot_eff:
						#	symbol(halos.Mv[i], log10(abs(f)), color="black", symbolName="plus", symbolsize="5pt")
						
				
					#select(0,1)
					if not plot_eff:
						symbol(halos.Mv[i], log10(Ws[i]/refWs[i]), color=color, symbolName=symbolName if (f) <= 1 else symbolName[:-5], symbolsize=symbolsize, linewidth="0.4pt")
					#select(0,0)
			print log10(refWs)-log10(Ws)
		
		#for j, (color, subscope) in enumerate(zip(self.colors, self.subscopes)[:]):
		if 0:
			color = "blue"
			subscope = self.subscopes[0] # nfw
			subscope.load()
			#import pdb
			#pdb.set_trace()
			density = subscope.subscope["dm_density"]
			M3 = density.enclosed_mass(r3)
			slope3 = density.logslope(r3)
			density = mab.gd.potential.Einasto.from_r3(M3, r3, slope3, alpha=0.22)
			@mab.parallelize.timed
			@mab.parallelize.parallelize()
			def do(i):
				halo = halos[i]
				#density = [halo.density_nfw, halo.density_hay, halo.density_ein][j]#, halo.density_moo] #, halo.density_bre]
				#print density, r3_nearest
				obj = density
				# not all density objects have analytic energy expressions
				if not hasattr(density, "enclosed_binding_energy"):
					obj = mab.gd.potential.ProfileNumerical1dWrapper(200, density)
				return -obj.enclosed_binding_energy(halo.rtid)
			Ws = array(do(range(len(halos)))) * to_erg
			refWs = halos.Wtid * to_erg
			print refWs, Ws
			if any(isnan(log10(Ws))):
				import pdb
				pdb.set_trace()
			converged = True
			#symbolsize = "2.2pt"
			for i in range(len(Ws)):
				dE = (Ws[i] - refWs[i])/2.
				f = -dE / (Esn*Nsn)			
				sign = 1
				if dE > 0:
					#symbol(halos.Mv[i], log10(dE), color=color, symbolName="circlesolid", symbolsize=symbolsize)
					pass
				else:
					sign = -1
					#symbol(halos.Mv[i], log10(-dE), color=color, symbolName="circle", symbolsize=symbolsize, linewidth="0.4pt")
					symbol(halos.Mv[i], log10(abs(f)), color=color, symbolName="circlesolid", symbolsize=symbolsize, linewidth="0.4pt")
				#if stripped[i]:
				#	symbol(halos.Mv[i], log10(dE*sign), color="black", symbolName="plus", symbolsize="5pt")
				select(0,1)
				symbol(halos.Mv[i], log10(Ws[i]/refWs[i]), color=color, symbolName="circlesolid" if (f) <= 1 else "circle", symbolsize=symbolsize, linewidth="0.4pt")
				select(0,0)
			print log10(refWs)-log10(Ws)
		if name == "fnx" and not plot_eff:
			#select(0,1)
			text("Need to lose energy", -8, 2.2)
			text("Need energy input",   -8, -1.5)
			#select(0,0)
		#legend(
		#ylim(-0.15, 0.15)
		#labels("Mv", "log10 abs(Ws - Wref)")
		#labels("M<sub>V</sub>", "log &epsilon;<sub>SN</sub>")
		#title("r3 = %.2f, r3_near = %.2f" % (r3, r3_nearest))
		#select(0, 1)
		#title(name[0].upper() + name[1:])
		Name = name[0].upper() + name[1:]
		if plot_eff:
			text(Name, -14.5, 1.9)
			ylim(-4.99, 2.5)
		else:
			text(Name, -14.5, 2.5)
			ylim(-2.5, 3.2)
		xlim(halos.Mv.min(), -5)
		#select(0, 0)
		#xlim(halos.Mv.min(), halos.Mv.max())
		draw()

		
		
		
class PlotAquariusHalosTidalRadius(object):
	def __init__(self, subhaloslist):
		self.subhaloslist = subhaloslist
		
	def run(self, argv, opts, scope):
		halos = self.subhaloslist.load()
		
		halos = halos[halos.rho0 != -1]
		halos = halos[halos.rho0_hay != -1]
		halos = halos[halos.rho0_hay_p != -1]
		try:
			mainhalo = scope["mainhalo"]
		except:
			mainhalo = None
		if mainhalo:
			halos = halos[halos.halo == "Aq-%s-2" % mainhalo]
		print `halos.halo`
		print (halos.halo == "Aq-B-2") | (halos.halo == "Aq-A-2")
		halos = halos[(halos.halo == "Aq-B-2") | (halos.halo == "Aq-A-2") | (halos.halo == "Aq-D-2")]
		#halos = halos[(halos.halo == "Aq-A-2")]
		
		r3 = scope["r3"]
		
		#box()
		document(size="20cm,30cm")
		mozaic(2,3,box)
		symbolsize="15pt"
		mask_disrupted = (halos.rt_hay/halos.rs_hay) < 1.
		#mask_badfit = halos.rt_hay/halos.rs_hay > 50
		#mask_badfit = halos.rconv > halos.rs_hay
		#mask_badfit = halos.rt_hay/halos.rs_hay < 1.1
		mask_badfit = log10(halos.Mtid) < 8.
		mask_disrupted = log10(halos.Mtid) == 8.
		mask_no_conv = halos.r2 < 1.
		mask_ok = ~(mask_disrupted | mask_badfit)
		
		if 0:
			print halos.rtid
			print halos.rs_hay
			print halos.rtid/halos.rs_hay
			issue = log10(halos.rtid/halos.rs_hay) < -3
			print halos[issue].halo
			print halos[issue].treepos
			
			
			import pdb
			pdb.set_trace()
	
		if 1:
			scatter(log10(halos.rtid/halos.rs_hay)[mask_ok], log10(halos.rt_hay/halos.rs_hay)[mask_ok], symbolsize=symbolsize)
			if any(mask_badfit):
				scatter(log10(halos.rtid/halos.rs_hay)[mask_badfit], log10(halos.rt_hay/halos.rs_hay)[mask_badfit], color="blue", symbolsize=symbolsize)
			if any(mask_disrupted):
				scatter(log10(halos.rtid/halos.rs_hay)[mask_disrupted], log10(halos.rt_hay/halos.rs_hay)[mask_disrupted], color="red", symbolsize=symbolsize)
		else:
			scatter(log10(halos.rtid/halos.rs_moo)[mask_ok], log10(halos.rt_moo/halos.rs_moo)[mask_ok], symbolsize=symbolsize)
			if any(mask_badfit):
				scatter(log10(halos.rtid/halos.rs_moo)[mask_badfit], log10(halos.rt_moo/halos.rs_moo)[mask_badfit], color="blue", symbolsize=symbolsize)
			if any(mask_disrupted):
				scatter(log10(halos.rtid/halos.rs_moo)[mask_disrupted], log10(halos.rt_moo/halos.rs_moo)[mask_disrupted], color="red", symbolsize=symbolsize)
			
		labels("log r<sub>tid</sub>/r<sub>s,hay</sub>", "log r<sub>t,hay</sub>/r<sub>s,hay</sub>")
		
		x = []
		y = []
		y2 = []
		if 1:
			for logrtid_hay in arange(-0.5, 1.5, 0.15):
				rtid_hay = 10**logrtid_hay * 1
				dens = mab.gd.potential.NFWCut(1., 1., rtid_hay)
				Mtotal = dens.enclosed_mass(inf)
				def f(args):
					r = 10**args[0]
					M = dens.enclosed_mass(r)
					diff = (M/Mtotal - 0.95)**2
					print r, M, M/Mtotal, diff
					return diff
				values = scipy.optimize.fmin(f, [1.])
				logr_tid = values[0]
				y.append(logrtid_hay)
				x.append(logr_tid)
				y2.append(logrtid_hay-logr_tid)
				graph(x, y, color="blue")
				#print logr
				
		
		
			select(0, 1)
			scatter(log10(halos.rtid/halos.rs_hay)[mask_ok], log10((halos.rt_hay/halos.rtid))[mask_ok], symbolsize=symbolsize)
			if any(mask_badfit):
				scatter(log10(halos.rtid/halos.rs_hay)[mask_badfit], log10((halos.rt_hay/halos.rtid))[mask_badfit], color="blue", symbolsize=symbolsize)
			if any(mask_disrupted):
				scatter(log10(halos.rtid/halos.rs_hay)[mask_disrupted], log10((halos.rt_hay/halos.rtid))[mask_disrupted], color="red", symbolsize=symbolsize)
			graph(x, y2, color="blue")
			labels("log r<sub>tid</sub>/r<sub>s,hay</sub>", "log r<sub>t,hay</sub>/r<sub>tid</sub>")
		select(1, 0)
		
		@vectorize
		def calc_mass(halo):
			density_halo = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
			#return density_halo.enclosed_mass(inf)
			return density_halo.enclosed_mass(halo.rtid)
		masses = calc_mass(halos)
		#scatter(log10(halos.Mtid), log10(halos.Mtid/masses))
		Mmin = min([min(halos.Mtid), min(masses)])
		Mmax = max([max(halos.Mtid), max(masses)])
		#line(log10(Mmin), log10(Mmin), log10(Mmax), log10(Mmax), color="red")
		print masses
		print halos.halo[masses<0], halos.treepos[masses<0]
		line(log10(Mmin), 0, log10(Mmax), 0, color="red")
		if 1:
			scatter(log10(halos.Mtid)[mask_ok], log10(masses/halos.Mtid)[mask_ok], symbolsize=symbolsize)
			if any(mask_badfit):
				scatter(log10(halos.Mtid)[mask_badfit], log10(masses/halos.Mtid)[mask_badfit], color="blue", symbolsize=symbolsize)
			if any(mask_disrupted):
				scatter(log10(halos.Mtid)[mask_disrupted], log10(masses/halos.Mtid)[mask_disrupted], color="red", symbolsize=symbolsize)

		# NFW
		@vectorize
		def calc_mass(halo):
			density_halo = mab.gd.potential.NFW(1, halo.rs)
			density_halo.rho0 = halo.rho0
			#density_halo.update()
			#return density_halo.enclosed_mass(inf)
			return density_halo.enclosed_mass(halo.rtid)
		masses = calc_mass(halos)
		scatter(log10(halos.Mtid), log10(masses/halos.Mtid), symbolsize=symbolsize, color="red")
		selection = log10(halos.Mtid) > 9.7
		print "carlos", halos.halo[selection], halos.treepos[selection]
		
		# Einasto
		@vectorize
		def calc_mass(halo):
			density_halo = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
			#return density_halo.enclosed_mass(inf)
			return density_halo.enclosed_mass(halo.rtid)
		masses = calc_mass(halos)
		scatter(log10(halos.Mtid), log10(masses/halos.Mtid), symbolsize=symbolsize, color="green")
		
		
				
		#print math.log(math.e)
		#dsa
		#Mtot = halos.rho0_hay * 4 * pi * (c**2 * (9*c*(-c**3+(c**3+2)*math.log(c)+1)+2*3**0.5 * pi * (2*c+1)*(c-1)**2))*(9*(c**3-1)**2)
		#print log(e)
		#dsa
		#c = halos.rt_hay / halos.rs_hay
		#Mtot = 4 * pi * halos.rs_hay**3 * halos.rho0_hay * (c**2 * (2 * sqrt (3) * (-1 + c)**2 * (1 + 2 * c) * pi + 9 * c * (1 - c**3 + (2 + c**3) * log(c)))) / (9 * (-1 + c**3)**2)
		#scatter(log10(Mtot), log10(masses), color="orange")


			
		labels("log M<sub>tid</sub>", "log M<sub>tid,analytic</sub>")
		
		
		
		select(1,1)
		
		if 0:
			to_erg = 1.9891e43
			
			
			@vectorize
			def calc_W(halo):
				density_halo = mab.gd.potential.NFW(1, halo.rs)
				density_halo.rho0 = halo.rho0
				#pot = mab.gd.potential.ProfileNumerical1dWrapper(200, density_halo)
				E1 = density_halo .enclosed_binding_energy(halo.rtid)# * to_erg
				return -E1
			Ws = calc_W(halos)
			#scatter(log10(halos.Wtid*to_erg)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
			scatter(log10(halos.Mtid)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="red")
			
			@vectorize
			def calc_W(halo):
				density_halo = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
				pot = mab.gd.potential.ProfileNumerical1dWrapper(200, density_halo)
				#E1 = pot.enclosed_binding_energy(halo.rtid)# * to_erg
				E1 = pot.enclosed_binding_energy(halo.rtid)# * to_erg
				#E1 = pot.enclosed_binding_energy(inf)# * to_erg
				return -E1
			Ws = calc_W(halos)
			Wmin = min([min(halos.Wtid), min(Ws)])
			Wmax = max([max(halos.Wtid), max(Ws)])
			#line(log10(Wmin), log10(Wmin), log10(Wmax), log10(Wmax), color="red")
			#line(log10(Wmin*to_erg), 0, log10(Wmax*to_erg), 0, color="red")
			#scatter(log10(halos.Wtid*to_erg)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize)
			scatter(log10(halos.Mtid)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="orange")
			#if any(mask_badfit):
			#	scatter(log10(halos.Wtid*to_erg)[mask_badfit], log10(Ws/halos.Wtid)[mask_badfit], color="blue", symbolsize=symbolsize)
			#if any(mask_disrupted):
			#	scatter(log10(halos.Wtid*to_erg)[mask_disrupted], log10(Ws/halos.Wtid)[mask_disrupted], color="red", symbolsize=symbolsize)
		
			@vectorize
			def calc_W(halo):
				density_halo = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
				pot = mab.gd.potential.ProfileNumerical1dWrapper(200, density_halo)
				#E1 = pot.enclosed_binding_energy(halo.rtid)# * to_erg
				E1 = pot.enclosed_binding_energy(halo.rtid)# * to_erg
				#E1 = pot.enclosed_binding_energy(inf)# * to_erg
				return -E1
			Ws = calc_W(halos)
			scatter(log10(halos.Mtid)[mask_ok], log10(Ws/halos.Wtid)[mask_ok], symbolsize=symbolsize, color="green")

		labels("log Wtid", "log Wanalytic")
		if 0:
			select(0, 1)
			scatter(log10(halos.Wtid)[~mask], log10(Ws/halos.Wtid)[~mask])
			scatter(log10(halos.Wtid)[mask], log10(Ws/halos.Wtid)[mask], color="blue")
			
			
		select(0,2)
		
		i = 0
		if 1:
			symbolsize = "5pt"
			for converged, halo in zip(mask_no_conv, halos):
				density_nfw = mab.gd.potential.NFW(halo.rho0, halo.rs)
				density_nfw.rho0 = halo.rho0
				density_hay = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
				density_ein = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
				density_moo = mab.gd.potential.Moore(halo.rho0_moo, halo.rs_moo, halo.rt_moo, halo.rd_moo)
				density_bre = mab.gd.potential.Breddels(halo.rho0_bre, halo.rs_bre, halo.alpha_bre, halo.rt_bre, halo.rd_bre)
				densities = [density_nfw, density_hay, density_ein, density_moo, density_bre]
				line(log10(halo.Mtid), density_hay.logslope(r3), log10(halo.Mtid), density_ein.logslope(r3), color="lightgrey", linewidth="0.3px")
				for color, density in zip(nicecolors[1:], densities):
					#symbol(log10(halo.rs), density.logslope(r3), color=color)
					#symbol(log10(halo.rtid), density.logslope(r3), color=color)
					#symbol(i, density.logslope(r3), color=color)
					symbol(log10(halo.Mtid), density.logslope(r3), color=color, symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
				symbol(log10(halo.Mtid), halo.g1000, color="grey", symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
				#symbol(log10(halo.Mtid), halo.g300, color="grey", symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
				
				i += 1
			#labels("index", "slope at 1kpc")
			labels("log Mtid", "slope at 1kpc")
			ylim(-7, 0)
		
			select(1,2)
			for converged, halo in zip(mask_no_conv, halos):
				density_hay = mab.gd.potential.NFWCut(halo.rho0_hay, halo.rs_hay, halo.rt_hay)
				density_nfw = mab.gd.potential.NFW(halo.rho0, halo.rs)
				density_nfw.rho0 = halo.rho0
				density_ein = mab.gd.potential.Einasto(halo.rho2, halo.r2, halo.alpha)
				density_moo = mab.gd.potential.Moore(halo.rho0_moo, halo.rs_moo, halo.rt_moo, halo.rd_moo)
				density_bre = mab.gd.potential.Breddels(halo.rho0_bre, halo.rs_bre, halo.alpha_bre, halo.rt_bre, halo.rd_bre)
				#densities = [density_nfw, density_hay, density_ein]
				#densities = [density_nfw, density_hay, density_ein]
				densities = [density_nfw, density_hay, density_ein, density_moo, density_bre]
				#densities = [density_nfw, density_ein, density_moo]
				for color, density in zip(nicecolors[1:], densities):
					#symbol(log10(halo.rs), density.logslope(r3), color=color)
					#symbol(log10(halo.rtid), density.logslope(r3), color=color)
					#symbol(i, density.logslope(r3), color=color)
					logM = log10(halo.Mtid)
					#refslope = halo.g300
					refslope = halo.g1000
					#print r3
					#refslope = density_ein.logslope(r3)
					if (logM < 7.5) and (logM > 7.4):
						symbol(log10(halo.Mtid), refslope-density.logslope(r3), color="cyan", symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
						print "issue", halo.halo, halo.treepos
					elif int(halo.treepos) == 36671:
						symbol(log10(halo.Mtid), refslope-density.logslope(r3), color="purple", symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
					else:
						symbol(log10(halo.Mtid), refslope-density.logslope(r3), color=color, symbolName="circlesolid" if converged else "plus", symbolsize=symbolsize)
					
					if int(halo.treepos) == 288:
						print halo.halo, halo.treepos, density.__class__, density.logslope(r3), refslope
				#if (density_hay.logslope(r3)-density_ein.logslope(r3)) > 1.7:
				#	print halo.halo, halo.treepos, log10(halo.Mtid), density_hay.logslope(r3)-density_ein.logslope(r3)
				i += 1
			#labels("index", "slope at 1kpc")
			
			labels("log Mtid", "delta slope at 1kpc")
			ylim(-3, 3)
		print "r3", r3, log10(r3)
		
		draw()
		

class PlotGrid(object):
	def __init__(self, grid, fillregion=None):
		self.grid = grid
		self.fillregion = fillregion
		
	def run(self, opts, args, scope):
		self.grid.load()
		document(size="12cm,10cm")
		page(fontsize="14pt")
		box()
		self.plot_grid(self.grid)
		draw()
		
	def plot_grid(self, grid):
		#colours = [Color.blue.blend(Color.white, a) for a in [0.2, 0.35, 0.5]][::-1]
		alphas = [0.45, 0.3, 0.15]
		#indexedimage(grid.grid.T, resize=grid.resize, colormap="whiterainbow")
		if self.fillregion:
			vfill(*self.fillregion, color="lightgrey")
		for i, level in list(enumerate(grid.levels))[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			ymin = grid.contours[:,i,0]
			ymax = grid.contours[:,i,1]
			#print len(grid.profile.u_centers), len(ymin)
			print len(grid.profile.u_centers), len(ymin), len(ymax)
			o2 = fillrange(grid.profile.u_centers, ymin, ymax, color=color)
			#graph(grid.profile.u_centers, ymin, color="black", linewidth="0.2pt")
			#graph(grid.profile.u_centers, ymax, color="black", linewidth="0.2pt")
		o1 = graph(grid.profile.u_centers, grid.median, color="green")
		labels(grid.profile.label, grid.label)
		(x1, y1), (x2, y2) = grid.resize
		xlim(x1, x2)
		ylim(y1, y2)
		return o1, o2

class EnclosedMass(ParameterSet):
	def __init__(self, parameterset, enclosed_mass_grid):
		ParameterSet.__init__(self, parameterset)
		self.enclosed_mass_grid = enclosed_mass_grid
		
	def plot_grid(self, grid):
		#colours = [Color.blue.blend(Color.white, a) for a in [0.2, 0.35, 0.5]][::-1]
		alphas = [0.45, 0.3, 0.15]
		#indexedimage(grid.grid.T, resize=grid.resize, colormap="whiterainbow")
		for i, level in list(enumerate(grid.levels))[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			ymin = grid.contours[:,i,0]
			ymax = grid.contours[:,i,1]
			#print len(grid.profile.u_centers), len(ymin)
			print len(grid.profile.u_centers), len(ymin), len(ymax)
			fillrange(grid.profile.u_centers, ymin, ymax, color=color)
		graph(grid.profile.u_centers, grid.median, color="green")
		labels(grid.profile.label, grid.label)
		(x1, y1), (x2, y2) = grid.resize
		xlim(x1, x2)
		ylim(y1, y2)
		
	def run(self, argv, opts, scope):
		self.enclosed_mass_grid.load()
		self.load()
		document(size="8cm,8cm")
		mozaic(1,1,box)
		select(0, 0)
		self.plot_grid(self.enclosed_mass_grid)
		if 0:
			symbolsize = "20pt"
			symbol(log10(1.8), log10(2.4*1e8), symbolsize=symbolsize) # g
			#symbol(log10(1.8), log10((2.4+1.1)*1e8))
			#symbol(log10(1.8), log10((2.4-0.7)*1e8))
			symbol(log10(0.375), log10(2.25e7), symbolsize=symbolsize) # wolf
			#symbol(log10(0.375), log10((2.25+0.16)*1e7))
			#symbol(log10(0.375), log10((2.25-0.16)*1e7))
			symbol(log10(0.3), log10(1.2*1e7), symbolsize=symbolsize) # strigari
			#symbol(log10(0.3), log10((1.2+0.11)*1e7))
			#symbol(log10(0.3), log10((1.2-0.37)*1e7))
		xlim(-1, 1)
		ylim(6, 9.)
		#box()
		#autolegend("prior", "using prior", "mass-conc. relation", location="right,bottom", edgespacing="5mm")
		draw()


class MassConcentration(ParameterSet):
	def __init__(self, parameterset, mass_concentration_grid):
		ParameterSet.__init__(self, parameterset)
		self.mass_concentration_grid = mass_concentration_grid
		
	def plot_grid(self, grid):
		#colours = [Color.blue.blend(Color.white, a) for a in [0.2, 0.35, 0.5]][::-1]
		alphas = [0.45, 0.3, 0.15]
		#indexedimage(grid.grid.T, resize=grid.resize, colormap="whiterainbow")
		for i, level in list(enumerate(grid.levels))[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			ymin = grid.contours[:,i,0]
			ymax = grid.contours[:,i,1]
			#print len(grid.profile.u_centers), len(ymin)
			print len(grid.profile.u_centers), len(ymin), len(ymax)
			fillrange(grid.profile.u_centers, ymin, ymax, color=color)
		graph(grid.profile.u_centers, grid.median, color="green")
		labels(grid.profile.label, grid.label)
		(x1, y1), (x2, y2) = grid.resize
		xlim(x1, x2)
		ylim(y1, y2)
		
	def run(self, argv, opts, scope):
		self.mass_concentration_grid.load()
		self.load()
		document(size="15cm,8cm")
		mozaic(2,1,box)
		#box()
		select(0, 0)
		self.plot_parameters(0, 1, scope, color="black")
		#self.plot_parameters(0, 1, scope, color="black", fill=False)
		resize = self.makeresize(0, 1)
		M200_grid = self.mass_concentration_grid.M200_grid
		concentration_grid = self.mass_concentration_grid.c_grid
		levels = [8.5,9.,9.5,10.,10.5,11.]
		clearautolegend()
		for level, linestyle, color in zip(levels, alllinestyles*3, (goodcolors*3)[1:]):
				color = "blue"
				if level == levels[0]:
						linestyle = "dot"
				else:
						linestyle = "normal"
				contour(log10(M200_grid).transpose(), levels=[level], resize=resize, color=color, linestyle=linestyle, addlegend=level==levels[0])
		levels = arange(10, 41, 5)
		#for level, linestyle, color in zip(levels, alllinestyles*3, (goodcolors*3)[1:]):
		for level in levels:
				color = "orange"
				if level == levels[0]:
						linestyle = "dash"
				else:
						linestyle = "normal"
				contour(concentration_grid.transpose(), levels=[level], resize=resize, color=color, linestyle=linestyle, addlegend=level==levels[0])
		autolegend("log M<sub>200</sub> = 8.5", "c=10", location="right,bottom", edgespacing="5mm")
		fill=False
		select(1, 0)
		p = probimage2d(self.mass_concentration_grid.pdf, 0, 1, resize=resize, colormap=None, color="red", drawcontourlines=True, premultiply=True, fill=fill)
		M200_min, M200_max = M200_grid.min(), M200_grid.max()
		n_M200 = 100
		logM200s = arange(n_M200) / (n_M200+1.0) * (log10(M200_max) - log10(M200_min)) + log10(M200_min)
		logM1kpcs = []
		rss = []
		#self.plot_parameters(0, 1, scope, color="blue")
		self.plot_parameters(0, 1, scope, color="black", addlabels=False, addlegend=False)
		#self.plot_parameters(0, 1, scope, prior_grid=self.mass_concentration_grid.pdf, color="lightgreen", fill=True, addlabels=False)
		self.plot_parameters(0, 1, scope, fill=False, prior_grid=self.mass_concentration_grid.pdf, color="green", linewidth="1.5pt")
		import mab
		for i in range(n_M200):
			M200 = 10**logM200s[i]
			c = self.mass_concentration_grid.c(M200)
			pot = mab.gd.potential.NFW(M200, 1.)
			rs = pot.r200/c
			rss.append(rs)
			pot = mab.gd.potential.NFW(M200, rs)
			M1kpc = pot.enclosed_mass(1.0)
			logM1kpcs.append(log10(M1kpc))
		graph(logM1kpcs, log10(rss), linewidth="3pt", linestyle="dot")
		(x1, y1), (x2, y2) = resize
		xlim(x1, x2)
		ylim(y1, y2)
		autolegend("prior", "using prior", "mass-conc. relation", location="right,bottom", edgespacing="5mm")
		draw()

class ParameterSetGd(ParameterSet):
	def __init__(self, parameterset, mass_enclosed_grid, logslope_grid, anisotropy_grid, mass_enclosed_grid2, logslope_grid2, anisotropy_grid2, vdisp_grid, plot_binned_data, plot_moments, kurtosis_grid, mass_measurements=[], bestfitscope=None):
		ParameterSet.__init__(self, parameterset)
		self.mass_enclosed_grid = mass_enclosed_grid
		self.logslope_grid = logslope_grid
		self.anisotropy_grid = anisotropy_grid
		self.mass_enclosed_grid2 = mass_enclosed_grid2
		self.logslope_grid2 = logslope_grid2
		self.anisotropy_grid2 = anisotropy_grid2
		self.vdisp_grid = vdisp_grid
		self.plot_binned_data = plot_binned_data
		self.plot_moments = plot_moments
		self.kurtosis_grid = kurtosis_grid
		self.mass_measurements = mass_measurements
		self.bestfitscope = bestfitscope
		
	def plot_grid(self, grid):
		#colours = [Color.blue.blend(Color.white, a) for a in [0.2, 0.35, 0.5]][::-1]
		alphas = [0.45, 0.3, 0.15]
		#indexedimage(grid.grid.T, resize=grid.resize, colormap="whiterainbow")
		for i, level in list(enumerate(grid.levels))[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			ymin = grid.contours[:,i,0]
			ymax = grid.contours[:,i,1]
			#print len(grid.profile.u_centers), len(ymin)
			print len(grid.profile.u_centers), len(ymin), len(ymax)
			o2 = fillrange(grid.profile.u_centers, ymin, ymax, color=color)
			#graph(grid.profile.u_centers, ymin, color="black", linewidth="0.2pt")
			#graph(grid.profile.u_centers, ymax, color="black", linewidth="0.2pt")
		o1 = graph(grid.profile.u_centers, grid.median, color="green")
		labels(grid.profile.label, grid.label)
		(x1, y1), (x2, y2) = grid.resize
		xlim(x1, x2)
		ylim(y1, y2)
		return o1, o2
		
	def run(self, argv, opts, scope):
		self.mass_enclosed_grid.load()
		self.logslope_grid.load()
		self.anisotropy_grid.load()
		self.mass_enclosed_grid2.load()
		self.logslope_grid2.load()
		self.anisotropy_grid2.load()
		self.vdisp_grid.load()
		self.plot_binned_data.load()
		self.plot_moments.load()
		self.kurtosis_grid.load()
		self.load()
		#self.parameterset.compute_weights()
		dim = self.parameterset.dimension
		self.galaxy = None
		try:
			self.galaxy = scope["galaxy"]
		except:
			pass
		if dim == 2:
			#w = 7# * dim
			#document("%fcm,%fcm" % (w,w))
			#mozaic(dim, dim, box)
			#document("25cm,25cm")
			extra = False
			#import pdb
			try:
				short = scope["short"]
			except:
				short = False
		
			
			if short:
				document("15cm,7cm")
				mozaic(2, 1, box)
			elif extra:
				document("25cm,25cm")
				mozaic(3, 3, box)
			else:
				document("15cm,20cm")
				mozaic(2, 3, box)
			if short:
				select(0, 0)
				self.plot_parameters(0, 1, scope)
			else:
				select(0, 2)
				self.plot_parameters(0, 1, scope)
				select(0, 1)
				self.plot_parameters(0, 0, scope)
				select(0, 0)
				self.plot_parameters(1, 1, scope)
			if short:
				select(1,0)
				vfill(0, 0.1, color="lightgrey")
				self.plot_grid(self.anisotropy_grid)
				if self.galaxy and hasattr(self.galaxy, "beta"): 
					hline(self.galaxy.beta, color="red", linestyle="dash")
				
			else:
				select(1,0)
				obj1, obj2 = self.plot_grid(self.mass_enclosed_grid)
				if self.mass_measurements:
					clearautolegend()
					x = [k[0] for k in self.mass_measurements]
					y = [k[1] for k in self.mass_measurements]
					labels = [k[2] for k in self.mass_measurements]
					sizes = "7pt 3pt 6pt".split()
					symbolnames = "square triangle circle".split()
					colors = "black purple red".split()
					
					i = 0
					for x, y, label in self.mass_measurements:
						#symbol((0.3), log10(1.2*1e7), symbolsize="7pt", symbolName="square", color="black") # strigari
						symbol(x, y, symbolsize=sizes[i], symbolName=symbolnames[i], color=colors[i])
						i += 1
					#symbol((0.375), log10(2.25e7), symbolsize="6pt", symbolName="circle", color="red") # wolf
					#autolegend("Strigari et al. 2008", "Wolf et al. 2010", location="right,bottom", edgespacing="5mm")
					autolegend(*labels, location="right, bottom", edgespacing="5mm")
				
				try:
					dm_density = scope["dm_density"]
					menc = dm_density.enclosed_mass(self.mass_enclosed_grid.profile.centers)
					obj0 = graph(self.mass_enclosed_grid.profile.centers, log10(menc), color="red", linestyle="dash")
				except:
					pass
				ylim(6, 8.5)
				#self.plot_grid(self.kurtosis_grid)
				#self.plot_moments.plot_kurtosis(self.plot_moments.binned_data[1], symbolsize="5pt", alpha=0.5)
				select(1,1)
				self.plot_grid(self.logslope_grid)
				try:
					dm_density = scope["dm_density"]
					slope = dm_density.logslope(self.logslope_grid.profile.centers)
					graph(self.logslope_grid.profile.centers, slope, color="red", linestyle="dash")
					legend(["line", "line", "squaresolid"], ["true value", "median", "confidence intervals"], [obj0, obj1, obj2], location="left, top")
				except:
					legend(["line", "squaresolid"], ["median", "confidence intervals"], [obj1, obj2], location="left, top")
				select(1,2)
				vfill(0, 0.1, color="lightgrey")
				#vline(0.05, linestyle="dash")
				self.plot_grid(self.anisotropy_grid)
				if self.galaxy and hasattr(self.galaxy, "beta"): 
					hline(self.galaxy.beta, color="red", linestyle="dash")
				if self.bestfitscope:
					self.bestfitscope.load()
					solution = self.bestfitscope.subscope["solution"]
					xborders_3d = self.bestfitscope.subscope["schwmodel"].storage_3d.xborders
					solution.load()
					moments3d = solution.moments3d_solution
					varvr = moments3d[1] 
					varvphi = moments3d[2]/2
					#self.varvtheta = moments3d[6]
					anisotropy  = 1 - (varvphi)/(varvr)
					graph(xborders_3d, anisotropy, color="orange")#, linestyle="dot")
					print xborders_3d, anisotropy
					
					print solution, dir(solution)
					#dsa
					
				if extra:
					select(2,0)
					self.plot_grid(self.mass_enclosed_grid2)
					xlim(-2, 1)
					ylim(4, 10)
					select(2,1)
					self.plot_grid(self.vdisp_grid)
					print "vdisp grid", self.vdisp_grid.grid.shape
					self.plot_binned_data.plot_data(symbolsize="5pt", alpha=0.5)
					select(2,2)
					self.plot_grid(self.anisotropy_grid2)
		if dim == 3:
			#w = 7# * dim
			#document("%fcm,%fcm" % (w,w))
			#mozaic(dim, dim, box)
			extra = False
			if extra:
				document("25cm,25cm")
				mozaic(3, 4, box)
			else:
				document("20cm,15cm")
				mozaic(3, 3, box)
			
			select(0, 2)
			self.plot_parameters(0, 1, scope)
			select(1, 2)
			self.plot_parameters(0, 2, scope)
			select(2, 2)
			self.plot_parameters(1, 2, scope)

			select(0, 1)
			self.plot_parameters(0, 0, scope)
			select(1, 1)
			self.plot_parameters(1, 1, scope)
			select(2, 1)
			self.plot_parameters(2, 2, scope)

			select(0,0)
			vfill(0, 0.1, color="lightgrey")
			self.plot_grid(self.anisotropy_grid)
			if self.galaxy and hasattr(self.galaxy, "beta"): 
				hline(self.galaxy.beta, color="red", linestyle="dash")
			#vline(0.0275)
			select(1,0)
			obj1, obj2 = self.plot_grid(self.mass_enclosed_grid)
			if self.mass_measurements:
				clearautolegend()
				x = [k[0] for k in self.mass_measurements]
				y = [k[1] for k in self.mass_measurements]
				labels = [k[2] for k in self.mass_measurements]
				sizes = "7pt 3pt 6pt".split()
				symbolnames = "square triangle circle".split()
				colors = "black purple red".split()
				
				i = 0
				for x, y, label in self.mass_measurements:
					#symbol((0.3), log10(1.2*1e7), symbolsize="7pt", symbolName="square", color="black") # strigari
					symbol(x, y, symbolsize=sizes[i], symbolName=symbolnames[i], color=colors[i])
					i += 1
				#symbol((0.375), log10(2.25e7), symbolsize="6pt", symbolName="circle", color="red") # wolf
				#autolegend("Strigari et al. 2008", "Wolf et al. 2010", location="right,bottom", edgespacing="5mm")
				autolegend(*labels, location="left, top", edgespacing="4mm")
			if 0:
				symbolsize="20pt"
				#symbol(log10(1.8), log10(2.4*1e8), symbolsize=symbolsize) # g
				#symbol(log10(1.8), log10((2.4+1.1)*1e8))
				#symbol(log10(1.8), log10((2.4-0.7)*1e8))
				clearautolegend()
				symbol((0.3), log10(1.2*1e7), symbolsize="7pt", symbolName="square", color="black") # strigari
				symbol((0.375), log10(2.25e7), symbolsize="6pt", symbolName="circle", color="red") # wolf
				#symbol(log10(0.375), log10((2.25+0.16)*1e7))
				#symbol(log10(0.375), log10((2.25-0.16)*1e7))
				autolegend("Strigari et al. 2008", "Wolf et al. 2010", location="right,bottom", edgespacing="5mm")
			try:
				dm_density = scope["dm_density"]
				menc = dm_density.enclosed_mass(self.mass_enclosed_grid.profile.centers)
				obj0 = graph(self.mass_enclosed_grid.profile.centers, log10(menc), color="red", linestyle="dash")
			except:
				pass
			ylim(6, 8.5)
			select(2,0)
			self.plot_grid(self.logslope_grid)
#<<<<<<< Updated upstream
			try:
				dm_density = scope["dm_density"]
				slope = dm_density.logslope(self.logslope_grid.profile.centers)
				graph(self.logslope_grid.profile.centers, slope, color="red", linestyle="dash")
				legend(["line", "line", "squaresolid"], ["true value", "median", "confidence intervals"], [obj0, obj1, obj2], location="left, top", edgespacing="5mm")
			except:
				legend(["line", "squaresolid"], ["median", "confidence intervals"], [obj1, obj2], location="left, top", edgespacing="5mm")
#=======
			#dm_density = scope["dm_density"]
			#slope = dm_density.logslope(self.logslope_grid.profile.centers)
			#graph(self.logslope_grid.profile.centers, slope, color="red")
#>>>>>>> Stashed changes
			#try:
			#exceot
			#select(1,2)
			if extra:
				select(0,3)
				self.plot_grid(self.vdisp_grid)
				self.plot_binned_data.plot_data(symbolsize="5pt", alpha=0.5)
				select(1,3)
				#self.plot_grid(self.anisotropy_grid2)
				self.plot_grid(self.kurtosis_grid)
				self.plot_moments.plot_kurtosis(self.plot_moments.binned_data[1], symbolsize="5pt", alpha=0.5)
				select(2,3)
				self.plot_grid(self.mass_enclosed_grid2)
				
				#symbol(log10(1.8), log10(2.4*1e8)) # g
				#symbol(log10(1.8), log10((2.4+1.1)*1e8))
				#symbol(log10(1.8), log10((2.4-0.7)*1e8))
				#symbol(log10(0.375), log10(2.25e7)) # wolf
				#symbol(log10(0.375), log10((2.25+0.16)*1e7))
				#symbol(log10(0.375), log10((2.25-0.16)*1e7))
				#symbol(log10(0.3), log10(1.2*1e7)) # strigari
				#symbol(log10(0.3), log10((1.2+0.11)*1e7))
				#symbol(log10(0.3), log10((1.2-0.37)*1e7))
				xlim(-1, 0.5)
				ylim(6, 9)
		#select(1,0)
		if 0:
			i, j = 0, 1
			x = [pv.values_org[i] for pv in self.parameterset.parameter_values]
			y = [pv.values_org[j] for pv in self.parameterset.parameter_values]
			#colors = [pv.weight for pv in self.parameterset.parameter_values]
			#scatter(x, y, colormap="whiteblack")
			for pv in self.parameterset.parameter_values:
				print "r", pv.left[i], pv.left[j], pv.right[i], pv.right[j]
				rectangle(pv.left[i], pv.left[j], pv.right[i], pv.right[j])
			scatter(x, y, color="red", symbolsize="20pt")
			grow(1.1)
		draw()

class PlotSchwLosvds(object):
	def __init__(self, storage_2d_losvd, dfgrid):
		self.storage_2d_losvd = storage_2d_losvd
		self.dfgrid = dfgrid
		
		
	def run(self, argv, opts, scope):
		self.storage_2d_losvd.load()
		losvds = self.storage_2d_losvd.losvds
		n_I1, n_I2 = self.dfgrid.n_I1, self.dfgrid.n_I2
		newshape = (n_I2, n_I1) + losvds.shape[1:]
		print newshape
		losvds = losvds.reshape(newshape)
		print losvds.shape
		scale = 2.2
		document("%scm,%scm" % (n_I1*scale, n_I2*scale))
		mozaic(n_I1, n_I2)
		for i1 in range(n_I1):
			for i2 in range(n_I2):
				select(i1,i2)
				border()
				spacer("1mm")
				losvd = losvds[i2,i1]
				print losvd.min(), losvd.max(), losvd.std(), losvd.shape
				indexedimage(losvd, colormap="heat2")
				if 0:
					x = arange(losvd.shape[0])
					Ny, Nx = losvd.shape
					for i in range(Nx)[::2]:
						y = losvd[:,i]
						if y.max() > 0:
							graph(x, y/y.max() * Ny, linewidth="0.5pt")
		#select(0, 0)
		#labels("r", "vlos")
		draw()
		

class PlotSchwDensity2d(object):
	def __init__(self, storage_2d_losvd, dfgrid):
		self.storage_2d_losvd = storage_2d_losvd
		self.dfgrid = dfgrid
		
		
	def run(self, args, opts, scope):
		print args
		if len(args) > 1:
			Lz_index = int(args[1])
			return_value = 2
		else:
			Lz_index = None
			return_value = 1
		self.storage_2d_losvd.load()
		masses = self.storage_2d_losvd.masses
		print self.storage_2d_losvd
		print self.storage_2d_losvd.masses.shape
		print self.storage_2d_losvd.losvds.shape
		n_I1, n_I2 = self.dfgrid.n_I1, self.dfgrid.n_I2
		print n_I1, n_I2, self.storage_2d_losvd.NLz
		newshape = (n_I2, n_I1) + masses.shape[1:]
		print newshape
		#dsa
		masses = masses.reshape(newshape)
		print masses.shape
		scale = 2.2
		#document("%scm,%scm" % (n_I1*scale, n_I2*scale))
		document("%scm,%scm" % ((self.storage_2d_losvd.NLz+1)*scale, n_I2*scale))
		mozaic(self.storage_2d_losvd.NLz+1, n_I2)
		#for i1 in range(n_I1):
		i1 = 12
		for Lz_index in range(self.storage_2d_losvd.NLz):
			for i2 in range(n_I2):
				#select(i1,i2)
				select(Lz_index,i2)
				border()
				spacer("1mm")
				if Lz_index is None:
					mass = sum(masses[i2,i1], axis=0)
				else:
					mass = masses[i2,i1, Lz_index]
				#print mass.shape
				mass = mass.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				#print losvd.min(), losvd.max(), losvd.std(), losvd.shape
				if 0:
					mass = log10(mass)
					mass -= mass.max()
					mass[mass<-2] = -2
				indexedimage(mass.T, colormap="heat2")
				if 0:
					x = arange(losvd.shape[0])
					Ny, Nx = losvd.shape
					for i in range(Nx)[::2]:
						y = losvd[:,i]
						if y.max() > 0:
							graph(x, y/y.max() * Ny, linewidth="0.5pt")
		if 1:
			for i2 in range(n_I2):
				#select(i1,i2)
				select(self.storage_2d_losvd.NLz,i2)
				border()
				spacer("1mm")
				mass = sum(masses[i2,i1], axis=0)
				mass = mass.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				indexedimage(mass.T, colormap="heat2")
		#select(0, 0)
		#labels("r", "vlos")
		draw()		
		return return_value

class PlotSchwVelocity2d(object):
	def __init__(self, storage_2d_losvd, dfgrid):
		self.storage_2d_losvd = storage_2d_losvd
		self.dfgrid = dfgrid
		
		
	def run(self, args, opts, scope):
		print args
		if len(args) > 1:
			Lz_index = int(args[1])
			return_value = 2
		else:
			Lz_index = None
			return_value = 1
		self.storage_2d_losvd.load()
		losvds = self.storage_2d_losvd.losvds
		masses = self.storage_2d_losvd.masses
		print losvds.shape
		n_I1, n_I2 = self.dfgrid.n_I1, self.dfgrid.n_I2
		print n_I1, n_I2, self.storage_2d_losvd.NLz
		newshape = (n_I2, n_I1) + losvds.shape[1:]
		print newshape
		#dsa
		losvds = losvds.reshape(newshape)
		newshape = (n_I2, n_I1) + masses.shape[1:]
		masses = masses.reshape(newshape)
		print losvds.shape
		scale = 2.2
		#document("%scm,%scm" % (n_I1*scale, n_I2*scale))
		document("%scm,%scm" % ((self.storage_2d_losvd.NLz)*scale*3, n_I2*scale))
		mozaic(self.storage_2d_losvd.NLz*3, n_I2)
		#for i1 in range(n_I1):
		i1 = 12
		for Lz_index in range(self.storage_2d_losvd.NLz):
			for i2 in range(n_I2):
				#select(i1,i2)
				select(Lz_index*3,i2)
				border()
				spacer("1mm")
				losvd = losvds[i2,i1, Lz_index]
				mass = masses[i2,i1, Lz_index]
				vs = self.storage_2d_losvd.vs
				#NA = self.storage_2d_losvds.
				meanv = mean((losvd.T*vs).T, axis=0)
				varv = mean((losvd.T*vs**2).T, axis=0)
				mask = mass>0
				meanv[mask] /= mass[mask]
				varv[mask] /= mass[mask]
				sigmav = sqrt(varv-meanv**2)
				#import pdb
				#pdb.set_trace()
				
				#means = 1
				#print mass.shape
				mass = mass.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				meanv = meanv.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				sigmav = sigmav.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				losvd = losvd.reshape((self.storage_2d_losvd.Nv, self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				#print losvd.min(), losvd.max(), losvd.std(), losvd.shape
				if 0:
					mass = log10(mass)
					mass -= mass.max()
					mass[mass<-2] = -2
				select(Lz_index*3,i2)
				border()
				spacer("1mm")
				indexedimage(mass.T, colormap="heat2")
				select(Lz_index*3+1,i2)
				border()
				spacer("1mm")
				indexedimage(meanv.T, colormap="heat2")
				select(Lz_index*3+2,i2)
				print meanv.max()
				print sigmav.max()
				border()
				spacer("1mm")
				indexedimage(sigmav.T, colormap="heat2")
				if 0:
					x = arange(losvd.shape[0])
					Ny = losvd.shape[0]
					#Ny, Nx = losvd.shape
					#for i in range(Nx)[::2]:
					#for i, j in [(
					for i in [10, 15,20]:
						for j in [10, 15,20]:
							y = losvd[:,i,j]
							if y.max() > 0:
								graph(x, y/y.max() * Ny, linewidth="0.5pt", color="green")
		if 0:
			for i2 in range(n_I2):
				#select(i1,i2)
				select(self.storage_2d_losvd.NLz,i2)
				border()
				spacer("1mm")
				mass = sum(masses[i2,i1], axis=0)
				mass = mass.reshape((self.storage_2d_losvd.NR, self.storage_2d_losvd.NR))
				indexedimage(mass, colormap="heat2")
		#select(0, 0)
		#labels("r", "vlos")
		draw()		
		return return_value

class PlotSchwOneLosvd(object):
	def __init__(self, storage_2d_losvd, dfgrid):
		self.storage_2d_losvd = storage_2d_losvd
		self.dfgrid = dfgrid
		
		
	def run(self, argv, opts, scope):
		i1 = int(argv[1])
		i2 = int(argv[2])
		self.storage_2d_losvd.load()
		losvds = self.storage_2d_losvd.losvds
		n_I1, n_I2 = self.dfgrid.n_I1, self.dfgrid.n_I2
		newshape = (n_I2, n_I1) + losvds.shape[1:]
		print newshape
		losvds = losvds.reshape(newshape)
		print losvds.shape
		#document("%scm,%scm" % (n_I1*2, n_I2*2))
		mozaic(2, 2, box)
		#ox()
		#border()
		#spacer("1mm")
		print i2, i1
		print losvds.shape
		losvd = losvds[i2,i1]
		#print losvd.min(), losvd.max(), losvd.std(), losvd.shape
		indexedimage(losvd, colormap="whiteblack")
		E = self.dfgrid.Es[i1]
		print E
		l = self.dfgrid.ls[i2]
		profile_model = self.dfgrid.profile_model
		L = l * profile_model.Lmax_at_E(E)
		print "l",l, L
		rmin, rmax, = 0, 10. #self.storage_2d_losvd.Rmax
		nr = 10000
		r = (arange(nr) + 0.5)/ (nr + 1.) * (rmax - rmin) + rmin
		Ekinr = 2*(E - profile_model.potentialr(r)) - L**2/r**2
		mask = Ekinr>=0 
		mass = 1 / sqrt(Ekinr) * r**1
		mass[Ekinr<0] = 0
		vr = sqrt(Ekinr)
		vr[Ekinr < 0] = 0
		select(1, 0)
		graph(r, vr * mass)
		vt = L/r
		vt[Ekinr < 0] = 0
		print vr
		print vt
		print r[Ekinr>0]
		#dsa
		
		r = r[mask]
		vr = vr[mask]
		weight = mass = mass[mask]
		vt = vt[mask]
		
		M = len(r)
		costheta = numpy.random.random(M) * 2 - 1
		phi = numpy.random.random(M) * 2 * pi
		eta = numpy.random.random(M) * 2 * pi
		theta = numpy.arccos(costheta)
		#sintheta = numpy.sqrt(1-costheta**2)
		sintheta = numpy.sin(theta)
		#print r.shape, sintheta.shape, phi.shape, len(dt) 
		x = r * sintheta * numpy.cos(phi)
		y = r * sintheta * numpy.sin(phi)
		R = sqrt(x**2+y**2)
		vlos = vr * costheta - vt * numpy.cos(eta) * sintheta
		
		select(0, 1)
		vmax = self.storage_2d_losvd.vmax
		Rmax = self.storage_2d_losvd.Rmax
		d, _, _ = numpy.histogram2d(R, vlos, bins=[30,30], range=[(0, Rmax), (-vmax, vmax)], normed=True)
		indexedimage(d.T, colormap="whiteblack")
		
		#histogram(vlos, R)
		select(1,1)
		d, _ = numpy.histogram(R, bins=100, range=(0, Rmax), weights=weight, normed=True, new=True)
		NR = 100
		R = (arange(NR) + 0.5) / NR * Rmax
		graph(R, d)
		
		#histogram(R, function=lambda x: x*weight, datamin=0, datamax=Rmax, bincount=100, normalize=True)
		mass = sum(losvd, axis=0)
		NR = self.storage_2d_losvd.NR
		R = (arange(NR) + 0.5) / NR * Rmax
		assert len(mass) == NR
		#print R, mass 
		graph(R, mass/mass.max() / 2, color="red")
		
		
		print r
		
		print l, L
		draw()		
				
class PlotSchwLosvdDFGridFit(object):
	def __init__(self, storage_2d_losvd, dfgrid, observation, vmean, filters=[]):
		self.storage_2d_losvd = storage_2d_losvd
		self.dfgrid = dfgrid
		self.observation = observation
		self.vmean = vmean
		self.filters = filters
		
		
	def run(self, argv, opts, scope):
		self.storage_2d_losvd.load()
		losvds = self.storage_2d_losvd.losvds
		n_I1, n_I2 = self.dfgrid.n_I1, self.dfgrid.n_I2
		newshape = (n_I2, n_I1) + losvds.shape[1:]
		print newshape
		losvds = losvds.reshape(newshape)
		print losvds.shape
		stars = self.observation.load()
		for filter in self.filters:
			stars = filter(stars)
		stars, aperture_index = self.storage_2d_losvd.stars_to_apertureindices(stars)
		velocity_indices = self.storage_2d_losvd.velocity_to_index(stars.vlos-self.vmean)
		document("35cm,15cm")
		sigma_v = 2.
		for i1 in range(n_I1):
			for i2 in range(n_I2):
				pass
				#losvds[i2,i1] /= (delta_v * delta_R)
				if 0:
					for ai in range(losvds.shape[-1]):
						if sum(sum(losvds[i2,i1,:,ai])) > 0:
							losvds[i2,i1,:,ai] /= sum(losvds[i2,i1,:,ai])
				else: 
					losvds[i2,i1] /= sum(losvds[i2,i1]) 
		
		if 0:
			document("%scm,%scm" % (n_I1*2, n_I2*2))
			mozaic(n_I1, n_I2)
			for i1 in range(n_I1):
				for i2 in range(n_I2):
					select(i1,i2)
					border()
					spacer("1mm")
					losvd = losvds[i2,i1]
					print losvd.min(), losvd.max(), losvd.std(), losvd.shape
					indexedimage(losvd, colormap="whiteblack")
					scatter(aperture_index, velocity_indices, color="red")
		else:
			mozaic(2,1,box)
			pgrid = zeros_like(losvds[:,:,0,0])
			print aperture_index
			print velocity_indices
			
			best_index = -1
			best_p = -1e8
			for i1 in range(n_I1):
				for i2 in range(n_I2):
					ps = losvds[i2,i1,velocity_indices,aperture_index]
					#print ps
					p = sum(log(ps))
					#print p,
					print sum(losvds[i2,i1])
					if p > best_p:
						best_index = (i2,i1)
						best_p = p
						print (i2, i1)
					pgrid[i2,i1] = p
			pgrid -= pgrid.max()
			print pgrid
			pgrid[pgrid<-100] = -100
			#indexedimage(exp(pgrid), colormap="whiteblack")
			probimage2d(exp(pgrid.T), 0, 1, range(n_I1), range(n_I2), colormap=None, drawcontourlines=True) 
			#"whiteblack")
			select(1,0)
			if 0:
				v = stars.vlos-self.vmean
				scatter(stars.rc/3600, v, symbolsize="20pt")
				ylim(-abs(v).max(), abs(v).max())
				print abs(v).max(), self.storage_2d_losvd.vmax
			else:
				#i1 = int(argv[1])
				#i2 = int(argv[2])
				i2, i1 = best_index
				losvd = losvds[i2,i1]
				print best_index
				#ps = losvds[i2,i1,velocity_indices,aperture_index]
				ps = losvd[velocity_indices,aperture_index]
				p = sum(log(ps))
				print
				
				print p

				indexedimage(losvd, colormap="whiteblack")
				scatter(aperture_index, velocity_indices, color="red", symbolsize="20pt")
		draw()


class PlotDF(object):
	def __init__(self, df):
		self.df = df
		
		
	def run(self, args, opts, scope):
		self.df.load()
		
		x = self.df.es
		y = self.df.fE
		x0 = scope["profile_model_2c"].potentialr(1e-5)
		
		page(fontsize="14pt")
		print x
		print y
		mask = y > 0
		#graph(log(x)[mask], log10(y)[mask])
		graph((-x/x0)[mask], log10(y)[mask] + 5)
		#labels("log10(-E) (binding energy)", "log10 f_E(E)")
		labels("-E/E<sub>0</sub>", "log f<sub>E</sub>(E)")
		xlim(0., 1.00001)
		draw()

class PlotDFSamples(object):
	def __init__(self, df_sampling, galaxy):
		self.df_sampling = df_sampling
		self.galaxy = galaxy
		
	def run(self, args, opts, scope):
		self.df_sampling.load()
		self.plot(scope, True)
		
	def plot(self, scope, dm=False):
		if not dm:
			samples = self.df_sampling.df_samples
		else:
			samples = self.df_sampling.df_samples_dm
			
		Nperbin = 1000/5
		N = self.df_sampling.N
		if dm:
			N /= 1000000
				
		galaxy = self.galaxy	
		mozaic(3,2,box)
		select(0,1)
		rmax = samples.r3d.max()
		rmax = 1000
		print samples.r3d
		c1 = cumhistogram(log10(samples.r3d), binwidth=0.01, color="red", datamin=-4, datamax=log10(rmax), normalize=True, differential=False)
		numpy.random.seed(1)
		if dm:
			r2 = galaxy.profile_model.dm_profile.sample_r(N=1000, rmax=rmax)
			#print r2
			#dsa
			print galaxy.profile_model.dm_profile
		else:
			r2 = galaxy.light_model.sample_r(N=N*10)
			print galaxy.light_model
		
		c2 = cumhistogram(log10(r2), binwidth=0.01, datamin=-4, datamax=log10(rmax), normalize=True, differential=False)
		print sum(c2.data)
		logr = c2.bins
		r = 10**logr
		dlogr = logr[1] - logr[0]
		if dm:
			scale = galaxy.profile_model.dm_profile.enclosed_mass(rmax)
			#graph(logr, galaxy.profile_model.dm_profile.densityr(r)*log(10)*r**3*4*pi/scale, color="green")
			graph(logr, galaxy.profile_model.dm_profile.enclosed_mass(r)/scale, color="green", linestyle="dot")
		else:
			print galaxy.profile_model.light_profile
			y = galaxy.profile_model.light_profile.densityr(r, M=1.)*log(10)*r**3*4*pi
			graph(logr, cumsum(y*dlogr), color="green")
			print sum(y*dlogr)
			
		labels("log r", "dM/M_tot")
		select(0, 0)
		diff = c2.data-c1.data
		print len(c2.bins), len(c2.data)
		#print diff
		hline(color="lightgrey")
		histogramline(c2.bins, diff, drawsides=False, fill=False, drawverticals=False)
		labels("log r", "deltaM")
		ylim(-0.01/1, 0.01/1)
		if 1:
			xs = []
			betas = []
			varvrs = []
			varvphis = []
			varvthetas = []
			for n, r, vr, vphi, vtheta in binrange(Nperbin, samples.r3d, samples.vr, samples.vphi, samples.vtheta):
				#print n, mean(vr), std(vr), std(vphi)
				varvphi = var(vphi)
				varvtheta = var(vtheta)
				varvr = var(vr)
				beta = 1 - (varvphi+varvtheta)/(2*varvr)
				#print mean(r), beta
				xs.append(log10(mean(r)))
				betas.append(beta)
				varvrs.append(varvr)
				varvphis.append(varvphi)
				varvthetas.append(varvtheta)
				
			select(1,1)
			graph(xs, sqrt(varvrs), color="red")
			graph(xs, sqrt(varvphis), color="green")
			graph(xs, sqrt(varvthetas), color="blue")
			labels("log r", "sigma_(r,t)")
			logrs = arange(min(xs), max(xs), 0.01)
			#if isinstance(galaxy, Galaxy_constant_anisotropy)
			if hasattr(galaxy, "beta") and not dm:
				varvrs = galaxy.jeans().sigmar(10**logrs)**2
				varvts = (1-galaxy.beta) * varvrs
				#print varvrs
				graph(logrs, sqrt(varvrs), color="black", linestyle="dash")
				graph(logrs, sqrt(varvts), color="black", linestyle="dash")
				
			ylim(0, 120)
			
			
			select(1,0)
			#print len(xs), len(betas)
			#graph(xs, betas, color="black")
			xs = array(xs)
			graph(10**xs, betas, color="black")
			xlim(0, 1.5)
			if hasattr(galaxy, "beta0") and not dm:
				hline(galaxy.beta0)
				hline(galaxy.betainf)
			ylim(-2,1)
			labels("log r", "anisotropy")
		
		if 1:
			Rs = []
			betas = []
			sigma_vlos = []
			Nperbin = 1000
			R = sqrt(samples.x**2 + samples.y**2)
			for n, R, vz in binrange(Nperbin, R, samples.vz):
				#print n, mean(vr), std(vr), std(vphi)
				varvz = var(vz)
				Rs.append(mean(R)) #log10(mean(R)))
				sigma_vlos.append(varvz**0.5)
			select(2,1)
			graph(Rs, sigma_vlos, color="red")
			ylim(0, 20)
			xlim(0, 1.5)
			labels("log r", "sigma_vlos")
		if 1:
			select(2,0)
			if hasattr(galaxy, "beta") and not dm:
				rs = arange(1e-3, 10.1, 0.1/100)
				varvrs = galaxy.jeans().sigmar(rs)**2
				varvts = (1-galaxy.beta) * varvrs
				#print varvrs
				graph(rs, sqrt(varvrs), color="black", linestyle="dash")
				graph(rs, sqrt(varvts), color="black", linestyle="dash")
			labels("r", "sigma")
			
		draw()
		
		
class BinVlos(object):
	def __init__(self, observation):
		self.observation = observation
		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		if 1:
			Rs = []
			betas = []
			sigma_vlos = []
			Nperbin = 40
			R = sqrt(stars.x**2 + stars.y**2)
			for n, R, vz in binrange(Nperbin, R, stars.vz):
				#print n, mean(vr), std(vr), std(vphi)
				varvz = var(vz)
				Rs.append(mean(R)) #log10(mean(R)))
				sigma_vlos.append(varvz**0.5)
			#select(2,1)
			graph(Rs, sigma_vlos, color="red")
			ylim(0, 20)
			xlim(0, 1.5)
		
		#box()
		draw()
		
class MassProfiles(object):
	def __init__(self, profiles, param_r, N):
		self.profiles = profiles
		self.param_r = param_r
		self.N = N
		
	def run(self, args, opts, scope):
		mozaic(2,1,box)
		
		x = mgrid[self.param_r.min:self.param_r.max:1j*self.N]
		r = self.param_r.f(x)
		for profile, color in zip(self.profiles, nicecolors):
			vcirc = profile.vcirc(r)
			graph(x, vcirc, color=color)
		ylim(0, 300)
		select(1,0)
		for profile, color in zip(self.profiles, nicecolors):
			menc = profile.enclosed_mass(r)
			graph(x, log10(menc), color=color)
			symbol(self.param_r.finv(profile.r200), log10(profile.M200), color=color, symbolsize="20pt")
		labels(self.param_r.label, "log M")
		draw()
			
		
		
		
		
		
		
		
	
	
	
class Photometry2d(object):
	def __init__(self, photometry):
		self.photometry = photometry
		
	def run(self, args, opts, scope):
		self.photometry.load()
		document(size="20cm,20cm")
		box()
		resize = (-1, -1), (1,1)
		image = log10(self.photometry.grid2d.T)
		image -= image.max()
		logmin = -3
		image[image<logmin] = logmin
		indexedimage(image, resize=resize)
		contour(image, levels=10, resize=resize)
		labels("x", "y")
		draw()
		
		

class Counts(object):
	def __init__(self, observation, aperture, member_filter, light_profile, model, foreground_model, galaxy_model, member_ratios):
		self.observation = observation
		self.aperture = aperture
		self.member_filter = member_filter
		self.light_profile = light_profile
		self.model = model
		self.foreground_model = foreground_model
		self.galaxy_model = galaxy_model
		self.member_ratios = member_ratios
		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		self.member_ratios.load()
		self.aperture.load()
		stars = stars.filter(lambda star: star.re < self.aperture.aperture_rborders[-1])
		bincount = len(self.aperture)
		nx = int(math.ceil(sqrt(bincount)))+1
		ny = int(math.ceil(float(bincount)/nx))+1
		if nx * ny >= bincount:
			ny += 1
		print bincount, nx, ny
		nx = 3
		ny = 3
		mozaic(nx, ny, box)
		allstars = stars
		members = []
		non_members = []
		counts = []
		members_predict = []
		non_members_predict = []
		N_all_member = len(self.member_filter.filter(allstars))
		N_all = len(allstars)
		N_all_non_member = N_all - N_all_member
		test = []
		
		mean_sigma = 2.
		v1 = self.member_filter.v - 3*self.member_filter.sigma
		v2 = self.member_filter.v + 3*self.member_filter.sigma
		
		if self.model:
			self.model.load()
			#self.model.model_wrapper.model.fraction_foreground = 0.35
		
		f_m_in = scipy.integrate.quad(lambda v: self.galaxy_model(v, mean_sigma), v1, v2)[0]
		f_n_in = scipy.integrate.quad(lambda v: self.foreground_model(v, mean_sigma), v1, v2)[0]
		f_m_out = scipy.integrate.quad(lambda v: self.galaxy_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.galaxy_model(v, mean_sigma), v2, inf)[0]
		f_n_out = scipy.integrate.quad(lambda v: self.foreground_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.foreground_model(v, mean_sigma), v2, inf)[0]
		
		print "=" * 70
		print f_m_in, f_n_in
		print f_m_out, f_n_out 
		print "=" * 70
		for p in self.galaxy_model.parameters:
			print p.name, p.get()
		print v1, v2
		#dsa
		M = numpy.matrix([[f_m_in, f_n_in], [f_m_out, f_n_out]]) #.T
		scaling = (1-0.33)
		scaling = 1.
		print M
		#M = numpy.matrix([[f_m_in, f_n_in], [f_m_out, f_n_out]])#$.T
		#self.light_profile.b *= 1/(1-0.33)
		#self.aperture.aperture_rborders_kpc[-1] -= 0.2
		def rho_m(R):
			return self.light_profile.densityR(R) / scipy.integrate.quad(lambda R: 2*pi*self.light_profile.densityR(R)*R, 0, self.aperture.aperture_rborders_kpc[-1])[0]
		def rho_n(R):
			return 1./scipy.integrate.quad(lambda R: 2 * pi * R, 0, self.aperture.aperture_rborders_kpc[-1]*scaling)[0]
		for i in range(bincount):
			#select(i%nx, i/nx)
			#filter = self.aperture.aperture_filter(i)
			#stars = allstars.filter(filter)
			Rmin = self.aperture.aperture_rborders[i]
			Rmax = self.aperture.aperture_rborders[i+1]
			stars = allstars.filter(lambda star: (star.re >= Rmin) & (star.re < Rmax))
			#histogram(stars.vlos_helio, datamin=-200, datamax=300, binwidth=2.)
			N = len(stars)
			N_member = len(self.member_filter.filter(stars))
			counts.append(N)
			Rmin = self.aperture.aperture_rborders_kpc[i]
			Rmax = self.aperture.aperture_rborders_kpc[i+1]
			
			#print self.member_filter
			#import pdb
			#pdb.set_trace()
			
			C_in = N_member
			C_out = N - C_in
			if 0:
				members.append(N_member)
				non_members.append(N-N_member)
			else:
				print "before:", N_member, N-N_member
				print M.I
				y = dot(M.I, array([C_in, C_out])).T
				print y
				y = array(y.flat)
				Nm = y[0]
				Nn = y[1]
				print "after:", Nm, Nn
				members.append(Nm)
				non_members.append(Nn)
				if Nm < 0:
					import pdb
					pdb.set_trace()
				
			
			#C_in = len(self.member_filter.filter(stars))
			#M = [[
			
			
			#print scipy.integrate.quad(lambda R: rho_n(R) * 2 * pi * R, 0, self.aperture.aperture_rborders_kpc[-1])[0]
			
			R1 = Rmin
			R2 = Rmax
			#R1 = 0
			#R2 = self.aperture.aperture_rborders_kpc[-1]
			#print dir(self.model.model_wrapper.model)
			#print self.model.model_wrapper.model.fraction_foreground
			#dsa
			#f = 0.264414637969
			#f = self.model.fitter.get_model().fraction_foreground
			#f = self.model.fitter.get_model().fraction_foreground
			model = self.model.fitter.get_model()
			#f = 1.
			ratio = abs(model.ratio)
			w1 =  ratio / (1. + ratio);
			w2 = 1 / (1. + ratio);
			w1 = w1*1.7# * 0.7
			f = w1
			#w1 = 0.84
			print "ratio,w1,w2", ratio, w1, w2, w1/w2
			#f = 0.4
			#f = 0.3
			#w1 = 0.14
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all) * 2 * pi * R, R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*(1-f)) * 2 * pi * R, R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*f) * 2 * pi * R, R1, R2)[0]
			#f_total = f_member + f_non_member 
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all)/ (rho_m(R)*N_all_member/N_all + rho_n(R)*N_all_non_member/N_all) , R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all)/ (rho_m(R)*N_all_member/N_all + rho_n(R)*N_all_non_member/N_all) , R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
			#print "ftotal", f_total
			
			if 1:
				mcounts = 0
				nmcounts = 0
				for star in stars:
					R = self.aperture.light_model.arcsec_to_kpc(star.re)
					f_member = rho_m(R)*(1-w1)
					f_non_member = rho_n(R)*(w1)
					f = f_member + f_non_member
					f_member = f_member / f
					f_non_member = f_non_member/ f
					mcounts += f_member
					nmcounts += f_non_member
				f_member = mcounts / N
				f_non_member = nmcounts / N
				nm = f_member * N
				nnm = f_non_member * N
				print ">>>>", f_member, f_non_member, nm, nnm, nm+nnm, N, w1
			elif 1:
				expected = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member + rho_m(R)*N_all_member) * 2 * pi * R, R1, R2)[0]
				bias = N/expected
				#dsa
				#* N_all_member / N_all
				nm = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member) * 2 * pi * R, R1, R2)[0] * bias
				nnm = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member) * 2 * pi * R, R1, R2)[0] * bias
			elif 0:
				f_member, f_non_member = self.member_ratios.ratios(R1, R2)
				nm = f_member * N
				nnm = f_non_member * N
			else:
				nm = N * f_member/f_total
				nnm = N * f_non_member/f_total
				
			print "DSADSADSA"
			print nm+nnm, N
			#print expected, N, bias
			print "\t", nm, N_member, nnm, N-N_member, nm+nnm, N
			non_members_predict.append(nnm)
			members_predict.append(nm)
			test.append(f_member)
			continue
			
			f = scipy.integrate.quad(lambda R: rho_n(R) * 2 * pi * R, Rmin, Rmax)[0] * N_all_member / N_all
			#dsa
			#f = 1#scipy.integrate.quad(lambda R: R, Rmin, Rmax)[0] / 
			#scipy.integrate.quad(lambda R: R, Rmin, self.aperture.aperture_rborders_kpc[-1])[0]
			p = N * f
			print f
			#p = 1./scipy.integrate.quad(lambda R: self.light_profile.densityR(R), Rmin, Rmax)[0]
			#p = scipy.integrate.quad(lambda R: R**1, Rmin, Rmax)[0]/scipy.integrate.quad(lambda R: self.light_profile.densityR(R)*R**1, Rmin, Rmax)[0] 
			#p = nnm
		select(nx-2, ny-1)
		clearautolegend()
		histogram(allstars.vlos_helio, datamin=-200, datamax=300, binwidth=2., normalize=True)
		if self.model:
			f = model
			x = arange(-200, 300, 0.25)
			clearautolegend()
			graph(x, [f(x_i, 0.*x_i) for x_i in x], color="red")
			#autolegend("data", "model")
			autolegend("model fit")
		labels("vlos, helio (km/s)", "p(v)")
			
		
		print "!", sum(test)
		#dsa
		members = array(members)
		non_members = array(non_members)
		counts = array(counts)
		select(nx-2, ny-2)
		clearautolegend()
		graph(non_members, color="green")
		graph(non_members_predict, color="green", linestyle="dash")
		autolegend("counts", "prediction")
		labels("bin number", "# non members")
		select(nx-2, ny-3)
		non_members_predict = array(non_members_predict)
		non_members = array(non_members)
		graph(non_members-non_members_predict, color="green")
		labels("bin number", "counts-prediction (non mem.)")
		hline(0)
		select(nx-3, ny-2)
		graph(members, color="green")
		graph(members_predict, color="green", linestyle="dash")
		labels("bin number", "# members")
		select(nx-3, ny-3)
		members_predict = array(members_predict)
		members = array(members)
		graph(members-members_predict, color="green")
		labels("bin number", "counts-prediction (members)")
		hline(0)
		print "non_members", sum(non_members_predict), sum(non_members)
		select(nx-1, ny-1)
		#y = non_members*1.#/counts
		#y /= sum(y)
		#print y
		print "nm", non_members
		graph(non_members/members, color="green")
		print non_members, members
		print non_members_predict
		#graph(members, color="green")
		#graph(non_members, color="red")
		#fraction_non_member_predict = array(fraction_non_member_predict)
		#fraction_non_member_predict /= sum(fraction_non_member_predict)
		members_predict = array(members_predict)# * 3
		non_members_predict = array(non_members_predict)
		graph(non_members_predict/members_predict, color="black")
		#graph(members_predict, color="black")
		#graph(non_members_predict, color="black")
		print "1)", sum(non_members_predict), N_all_non_member, sum(non_members)
		print "2)", sum(members_predict), N_all_member, sum(members)
		select(nx-3, ny-1)
		ylim(-1, 1)
		hline(0)
		graph(non_members/members - non_members_predict/members_predict, color="black")
		
		select(nx-1, ny-2)
		
		Rkpc = arange(0, 2, 0.01)
		Rarcsec = self.aperture.light_model.kpc_to_arcsec(Rkpc)
		#graph(Rkpc, log10(rho_m(Rkpc)/rho_n(Rkpc)*f_m_in/f_n_in ))
		ratio = rho_m(Rkpc)/rho_n(Rkpc)*f_m_in/f_n_in 
		graph(Rarcsec, log10(ratio))
		index = argmin(abs(Rarcsec-3400))
		print "min", Rarcsec[index], ratio[index]
		vline(3400)
		index = argmin(abs(ratio-0.5))
		print "min", Rarcsec[index], ratio[index]
		
		#print f_m_in, f_n_in
		#print f_m_in+f_n_in
		#dsa
		
		
		
		draw()
		
		
class Counts2(object):
	def __init__(self, observation, aperture, member_filter, light_profile, simplemodel):
		self.observation = observation
		self.aperture = aperture
		self.member_filter = member_filter
		self.light_profile = light_profile
		self.simplemodel = simplemodel
		
	def run(self, args, opts, scope):
		stars = self.observation.load()
		#import pdb
		#pdb.set_trace()
		self.simplemodel.load()
		self.aperture.load()
		stars = stars.filter(lambda star: star.re < self.aperture.aperture_rborders[-1])
		bincount = len(self.aperture)
		nx = int(math.ceil(sqrt(bincount)))+1
		ny = int(math.ceil(float(bincount)/nx))+1
		if nx * ny >= bincount:
			ny += 1
		print bincount, nx, ny
		nx = 3
		ny = 2
		mozaic(nx, ny, box)
		allstars = stars
		members = []
		non_members = []
		counts = []
		members_predict = []
		non_members_predict = []
		N_all_member = len(self.member_filter.filter(allstars))
		N_all = len(allstars)
		N_all_non_member = N_all - N_all_member
		test = []
		
		mean_sigma = 2.
		v1 = self.member_filter.v - 3*self.member_filter.sigma
		v2 = self.member_filter.v + 3*self.member_filter.sigma
		
		
		f_m_in = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v1, v2)[0]
		f_n_in = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v1, v2)[0]
		f_m_out = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v2, inf)[0]
		f_n_out = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v2, inf)[0]
		
		print "=" * 70
		print f_m_in, f_n_in
		print f_m_out, f_n_out 
		print "=" * 70
		#for p in self.galaxy_model.parameters:
		#	print p.name, p.get()
		print v1, v2
		#dsa
		M = numpy.matrix([[f_m_in, f_n_in], [f_m_out, f_n_out]]) #.T
		scaling = (1-0.33)
		scaling = 1.
		print M
		#M = numpy.matrix([[f_m_in, f_n_in], [f_m_out, f_n_out]])#$.T
		#self.light_profile.b *= 1/(1-0.33)
		#self.aperture.aperture_rborders_kpc[-1] -= 0.2
		def rho_m(R):
			return self.light_profile.densityR(R) / scipy.integrate.quad(lambda R: 2*pi*self.light_profile.densityR(R)*R, 0, self.aperture.aperture_rborders_kpc[-1])[0]
		def rho_n(R):
			return 1./scipy.integrate.quad(lambda R: 2 * pi * R, 0, self.aperture.aperture_rborders_kpc[-1]*scaling)[0]
		for i in range(bincount):
			#select(i%nx, i/nx)
			#filter = self.aperture.aperture_filter(i)
			#stars = allstars.filter(filter)
			Rmin = self.aperture.aperture_rborders[i]
			Rmax = self.aperture.aperture_rborders[i+1]
			stars = allstars.filter(lambda star: (star.re >= Rmin) & (star.re < Rmax))
			#histogram(stars.vlos_helio, datamin=-200, datamax=300, binwidth=2.)
			N = len(stars)
			N_member = len(self.member_filter.filter(stars))
			counts.append(N)
			Rmin = self.aperture.aperture_rborders_kpc[i]
			Rmax = self.aperture.aperture_rborders_kpc[i+1]
			
			#print self.member_filter
			#import pdb
			#pdb.set_trace()
			
			C_in = N_member
			C_out = N - C_in
			if 0:
				members.append(N_member)
				non_members.append(N-N_member)
			else:
				print "before:", N_member, N-N_member
				print M.I
				y = dot(M.I, array([C_in, C_out])).T
				print y
				y = array(y.flat)
				Nm = y[0]
				Nn = y[1]
				print "after:", Nm, Nn
				members.append(Nm)
				non_members.append(Nn)
				if Nm < 0:
					import pdb
					#pdb.set_trace()
				
			
			#C_in = len(self.member_filter.filter(stars))
			#M = [[
			
			
			#print scipy.integrate.quad(lambda R: rho_n(R) * 2 * pi * R, 0, self.aperture.aperture_rborders_kpc[-1])[0]
			
			R1 = Rmin
			R2 = Rmax
			#R1 = 0
			#R2 = self.aperture.aperture_rborders_kpc[-1]
			#print dir(self.model.model_wrapper.model)
			#print self.model.model_wrapper.model.fraction_foreground
			#dsa
			#f = 0.264414637969
			#f = self.model.fitter.get_model().fraction_foreground
			#f = self.model.fitter.get_model().fraction_foreground
			if 0:
				model = self.model.fitter.get_model()
				#f = 1.
				ratio = abs(model.ratio)
				w1 =  ratio / (1. + ratio);
				w2 = 1 / (1. + ratio);
				w1 = w1*1.7# * 0.7
				f = w1
				#w1 = 0.84
				print "ratio,w1,w2", ratio, w1, w2, w1/w2
				#f = 0.4
			#f = 0.3
			#w1 = 0.14
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all) * 2 * pi * R, R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*(1-f)) * 2 * pi * R, R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*f) * 2 * pi * R, R1, R2)[0]
			#f_total = f_member + f_non_member 
			#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all)/ (rho_m(R)*N_all_member/N_all + rho_n(R)*N_all_non_member/N_all) , R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all)/ (rho_m(R)*N_all_member/N_all + rho_n(R)*N_all_non_member/N_all) , R1, R2)[0]
			#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
			#print "ftotal", f_total
			
			if 1:
				mcounts = 0
				nmcounts = 0
				for star in stars:
					R = self.aperture.light_model.arcsec_to_kpc(star.re)
					#f_member = rho_m(R)*(1-w1)
					#f_non_member = rho_n(R)*(w1)
					if 1:
						counts1 = 0
						f_member = 0
						f_non_member = 0
						for j in range(2):
							if (1<<j) & int(star.catalogue_mask):
								f_member     += self.simplemodel.memberRe(star.re,catalogue=j)
								f_non_member += self.simplemodel.non_memberRe(star.re,catalogue=j)
								counts1 += 1
						f_member /= counts1
						f_non_member /= counts1
						f = 1.
					else:
						f_member = self.simplemodel.memberRetot(star.re, star.catalogue_mask)
						f_non_member = 1-f_member
						f = 1
					f_member = f_member / f
					f_non_member = f_non_member/ f
					mcounts += f_member
					nmcounts += f_non_member
				f_member = mcounts / N
				f_non_member = nmcounts / N
				nm = f_member * N
				nnm = f_non_member * N
				print ">>>>", f_member, f_non_member, nm, nnm, nm+nnm, N
			elif 1:
				expected = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member + rho_m(R)*N_all_member) * 2 * pi * R, R1, R2)[0]
				bias = N/expected
				#dsa
				#* N_all_member / N_all
				nm = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member) * 2 * pi * R, R1, R2)[0] * bias
				nnm = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member) * 2 * pi * R, R1, R2)[0] * bias
			elif 0:
				f_member, f_non_member = self.member_ratios.ratios(R1, R2)
				nm = f_member * N
				nnm = f_non_member * N
			else:
				nm = N * f_member/f_total
				nnm = N * f_non_member/f_total
				
			print "DSADSADSA"
			print nm+nnm, N
			#print expected, N, bias
			print "\t", nm, N_member, nnm, N-N_member, nm+nnm, N
			non_members_predict.append(nnm)
			members_predict.append(nm)
			test.append(f_member)
			continue
			
			f = scipy.integrate.quad(lambda R: rho_n(R) * 2 * pi * R, Rmin, Rmax)[0] * N_all_member / N_all
			#dsa
			#f = 1#scipy.integrate.quad(lambda R: R, Rmin, Rmax)[0] / 
			#scipy.integrate.quad(lambda R: R, Rmin, self.aperture.aperture_rborders_kpc[-1])[0]
			p = N * f
			print f
			#p = 1./scipy.integrate.quad(lambda R: self.light_profile.densityR(R), Rmin, Rmax)[0]
			#p = scipy.integrate.quad(lambda R: R**1, Rmin, Rmax)[0]/scipy.integrate.quad(lambda R: self.light_profile.densityR(R)*R**1, Rmin, Rmax)[0] 
			#p = nnm
		#select(nx-2, ny-1)
		select(2,0)
		clearautolegend()
		histogram(allstars.vlos_helio, datamin=-200, datamax=300, binwidth=2., normalize=True)
		vs = arange(-200, 300, 0.5)
		distribution = vs * 0
		if 1:
			for star in allstars:
				distribution1 = distribution * 0
				counts = 0
				for j in range(2):
					if (1<<j) & int(star.catalogue_mask):
						#$pm = self.simplemodel.memberRe(star.re, catalogue=j)
						#pnm = 1-pm
						distribution1 += exp(array([self.simplemodel.logL(star.re, v, star.e_vlos, j) for v in vs]))
						print sum(distribution1)
						counts += 1
				distribution += distribution1/counts
		distribution /= len(allstars)
		graph(vs, distribution, color="red")
					
		if 0:#self.model:
			f = model
			x = arange(-200, 300, 0.25)
			clearautolegend()
			graph(x, [f(x_i, 0.*x_i) for x_i in x], color="red")
			#autolegend("data", "model")
			autolegend("model fit")
		labels("vlos, helio (km/s)", "p(v)")
			
		
		print "!", sum(test)
		#dsa
		members = array(members)
		non_members = array(non_members)
		counts = array(counts)
		select(1, 1)
		clearautolegend()
		graph(non_members, color="green")
		graph(non_members_predict, color="green", linestyle="dash")
		autolegend("counts", "prediction")
		labels("bin number", "# non members")
		select(1, 0)
		non_members_predict = array(non_members_predict)
		non_members = array(non_members)
		graph(non_members-non_members_predict, color="green")
		labels("bin number", "counts-prediction (non mem.)")
		ylim(-20, 20)
		hline(0)
		if 1:
			select(0, 1)
			graph(members, color="green")
			graph(members_predict, color="green", linestyle="dash")
			labels("bin number", "# members")
		select(0, 0)
		members_predict = array(members_predict)
		members = array(members)
		graph(members-members_predict, color="green")
		labels("bin number", "counts-prediction (members)")
		ylim(-20, 20)
		hline(0)
		if 0:
			print "non_members", sum(non_members_predict), sum(non_members)
			select(nx-1, ny-1)
			#y = non_members*1.#/counts
			#y /= sum(y)
			#print y
			print "nm", non_members
			graph(non_members/members, color="green")
			print non_members, members
			print non_members_predict
			#graph(members, color="green")
			#graph(non_members, color="red")
			#fraction_non_member_predict = array(fraction_non_member_predict)
			#fraction_non_member_predict /= sum(fraction_non_member_predict)
			members_predict = array(members_predict)# * 3
			non_members_predict = array(non_members_predict)
			graph(non_members_predict/members_predict, color="black")
			#graph(members_predict, color="black")
			#graph(non_members_predict, color="black")
			print "1)", sum(non_members_predict), N_all_non_member, sum(non_members)
			print "2)", sum(members_predict), N_all_member, sum(members)
		if 0:
			select(nx-3, ny-1)
			ylim(-1, 1)
			hline(0)
			graph(non_members/members - non_members_predict/members_predict, color="black")
		
		#if 0:
		select(2,1)
		
		Rkpc = arange(0, 2, 0.01)
		Rarcsec = self.aperture.light_model.kpc_to_arcsec(Rkpc)
		#graph(Rkpc, log10(rho_m(Rkpc)/rho_n(Rkpc)*f_m_in/f_n_in ))
		if 0:
			ratio = rho_m(Rkpc)/rho_n(Rkpc)*f_m_in/f_n_in 
			graph(Rarcsec, log10(ratio))
		
		#for color, i in zip(nicecolors[1:], range(self.simplemodel.catalogues)):
		weight = f_m_in/f_n_in # additional weight, since we are interested in the '+/- 3 sigma region
		#for color, i in zip(nicecolors[1:], range(1)):
		for color, i in zip(nicecolors[1:], range(self.simplemodel.catalogues)):
			f_member = self.simplemodel.memberRe(Rarcsec, weight=weight, catalogue=i) # i)
			f_non_member = self.simplemodel.non_memberRe(Rarcsec, weight=weight, catalogue=i)#, i)
			#graph(Rarcsec, log10(f_member), color=color)
			#graph(Rarcsec, log10(f_non_member), color=color, linestyle="dash")
			ratio = f_member/f_non_member
			graph(Rkpc, log10(ratio), color=color)
			index = argmin(abs(ratio-3))
			print "min", Rarcsec[index], ratio[index]
			#vline(Rarcsec[index], color=color)
			vline(Rkpc[index], color=color)
		labels("R (kpc)", "ratio of 3 sigma members")
		
		#vline(3400)
		#ylim(-1, 2)

		if 0:
			index = argmin(abs(Rarcsec-3400))
			print "min", Rarcsec[index], ratio[index]
			index = argmin(abs(ratio-0.5))
			print "min", Rarcsec[index], ratio[index]
		
		#print f_m_in, f_n_in
		#print f_m_in+f_n_in
		#dsa
		
		
		
		draw()
		
		
		
		
class MemberDensities(object):
	def __init__(self, observation, member_filter, galatic_distribution):
		self.observation = observation
		self.member_filter = member_filter
		self.galatic_distribution = galatic_distribution
		
	def run(self, args, opts, scope):
		allstars = self.observation.load()
		
		members = self.member_filter.filter(allstars)
		non_members = self.member_filter.filter(allstars, invert=True)
		print scipy.integrate.quad(lambda v: self.galatic_distribution(v, 0), -inf, inf)
		dsa
		
		
		#non_members = allstars.filter(lambda star: star not in members)
		#box()
		mozaic(2,2,box)
		select(0,0)
		h1 = histogram(members.rc, datamin=0, datamax=8000, bincount=40)
		h2 = histogram(non_members.rc, datamin=0, datamax=8000, bincount=40, color="red")
		vline(3400)
		select(1,0)
		mask = (h2.data > 0) & (h1.data > 0)
		print mask
		#pr
		c = h1.data[mask]/h2.data[mask]
		print log10(c), (c)
		#scipy.stats.gaussian_kde(members.rc, 
		graph(h1.bins[mask], log10(c))
		vline(3400)
		
		bla = scipy.stats.gaussian_kde(members.rc)
		rs = arange(0, 8000)
		y = bla(rs)
		select(0, 1)
		graph(rs, y * len(members))
		
		bla = scipy.stats.gaussian_kde(non_members.rc)
		rs = arange(0, 8000)
		y_nm = bla(rs)
		graph(rs, y_nm * len(non_members))
		
		y = y / sum(y) * len(members)
		y_nm = y_nm / sum(y_nm) * len(non_members)
		y = y[:len(y)/2]
		y_nm = y_nm[:len(y)]
		select(1,1)
		ratio = y/y_nm
		graph(rs[:len(y)], log10(ratio))
		mask = ratio < 0.333
		#graph(rs[:len(y)][mask], log10(ratio)[mask], color="red")
		distance = abs(ratio-1./3)
		index = argmin(distance)
		Rcut = rs[index]
		vline(Rcut)
		print "Rcut", Rcut
		print ratio[3400]
		vline(3400)
		
		
		
		#import pdb
		#pdb.set_trace()
		
		
		draw()
		#print len(allstars)
		#print len(members)
		#print len(non_members)
		#print len(members)+len(non_members)
		
		
		
class BestFits(object):
	def __init__(self, bestfit1, bestfit2):
		self.bestfit1 = bestfit1
		self.bestfit2 = bestfit2
		
		
	def run(self, opts, args, scope):
		#import pdb
		#pdb.set_trace()
		document(size="30cm,15cm")
		page(fontsize="12pt")
		mozaic(4,2,box)
		self.bestfit1.load()
		losvd_per_bin1 = self.bestfit1.subscope["losvd_per_bin"]
		losvd_per_bin1.load()
		object1 = self.do(losvd_per_bin1, scope=self.bestfit1.subscope, labels=True, color="red", linestyle="dash")

		self.bestfit2.load()
		losvd_per_bin2 = self.bestfit2.subscope["losvd_per_bin"]
		losvd_per_bin2.load()
		object2 = self.do(losvd_per_bin2, scope=self.bestfit2.subscope, labels=True, color="blue", linestyle="dot")
		
		bins = losvd_per_bin2.losvd_per_bin.shape[1]
		print "bins", bins
		print len(losvd_per_bin1.ks_p)
		for j in range(2):
			for i in range(4):
				index = i + j * 4
				if index < bins:
					select(i, 1-j)
					#if labels:
					#	title("bin #%d" % index)
					#graph(vs, losvd[:,index], **kwargs)
					histogram(losvd_per_bin2.stars_group[index].vlos, datamin=-50, datamax=50, binwidth=2., normalize=True)
					labels("v<sub>los,sys</sub> (km/s)", "p(v)")
					D, p = losvd_per_bin1.ks_D[index], losvd_per_bin1.ks_p[index]
					#legend(["line", "line"], ["cored (D=%.3f, p=%.3f)" % (D, p), "NFW (D=%.3f, p=%.3f)" % (D, p)], [object1, object2], linelength="0cm", edgespacing="2mm", fontsize="8pt", location="left, right")
					fontsize = "10pt"
					legend(["line"], ["alpha=0\nD=%.3f\np=%.4f" % (D, p)], [object1], linelength="0cm", edgespacing="2mm", fontsize=fontsize, location="left, top")
					D, p = losvd_per_bin2.ks_D[index], losvd_per_bin2.ks_p[index]
					legend(["line"], ["NFW\nD=%.3f\np=%.4f" % (D, p)], [object2], linelength="0cm", edgespacing="2mm", fontsize=fontsize, location="right, top")
		draw()
		
	def do(self, losvd_per_bin, scope, labels=False, **kwargs):
		losvd = losvd_per_bin.losvd_per_bin
		bins = losvd.shape[1]
		vs = losvd_per_bin.storage_losvd.vcenters
		aperture = scope["aperture_m2"]
		aperture.load()
		
		for j in range(2):
			for i in range(4):
				index = i + j * 4
				select(i, 1-j)
				if index < bins:
					if labels:
						#title("bin #%d" % (index+1))
						rs = aperture.aperture_rborders_kpc
						title("R = %.2f-%.2f kpc" % (rs[index], rs[index+1]))
					g = graph(vs, losvd[:,index], **kwargs)
					ylim(0, 0.065)
				else:
					current.container.context.color = "white"
				
		return g
				
		
		
class DMMassProfiles(object):
	def __init__(self, subscopes, contexts):
		self.subscopes = subscopes
		self.contexts = contexts
		
	def run(self, args, opts, scope_):
		document(size="8cm,16cm")
		mozaic(1,4,box)
		#first = True
		logrs = arange(-1.5, 1.5 +0.01, 0.1/4)
		rs = 10**logrs
		
		logrborders = arange(-1.5-0.1/4/2, 1.5+0.1/4/2 +0.01, 0.1/4)
		assert len(logrborders) == (1+len(rs))
		
		scales = []
		light_profile = profile = scope_["light_profile"]
		light_model = scope_["light_model"]
		scales.append(profile.scale)
		Mtot = profile.enclosed_mass(inf)
		def f(r):
			mass = profile.enclosed_mass(r)[0]/Mtot
			return (mass-0.5)**2
		#import pdb
		#pdb.set_trace()
		masses = array([profile.cumdensityR(0, R)/Mtot for R in rs])
		rhalf = scipy.optimize.fmin(f, profile.scale*10)[0]
		print "rhalf 2", rhalf, profile.scale
		scales.append(rhalf)
		obs = scope_["observation_clean"]
		filters = scope_["filters_aperture"]
		stars = obs.load()
		for filter in filters:
			stars = filter(stars)
		rc = light_model.arcsec_to_kpc(stars.rc )
		rc_max = max(rc)
		
		def f(r):
			r = abs(r)
			s = profile.logslope(r)[0]
			return (s+3.0)**2
		#import pdb
		#pdb.set_trace()
		
		Mmax = profile.cumdensityR(0, rc_max)
		Mmax_total = profile.cumdensityR(0, 1e3)
		#masses = array([profile.cumdensityR(0, R)/Mmax for R in rs])
		masses = array([profile.cumdensityR(0, R)/Mmax_total for R in rs])
		r3 = abs(scipy.optimize.fmin(f, profile.scale*10)[0])
		print "r3", r3

		def f(r):
			mass = profile.cumdensityR(0, r, 1.)
			return (mass-0.5)**2
		Rhalf = scipy.optimize.fmin(f, profile.scale*10)[0]
		
		name = scope_["object_name"]
		print(name)
		with open(name + "_log10r.txt", "w") as f:
		  f.write(",".join(map(str, logrs.tolist())))
		
		
		mass_profiles = []
		for contexts, subscope in zip(self.contexts, self.subscopes):
			print subscope.includes
			subscope.load()
			#scope = scope.scope
			profile = subscope.subscope["dm_profile"]
			select(0, 1)
			graph(logrs, profile.logslope(rs), **contexts)

			select(0, 2)
			graph(logrs, log10(profile.enclosed_mass(rs)), **contexts)
			mass = log10(profile.enclosed_mass(rs))
			mass_profiles.append(mass)
			pname = subscope.subscope["model_set"]
			with open(name + "_" + pname +"_log10M.txt", "w") as f:
			  f.write(",".join(map(str, mass.tolist())))
		select(0, 2)
		labels("log r/kpc", "log M(r) / M<sub>sun</sub>")
		light_profile = subscope.subscope["light_profile"]
		graph(logrs, log10(light_profile.enclosed_mass(rs)), color="black")

		mass = log10(light_profile.enclosed_mass(rs))
		pname = "light"
		with open(name + "_" + pname +"_log10M.txt", "w") as f:
		  f.write(",".join(map(str, mass.tolist())))
		raise "dsa"

		vline(log10(r3), color="red", linestyle="dash")
		vline(log10(rhalf), color="black", linestyle="normal")
		ylim(3, 10.)
		#for scale in scales:
		#	vline(log10(scale))
		#histogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=30, normalize=True)
		#cumhistogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=100, normalize=True)

		select(0, 1)
		graph(logrs, profile.logslope(rs), **contexts)
		labels("log r/kpc", "dlog&rho;(r)/dlogr")
		#for scale in scales:
		#	vline(log10(scale))
		#histogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=30, normalize=True)
		#cumhistogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=100, normalize=True)
		#histogramline(logrborders, masses, color="red")
		graph(logrs, light_profile.logslope(rs), color="black")
		#for x in [-1.0, -1.2, -1.4, -3.]:
		#	hline(x)
		vline(log10(r3), color="red", linestyle="dash")
		vline(log10(rhalf), color="black", linestyle="normal")
		ylim(-4, 1)

		
		if 0:
			select(0, 1)
			labels("log r/kpc", "log M(r) / M<sub>nfw</sub>(r)")
			
			scope1 = self.subscopes[0]
			profile = scope1.subscope["dm_profile"]
			logM1 = log10(profile.enclosed_mass(rs))
			for contexts, subscope in zip(self.contexts, self.subscopes)[1:]:
				profile = subscope.subscope["dm_profile"]
				graph(logrs, log10(profile.enclosed_mass(rs)) - logM1 , **contexts)
				hline(0)
			#for scale in scales:
			#	vline(log10(scale))
			#histogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=30, normalize=True)
			#cumhistogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=100, normalize=True)
			#histogramline(logrborders, masses, color="red")
			line_r3 = vline(log10(r3), color="red", linestyle="dash")
			line_rhalf = vline(log10(rhalf), color="black", linestyle="normal")
			legend(["line"] * 2, ["r<sub>-3</sub>", "r<sub>1/2</sub>"], [line_r3, line_rhalf], location='right, top', edgespacing="4mm")
			#ylim(-0.75, 0.75)
			ylim(-1.0, 1.)
		
		if 1:
			select(0, 3)
			name = scope_["object_name"]
			print "slopes of ", name,"",
			slopes = []
			logr1 = log10(r3)-0.75
			logr2 = log10(r3)+0.75
			maxdiff = None
			for i in range(len(mass_profiles)):
				for j in range(i+1, len(mass_profiles)):
					if maxdiff is None:
						maxdiff = abs(mass_profiles[i]-mass_profiles[j])
					else:
						maxdiff = maximum(maxdiff, mass_profiles[i]-mass_profiles[j])
			
					
			i1 = argmin(abs(logrs-logr1))
			i2 = argmin(abs(logrs-logr2))
			maxdiff1 = maxdiff * 1.
			maxdiff2 = maxdiff * 1.
			i3 = argmin(abs(logrs-log10(r3)))
			maxdiff2[:i3] = 1e4
			maxdiff1[i3:] = 1e4
			i1 = argmin(abs(maxdiff1-0.05))
			i2 = argmin(abs(maxdiff2-0.05))
			print i1, i2
			for contexts, subscope in zip(self.contexts, self.subscopes):
				profile = subscope.subscope["dm_profile"]
				#rhoavg = profile.enclosed_mass(rs) * 3 / (4*pi*rs**3)
				#gamma  = 3 * (1-profile.densityr(rs)/rhoavg)
				dmdr = 4 * pi * rs**2 * profile.densityr(rs)
				M = profile.enclosed_mass(rs)
				dlogmdlogr = rs/M * dmdr
				#print gamma
				#graph(logrs, gamma, **contexts)
				x = dlogmdlogr[i1:i2]
				#print mean(x),
				slopes.append(mean(x))
				select(0, 3)
				graph(logrs, dlogmdlogr, **contexts)
				#select(0, 4)
				#y = log10(profile.densityr(rs))
				#graph(logrs, y, **contexts)
			def f(params):
				a, b = params
				#print "param:", a, b
				y = log10(a*(10**(logrs[i1:i2]-log10(r3)))**b)
				diffs = []
				for mass in mass_profiles:
					diffs.extend(y-mass[i1:i2])
				return array(diffs)
			
			res = scipy.optimize.leastsq(f, [7, 1.4], full_output=True)
			params = res[0]
			jac = res[1]
			covariance = numpy.matrix(jac).I
			a, b = log10(params[0]), params[1]
			print "%s: fit mass = %e r ^ %.2f" % (name, 10**a , b)
			y = log10(10**a *(10**(logrs-log10(r3)))**b)
			select(0, 2)
			graph(logrs, y, linestyle="dash")
			#print covariance
			#print res
			#import pdb
			#pdb.set_trace()
				#def f(x):
					
			select(0, 2)
			vline(logrs[i1])
			vline(logrs[i2])
			select(0, 3)
			vline(logrs[i1])
			vline(logrs[i2])
			print ", mean: ", mean(slopes)
			
			#labels("log r/kpc", "gamma")
			#labels("log r/kpc", "dlogM/dlogr")
			select(0, 3)
			labels("log r/kpc", "dlogM/dlogr")
			vline(log10(r3), color="red", linestyle="dash")
			vline(log10(rhalf), color="black", linestyle="normal")
			if 0:
				select(0, 4)
				labels("log r/kpc", "rho(r)")
				vline(log10(r3), color="red", linestyle="dash")
				vline(log10(rhalf), color="black", linestyle="normal")


		if 1:
			select(0,0)
			vline(log10(r3), color="red", linestyle="dash")
			vline(log10(rhalf), color="black", linestyle="normal")
			#vline(log10(Rhalf), color="blue", linestyle="normal")
			cumhistogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=100, normalize=True)
			histogramline(logrborders, masses, color="red")
			#histogramline(logrborders, masses*Mmax/Mmax_total, color="blue")
			labels("log R/kpc", "f(&lt;R)")
			ylim(0, 1.1)
		#select(0, 0)
		#labels("log r/kpc", "dlog&rho;(r)/dlogr - nfw")
		#for scale in scales:
		#	vline(log10(scale))
		#graph(logrs, light_profile.logslope(rs) - slope1)
		#histogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=30, normalize=True)
		#cumhistogram(log10(rc), datamin=-1.5, datamax=1.5, bincount=100, normalize=True)
		#histogramline(logrborders, masses, color="red")
		#ylim(-3, 2)
			
		draw()
		
class DMMassProfilesAmina(object):
	def __init__(self, subscopes, contexts, mass_measurements=[], extra_contexts=[]):
		self.subscopes = subscopes
		self.contexts = contexts
		#self.extra_profiles = extra_profiles
		self.mass_measurements = mass_measurements
		self.extra_contexts = extra_contexts
		
	def run(self, args, opts, scope_):
		document(size="10cm,10cm")
		#mozaic(1,3,box)
		#first = True
		logrs = arange(-1.5, 1.5 +0.01, 0.1/4)
		rs = 10**logrs
		
		logrborders = arange(-1.5-0.1/4/2, 1.5+0.1/4/2 +0.01, 0.1/4)
		assert len(logrborders) == (1+len(rs))
		
		scales = []
		light_profile = profile = scope_["light_profile"]
		light_model = scope_["light_model"]
		scales.append(profile.scale)
		Mtot = profile.enclosed_mass(inf)
		def f(r):
			mass = profile.enclosed_mass(r)[0]/Mtot
			return (mass-0.5)**2
		#import pdb
		#pdb.set_trace()
		masses = array([profile.cumdensityR(0, R)/Mtot for R in rs])
		rhalf = scipy.optimize.fmin(f, profile.scale*10)[0]
		print "rhalf 2", rhalf, profile.scale
		scales.append(rhalf)
		obs = scope_["observation_clean"]
		filters = scope_["filters_aperture"]
		stars = obs.load()
		for filter in filters:
			stars = filter(stars)
		rc = light_model.arcsec_to_kpc(stars.rc )
		rc_max = max(rc)
		
		def f(r):
			r = abs(r)
			s = profile.logslope(r)[0]
			return (s+3.0)**2
		#import pdb
		#pdb.set_trace()
		
		Mmax = profile.cumdensityR(0, rc_max)
		masses = array([profile.cumdensityR(0, R)/Mmax for R in rs])
		r3 = abs(scipy.optimize.fmin(f, profile.scale*10)[0])
		print "r3", r3
		
		
		
		names = []
		for contexts, subscope in zip(self.contexts, self.subscopes):
			print subscope.includes
			subscope.load()
			#scope = scope.scope
			profile = subscope.subscope["dm_profile"]
			#select(0, 1)
			#graph(logrs, profile.logslope(rs), **contexts)
			names.append(subscope.subscope["dm_name"])
			#select(0, 2)
			graph(logrs, log10(profile.enclosed_mass(rs)), **contexts)
		#mass_14 = rs**1.4 / 0.3**1.4 * 1.1e7
		#g = graph(logrs, log10(mass_14), color="orange") #, linestyle="dash")
		#names.append("r<sup>1.4</sup>")
		autolegend(*names, location="left, top", edgespacing="3mm")
		
		clearautolegend()
		names = []
		for contexts, mass_measurements in zip(self.extra_contexts, self.mass_measurements):
			r, logm, logmpos, logmneg, name = mass_measurements
			names.append(name)
			symbol(log10(r), logm, **contexts) #symbolName="plus")
			errpos = logmpos - logm
			errneg = logmpos - logm
			#errorbars([log10(r)], [logm], ypos=[errpos], yneg=[errneg], **contexts)
			#graph(logrs, log10(profile.enclosed_mass(rs)), **contexts)

		if names:
			autolegend(*names, location="right, bottom", edgespacing="3mm")
		#select(0, 2)
		labels("log r/kpc", "log M(r) / M<sub>sun</sub>")
		light_profile = subscope.subscope["light_profile"]
		graph(logrs, log10(light_profile.enclosed_mass(rs)), color="black")
		#vline(log10(r3), color="red", linestyle="dash")
		#vline(log10(rhalf), color="black", linestyle="normal")
		draw()
		
		
class PdfCombined(ParameterSet):
	def __init__(self, obs_groups):
		self.obs_groups = obs_groups
		
		
	def run(self, args, opts, scope):
		grids = []
		grid = None
		for obs_group in range(self.obs_groups):
			subscope = scope.clone()
			subscope["obs_group"] = obs_group
			parameterset = subscope["command_plot_parameter_set"]
			parameterset.load()
			grids.append(parameterset.smooth)
			if grid is None:
				grid = parameterset.smooth ** 0.2
			else:
				grid = grid * parameterset.smooth ** 0.8
				pass
		
		box()
		#self.smooth = grid
		grid = grid
		self.parameterset = parameterset.parameterset
		#self.plot_parameters(0, 1, scope)
		i, j = 0, 1
		sigmas = 3
		color = "black"
		probimage2d(grid, i, j, resize=self.makeresize(i,j), drawcontourlines=True, fill=False, colormap=None, premultiply=True, color=color, sigmas=sigmas)
		scope["obs_bias"] = ""
		parameterset = scope["command_plot_parameter_set"]
		
		parameterset.load()
		parameterset.plot_parameters(0, 1, scope, color="red")
		draw()
		
		
		

class PlotSummary(object):
	def __init__(self, observations, binned_data, light_model):
		self.observations = observations
		self.binned_data = binned_data
		self.light_model = light_model
		
	def run(self, opts, args, scope):
		self.binned_data.load()
		stars = self.observations.load()
		#box()
		page(fontsize="15pt")
		mozaic(1,2, box)
		
		select(0, 1)
		#graph(
		self.binned_data.aperture.load()
		kpc = self.binned_data.aperture.aperture_rcenters_kpc
		print dir(self.binned_data)
		print (self.binned_data.moments.shape)
		scatter(kpc, self.binned_data.moments[2,:]**0.5, symbolsize="25pt")
		xlim(0, 1.5)
		ylim(0, 15)
		labels("afstand tot het centrum in kiloparsec", "snelheidsdispersie")
		
		select(0, 0)
		xlim(0, 1.5)
		kpc = self.light_model.arcsec_to_kpc(stars.rc)
		scatter(kpc, stars.vlos)
		xlim(0, 1.5)
		labels("afstand tot het centrum in kiloparsec", "snelheid van ster")
		ylim(-40, 40)
		
		draw()
		
		
		
		
		
class BesanconTest(object):
	def __init__(self, data):
		self.data = data
		
	def run(self, args, opts, scope):
		data = self.data.load()
		
		print data.l
		print data.b
		print len(data.l)
		vsplit(box)
		#box()
		scatter(data.l, data.b)
		symbol(mean([min(data.l), max(data.l)]), mean([min(data.b), max(data.b)]), color="red")
		width = max(data.l) - min(data.l)
		height = max(data.b) - min(data.b)
		ml = (min(data.l) + max(data.l))/2
		mb = (min(data.b) + max(data.b))/2
		print ml, mb
		symbol(ml, mb, color="green")
		vline(ml)
		hline(mb)
		#def count(
		left = ml - width*0.15
		vline(left)
		right = ml + width*0.15
		vline(right)
		print sum(data.l < left), sum(data.l > right)

		top = mb - height*0.15
		hline(top)
		bottom = mb + height*0.15
		hline(bottom)
		print sum(data.b < bottom), sum(data.b > top)

		symbol(min(data.l), min(data.b), color="red")
		symbol(min(data.l), max(data.b), color="red")
		symbol(max(data.l), min(data.b), color="red")
		symbol(max(data.l), max(data.b), color="red")
		
		a = histogram2d(data.l, data.b, [5,5])
		d = a[0]
		print a
		select(1)
		indexedimage(d, colormap="whiteblack")
		
		select(0)
		for i in a[1]:
			vline(i, color="green", linestyle="dash")
		for i in a[2]:
			hline(i, color="green", linestyle="dash")
		
		grow(1.1)
		
		draw()
		
class PlotGovernatoSlope(object):
	def __init__(self, data):
		self.data = data
	
	def run(self, opts, args, scope):
		document(size="10cm,15cm")
		#$box()
		hsplit(box)
		
		masses = mass = 10**arange(log10(1e3/2), log10(2e10), 0.2)
		masses_vir = m_vir = 1e10*(masses/2e7)**0.5
		masses = mass = 10**arange(log10(4*1e8), 12.2, 0.1)
		alpha = -0.5 + 0.35 * log10(mass/1e8)
		graph(log10(mass), alpha, color="red")
		if 0:
			vline(6)
			vline(7)
		if 0:
			hline(-0.5 + 0.35 * log10(1e6/1e8))
			hline(-0.5 + 0.35 * log10(1e7/1e8))
		for d in self.data:
			name, mass, value, sigma = d
			symbol(log10(mass), value)
			line(log10(mass), value-sigma, log10(mass), value+sigma, color="blue")
			
		alphas_stellar = []
		import mab
		for m in masses:
			m_vir = 1e10*(m/2e7)**0.5
			m_vir = m
			c_target = mab.utils.mass_concentration.c_maccio(m_vir, h=0.72)
			def f(logrs):
				rs = 10**logrs[0]
				nfw = mab.gd.potential.NFW(m_vir, rs)
				#print ">", rs, nfw.c, c_target
				return ((c_target - nfw.c)/c_target)**2
			logrs = scipy.optimize.fmin(f, 0)
			print logrs
			rs = 10**logrs[0]
			print m, m_vir, rs, 
			nfw = mab.gd.potential.NFW(m_vir, rs)
			alpha = nfw.logslope(0.5)
			print nfw.c, c_target, alpha
			alphas_stellar.append(alpha)
		graph(log10(masses), alphas_stellar)
		
		select(1)
		alphas_vir = []
		for m_vir in masses_vir:
			c_target = mab.utils.mass_concentration.c_maccio(m_vir, h=0.72)
			def f(logrs):
				rs = 10**logrs[0]
				nfw = mab.gd.potential.NFW(m_vir, rs)
				return ((c_target - nfw.c)/c_target)**2
			logrs = scipy.optimize.fmin(f, 0)
			rs = 10**logrs[0]
			nfw = mab.gd.potential.NFW(m_vir, rs)
			alpha = nfw.logslope(0.5)
			alphas_vir.append(alpha)
		graph(log10(masses_vir), alphas_vir)
		masses_stellar = masses_vir**2/1e10**2 * 2e7
		alphas = -0.5 + 0.35 * log10(masses_stellar/1e8)
		graph(log10(masses_vir), alphas, color="red")
		ylim(-2.2, 0.7)
		draw()
		
		
		
		
class PlotMassProfileDfs(object):
	def __init__(self, solution):
		self.solution = solution
		
	def run(self, args, opts, scope):
		#box()
		
		
		try:
			vdisp = scope["vdisp"]
		except:
			vdisp = False
		try:
			anisotropy = scope["anisotropy"]
		except:
			anisotropy = False
		try:
			light = scope["light"]
		except:
			light = False
		try:
			plotdf = scope["plotdf"]
		except:
			plotdf = False
		try:
			plot_density = scope["density"]
		except:
			plot_density = False
		if vdisp or light or plotdf or anisotropy:
			document(size="10cm,10.5cm")
			page(fontsize="17pt")
			box()
			#if not plot_density:
			if 1:
				current.container.drawOutsideRight = True #scope["doleft"]
				current.container.drawOutsideTop = True #scope["doleft"]
				current.container.drawOutsideLeft = True #scope["doleft"]
				current.container.drawOutsideBottom = True #scope["dobottom"]
				current.container.leftaxes[0].drawLabels = scope["doleft"]
				current.container.bottomaxes[0].drawLabels = scope["dobottom"]
		else:
			document(size="30cm,25cm")
			page(fontsize="18pt")
			mozaic(3,3,box)
			#draw()
			#document(size="30cm,45cm")
			#page(fontsize="18pt")
			#mozaic(3,6,box)
		self.solution.load()
		df = self.solution.orbitweights
		nE, nL = 20, 8
		df = df.reshape((nL, nE))
		df  /= sum(df)
		
		
		if 0:
			import Watershed
			ws = Watershed.Watershed()
			
		if 1:
			from scipy import ndimage
			df = (df-df.min())/(df.max()-df.min()) * 255
			#df = 255-df
			debug = False
			if debug:
				select(0, 0)
				indexedimage(df)
			if 0:
				markers = (df * 0).astype(int8)
				markers[1,8] = 1
				markers[7,8] = 2
				
				select(1, 0)
				markers[df<10] = -1
				
				image = df
				Nx, Ny = markers.shape
				maxiters = 1000
				iteration = 0

				moved = True
				while moved:
					moved = False
					for i in range(Nx)[::-1]:
						#for j in range(Ny):
						for j in range(Ny)[::-1]:
							value = image[i,j]
							marker = markers[i,j]
							print i, j, marker, value
							if markers[i,j] > 0:
								if i > 0: # left
									if (markers[i-1,j] == 0) and (image[i-1,j] > value):
										markers[i-1,j] = marker
										markers[i,j] = 0
										print "left"
										moved = True
								if i < Nx-1: # right
									if (markers[i+1,j] == 0) and (image[i+1,j] > value):
										markers[i+1,j] = marker
										markers[i,j] = 0
										print "right"
										moved = True
								if j > 0: # above
									if (markers[i,j-1] == 0) and (image[i,j-1] > value):
										markers[i,j-1] = marker
										markers[i,j] = 0
										print "above"
										moved = True
								if j < Ny -1: # below
									if (markers[i,j+1] == 0) and (image[i,j+1] > value):
										markers[i,j+1] = marker
										markers[i,j] = 0
										print "below"
										moved = True
				
				moved = True
				while (sum(markers == 0) > 0) and moved: # if contains zeros.. continue
					print "#" * 70
					print markers
					for i in range(Nx)[::-1]:
						for j in range(Ny)[::-1]:
							value = image[i,j]
							marker = markers[i,j]
							print i, j, marker, value
							if markers[i,j] > 0:
								if i > 0: # left
									if (markers[i-1,j] == 0) and (image[i-1,j] < value):
										markers[i-1,j] = marker
										print "left"
										moved = True
								if i < Nx-1: # right
									if (markers[i+1,j] == 0) and (image[i+1,j] < value):
										markers[i+1,j] = marker
										print "right"
										moved = True
								if j > 0: # above
									if (markers[i,j-1] == 0) and (image[i,j-1] < value):
										markers[i,j-1] = marker
										print "above"
										moved = True
								if j < Ny -1: # below
									if (markers[i,j+1] == 0) and (image[i,j+1] < value):
										markers[i,j+1] = marker
										print "below"
										moved = True
			if 1:
				import mab.watershed as ws
				markers1 = ws.markmax(df, 1, 8)
				ws.watershed(df, markers1)
				markers2 = ws.markmax(df, 7, 8)
				ws.watershed(df, markers2)
				#iteration += 1
				#if iteration > maxiters:
				#	break
			if debug:
				indexedimage(markers*1.0)
			#df = 255-df
			#print "markers"
			#print markers
			#ws = ndimage.watershed_ift(df.astype(uint8), markers)
			#ws = ndimage.watershed_ift(df.astype(uint8), ws)
			#ws[df < 2] = -1
			if debug:
				select(2,0)
				indexedimage(markers*1.0)
				#print ws
			
				draw()
				sys.exit(0)
		if 0:
			import cv,cv2
			import Image
			df = (df-df.min())/(df.max()-df.min()) * 255
			pilImage = Image.fromarray(df)
			#pilImage.show()
			#graymat = cv.fromarray(df)
			#gray = graymat.to_pil_image()
			#print gray
			#import pyopencv
			#mat = cv2.from_pil_image(image)
			df = df.astype(float32)
			#h,w = df.shape
			#print w, h
			#df2 = cv.CreateMat(h, w, cv.CV_32FC3)
			#df0 = cv.fromarray(df)
			#color = cv.CvtColor(df0, df2, cv.CV_GRAY2BGR)
			#gray = cv2.cvtColor(df)
			#print gray
			#print df0
			df = df.astype(uint8)
			print df
			markers = (df * 0).astype(int32)
			#ret,thresh = cv2.threshold(df,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
			dfcolor = cv2.cvtColor(df,cv2.COLOR_GRAY2BGR)
			ret,thresh = cv2.threshold(df,5,255,cv2.THRESH_BINARY)#+cv2.THRESH_OTSU)
			
			print thresh
			print "*"*70
			fg = cv2.erode(thresh,None,iterations = 1)
			indexedimage(df*1.)
			print fg
			print "="*70
			bgt = cv2.dilate(thresh,None,iterations = 1)
			select(1,0)
			indexedimage(bgt*1.)
			print bgt
			print "="*70
			ret,bg = cv2.threshold(bgt,1,128,1)
			select(2,0)
			indexedimage(bg*1.)
			print bg
			print "*"*70
			marker = cv2.add(fg,bg)
			select(0,1)
			indexedimage(marker*1.)
			print marker
			marker32 = numpy.int32(marker)
			cv2.watershed(dfcolor,marker32)
			m = cv2.convertScaleAbs(marker32)
			select(1,1)
			indexedimage(m*1.)
			print "#" * 10, "m"
			print marker32
			print m
			
			ret,thresh = cv2.threshold(m,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
			res = cv2.bitwise_and(dfcolor,dfcolor,mask = thresh)
			print res.shape
			res = res * 1.0
			#dsa
			print res[:,:,0]
			select(2,1)
			indexedimage(res[:,:,0])
			draw()
			#print res
			
			if 0:
				cv2.watershed(dfcolor, markers)
				print markers
				sys.exit(0)
				print ret
				print thresh
				cv2.imshow("test", ret)
				cv2.waitKey(0)
				#sys.exit(0)
		if 0:
			df = (df-df.min())/(df.max()-df.min())
			#im = df.T * 255 
			if 1:
				import pymf
				import pdb
				N = 3
				colormap = "whiteblack"
				for y in range(5):
					nmf_mdl = pymf.NMF(df, num_bases=N, niter=5500)
					nmf_mdl.initialization()
					nmf_mdl.factorize()
					#pdb.set_trace()
					#box()
					print nmf_mdl.W.shape
					print nmf_mdl.H.shape
					if 0:
						select(0, 0)
						indexedimage(df)
						select(1,0)
						indexedimage(rec)
					
					#mozaic(3,3,box)
					#for i in range(3):
					for x in range(3):
						select(x, y)
						index = x# + y * 3
						if index < N:
							rec = dot(nmf_mdl.W[:,index:index+1], nmf_mdl.H[index:index+1,:])
							indexedimage(rec, colormap=colormap)
				rec = dot(nmf_mdl.W[:,:], nmf_mdl.H[:,:])
				select(2,5)
				indexedimage(rec, colormap=colormap)
				select(1,5)
				indexedimage(df, colormap=colormap)
				draw()


				
				
				sys.exit(0)
		if 0:
			import PIL
			import pdb
			im1 = PIL.Image.frombuffer("F", df.shape[::-1], (df*255).astype(float32))
			#im = PIL.Image.merge("RGB", (im1, im1, im1))
			im = im1.convert("RGB")
			im.save("df-teresa.png")
			#pdb.set_trace()
			sys.exit(0)
			
		global model
			
		if 0:
			dfL = df.sum(axis=1)
			dfL /= sum(dfL)
			graph(dfL)
			import emcee
			if 1:
				Ng = 2
				Ls = arange(len(dfL)) * 1.
				def lnprob(x, ivar=None):
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

				#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
				#sampler.run_mcmc(p0, 10000)
				#u = [1, 1, 1, 1, 1.2, 1, 1, 3, 1, 1, 5, 1]
				u =  [1, 1, 1, 1, 5, 1]
				bounds = None
				x = scipy.optimize.fmin_l_bfgs_b(lnprob, u, None, bounds=bounds, approx_grad=True, iprint=1,factr=1e-2,maxfun=200000)[0]
				print u
				
				
				#sdsa
			print model
			print dfL
			graph(model, color="red")
			
			for g in range(Ng):
				w, mu, sigma = x[g*3:g*3+3]
				w = (arctan(w)*2/pi+1)*50 # between 0 and 5
				sigma = exp(sigma)
				modelpart = w * kaplot.gaussian(Ls, mu, sigma)
				graph(modelpart, color="blue")
			
			draw()
			
			sys.exit(0)
			
		print df.shape
		ls = arange(nL)
		es = arange(nE)
		#ls, es = numpy.meshgrid(ls, es)
		es, ls = numpy.meshgrid(es, ls)
		df = df / sum(df)
		dforg = 1 * df
		self.solution.storage_2d_m0.load()
		self.solution.storage_2d_m2.load()
		self.solution.storage_3d.load()
		
		def do(df, **kwargs):
			
			if 0:
				import emcee
				
				def mvg(x, y, mux, muy, sigmax, sigmay, rho):
					u = (x-mux)**2/sigmax**2 + (y-muy)**2/sigmay**2 - 2 * rho*(x-mux)*(y-muy)/(sigmax*sigmay)
					return 1./(2*pi*sigmax*sigmay*sqrt(1-rho**2))*exp(-1./(2*(1-rho**2))*u)
				def lnprob(x, ivar):
					fom = 0
					mux, muy, sigmax, sigmay, rho = x
					rho = arctan(rho)*2/pi
					sigmax = exp(sigmax)
					sigmay = exp(sigmay)
					print "%10.4f %10.4f %10.4f %10.4f %10.4f" %(mux, muy, sigmax, sigmay, rho)
					for e in range(nE):
						for l in range(nL):
							fom += df[l,e] * log(mvg(e, l, mux, muy, sigmax, sigmay, rho))
					print fom
					if isnan(fom):
						sys.exit(0)
					return fom
					#print x
					#print ivar
					#dsa
					#print x
					#return -0.5
					
				n = 1
				ndim, nwalkers = n*5, 10
				ivar = 1. / numpy.random.rand(ndim)
				p0 = [numpy.random.rand(ndim) for i in range(nwalkers)]

				sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
				sampler.run_mcmc(p0, 1000)
				
				exit(0)
			
			
			#print ls.shape
			#print ls>0
			#df[(ls+es*0.5)>7] = 0
			#df[(ls+es*0.5)>7] = 0
			#indexedimage(df, colormap="whiteblack")
			
			orbitweights = df.reshape(-1)
			#orbitweights /= sum(orbitweights)
			#moments_m0 = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
			moments_m0 = tensordot(orbitweights, self.solution.storage_3d.moments3d, axes=[(0,), (0,)])
			moments_m2 = tensordot(orbitweights, self.solution.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
			moments_m0_projected = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
			cumu_mass = cumsum(moments_m0[0])/sum(moments_m0[0])
			cumu_mass_projected =  cumsum(moments_m0_projected[0])/sum(moments_m0_projected[0])
			
			rcenters_kpc_projected = Rs = self.solution.storage_2d_m0.projection.gridR.aperture_rcenters_kpc
			rcenters_kpc_3d = Rs = self.solution.storage_3d.rs
			
			dr = rcenters_kpc_projected[1] - rcenters_kpc_projected[0]
			density2d = moments_m0_projected[0]/sum(moments_m0_projected[0])/(2*pi*rcenters_kpc_projected*dr)
			density2d = moments_m0_projected[0]/(2*pi*rcenters_kpc_projected*dr)
			
			
			
			
			Rs = self.solution.storage_3d.x
			index = argmin(abs(cumu_mass-0.5))
			Rhalf = Rs[index]
			print "r_1/2 = %s" % Rs[index], "index =", index
			rcenters_kpc = Rs = self.solution.storage_2d_m2.projection.gridR.aperture_rcenters_kpc
			self.solution.binned_data_m2.load()
			sigma = self.solution.binned_data_m2.moments[2]**0.5
			kappa = self.solution.binned_data_m2.moments[4]/self.solution.binned_data_m2.moments[2]**2
			errors = sqrt(self.solution.binned_data_m2.e_moments[2]**2/(2*sigma)**2)
			kappa_errors = sqrt(self.solution.binned_data_m2.e_moments[2]**2/(2*sigma)**2)

			symbolsize = "6pt"
			fontSize = "6pt"

			#print moments_m0.shape
			mask = moments_m2[0] > 0
			if anisotropy:
				varvr = moments_m0[self.solution.storage_3d.v20index]
				varvphi = moments_m0[self.solution.storage_3d.v02index]/2
				#varvtheta = moments[6]
				mask = moments_m0[0] > 0
				betas = varvphi * 0 
				betas[mask] = 1 - (varvphi[mask])/(varvr[mask])
				graph(rcenters_kpc_3d, betas, **kwargs)
				hline(0)
				ylim(-3, 1)
				xlim(0, 1.25)
			elif plotdf:
				pass
			elif light:
				if 1:
					try:
						kinematic_bootstrap_band = scope["kinematic_bootstrap_band"]
					except:
						kinematic_bootstrap_band = None
					if kinematic_bootstrap_band:
						path = scope["model_id_path"]
						for i in range(150):
							filename = os.path.join(path, "massr-" + kwargs["color"] +"-" +str(i) + ".txt.npy")
							y = numpy.load(filename)
							mask = ~isnan(y)
							#print y
							#print mask
							if any(~isnan(y)):
								#graph(rcenters_kpc_3d[mask], y[mask], alpha=0.2, **kwargs)
								kwargs2 = dict(kwargs)
								del kwargs2["color"]
								color = kaplot.utils.getColor(kwargs["color"])
								color = color.blend(kaplot.utils.getColor("white"), 0.4)
								graph(rcenters_kpc_3d[mask], y[mask], color=color, **kwargs2)
					if seed is not None:
						if 1:
							path = scope["model_id_path"]
							filename = os.path.join(path, "massr-" + kwargs["color"] +"-" +str(seed) + ".txt")
							print "#" * 10, seed, filename
							numpy.save(filename, cumu_mass)
						#with open(filename, "a") as f:
						#print >>f, "%d, %f, %f, %f" % (seed, Rred, Rblue, Rred/Rblue)
						#logger.info("write radii to %s" % filename)
				if not plot_density:
					graph(rcenters_kpc_3d, cumu_mass, linewidth="2pt", **kwargs)
				print "r_1/2", scope["rhalf"], "r3", scope["r3"]
				#import pdb
				#pdb.set_trace()
				#graph(rcenters_kpc_projected, cumu_mass_projected, linewidth="2pt", **kwargs)
				light_model = scope["light_model"]
				#print dir(light_model)
				#print light_model.light_profile.b
				#print "!!!!!!!!!!df normalization", sum(df)
				# total mass of the profile should de equal to mass of df
				fit = mab.gd.potential.Plummer(sum(df), 1.0)
				#fit = mab.gd.potential.Plummer(1.0, 1.0)
				def f(params):
					b = params[0]
					fit.b = b
					#print b
					skip = 5
					#return sum((log10(fit.densityR(rcenters_kpc_projected[skip])) - log10(density2d[skip]))**2)
					return sum(((fit.densityR(rcenters_kpc_projected[skip])) - (density2d[skip]))**2)
				b, = scipy.optimize.fmin(f, 1.)
				fit.b = b
				print "b =", b
				
				if 0:
					fitexp = mab.gd.potential.ProjectedExponential(sum(df), 1.0)
					def f(params):
						#b = params[0]
						fitexp.scale = abs(params[0])
						fitexp.rho0 = 1.
						#sprint self.scale
						Mcurrent = fitexp.enclosed_mass(inf)
						fitexp.rho0 = sum(df)/Mcurrent
						
						#Mtot = fitexp.enclosed_mass(inf)
						#fitexp.rho0 = fitexp.rho0/Mtot
						
						#fitexp.r = b
						#print b
						skip = 5
						#return sum((log10(fit.densityR(rcenters_kpc_projected[skip])) - log10(density2d[skip]))**2)
						return sum(((fitexp.densityR(rcenters_kpc_projected[skip])) - (density2d[skip]))**2)
					scale, = scipy.optimize.fmin(f, 1.)
					fitexp.scale = abs(scale)
					fitexp.rho0 = 1.
					Mtot = fitexp.enclosed_mass(inf)
					fitexp.rho0 = sum(df)/Mtot
					print "b =", b
				
				#dsa
				#color = kwargs["
				kwargs2 = dict(kwargs)
				color = kaplot.utils.getColor(kwargs["color"])
				kwargs2["color"] = color.blend(kaplot.utils.getColor("white"), 0.2)
				if plot_density:
					#graph(rcenters_kpc_projected, fit.densityR(rcenters_kpc_projected), linewidth="3pt", **kwargs2)
					#graph(rcenters_kpc_projected, density2d, linewidth="1pt", **kwargs)
					graph(log10(rcenters_kpc_projected), log10(fit.densityR(rcenters_kpc_projected)), linewidth="3pt", **kwargs2)
					#kwargs2["color"] = color.blend(kaplot.utils.getColor("white"), 0.4)
					#graph(log10(rcenters_kpc_projected), log10(fitexp.densityR(rcenters_kpc_projected)), linewidth="1.5pt", **kwargs2)
					y = log10(density2d)
					mask = density2d > 0
					graph(log10(rcenters_kpc_projected[mask]), y[mask], linewidth="1pt", **kwargs)
				else:
					hline(0.5)
				
			elif vdisp:
				scatter(rcenters_kpc, sigma, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, sigma, yerr=errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				import pdb
				#pdb.set_trace()
				graph(Rs[mask], ((moments_m2[2]/moments_m2[0])**0.5)[mask], **kwargs)
				ylim(0, 15)
				xlim(0, 1.25)
			else:
				select(1,1)
				graph(cumu_mass, **kwargs)
				hline(0.5)
				select(1, 0)
				scatter(rcenters_kpc, sigma, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, sigma, yerr=errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				mask = moments_m2[0] > 0
				graph(Rs[mask], ((moments_m2[2]/moments_m2[0])**0.5)[mask], **kwargs)
				ylim(0, 15)
				xlim(0, 1.5)
				select(2,0)
				graph(Rs[mask], ((moments_m2[4]/moments_m2[0])/(moments_m2[2]/moments_m2[0])**2)[mask], **kwargs)
				scatter(rcenters_kpc, kappa, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, kappa, yerr=kappa_errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				xlim(0, 1.5)
				ylim(1, 5)
				
				select(2,1)
				
				rdeg = rcenters_kpc_projected/1.4
				mask = moments_m0_projected[0] > 0
				graph(log10(rdeg[mask]), log10(moments_m0_projected[0]/rdeg)[mask], **kwargs)
				#graph(rcenters_kpc_projected, (moments_m0_projected[0]/rcenters_kpc_projected), **kwargs)
				#xlim(0, 0.7)
				
				select(0, 1)
				
				
				graph(moments_m0[0], **kwargs)
			return Rhalf, moments_m0_projected[0]
		
		try:
			seed = scope["kinematic_bootstrap_seed"]
		except:
			seed = None
		try:
			bootstrap_band = scope["kinematic_bootstrap_band"]
		except:
			bootstrap_band = False
		
		colormap = "whiteblack"
		#mask = (ls > 2) | (es < 5)
		weight = 1./(1.+exp((ls-1.5)*2.5))
		#weight = weight * 1./(1.+exp((es-10.2)*2.5))
		#weight = 1./(1.+exp(((ls+es)-10)*2.5))
		dfbetween = df * 1.
		dfbetween[~((markers1 == markers2))] = 0
		print sum(dfbetween)
		print ~((markers1 == markers2))
		#import pdb
		#pdb.set_trace()
		if 1:
			#mask = (ls > 2)
			mask1 = mask = markers1 == 1
			#mask = mask & (
			df[~mask] = 0
			#df[(markers1==markers2)] /= 2
			df[(markers1==markers2)] = 0
			weights = (df/df.max()) > 0.1
			df1 = df * 1.
			import pdb
			#pdb.set_trace()
			#df1[ls<4] += dfbetween[ls<4]
			df = df1
			#weights[(markers1==markers2)] /= 2
		else:
			df = df * weight
		print "red/MR component"
		Rred, mass_proj_red = do(df1, color="red", linestyle="dash")
		if not (vdisp or light or plotdf or anisotropy):
			select(0,2)
			#indexedimage(df, colormap=colormap)
			indexedimage(df, colormap="whiteblue")
			#contourfill(df/df.max(), 0.5, 1., color="blue")
		df = dforg * 1.
		#mask = (ls > 5)
		if 1:
			#df[~mask] = 0
			mask = markers2 == 1
			df[~mask] = 0
			#df[(markers1==markers2)] /= 2
			df[(markers1==markers2)] = 0
			df2 = df * 1.
			df2 = df2# + dfbetween
			df2[ls>2] += dfbetween[ls>2]
			df = df2
			weights = markers2 * 1.
			#weights[(markers1==markers2)] /= 2
		else:
			df = df * (1-weight)
		print "blue/MP component"
		Rblue, mass_proj_blue = do(df, color="royalblue", linestyle="dot")
		print "ratio", Rred/Rblue
		print "rhalf", Rred, Rblue
		
		print "MASSES: in MR, between, MP", sum(df1), sum(dfbetween), sum(df2)
		
		if 0:
			select(2,2)
			ratios = mass_proj_blue/mass_proj_red
			graph(ratios[mass_proj_red>0][:50])
			#select(2,1)
			import pdb
			#pdb.set_trace()
			rs = arange(1e-2, 1.5, 0.01)
			light_profile = scope["light_profile"]
			graph(log10(rs/1.4), log10(light_profile.densityR(rs, M=1.)*400/2/1.0))
			vline(log10(2.3e-2))
			vline(log10(7.5e-2))
			#import pdb
			#pdb.set_trace()
		#filename = scope["modeldir"]
		if seed is not None:
			path = scope["model_id_path"]
			filename = os.path.join(path, "light.txt")
			with open(filename, "a") as f:
				print >>f, "%d, %f, %f, %f" % (seed, Rred, Rblue, Rred/Rblue)
				logger.info("write radii to %s" % filename)
		#sys.exit(0)
		#print filename
		#dsa
		if not (vdisp or light or plotdf or anisotropy):
			select(1,2)
			indexedimage(df, colormap="whitered")
			#sel
			#contour(df/df.max(), [0.5], fill=True, color="red", alpha=1.)

				
			#contour(weights, [0.6], color="red")
		do(dforg, color="black")
		if not (vdisp or light or anisotropy):
			if not plotdf:
				select(0, 0)
			indexedimage(dforg, colormap=colormap)
			Ny, Nx = df.shape
			def maskcontour(mask, **kwargs):
				for x in range(Nx):
					for y in range(Ny):
						if mask[y,x] == 1:
							if (x > 0) and (mask[y,x-1] == 0):
								line(x-0.5, y-0.5, x-0.5, y+0.5, **kwargs)
							if (x < (Nx-1)) and (mask[y,x+1] == 0):
								line(x+0.5, y-0.5, x+0.5, y+0.5, **kwargs)
							if (y > 0) and (mask[y-1,x] == 0):
								line(x-0.5, y-0.5, x+0.5, y-0.5, **kwargs)
							if (y < (Ny-1)) and (mask[y+1,x] == 0):
								line(x-0.5, y+0.5, x+0.5, y+0.5, **kwargs)
			mask = (df1/df1.max()) > 0.05
			maskcontour(mask, color="red", linestyle="dash", alpha=0.5)
			mask = (df2/df2.max()) > 0.05
			maskcontour(mask, color="royalblue", linestyle="dot", alpha=0.5)
			
			#mask = (dfbetween/dfbetween.max()) > 0#0.05
			#maskcontour(mask, color="orange", linestyle="dot", alpha=0.5)
			
			
			
		if vdisp:
			labels("R (kpc)" if scope["dobottom"] else "", "&sigma;<sub>los</sub> (km/s)" if scope["doleft"] else "")
		if light:
			if plot_density:
				current.container.xaxis.interval = 1.0
				labels("log R (kpc)" if scope["dobottom"] else "", "log I(R)/I<sub>total</sub>" if scope["doleft"] else "") # (L<sub>&#x2609;</sub>/kpc<sup>2</sup>)")
				current.container.xaxis.subticks = 3
				current.container.yaxis.subticks = 3
				xlim(0, 1.25)
				ylim(-4.2, 0.8)
				xlim(-2, 0.5)
			else:
				labels("r (kpc)" if scope["dobottom"] else "", "M(&lt;r)/M<sub>total</sub>" if scope["doleft"] else "")
				current.container.xaxis.subticks = 4
				current.container.yaxis.subticks = 3
			#labels("R (kpc)" if scope["dobottom"] else "", "const * I(R) (L<sub>&#x2609;</sub>/kpc<sup>2</sup>)" if scope["doleft"] else "")
			r3 = scope["r3"]
			print "r3", r3
			#vline(r3, color="black")
			#vline(0.3, color="black")
			if not plot_density:
				try:
					label = "ratio=%.2f" % (Rred/Rblue)
					#label += "\n" + scope["label"]
				except:
					label = "ratio=%.2f" % (Rred/Rblue)
				text(label, 1, 0.2)
	
	
		#import pdb
		#pdb.set_trace()
			
		if plotdf:
			labels("energy index" if scope["dobottom"] else "", "angular momentum index" if scope["doleft"] else "")
			current.container.xaxis.subticks = 4
			current.container.yaxis.subticks = 1
		
		if anisotropy:
			labels("r (kpc)" if scope["dobottom"] else "", "vel. anisotropy" if scope["doleft"] else "")
			pass
				
			
			
			
		#xlim(-2, 1)
		
		#if not plotdf:
		#select(2,2)
		#indexedimage(weight, colormap=colormap)
		draw()
		
class CommandCalcScaleRadii(object):
	def run(self, args, opts, scope):
		path = scope["model_id_path"]
		filename = os.path.join(path, "light.txt")
		with open(filename, "r") as f:
			lines = f.readlines()
		values = [eval(k) for k in lines]
		ratios = [k[3] for k in values]
		m = median(ratios)
		s = std(ratios)
		ratios = array(ratios)
		mask = (ratios > (m-3*s)) & (ratios < (m+3*s))
		print mask
		print mask.shape, ratios.shape
		ratios = ratios[mask]
		print mean(ratios), std(ratios), median(ratios)
		
class CommandPlotScaleRadii(object):
	def run(self, args, opts, scope):
		path = scope["model_id_path"]
		filename = os.path.join(path, "light.txt")
		with open(filename, "r") as f:
			lines = f.readlines()
		values = [eval(k) for k in lines]
		ratios = [k[3] for k in values]
		
		m = median(ratios)
		s = std(ratios)
		ratios = array(ratios)
		mask = (ratios > (m-3*s)) & (ratios < (m+3*s))
		#print mask
		#print mask.shape, ratios.shape
		ratios = ratios[mask]
		print mean(ratios), std(ratios), median(ratios)
		
		
		box()
		histogram(ratios, datamin=0.6, datamax=1.2, binwidth=0.05)
		clearautolegend()
		vline(0.84, color="red")
		#print ratios, max(ratios)
		vline(mean(ratios), color="green")
		vline(median(ratios), color="blue")
		autolegend("best fit", "mean", "median")
		labels("ratio", "N")
		grow(top=1.1)
		draw()
		#print mean(ratios), std(ratios), median(ratios)
		
		
		
class PlotMassProfileDfsNMF(object):
	def __init__(self, solution):
		self.solution = solution
		
	def run(self, args, opts, scope):
		document(size="30cm,25cm")
		page(fontsize="18pt")
		mozaic(4,6,box)
		self.solution.load()
		self.solution.storage_2d_m0.load()
		self.solution.storage_2d_m2.load()
		self.solution.storage_3d.load()
		df = self.solution.orbitweights
		nE, nL = 20, 8
		df = df.reshape((nL, nE))
		if 1:
			df = (df-df.min())/(df.max()-df.min())
			#im = df.T * 255 
			if 1:
				import pymf
				import pdb
				N = 3
				colormap = "whiteblack"
				for y in range(5):
					class MyNMF(pymf.NMF):
						def frobenius_norm(self):
							""" Frobenius norm (||data - WH||) for a data matrix and a low rank
							approximation given by WH
							
							Returns:
								frobenius norm: F = ||data - WH||
							"""		
						
							#err = sqrt( sum((self.data[:,:]/dot(self.W, self.H))**2 ))
							err = sqrt( sum(log10(self.data[:,:])-log10(dot(self.W, self.H))**2 ))
							return err
						
					#nmf_mdl = MyNMF(df, num_bases=N, niter=500)
					nmf_mdl = pymf.NMF(df, num_bases=N, niter=5000)
					nmf_mdl.initialization()
					nmf_mdl.factorize()
					#pdb.set_trace()
					#box()
					print nmf_mdl.W.shape
					print nmf_mdl.H.shape
					if 0:
						select(0, 0)
						indexedimage(df)
						select(1,0)
						indexedimage(rec)
					
					#mozaic(3,3,box)
					#for i in range(3):
					recs = [dot(nmf_mdl.W[:,index:index+1], nmf_mdl.H[index:index+1,:]) for index in range(N)]
					Ls = [recs[k].sum(axis=1) for k in range(N)]
					avgL = [sum(Ls[k] * range(8)) / sum(Ls[k]) for k in range(N)]
					indices = argsort(avgL)
					for x in range(N):
						select(x, y)
						index = indices[x]#x# + y * 3
						if index < N:
							rec = dot(nmf_mdl.W[:,index:index+1], nmf_mdl.H[index:index+1,:])
							
							indexedimage(rec, colormap=colormap)
							print rec.min(), rec.max()
					rec = dot(nmf_mdl.W[:,:], nmf_mdl.H[:,:])
					select(N,y)
					indexedimage(rec, colormap=colormap)
					
					def calc_rhalf(df):
						orbitweights = df.reshape(-1)
						orbitweights /= sum(orbitweights)
						#moments_m0 = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
						moments_m0 = tensordot(orbitweights, self.solution.storage_3d.moments3d, axes=[(0,), (0,)])
						moments_m2 = tensordot(orbitweights, self.solution.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
						moments_m0_projected = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
						cumu_mass = cumsum(moments_m0[0])/sum(moments_m0[0])
						cumu_mass_projected =  cumsum(moments_m0_projected[0])/sum(moments_m0_projected[0])
						
						rcenters_kpc_projected = Rs = self.solution.storage_2d_m0.projection.gridR.aperture_rcenters_kpc
						rcenters_kpc_3d = Rs = self.solution.storage_3d.rs
						
						Rs = self.solution.storage_3d.x
						index = argmin(abs(cumu_mass-0.5))
						Rhalf = Rs[index]
						print "R_1/2 = %s" % Rhalf, "index =", index
						return Rhalf
					
					Rhalf_MR = calc_rhalf(recs[0])
					Rhalf_MP = calc_rhalf(sum(rec for rec in recs[1:]))
					ratio = Rhalf_MR/Rhalf_MP
					print ratio
					select(0, y)
					text("%.2f" % ratio, 10, 3, alighn="center,center")
					
					#df = rec
					
				select(1,5)
				indexedimage(df, colormap=colormap)
				draw()
				
				


				
				
				sys.exit(0)

			
		print df.shape
		ls = arange(nL)
		es = arange(nE)
		#ls, es = numpy.meshgrid(ls, es)
		es, ls = numpy.meshgrid(es, ls)
		
		dforg = 1 * df
		self.solution.storage_2d_m0.load()
		self.solution.storage_2d_m2.load()
		self.solution.storage_3d.load()
		
		def do(df, **kwargs):
			
			if 0:
				import emcee
				
				def mvg(x, y, mux, muy, sigmax, sigmay, rho):
					u = (x-mux)**2/sigmax**2 + (y-muy)**2/sigmay**2 - 2 * rho*(x-mux)*(y-muy)/(sigmax*sigmay)
					return 1./(2*pi*sigmax*sigmay*sqrt(1-rho**2))*exp(-1./(2*(1-rho**2))*u)
				def lnprob(x, ivar):
					fom = 0
					mux, muy, sigmax, sigmay, rho = x
					rho = arctan(rho)*2/pi
					sigmax = exp(sigmax)
					sigmay = exp(sigmay)
					print "%10.4f %10.4f %10.4f %10.4f %10.4f" %(mux, muy, sigmax, sigmay, rho)
					for e in range(nE):
						for l in range(nL):
							fom += df[l,e] * log(mvg(e, l, mux, muy, sigmax, sigmay, rho))
					print fom
					if isnan(fom):
						sys.exit(0)
					return fom
					#print x
					#print ivar
					#dsa
					#print x
					#return -0.5
					
				n = 1
				ndim, nwalkers = n*5, 10
				ivar = 1. / numpy.random.rand(ndim)
				p0 = [numpy.random.rand(ndim) for i in range(nwalkers)]

				sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
				sampler.run_mcmc(p0, 1000)
				
				exit(0)
			
			
			#print ls.shape
			#print ls>0
			#df[(ls+es*0.5)>7] = 0
			#df[(ls+es*0.5)>7] = 0
			#indexedimage(df, colormap="whiteblack")
			
			orbitweights = df.reshape(-1)
			orbitweights /= sum(orbitweights)
			#moments_m0 = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
			moments_m0 = tensordot(orbitweights, self.solution.storage_3d.moments3d, axes=[(0,), (0,)])
			moments_m2 = tensordot(orbitweights, self.solution.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
			moments_m0_projected = tensordot(orbitweights, self.solution.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
			cumu_mass = cumsum(moments_m0[0])/sum(moments_m0[0])
			cumu_mass_projected =  cumsum(moments_m0_projected[0])/sum(moments_m0_projected[0])
			
			rcenters_kpc_projected = Rs = self.solution.storage_2d_m0.projection.gridR.aperture_rcenters_kpc
			rcenters_kpc_3d = Rs = self.solution.storage_3d.rs
			
			Rs = self.solution.storage_3d.x
			index = argmin(abs(cumu_mass-0.5))
			Rhalf = Rs[index]
			print "R_1/2 = %s" % Rs[index], "index =", index
			rcenters_kpc = Rs = self.solution.storage_2d_m2.projection.gridR.aperture_rcenters_kpc
			self.solution.binned_data_m2.load()
			sigma = self.solution.binned_data_m2.moments[2]**0.5
			kappa = self.solution.binned_data_m2.moments[4]/self.solution.binned_data_m2.moments[2]**2
			errors = sqrt(self.solution.binned_data_m2.e_moments[2]**2/(2*sigma)**2)
			kappa_errors = sqrt(self.solution.binned_data_m2.e_moments[2]**2/(2*sigma)**2)

			symbolsize = "6pt"
			fontSize = "6pt"

			#print moments_m0.shape
			if light:
				graph(rcenters_kpc_3d, cumu_mass, **kwargs)
				graph(rcenters_kpc_projected, cumu_mass_projected, linestyle="dash", **kwargs)
				hline(0.5)
			elif vdisp:
				scatter(rcenters_kpc, sigma, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, sigma, yerr=errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				graph(Rs, (moments_m2[2]/moments_m2[0])**0.5, **kwargs)
				ylim(0, 15)
				xlim(0, 1.5)
			else:
				select(1,1)
				graph(cumu_mass, **kwargs)
				hline(0.5)
				select(1, 0)
				scatter(rcenters_kpc, sigma, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, sigma, yerr=errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				graph(Rs, (moments_m2[2]/moments_m2[0])**0.5, **kwargs)
				ylim(0, 15)
				xlim(0, 1.5)
				select(2,0)
				graph(Rs, ((moments_m2[4]/moments_m2[0])/(moments_m2[2]/moments_m2[0])**2)**1.0, **kwargs)
				scatter(rcenters_kpc, kappa, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
				errorbars(rcenters_kpc, kappa, yerr=kappa_errors, symbolsize=symbolsize, fontSize=fontSize, size="1mm", **kwargs)
				xlim(0, 1.5)
				ylim(1, 5)
				
				select(0, 1)
				
				
				graph(moments_m0[0], **kwargs)
			return Rhalf
		
		colormap = "whiteblack"
		#mask = (ls > 2) | (es < 5)
		weight = 1./(1.+exp((ls-1.5)*2.5))
		#weight = weight * 1./(1.+exp((es-10.2)*2.5))
		#weight = 1./(1.+exp(((ls+es)-10)*2.5))
		if 1:
			#mask = (ls > 2)
			mask = markers == 1
			df[~mask] = 0
		else:
			df = df * weight
		Rred = do(df, color="red")
		if not (vdisp or light):
			select(0,2)
			indexedimage(df, colormap=colormap)
		df = dforg * 1.
		#mask = (ls > 5)
		if 1:
			#df[~mask] = 0
			mask = markers == 2
			df[~mask] = 0
		else:
			df = df * (1-weight)
		Rblue = do(df, color="blue")
		print "ratio", Rred/Rblue
		if not (vdisp or light):
			select(1,2)
			indexedimage(df, colormap=colormap)
		do(dforg, color="black")
		if not (vdisp or light):
			select(0, 0)
			indexedimage(dforg, colormap=colormap)
		if vdisp:
			labels("R (kpc)", "&sigma;<sub>los</sub> (km/s)")
		if light:
			labels("R (kpc)", "M(&lt;R)/M<sub>total</sub>")
			
		select(2,2)
		indexedimage(weight, colormap=colormap)
		draw()
		
		

class PlotBindingEnergyJeans(object):
	def __init__(self, jeans):
		self.jeans = jeans
		
	def run(self, args, opts, scope):
		mozaic(2,2,box)
		#box()
		#rs = arange(0.001, 2, 0.1/4)
		logrs = arange(-2, 2, 0.1)
		rs = 10**logrs
		
		select(0, 0)
		labels("r / kpc", "W")
		dm_potential = scope["dm_profile"]
		#for potential, color in zip(self.potentials, nicecolors):
			#pot = potential.potentialr(rs)
		bind = dm_potential.enclosed_binding_energy(rs)
		print dm_potential.rs
		def K(s):
			def f(x):
				return 2*pi* x**2*dm_potential.densityr(x)*(3.-2*self.jeans.beta0)*self.jeans.sigmar(x)**2
			I, err = scipy.integrate.quad(f, 0, s)
			return I
		K = vectorize(K)
		#print K(1.)
		#dsa
		print rs
		print bind
		color = "black"
		kin = K(rs)
		graph(log10(rs), log10(-bind), color=color)
		graph(log10(rs), log10(kin), color="red")
		graph(log10(rs), log10(2*kin), color="red", linestyle="dash")
		select(0, 1)
		graph(log10(rs), log10(2*kin/abs(bind)), color="red")
		labels("log rs", "2K/-W")
		select(1, 1)
		graph(log10(rs), log10(2*kin/abs(bind)), color="red")
		print 2*kin/abs(bind)
		draw()
		#print args
		
		#jeans.


class PlotBindingEnergy(object):
	def __init__(self, potentials):
		self.potentials = potentials
		
		
	def run(self, args, opts, scope):
		
		mozaic(3,3,box)
		rs = arange(0.001, 2, 0.1/4)
		logrs = arange(-2, 2, 0.1)
		
		for potential, color in zip(self.potentials, nicecolors):
			print potential
			print potential.rs
			print potential.densityr(0), potential.densityr(potential.rs)
			print potential.enclosed_binding_energy(1.)
			
		to_erg = 1.9891e43
		#dsa
		select(0, 0)
		labels("r / kpc", "W")
		for potential, color in zip(self.potentials, nicecolors):
			#pot = potential.potentialr(rs)
			bind = potential.enclosed_binding_energy(rs) * to_erg
			graph(rs, log10(-bind), color=color)
		select(0, 1)
		labels("log r / kpc", "W")
		for potential, color in zip(self.potentials, nicecolors):
			#pot = potential.potentialr(rs)
			bind = potential.enclosed_binding_energy(10**logrs) * to_erg
			graph(logrs, log10(-bind), color=color)
			
		select(1, 0)
		for potential, color in zip(self.potentials, nicecolors):
			mass = potential.enclosed_mass(rs)
			graph((rs), log10(mass), color=color)
		select(1, 1)
		for potential, color in zip(self.potentials, nicecolors):
			mass = potential.enclosed_mass(10**logrs)
			graph(logrs, log10(mass), color=color)
		labels("log r / kpc", "delta W")
		
		select(2, 0)
		for potential, color in zip(self.potentials, nicecolors):
			density = potential.densityr(rs)
			graph((rs), log10(density), color=color)
		select(2, 1)
		for potential, color in zip(self.potentials, nicecolors):
			density = potential.densityr(10**logrs)
			graph(logrs, log10(density), color=color)
		select(0, 2)
		ref = self.potentials[0]
		bind_ref = ref.enclosed_binding_energy(10**logrs) * to_erg
		for potential, color in zip(self.potentials, nicecolors)[1:]:
			#pot = potential.potentialr(rs)
			bind = potential.enclosed_binding_energy(10**logrs) * to_erg
			graph(logrs, log10(-bind) - log10(-bind_ref), color=color)
		ylim(-1.5, 1.5)
		labels("log r / kpc", "delta W")
		
			
		draw()
			
		
		
				
				
class HaloEnergy(object):
	def __init__(self, filename, potential):
		self.filename = filename
		self.potential = potential
		
		
	def run(self, args, opts, scope):
		halo = mab.asciifile.readsimple(self.filename)
		
		rmin, rmax = min(halo.r), max(halo.r)
		rmax = 10**3
		rs = 10**arange(log10(rmin), log10(rmax), 0.05)

		#print halo.r
		#box()
		mozaic(3,2,box)
		select(0, 1)
		graph(log10(halo.r), log10(-halo.pot))
		graph(log10(rs), log10(-self.potential.potentialr(rs)), color="red")
		
		select(0, 0)
		graph(log10(halo.r), (halo.pot))
		graph(log10(rs), (self.potential.potentialr(rs)), color="red")
		#graph(log10(halo.r), (-self.potential.potentialr(rs)), color="red")
		
		select(1, 0)
		#graph(log10(halo.r), log10(-cumsum(halo.pot)))
		graph(log10(rs), log10(-self.potential.enclosed_binding_energy(rs)), color="red")
		#graph(log10(halo.r), log10(-cumsum(halo.pot)*2.87e5/0.73))
		G = self.potential.G
		#bla = (arange(len(halo.r))+1)/halo.r
		bla = (arange(len(halo.r))+1)/halo.r
		print bla
		graph(log10(halo.r), log10(cumsum(bla)*G*(2.87e5/0.73)**2))
		grow(bottom=0.5)
		
		
		select(2,0)
		#graph(log10(halo.r), log10(cumsum(bla)*G*(2.87e5/0.73)**2) - log10(-self.potential.enclosed_binding_energy(halo.r)))
		bla = ones(len(halo.r)) #arange(len(halo.r))+1
		graph(log10(halo.r), log10(cumsum(bla)*(2.87e5/0.73)) )
		graph(log10(rs), log10(self.potential.enclosed_mass(rs)), color="red")
		labels("log rs", "M")
		
		select(2,1)
		graph(log10(halo.r), log10(cumsum(bla)*(2.87e5/0.73))  - log10(self.potential.enclosed_mass(halo.r)), color="red")
		labels("log rs", "M")
		
		#graph(log10(rs), log10(-self.potential.enclosed_binding_energy(rs)), color="red")
		#scatter(log10(halo.r), log10(-cumsum(halo.pot)) - log10(-self.potential.enclosed_binding_energy(halo.r)), color="red")
		
		if 0:
			p = mab.gd.potential.Plummer(1., 1.)
			N = 1000000
			rp = p.sample_r(N=N)
			select(2,0)
			indices = argsort(rp)
			#rp = sort
			rp = rp[indices]
			pot = p.potentialr(rp)
			scatter(log10(rp), log10(-cumsum(pot)))
			
			rs = 10**arange(log10(rmin), log10(rmax), 0.05)
			scatter(log10(rs), log10(-p.enclosed_binding_energy(rs)*N), color="red")
			
		
		draw()
		
		
class Halos(object):
	def __init__(self, halos):
		self.halos = halos
		
	def run(self, args, opts, scope):
		halos = self.halos.load()
		
		#box()
		mozaic(2,2,box)
		select(0, 0)
		scatter(log10(halos.rho0_hay), log10(halos.rt_hay))
		select(1, 0)
		scatter(log10(halos.rs_hay), log10(halos.rt_hay))
		select(0, 1)
		scatter(halos.Mv, log10(halos.rt_hay))
		select(1, 1)
		scatter(log10(halos.rho0_hay), log10(halos.rt_hay))
		draw()
		
		
class GuisMrMp(object):
	def __init__(self, data_mp, data_mr, binned_data_m2=None, Rmax=None):
		self.data_mp = data_mp
		self.data_mr = data_mr
		self.binned_data_m2 = binned_data_m2
		self.Rmax = Rmax

	def run(self, args, opts, scope):
		self.data_mp.load()
		self.data_mr.load()
		
		document(size="10cm,10.5cm")
		page(fontsize="17pt")
		box()
		
		symbolsize = "6pt"
		
		self.binned_data_m2.load()
		sigma = self.binned_data_m2.moments[2]**0.5
		print "sigma", sigma
		#kappa = self.binned_data_m2.moments[4]/self.binned_data_m2.moments[2]**2
		errors = sqrt(self.binned_data_m2.e_moments[2]**2/(2*sigma)**2)
		self.binned_data_m2.aperture.load()
		print dir(self.binned_data_m2.aperture)
		x = self.binned_data_m2.aperture.aperture_rcenters_kpc #/60**2
		print "Re", self.data_mp.data.Re
		x = self.data_mp.data.Re[5:]
		color = "cyan"
		if 0:
			scatter(x, sigma, color=color, symbolsize=symbolsize,  symbolName="squaresolid")
			errorbars(x, sigma, yerr=errors, color=color)
		current.container.xaxis.interval = 0.5
		
		y = self.data_mp.data.vdisp
		yerr = (self.data_mp.data.error_low+self.data_mp.data.error_high)/2
		
		if 0:
			scatter(self.data_mp.data.Re, y * 1, color="orange", symbolName="squaresolid", symbolsize=symbolsize)
			errorbars(self.data_mp.data.Re, y * 1, yerr=yerr*1, color="orange") #, symbolName="squaresolid", symbolsize=symbolsize

		print len(self.data_mp.data.error_low), len(y), len(yerr)
		print errors
		print yerr[5:]
		print y[5:]
		print sigma
		if 0:
			yerr[5:] = sqrt((yerr[5:]**2 + errors**2)/4)
			y[5:] = sqrt((y[5:]**2 + sigma**2)/2)
		else:
			y[5:] = sigma
			yerr[5:] = errors
		scatter(self.data_mp.data.Re, y, color="royalblue", symbolName="squaresolid", symbolsize=symbolsize)
		errorbars(self.data_mp.data.Re, y, yerr=yerr, color="royalblue") #, symbolName="squaresolid", symbolsize=symbolsize
		
		
		scatter(self.data_mr.data.Re, self.data_mr.data.vdisp, color="red", symbolName="squaresolid", symbolsize=symbolsize)
		errorbars(self.data_mr.data.Re, self.data_mr.data.vdisp, yerr=self.data_mr.data.error_low, color="red") #, symbolName="squaresolid", symbolsize=symbolsize
		#scatter(rcenters_kpc, kappa, symbolName="squaresolid", symbolsize=symbolsize, fontSize=fontSize, **kwargs)
		
		if 0:
			autolegend("addin", "original", "mean")
		
		if self.Rmax:
			xlim(0, self.Rmax)
		#xlim(0, 2)
		#labels("R<sub>e</sub> (kpc), R (kpc)", "&sigma;<sub>los</sub> (km/s)")
		labels("R<sub>e</sub> (kpc)", "&sigma;<sub>los</sub> (km/s)")
		if 1:
			current.container.drawOutsideRight = True #scope["doleft"]
			current.container.drawOutsideTop = True #scope["doleft"]
			current.container.drawOutsideLeft = True #scope["doleft"]
			current.container.drawOutsideBottom = True #scope["dobottom"]
			current.container.leftaxes[0].drawLabels = True
			current.container.bottomaxes[0].drawLabels = True
		import pdb
		#pdb.set_trace()
		ylim(0, 15)
		draw()
		
		


class Posterior(object):
	def __init__(self, parameterset):
		self.parameterset = parameterset
		self.parameterset.load()
		grid = self.parameterset.probability_grid
		self.p = (grid)
		scale = 6000/200.
		scale = 1.
		self.p = exp(self.p*scale)
		self.p[isnan(grid)] = 0
		self.smooth = self.p
		index = argmax(self.smooth.flat)
		index = unravel_index(index, self.smooth.shape)
		#scales = parametersetcollect.probability_range
		#print scales
		#print index
		parameter_range_list = self.parameterset.parameter_range_list
		
		self.max_L = [parameter_range_list[k].min + (parameter_range_list[k].max - parameter_range_list[k].min) * (index[k]+0.5)/(self.smooth.shape[k]-1.) for k in range(len(index))]

		print "max Likelihood", self.max_L
		#self.smooth = self.p #
		if 1:
			self.smooth = self.p * 1.
			#self.smooth = scipy.ndimage.gaussian_filter(self.p, 1.5)
			#self.smooth = scipy.ndimage.gaussian_filter(self.p, 0.5)
			self.smooth_param = 1.
			print "self.smooth_param =", self.smooth_param
			self.smooth = scipy.ndimage.gaussian_filter(self.p, self.smooth_param)
			#self.smooth = self.p
		self.smooth = self.smooth / self.smooth.flat[np.argmax(self.smooth.flat)] 
		self.logp = np.log(self.smooth)
		self.xlist = np.linspace(parameter_range_list[0].min, parameter_range_list[0].max, self.p.shape[0])
		self.ylist = np.linspace(parameter_range_list[1].min, parameter_range_list[1].max, self.p.shape[1])
		self.logf = scipy.interpolate.interp2d(self.xlist, self.ylist, self.logp.T)#, fill_value=-np.inf)
		print "max", self.logp.flat[np.argmax(self.smooth.flat)] 
		print np.argmax(self.smooth.flat)
		#import pdb
		#pdb.set_trace()
		print self.logL(self.max_L[0], self.max_L[1])
		#print self.logL(self.max_L[0]+0.01, self.max_L[1]+0.01)
		#dsa
			
	def makeresize(self, i, j):
		return (self.parameterset.parameter_range_list[i].min, self.parameterset.parameter_range_list[j].min), (self.parameterset.parameter_range_list[i].max, self.parameterset.parameter_range_list[j].max)

	def logL(self, x, y):
		return self.logf(x,y)
		
class TransformHaloGeneral(object):
	def __init__(self, param_alpha, param_beta, param_gamma, param_logrs, subhalo_filename, posterior):
		self.param_alpha = param_alpha
		self.param_beta = param_beta
		self.param_gamma = param_gamma
		self.param_logrs = param_logrs
		self.subhalo_filename = subhalo_filename
		self.posterior = posterior
		
	def run(self, argv, opts, scope):
		name = scope["object_name"]
		r3 = scope["r3"]
		#Mv = scope["Mv"]
		Mvmin = scope["Mvmin"]
		Mvmax = scope["Mvmax"]

		halos = np.genfromtxt(self.subhalo_filename, names=True, dtype=None)
		mask = (halos["Mv"]  > Mvmin) & (halos["Mv"]  < Mvmax)#  & (halos["rho0"] > 0)
		halos = halos[mask]
		fn = "sub-info-target-2-list-2.dat"
		
		import emcee
		
		
		
		p = self.posterior.p * 0
		i = 0
		for x in self.posterior.xlist:
			j = 0
			for y in self.posterior.ylist:
				p[i,j] = self.posterior.logL(x, y)
				j += 1
			i += 1
		kwargs = {}
		i, j = 0, 1
		color = "red"
		p = np.exp(p)
		#print p
		#print self.posterior.logp
		#print self.posterior.xlist
		#print self.posterior.ylist
		if 0:
			box()
			obj = probimage2d(p, i, j, resize=self.posterior.makeresize(i,j), drawcontourlines=True, fill=False, colormap=None, premultiply=True, color=color, sigmas=2, **kwargs) #(0, 1), (0,1))
			draw()

		rmax_name = "2r3"
		rmax = r3 * 2
		parameters = self.param_alpha, self.param_beta, self.param_gamma
		for halo_index in range(3,len(halos)):
			target_mass_rmax_name = ("m" +name +"(%s)" % rmax_name).replace("(", "").replace(")", "")
			target_dens_rmax_name = ("d" +name +"(%s)" % rmax_name).replace("(", "").replace(")", "")
			#print halos.dtype
			target_dens_rmax = halos[target_dens_rmax_name][halo_index]
			target_mass_rmax = halos[target_mass_rmax_name][halo_index]
			#print target_dens_rmax, target_mass_rmax
			def fitall(x):
				#print logrs
				logmalpha,logmbeta,logmgamma,logrs, logscaling = x
				alpha = -10**logmalpha
				beta = -10**logmbeta
				gamma = 10**logmgamma
				rs = 10**logrs
				scaling = 10**logscaling
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=1.)
				#dsa
				profile.rho0 = scaling * target_dens_rmax / profile.densityr(rmax)
				#print profile.enclosed_mass(20)
				#print profile.enclosed_mass2(20)
				#print "dens", target_dens_rmax, profile.densityr(rmax)
				mass_rmax = profile.enclosed_mass(rmax)
				#print mass_rmax, mass_rmax
				M3 = profile.enclosed_mass(r3)
				slope3 = profile.logslope(r3)
				print "opt: ", alpha,beta, gamma, rs, scaling, slope3, np.log10(M3), np.log10(mass_rmax)
				#return (np.log10(mass_rmax) - np.log10(target_mass_rmax))**2 -self.posterior.logL(M3, slope3)[0]
				return (np.log10(mass_rmax) - np.log10(target_mass_rmax))**2 + (self.posterior.max_L[0] - np.log10(M3))**2 + (self.posterior.max_L[1] - slope3)**2 + (gamma - 2.0)**2 * 1e-4 + (alpha/beta)*1e-3 # + (beta - -3)**2 * 1e-2 
			
			startparams = scipy.optimize.fmin(fitall, [np.log10(1),np.log10(3), np.log10(2), 0, 0], disp=True)
			print "max L", self.posterior.max_L
			print "best fit", startparams
			logmalpha,logmbeta,logmgamma,logrs, logscaling = startparams
			alpha = -10**logmalpha
			beta = -10**logmbeta
			gamma = 10**logmgamma
			bestrs = rs = 10**logrs
			print "alpha,beta,gamma", alpha, beta, gamma
			scaling = 10**logscaling
			profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=-abs(alpha), beta=-abs(beta), gamma=abs(gamma), rs=rs, rho0=1.)
			profile.rho0 = scaling* target_dens_rmax / profile.densityr(rmax)
			mass_rmax = profile.enclosed_mass(rmax)
			diff = np.log10(mass_rmax) - np.log10(target_mass_rmax)
			#3if abs(diff) < 0.1:
			print "best fit mass(rmax)", np.log10(mass_rmax) , "target", np.log10(target_mass_rmax)
			M3 = profile.enclosed_mass(r3)
			slope3 = profile.logslope(r3)
			print "best fit m3", np.log10(M3), "target", self.posterior.max_L[0]
			print "best fit slope3", slope3 , "target", self.posterior.max_L[1]
			print -self.posterior.logL(np.log10(M3), slope3)[0]
			print -self.posterior.logL(np.log10(M3)+0.01, slope3)[0]
			print -self.posterior.logL(np.log10(M3), slope3+0.01)[0]
			print -self.posterior.logL(self.posterior.max_L[0], self.posterior.max_L[1])[0]
			print self.posterior.max_L[0], self.posterior.max_L[1]
			print 
			print self.posterior.makeresize(i,j)
			#sys.exit(0)
			
			def logL(x):
				logmalpha, logmbeta, logmgamma = x
				alpha = -10**logmalpha
				beta = -10**logmbeta
				gamma = 10**logmgamma
				def f(logrs):
					#print logrs
					rs = 10**logrs[0]
					#print alpha, beta, gamma, rs
					profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=1.)
					profile.rho0 = scaling * target_dens_rmax / profile.densityr(rmax)
					#print "dens", target_dens_rmax, profile.densityr(rmax)
					mass_rmax = profile.enclosed_mass(rmax)
					#print mass_rmax, target_mass_rmax, rmax, scaling
					return (np.log10(mass_rmax) - np.log10(target_mass_rmax))**2
				logrs, = scipy.optimize.fmin(f, np.log10(bestrs), disp=False)
				rs = 10**logrs
				#print "best rs", rs
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=1.)
				profile.rho0 = scaling * target_dens_rmax / profile.densityr(rmax)
				mass_rmax = profile.enclosed_mass(rmax)
				diff = np.log10(mass_rmax) - np.log10(target_mass_rmax)
				#3if abs(diff) < 0.1:
				#print "target mass", diff, "has", np.log10(mass_rmax) , "target", np.log10(target_mass_rmax)
				
				M3 = profile.enclosed_mass(r3)
				slope3 = profile.logslope(r3)
				logp = self.posterior.logL(np.log10(M3), slope3)[0] - (np.log10(mass_rmax) - np.log10(target_mass_rmax))**2  - (gamma - 2.0)**2 * 1e-2 
				print ">", alpha, beta, gamma, np.log10(M3), M3, slope3, mass_rmax, r3, logp
				#profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=rho0)
				to_erg = 1.9891e43
				try:
					E0 = profile.enclosed_binding_energy(rmax) * to_erg
				except:
					print "issue with", alpha, beta, gamma, rs
					raise
							
				if np.isnan(logp):
					logp = -np.inf
				return logp, [rs, profile.rho0, scaling, rmax, M3, slope3, E0, mass_rmax, target_mass_rmax]
			
			if 0:
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=-1.999, beta=-4, gamma=1, rs=1e-3, rho0=1.)
				profile.rho0 = target_dens_rmax / profile.densityr(rmax)
				mass_rmax_1 = profile.enclosed_mass(rmax)

				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=0, beta=0, gamma=1, rs=1e9, rho0=1.)
				profile.rho0 = target_dens_rmax / profile.densityr(rmax)
				mass_rmax_2 = profile.enclosed_mass(rmax)
				print "mass min/max", np.log10(mass_rmax_2), np.log10(mass_rmax_1)
				
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=0, beta=0, gamma=1, rs=1e9, rho0=1.)
				profile.rho0 = target_dens_rmax / profile.densityr(rmax)
				mass_rmax_2 = profile.enclosed_mass(rmax)
				print "mass min/max", np.log10(mass_rmax_2), np.log10(mass_rmax_1)

			print "logL", logL((alpha, beta, gamma))
			#continue
			#sys.exit(0)
			ndim, nwalkers = 3, 10
			#f
			#ivar = 1. / np.random.rand(ndim)
			#p0 = [np.random.rand(ndim) for i in range(nwalkers)]
			startpoints = []
			#alpha, beta, gamma = startparams[:3]
			#alpha = -abs(alpha)
			#beta = -abs(beta)
			#gamma = abs(gamma)
			for i in range(nwalkers):
				#startpoint = [p.min + np.random.random() * (p.max - p.min) for p in parameters]
				#startpoint = np.array(bestfit[:3]) + np.array([0.05*np.random.random() * (p.max - p.min) for p in parameters])
				startpoint = np.array(startparams[:3]) + np.array([0.001*(2*np.random.random()-1.) for p in parameters])
				startpoints.append(startpoint)

			sampler = emcee.EnsembleSampler(nwalkers, ndim, logL)
			pos, prob, state  = sampler.run_mcmc(startpoints, 100)[:3]
			sampler.reset()
			sampler.run_mcmc(pos, 1000)	
			fname = name + "_chain.npy"
			np.save(fname, sampler.chain)
			fname = name + "_lnprob.npy"
			np.save(fname, sampler.lnprobability)
			fname = name + "_blobs.npy"
			np.save(fname, sampler.blobs)
			print "wrote chain"
			sys.exit(0)
							
class PlotTransformHaloChain(object):
	def __init__(self, param_alpha, param_beta, param_gamma, param_logrs, subhalo_filename, posterior):
		self.param_alpha = param_alpha
		self.param_beta = param_beta
		self.param_gamma = param_gamma
		self.param_logrs = param_logrs
		self.subhalo_filename = subhalo_filename
		self.posterior = posterior
		
	def run(self, argv, opts, scope):
		Mv = scope["Mv"]
		
		name = scope["object_name"]
		r3 = scope["r3"]
		mozaic(3,3, box)
		fname = name + "_chain.npy"
		chain = np.load(fname)
		fname = name + "_lnprob.npy"
		lnprob = np.load(fname)
		fname = name + "_blobs.npy"
		blobs = np.load(fname)
		blobs = blobs.reshape(-1, blobs.shape[-1])
		#print blobs.shape, blobs[0]
		rs, rho0, scaling, rmax, M3, slope3, Es, mass_rmax, target_mass_rmax = blobs.T
		logmalpha = (chain[:,:,0].reshape(-1))
		logmbeta = (chain[:,:,1].reshape(-1))
		logmgamma = (chain[:,:,2].reshape(-1))
		alpha = -10**logmalpha
		beta = -10**logmbeta
		gamma = 10**logmgamma
		
		
		index = np.argmax(lnprob.flat)
		alpha_best = alpha[index]
		beta_best =  beta[index]
		gamma_best  = gamma[index]
		
		select(0, 2)
		bincount=25
		histogram(alpha, bincount=bincount)
		vline(alpha_best, color="red")
		select(1, 1)
		histogram(beta, bincount=bincount)
		vline(beta_best, color="red")
		select(2, 0)
		histogram(gamma, bincount=bincount)
		vline(gamma_best, color="red")
			
		select(0, 1)
		scatter(alpha, beta)
		symbol(alpha_best, beta_best, color="red")
		select(0, 0)
		scatter(alpha, gamma)
		symbol(alpha_best, gamma_best, color="red")
		select(1, 0)
		scatter(beta, gamma)
		symbol(beta_best, gamma_best, color="red")

		select(1, 2)
		scatter(np.log10(M3), slope3)
		xlim(7.5, 9.)
		ylim(-3, -0.3)
		
		select(1,2)
		
		if 0:
			Es = np.zeros_like(alpha)
			#for i in range(len(Es)):
			@np.vectorize
			def E(alpha, beta, gamma, rs, rho0, rmax):
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=rho0)
				#pot = mab.gd.potential.ProfileNumerical1dWrapper(200, profile)
				to_erg = 1.9891e43
				#E0 = pot.enclosed_binding_energy(rmax) * to_erg
				E0 = profile.enclosed_binding_energy(rmax) * to_erg
				return E0
			Es = E(alpha, beta, gamma, rs, rho0 * scaling, rmax)
		select(2,1)
		
		Mstar = 10**(-0.4*(Mv-4.84))
		Esn = 1e51
		Mmean = 0.4
		fsn = 0.0037
		Nsn = Mstar/Mmean * fsn
		Esn_total = Nsn * Esn
		
		eff = abs(Es / Esn_total)
		log_eff = np.log10(eff)
		histogram(log_eff, bincount=bincount)
		select(2,2)
		scatter(log_eff, (rs))
		
		
		i = 10
		profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha[i], beta=beta[i], gamma=gamma[i], rs=rs[i], rho0=rho0[i])
		to_erg = 1.9891e43
		E0 = profile.enclosed_binding_energy(rmax[i]) * to_erg
		print E0, Es[i]
		#dsa
		
		if 0:
			def f(logrs):
				#print logrs
				rs = 10**logrs[0]
				#print alpha, beta, gamma, rs
				profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=1.)
				profile.rho0 = scaling * target_dens_rmax / profile.densityr(rmax)
				#print "dens", target_dens_rmax, profile.densityr(rmax)
				mass_rmax = profile.enclosed_mass(rmax)
				#print mass_rmax, target_mass_rmax, rmax, scaling
				return (np.log10(mass_rmax) - np.log10(target_mass_rmax))**2
			logrs, = scipy.optimize.fmin(f, np.log10(bestrs), disp=False)
			rs = 10**logrs
			#print "best rs", rs
			profile = mab.gd.potential.TwoSlopeDensityWrapper(alpha=alpha, beta=beta, gamma=gamma, rs=rs, rho0=1.)
			profile.rho0 = scaling * target_dens_rmax / profile.densityr(rmax)
			mass_rmax = profile.enclosed_mass(rmax)
			diff = np.log10(mass_rmax) - np.log10(target_mass_rmax)
			
			
		#import atpy
		import asciitable
		table = dict(logeff=log_eff, alpha=alpha, beta=beta, gamma=gamma, rs=rs, logM3=np.log10(M3), slope3=slope3)
		asciitable.write(table=table, output="halo.asc", delimiter=",")
		draw()
	

		
	
