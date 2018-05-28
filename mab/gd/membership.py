# -*- coding: utf-8 -*-
from numpy import *
from kaplot import *
import os
import mab.cvsfile
import mab.gd.logging as logging
import scipy

logger = logging.getLogger("gd.membership")

def gaussian_(x, mu, sigma):
	return 1/(sigma*sqrt(2*pi)) * exp(-(mu-x)**2/(2*sigma**2))

class MemberRatiosZero(object):
	def load(self):
		pass
	def ratios(self, R1, R2):
		return 1.0, 0.0
	def ratios1(self, R):
		return 1.

class MemberRatios(object):
	def __init__(self, model, foreground_model, galaxy_model, light_profile, aperture):
		self.model = model
		self.foreground_model = foreground_model
		self.galaxy_model = galaxy_model
		self.light_profile = light_profile
		self.aperture = aperture
		#self.Rmax = Rmax
		
	def load(self):
		self.model.load()
		self.aperture.load()
		self.Rmax = self.aperture.getRmax()
		
	def ratios(self, R1, R2):
		model = self.model.fitter.get_model()
		ratio = abs(model.ratio)
		w1 =  ratio / (1. + ratio);
		w2 = 1 / (1. + ratio);
		f = w1
		
		def rho_m(R):
			return self.light_profile.densityR(R) / scipy.integrate.quad(lambda R: 2*pi*self.light_profile.densityR(R)*R, 0, self.Rmax)[0]
		def rho_n(R):
			return 1./scipy.integrate.quad(lambda R: 2 * pi * R, 0, self.Rmax)[0]
		
		#print "ratio,w1,w2", ratio, w1, w2, w1/w2
		#f = 0.4
		#f = 0.3
		#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all) * 2 * pi * R, R1, R2)[0]
		#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
		f_member = scipy.integrate.quad(lambda R: (rho_m(R)*(1-f)) * 2 * pi * R, R1, R2)[0]
		f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*f) * 2 * pi * R, R1, R2)[0]
		f_total = f_member + f_non_member 
		return f_member/f_total, f_non_member/f_total

	def ratios1(self, R):
		model = self.model.fitter.get_model()
		ratio = abs(model.ratio)
		w1 =  ratio / (1. + ratio);
		w2 = 1 / (1. + ratio);
		f = w1
		
		def rho_m(R):
			return self.light_profile.densityR(R) / scipy.integrate.quad(lambda R: 2*pi*self.light_profile.densityR(R)*R, 0, self.Rmax)[0]
		def rho_n(R):
			return 1./scipy.integrate.quad(lambda R: 2 * pi * R, 0, self.Rmax)[0]
		
		#print "ratio,w1,w2", ratio, w1, w2, w1/w2
		#f = 0.4
		#f = 0.3
		#f_member = scipy.integrate.quad(lambda R: (rho_m(R)*N_all_member/N_all) * 2 * pi * R, R1, R2)[0]
		#f_non_member = scipy.integrate.quad(lambda R: (rho_n(R)*N_all_non_member/N_all) * 2 * pi * R, R1, R2)[0]
		f_member =(rho_m(R)*(1-f)) # * 2 * pi * R
		f_non_member = (rho_n(R)*f)# * 2 * pi * R
		f_total = f_member + f_non_member 
		return f_member/f_non_member


class MembershipSigmaClipping(object):
	def __init__(self, input, output, membership_test, aperture, sigmafilter):
		self.input = input
		self.output = output
		self.membership_test = membership_test
		self.aperture = aperture
		self.sigmafilter = sigmafilter
		self.loaded = False
		
		
	def filter_aperture(self, stars, i, Nsigma=3.0):
			if not self.loaded:
				self.aperture.load()
				self.loaded = True
			aperture_stars = stars.filter(self.aperture.aperture_filter(i))
			N_total = len(aperture_stars)
			print N_total
			stars = self.sigmafilter.filter(aperture_stars)
			N_members = len(stars)
			fraction_members = N_members * 1. / N_total
			#mean_sigma, p_members = self.membership_test.calc_p_member(stars.vlos, stars.e_vlos, fraction_members)
			mean_sigma, p_members = self.membership_test.calc_p_member(aperture_stars.vlos, aperture_stars.e_vlos, fraction_members)
			print "mean_sigma", mean_sigma, "|", fraction_members
			#raise "dsa"
			
			N = 3.
			def newsigmatest(star):
				#return (star.vlos > -3*mean_sigma) & (star.vlos < 3*mean_sigma)
				return (star.vlos > -Nsigma*mean_sigma) & (star.vlos < Nsigma*mean_sigma)
			return stars.filter(newsigmatest)
			
			
		
	def run(self, args, opts, scope):
		self.input.load()
		allstars  = self.input.stars
		self.aperture.load()
		indices = range(len(self.aperture))
		print "aperture indices", indices
		#box()
		document(size="35cm,25cm")
		mozaic(6,3,box)
		
		#@mab.parallelize.parallelize(cores=opts.cores, info=scope["progressbar"])
		def do(i):
			#print allstars, self.aperture
			aperture_stars = allstars.filter(self.aperture.aperture_filter(i))
			N_total = len(aperture_stars)
			print N_total
			stars = self.sigmafilter.filter(aperture_stars)
			N_members = len(stars)
			fraction_members = N_members * 1. / N_total
			mean_sigma, p_members = self.membership_test.calc_p_member(stars.vlos, stars.e_vlos, fraction_members)
			
			print "mean_sigma", mean_sigma
			moststars = stars
			sigmas = []
			k4s = []
			k3s = []
			Ns = arange(2.5, 4, 0.1)
			
			for N in Ns: #[2.5, 3.0, 3.5, 4.0]:
				def newsigmatest(star):
					#return (star.vlos > -3*mean_sigma) & (star.vlos < 3*mean_sigma)
					return (star.vlos > -N*mean_sigma) & (star.vlos < N*mean_sigma)
				stars = moststars.filter(newsigmatest)
				sigmas.append(stars.vlos.std())
				k4s.append(mean(stars.vlos**4)/mean(stars.vlos**2)**2)
				k3s.append(mean(stars.vlos**3)/mean(stars.vlos**2)**(3./2))
			print "s2", sigmas
			print "k4", k4s
			print i, mean(k3s)
			
			select((i%2)*3+0, i/2)
			graph(Ns, sigmas)
			def mom(k, N):
				I, err = scipy.integrate.quad(lambda x: 2*gaussian(x, 0, mean_sigma)*x**k, 0, N*mean_sigma)
				return I
			graph(Ns, [mom(2, N)**0.5 for N in Ns], color="red")
			ylim(7, 12)
			
			select((i%2)*3+1, i/2)
			graph(Ns, k4s)
			graph(Ns, [mom(4, N)/mom(2, N)**2 for N in Ns], color="red")
			ylim(2,4)
			
			select((i%2)*3+2, i/2)
			histogram(aperture_stars.vlos, datamin=-50, datamax=50, binwidth=2., normalize=True)
			x = arange(-50, 50, 0.5)
			graph(x, gaussian(x, 0, mean_sigma), color="red")
			
			N = 3.
			def newsigmatest(star):
				#return (star.vlos > -3*mean_sigma) & (star.vlos < 3*mean_sigma)
				return (star.vlos > -N*mean_sigma) & (star.vlos < N*mean_sigma)
			stars = moststars.filter(newsigmatest)
			sigma = stars.vlos.std()
			graph(x, gaussian(x, stars.vlos.mean(), sigma), color="blue")
			
			
			print "===="
		
		
		for i in indices:
			do(i)
		draw()
		

class ModelWrapper(object):
	def __init__(self, model, objects, parameters, bounds, observation, filters=[]):
		self.model = model
		self.objects = objects
		self.parameters = parameters
		self.bounds = bounds
		self.observation = observation
		self.filters = filters
		
	def load(self):
		stars = self.observation.load()
		for filter in self.filters:
			stars = filter(stars)
		self.stars = stars
		
	def logp(self):
		return sum([log(self.model(star.vlos_helio, star.e_vlos)) for star in self.stars])
		
	def __len__(self):
		return len(self.parameters)
		
	def __getitem__(self, index):
		return getattr(self.objects[index], self.parameters[index])
	
	def __setitem__(self, index, value):
		print "set", self.parameters[index], "to", value
		return setattr(self.objects[index], self.parameters[index], value)

class SingleGaussian(object):
	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma
		
	def __call__(self, v, e_v):
		return gaussian_(v, self.mu, sqrt(self.sigma**2+e_v**2))
		
		
class WeightedModel(object):
	def __init__(self, fraction_foreground, foreground, galaxy):
		self.fraction_foreground = fraction_foreground
		self.foreground = foreground
		self.galaxy = galaxy
		
	def __call__(self, v, e_v):
		f1 = self.fraction_foreground
		f2 = 1. - f1
		#print f1 * self.foreground(v, e_v),  f2 * self.galaxy(v, e_v)
		return f1 * self.foreground(v, e_v) + f2 * self.galaxy(v, e_v)
		
class ModelFit(object):
	def __init__(self, model_wrapper, filename):
		self.model_wrapper = model_wrapper
		self.filename = filename
		
	def load(self):
		logger.info("loading fit parameters from: " + self.filename)
		values = numpy.load(self.filename)
		self.set(values)
	
	def save(self):
		logger.info("saving fit parameters to: " + self.filename)
		values = numpy.save(self.filename, self.x)

	def set(self, values):
		for i, value in enumerate(values):
			print i, value
			self.model_wrapper[i] = value
		
	def run(self, args, opts, scope):
		self.model_wrapper.load()
		bounds = self.model_wrapper.bounds
		x = [self.model_wrapper[k] for k in range(len(self.model_wrapper))]
		def f(values):
			self.set(values)
			y = -self.model_wrapper.logp()
			print y
			return y
		x = scipy.optimize.fmin_l_bfgs_b(f, x, None, bounds=bounds, approx_grad=True, iprint=1,factr=1e6,maxfun=200000)[0]
		print x
		self.x = x
		#self.x = [0.05]
		self.save()
		
		
class GalaticDistributionDoubleGaussian(object):
	def __init__(self, v1, sigma1, v2, sigma2, ratio):
		self.v1 = v1
		self.v2 = v2
		self.sigma1 = sigma1
		self.sigma2 = sigma2
		self.ratio = ratio
		self.f1 = 0.5
		self.f2 = 0.5
		
	def __call__(self, v, e_v):
		#sigma_eff = self.sigma1 + self.ratio * self.sigma2
		#k1 = 1./2 #/ (sqrt(2*pi) * sigma_eff)
		#k2 = 1./2 #self.ratio * k1
		#k1 = 1
		#k2 = self.sigma2/self.sigma1
		#return (k1 * gaussian_(v, self.v1, sqrt(self.sigma1**2+e_v**2)) + k2 * gaussian_(v, self.v2, sqrt(self.sigma2**2+e_v**2))) / (k1+k2)
		#self.f1 = self.ratio
		#self.f2 = 1 - self.ratio
		return self.f1 * gaussian_(v, self.v1, sqrt(self.sigma1**2+e_v**2)) + self.f2 * gaussian_(v, self.v2, sqrt(self.sigma2**2+e_v**2))
		
	def pdf(self, v, e_v=0.):
		return self(v, e_v=e_v)
		
	def rvs(self):
		return self.sample()
		
	def sample(self):
		#self.f1 = self.ratio
		#self.f2 = 1 - self.ratio
		if random.random() < self.f1:
			return random.normal(self.v1, self.sigma1)
		else:
			return random.normal(self.v2, self.sigma2)


class MembershipTest(object):
	def __init__(self, foreground):
		self.foreground = foreground
		
	#def calc_p_member(self, all_stars, member_stars, non_member_stars):
	def calc_p_member(self, vs, e_vs, fraction_members):
		#N_member = len(member_stars)
		#N_non_member = len(non_member_stars)
		#N_total = N_member + N_non_member
		#assert len(all_stars) == N_total
		f_member = fraction_members#N_member*1./N_total
		f_non_member = 1.-fraction_members#N_member*1./N_total
		
		mean_v = 0.
		def logp_j(v, e_v, mean_v, sigma):
			return log(fraction_members * gaussian_(v, mean_v, sqrt(sigma**2+e_v**2)) +  (1-fraction_members) * self.foreground(v, e_v))
			#ratio = 1.
			#return log((1-fraction_members) * self.foreground(v-mean_v, e_v))
			#return log(gaussian(v, mean_v, sqrt(sigma**2+e_v**2)))
		def logp(mean_v, sigma):
			#return sum([logp_j(star.vlos, star.e_vlos, mean_v, sigma) for star in stars])
			return sum([logp_j(v, e_v, mean_v, sigma) for v, e_v in zip(vs, e_vs)])
			
		logsigmas = arange(log(1.), log(100), 0.05)
		sigmas = 10**logsigmas
		#sigma = 10.
		logps = array([logp(mean_v, sigma) for sigma in sigmas])
		
		ps = exp(logps-logps.max())
		mean_sigma = sum(sigmas*ps)/sum(ps)
		#print "mean_sigma", mean_sigma
		if 0:
			clear()
			mozaic(2,1,box)
			#box()
			graph(logsigmas, ps)
			select(1,0)
			histogram(vs)
			draw()
		
		
		def p_member_j(v, e_v, mean_v, sigma):
			p_no_member = (1-fraction_members) * self.foreground(v, e_v)
			p_is_member = fraction_members * gaussian_(v, mean_v, sqrt(sigma**2+e_v**2))
			p_total =  p_is_member + p_no_member
			return p_is_member / (p_total)
		#p_members = [p_member_j(star.vlos, star.e_vlos, mean_v, mean_sigma) for star in stars]
		p_members = array([p_member_j(v, e_v, mean_v, mean_sigma) for v, e_v in zip(vs, e_vs)])
		if isnan(mean_sigma):
			import pdb; pdb.set_trace()
		return mean_sigma, p_members
		#for star, p in zip(stars, p_members):
		#	star.attributes.append("p_member")
		#	star.p_member = p
		#	outputstars.append(star)
		#p_members = array(p_members)

class MemberSingleGaussian(object):
	def __init__(self, v, sigma):
		self.v = v
		self.sigma = sigma
		
	def __call__(self, v, e_v):
		sigma = sqrt(self.sigma**2 + e_v**2)
		return gaussian_(v, self.v, sigma)

class MembershipTestVelocity(object):
	def __init__(self, modelpath, I0_background, I0_galaxy, light_model, aperture, observation, member_distribution, galatic_distribution):
		self.I0_background = I0_background
		self.I0_galaxy = I0_galaxy
		self.observation = observation
		self.member_distribution = member_distribution
		self.galatic_distribution = galatic_distribution
		self.light_model = light_model
		#self.output = output
		#self.filename = os.path.join(modelpath, self.output)
	
	def load(self):
		self.observation.load()
		
	def density_ratio(self, re):
		arcsec_per_kpc = self.light_model.arcsec_to_kpc(1.)
		re_kpc = self.light_model.arcsec_to_kpc(re)
		
		#re_kpc = self.light_model.arcsec_to_kpc(3600*0.8)
		#re_kpc = 2
		density_stars_per_kpcsq = self.light_model.densityR(re_kpc, self.I0_galaxy/arcsec_per_kpc**2)
		density_stars_per_arcsecsq = density_stars_per_kpcsq*arcsec_per_kpc**2
		return density_stars_per_arcsecsq  / (density_stars_per_arcsecsq + self.I0_background)
		
	def logp_sigma(self, stars, sigma, v, N_galaxy, N_mw):
		self.member_distribution.sigma = sigma
		self.member_distribution.v = v
		logp = 0
		N_tot = N_galaxy + N_mw
		f_galaxy = N_galaxy / float(N_tot)
		f_mw = N_mw / float(N_tot)
		f_mw *= 1.0
		tot = f_mw + f_galaxy
		f_mw = f_mw / tot
		f_galaxy = f_galaxy / tot
		
		for star in stars:
			#N_galaxy = self.density_ratio(star.re)
			#N_galactic = 1 - N_galaxy
			
			prv_no_member = self.galatic_distribution(star.vlos_helio, star.e_vlos) * f_mw
			prv_member = self.member_distribution(star.vlos_helio, star.e_vlos) * f_galaxy
			prv = prv_no_member + prv_member
			#prv = prv_member
			logp += log(prv)
		return logp
			
	def p_member(self, star, v_mean, v_sigma, N_galaxy, N_mw):
		self.member_distribution.sigma = v_sigma
		self.member_distribution.v = v_mean
		
		#N_galaxy = self.density_ratio(star.rc)
		#N_galactic = 1 - N_galaxy
		N_tot = N_galaxy + N_mw
		f_galaxy = N_galaxy / float(N_tot)
		f_mw = N_mw / float(N_tot)
		
		prv_no_member = self.galatic_distribution(star.vlos_helio, star.e_vlos) * f_mw
		prv_member = self.member_distribution(star.vlos_helio, star.e_vlos) * f_galaxy
		p_member_prior = 0.5
		p_no_member_prior = 1-p_member_prior
		p_member = prv_member * p_member_prior / (prv_member * p_member_prior + prv_no_member * p_no_member_prior)
		return p_member
		
		
	def test(self):
		stars = self.observation.stars
		#print len(stars)
		for star in stars:
			N_galaxy = self.density_ratio(star.re)
			N_galactic = 1 - N_galaxy
			
			prv_no_member = self.galatic_distribution(star.vlos_helio, star.e_vlos, N_galactic)
			prv_member = self.member_distribution(star.vlos_helio, star.e_vlos, N_galaxy)
			p_member_prior = 0.5
			p_no_member_prior = 1-p_member_prior
			p_member = prv_member * p_member_prior / (prv_member * p_member_prior + prv_no_member * p_no_member_prior)
			print "p member: ", p_member
			star.p_member = p_member
			while "p_member" in star.attributes: 
				star.attributes.remove("p_member")
			if "p_member" not in star.attributes: 
				star.attributes.append("p_member")
			

	def save(self):
		logger.info("writing to: %s" % self.observation.filename)
		mab.cvsfile.writecsv(self.observation.filename, self.observation.stars)