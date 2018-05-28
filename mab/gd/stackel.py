# -*- coding: utf-8 -*-
from numpy import *
import scipy.integrate
import mab.constants
import mab.astrounits
import numpy.linalg
import scipy.optimize

G = mab.constants.G.asNumber(mab.astrounits.KM**2/mab.astrounits.S**2 *mab.astrounits.KPC/mab.astrounits.MSOL)


class FitNFW(object):
	def __init__(self, profile, q):
		self.profile = profile
		self.q = q
		self.diff([2.35261954 , 0.01283445,  0.58450231])
		self.diff([1.24864126,  0.06140263, 0.98544004])
		self.x0 = [0.6991764 ,  0.078862  ,  0.40928223]
		self.x0 = [0.86221903585327619, 0.080825486706676303, 0.41423125273733458]
		self.x0 = [4.467375291643692, -0.074906181118484363, 0.35106409774745989]
		self.x0 = [1.4662257044471669, -0.25709520500299526, 0.45801066538846114]
		self.x0 = [1.4662279597173182, 0.18069066516877216, 0.69394152664209907]
		self.x0 = [2.5991020787284618, 0.0019461299880039896, 0.80257889373091051]
		self.x0 = [8.3144059397237093, -0.24695845951746026, 0.88798258782965778]
		self.x0 = [3.4476466244544652, -0.027790092809930906, 0.78207811625999668]
		self.x0 = [1.2173834016119298, 0.244754406292084, 0.69298367713825393]
		self.x0 = [0.59857743233259042, 0.50686261296026014, 0.6063393742725729]

		self.x0 = [0.49841435041817839, 0.6265050544039914, 0.52570042453200916]

		self.x0 = [0.63210052987625132, 0.35026282216154092, 0.83774155260668171]
		self.x0 = [0.59857920538488785, 0.50686186544067791, 0.60633871570758668]

		self.x0 = [0.10000000000000001, 0.54049196791175746, 0.9646311952304778]
		#self.x0 = [0.10000000000000001, 0.1139315593327422, 0.9696177515197355]
		self.x0 = [0.080000000000000001, 0.2139315593327422, 0.9696177515197355]

		
		print "diff", self.diff(self.x0)
		
	def load(self):
		pass
		
	def run(self, args, opts, scope):
		self.f()
		
	def f(self):
		#x0 = [1., 1., self.q]
		#x0 = [0.78460624,  0.00485917,  0.4]
		print self.diff(self.x0)
		bounds = [(0.1, 1e6), (-3, 3), (0.1, 0.99)]
		res = scipy.optimize.fmin_l_bfgs_b(self.diff, self.x0, bounds=bounds, approx_grad=True)
		print "self.x0 = %r" % (list(res[0]),)
		print res[0]
		
	def diff(self, x):
		self.rho0, self.rs, g = x
		#self.profile.gamma = -g
		self.q = g
		self.rs = 10**self.rs
		#print x
		n_u = 20
		u = arange(n_u) / (n_u - 1.) * 0.75 - 0.75
		r = 10**u
		print r.min(), r.max()
		#print r
		rho_nfw = self.density_dir(r, pi/20.)
		rho_fit = self.profile.density_dir(r, pi/20.)
		f = log10
		diff1 = sum((f(rho_nfw) - f(rho_fit)) ** 2)
		
		rho_nfw = self.density_dir(r, pi/2-pi/20.)
		rho_fit = self.profile.density_dir(r, pi/2 - pi/20.)
		diff2 = sum((f(rho_nfw) - f(rho_fit)) ** 2)
		
		diff = diff1+diff2
		print "d", diff
		return diff
		
	def density_dir(self, r, theta):
		x = r * cos(theta)
		y = r * sin(theta)
		return self.rho_xy(x, y)
	
	def rho_xy(self, x, y):
		#rho0, rs = 3., 1.
		m = sqrt(x**2+(y/self.q)**2)
		rho = self.rho0 * (m/self.rs)**-1 * (1+m/self.rs)**-2
		return rho
		
class FitEinasto(FitNFW):
	#def __init__(self, profile, q):
	#	self.profile = profile
	def rho_xy(self, x, y):
		#rho0, rs = 3., 1.
		self.alpha = 0.57
		m = sqrt(x**2+(y/self.q)**2)
		#rho = self.rho0 * (m/self.rs)**-1 * (1+m/self.rs)**-2
		x = m/self.rs
		return self.rho0 * exp(-2./self.alpha * (x**self.alpha-1))
		
		#return rho
		
		

class StackelProfileAxiOblate(object):
	"""
	Prolate spheroidal coordinates are used (lambda, nu) to describe oblate mass models
		-gamma < nu <= -alpha <= lambda < inf
		
	""" 
	def __init__(self, alpha, gamma):
		self.alpha = alpha
		self.gamma = gamma
		self.a1 = sqrt(-self.alpha)
		self.a2 = sqrt(-self.gamma)
		
	def uniform_to_nu(self, uniform):
		self.nu_min = -self.gamma
		self.nu_max = -self.alpha
		return uniform*(self.nu_max-self.nu_min) + self.nu_min
		
	def uniform_to_lambda(self, uniform, lambda_max):
		lambda_min = -self.alpha
		return uniform*(lambda_max-lambda_min) + lambda_min
	
	def xy_to_conf(self, x, y):
		#u = sqrt((y**2 + self.alpha - self.gamma)**2)
		#la = 1./.2 * (y**2 - self.alpha + u - self.gamma)
		#nu = 1./2. * (x**2 - self.alpha - u - self.gamma)
		#3c1 = -self.alpha - self.gamma + x**2 + y**2
		#c2 = self.alpha * self.gamma - self.gamma * x**2 - self.alpha * y**2
		#print "cc", c1, c2
		##print c1 + sqrt(c1**2 - 4*c2)
		#la = 1./2. *(c1 + sqrt(c1**2 - 4*c2))
		#nu = 1./2. *(c1 - sqrt(c1**2 - 4*c2))
		a1 = x**2 + y**2 -self.alpha - self.gamma
		a2 = x**4 + (y**2 + self.alpha - self.gamma) ** 2 + 2 * x**2*(y**2-self.alpha + self.gamma)
		#print "a1,a2", a1, a2
		la = 0.5 * (a1 + sqrt(a2))
		nu = 0.5 * (a1 - sqrt(a2))
		return nu, la
	
	def conf_to_xy(self, nu, la):
		x = sqrt(((la + self.alpha) * (nu + self.alpha))/(self.alpha - self.gamma))
		y = sqrt(((la + self.gamma) * (nu + self.gamma))/(self.gamma - self.alpha))
		return x, y
	
	def _Psi_tau(self, tau):
		P, err = scipy.integrate.quad(self.psi_tau, -self.gamma, tau)
		return P
	def Psi_tau(self, tau):
		return vectorize(self._Psi_tau)(tau)
	
	def rho_conf(self, nu, la):
		g_la = (la + self.alpha) / (la - nu)
		g_nu = (nu + self.alpha) / (nu - la) 
		return g_la**2 * self.psi_tau(la) + g_nu**2 * self.psi_tau(nu) + 2 * g_la * g_nu * (self.Psi_tau(la) - self.Psi_tau(nu))/(la - nu)
	
	def rho_xy(self, x, y):
		nus, las = self.xy_to_conf(x, y)
		return self.rho_conf(nus, las)
	
	def density_dir(self, r, theta):
		x = r * cos(theta)
		y = r * sin(theta)
		return self.rho_xy(x, y)
	
	
	def _G_tau(self, tau):
		#C = 2 * pi * G / sqrt(abs(tau+self.alpha))
		C = 2 * pi * G / sqrt(abs(tau+self.gamma))
		def f(sigma):
			#return sqrt(abs(sigma+self.alpha))/(2*(sigma+self.gamma)) * self.Psi_tau(sigma)
			return (sigma+self.alpha)/(2*(sigma+self.gamma)**(3./2)) * self.Psi_tau(sigma)
		#I, err = scipy.integrate.quad(f, -self.alpha, tau)
		I, err = scipy.integrate.quad(f, -self.gamma, tau)
		U = - 2 * pi * G * self.psiinf + C * I 
		return -U
		
	def _G_tau_prime(self, tau):
		#C = 2 * pi * G / sqrt(abs(tau+self.alpha))
		C = 2 * pi * G / sqrt(abs(tau+self.gamma))
		C1 = 2 * pi * G / sqrt(abs(tau+self.gamma))
		dtau = 1e-8
		C2 = 2 * pi * G / sqrt(abs(tau+dtau+self.gamma))
		Cprime = -1 * pi * G / abs(tau+self.gamma)**(3./2) * sign(tau+self.gamma)
		#print "Cprime", Cprime, (C2-C1)/dtau
		def f(sigma):
			#return sqrt(abs(sigma+self.alpha))/(2*(sigma+self.gamma)) * self.Psi_tau(sigma)
			return (sigma+self.alpha)/(2*(sigma+self.gamma)**(3./2)) * self.Psi_tau(sigma)
		I, err = scipy.integrate.quad(f, -self.gamma, tau)
		#I = self._G_tau(tau)/
		#I1, err = scipy.integrate.quad(f, -self.gamma, tau)
		#I2, err = scipy.integrate.quad(f, -self.gamma, tau+dtau)
		#print "G'", (I2-I1) / dtau, f(tau)
		Uprime = C * f(tau) + I * Cprime
		return -Uprime 
		
	def G_tau(self, tau):
		return vectorize(self._G_tau)(tau)
	def G_tau_prime(self, tau):
		return vectorize(self._G_tau_prime)(tau)
	
	def potential_conf(self, nu, la):
		#a = (la + self.alpha) * self.G_tau(la) - (nu + self.alpha) * self.G_tau(nu)
		a = (la + self.gamma) * self.G_tau(la) - (nu + self.gamma) * self.G_tau(nu)
		return -a/ (la - nu)
		#g_nu * self.U(nu) + g_la * self.U(la)
		#return (self.f(la) - self.f(nu)) / (la - nu)
	
	def dphi_dla(self, nu, la):
		a = (la + self.gamma) * self.G_tau(la) - (nu + self.gamma) * self.G_tau(nu)
		print a, (la - nu)**2, (self.G_tau(la)+(self.gamma+la) * self.G_tau_prime(la)), (la-nu)
		return a/ (la - nu)**2 - (self.G_tau(la)+(self.gamma+la) * self.G_tau_prime(la))/(la-nu) 
	
	def dphi_dnu(self, nu, la):
		a = (la + self.gamma) * self.G_tau(la) - (nu + self.gamma) * self.G_tau(nu)
		return -a/ (la - nu)**2 - (-self.G_tau(nu)-(self.gamma+nu) * self.G_tau_prime(nu))/(la-nu) 
	
	def dphi_dx(self, x, y):
		nu, la = self.xy_to_conf(x, y)
		Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))
		Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))
		a = x / (2*(la+self.alpha)*Psq)
		b = x / (2*(nu+self.alpha)*Qsq)
		c = y / (2*(la+self.gamma)*Psq)
		d = y / (2*(nu+self.gamma)*Qsq)
		return a * self.dphi_dla(nu, la) + b * self.dphi_dnu(nu, la)
		
	def dphi_dy(self, x, y):
		nu, la = self.xy_to_conf(x, y)
		Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))
		Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))
		a = x / (2*(la+self.alpha)*Psq)
		b = x / (2*(nu+self.alpha)*Qsq)
		c = y / (2*(la+self.gamma)*Psq)
		d = y / (2*(nu+self.gamma)*Qsq)
		return c * self.dphi_dla(nu, la) + d * self.dphi_dnu(nu, la)
	
	
	def potential_xy(self, x, y):
		nus, las = self.xy_to_conf(x, y)
		return self.potential_conf(nus, las)
	
	def C_tau_mod(self, nu, la):
		c1 = (la+self.alpha) * (nu+self.alpha) * (self.gamma + self.alpha)
		c2 = (la + self.gamma) * (nu + self.gamma)# * (-self.gamma + self.gamma) last term is zero
		return (la-self.alpha) / (2*sqrt(2) * sqrt(c1) * sqrt(c2))
	
	def N_tau(self, tau, E, i2, i3):
		return E - i2 / (tau+self.alpha) - i3 / (tau + self.gamma) + self.G_tau(tau)
	
	def n_tau(self, tau, E, i2, i3):
		return E - i2 / (tau+self.alpha) - i3 / (tau + self.gamma) + self.G_tau(tau)
	
	def p_sq_tau(self, tau, E, i2, i3):
		a = E - i2 / (tau+self.alpha) - i3 / (tau+self.gamma) + self.G_tau(tau)
		return a / (2*(tau+self.alpha))
	
	def p_tau(self, tau, E, i2, i3):
		return sqrt(self.p_sq_tau(tau, E, i2, i3))
	
	def v_conf(self, nu, la, E, i2, i3):
		Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))
		Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))
		P = sqrt(Psq)
		Q = sqrt(Qsq)
		vla = self.p_sq_tau(la, E, i2, i3)**0.5 / P
		vnu = self.p_sq_tau(nu, E, i2, i3)**0.5 / Q
		return vnu, vla
	
	def v_xy(self, x, y, E, i2, i3):
		nu, la = self.xy_to_conf(x, y) 
		Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))# * 1.2
		Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))# * 1.2
		P = sqrt(Psq)
		Q = sqrt(Qsq)
		a = x / (2*(la+self.alpha)*P)
		b = x / (2*(nu+self.alpha)*Q)
		c = y / (2*(la+self.gamma)*P)
		d = y / (2*(nu+self.gamma)*Q)
		#vx = vla * a / P + vnu * b / Q#* Qsq
		vnu, vla = self.v_conf(nu, la, E, i2, i3)
		vx = vla * a + vnu * b
		vy = vla * c + vnu * d
		return vx, vy
		
	
	def density_orbit_conf(self, nu, la, E, i2, i3):
		return vectorize(self._density_orbit_conf)(nu, la, E, i2, i3)
	def _density_orbit_conf(self, nu, la, E, i2, i3):
		if 1:
			#print nu, la
			dE = 1e-6*abs(E)
			di2 = dE
			di3 = dE
			try:
				if 0:
					m1 = [	(self.p_tau(nu, E+dE, i2    , i3    )  - self.p_tau(nu, E, i2, i3))/dE,
							(self.p_tau(nu, E   , i2+di2, i3    )  - self.p_tau(nu, E, i2, i3))/di2,
							(self.p_tau(nu, E   , i2    , i3+di3)  - self.p_tau(nu, E, i2, i3))/di3]
					m2 = [	(self.p_tau(la, E+dE, i2    , i3    )  - self.p_tau(la, E, i2, i3))/dE,
							(self.p_tau(la, E   , i2+di2, i3    )  - self.p_tau(la, E, i2, i3))/di2,
							(self.p_tau(la, E   , i2    , i3+di3)  - self.p_tau(la, E, i2, i3))/di3]
					m3 = [0, sqrt(2*i2), 0 ]
				else:
					#a = E - i2 / (tau+self.alpha) - i3 / (tau+self.gamma) + self.G_tau(tau)
					#return a / (2*(tau+self.alpha))
					p_nu = self.p_tau(nu, E, i2, i3)
					p_la = self.p_tau(la, E, i2, i3)
					
					c = (2*(nu+self.alpha))
					m1 = [1./p_nu / c, -1./p_nu/c/(nu+self.alpha), -1./p_nu/c/ (nu+self.gamma)]
					c = (2*(la+self.alpha))
					m2 = [1./p_la / c, -1./p_la/c/(la+self.alpha), -1./p_la/c/ (la+self.gamma)]
					m3 = [0, 1/sqrt(2*i2), 0]
				mat = array([m1, m3, m2])
				Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))# * 1.2
				Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))# * 1.2
				P = sqrt(Psq)
				Q = sqrt(Qsq)
				J = abs(numpy.linalg.det(mat))
				return J /(P * Q)
			except ZeroDivisionError:
				return 0
		else:
			if 0:
				w1 = sqrt( (e-i2/(la+self.alpha) - i3/(la+self.gamma) + self.G_tau(la)) / (la+self.gamma) )
				w2 = sqrt( (e-i2/(nu+self.alpha) - i3/(nu+self.gamma) + self.G_tau(nu)) / (nu+self.gamma) )
				o = (4 * sqrt(2) * (i3*(la+self.alpha) + (la+self.gamma)*(i2-E*(la+self.alpha)) - (la+self.alpha)*(la+self.gamma)*self.G_tau(la)) *\
						(i3*(nu+self.alpha) + (nu+self.gamma) * (i3-E*(nu+self.alpha)) - (nu+self.alpha) * (nu+self.gamma)*self.G_tau(nu)))
				print w1
				print w2
				print o
				#return -sqrt(i3) * (la + self.gamma) * (la - nu) * (nu + self.gamma) * w1 * w2 / o
			
			
			Psq = (la - nu) / (4*(la+self.alpha) * (la + self.gamma))# * 1.2
			Qsq = (nu - la) / (4*(nu+self.alpha) * (nu + self.gamma))# * 1.2
			P = sqrt(Psq)
			Q = sqrt(Qsq)
			a = sqrt(i3) * sqrt(la + self.gamma) * (la -nu) * sqrt(nu + self.gamma)
			b = sqrt(self.n_tau(la, E, i2, i3) * (la + self.alpha) * (la + self.gamma)) * \
				sqrt(self.n_tau(nu, E, i2, i3) * (nu + self.alpha) * (nu + self.gamma) *\
				32 * (la + self.alpha) * (la + self.gamma) * (nu + self.alpha) * (nu + self.gamma))
			print a, b
			return (a/b) / (P*Q)
		#return (self.n_tau(nu, E, i2, i3) > 0) * 1
		p_sq_nu = self.p_sq_tau(nu, E, i2, i3)
		p_sq_la = self.p_sq_tau(la, E, i2, i3)
		return ((p_sq_nu>0) & (p_sq_la>0)) * 1.
		return sqrt(-a/b)
	
	def V_eff_tau(self, tau, i2, i3):
		return i2/(tau+self.alpha) + i3 / (tau + self.gamma) - self.G_tau(tau)
		
		#return self.C_tau(la, nu) / sqrt(-self.N_tau(la) * self.N_tau(nu)) * sqrt(-self.N_tau(self.gamma))
		# avoid div by zero, nultiply original eq by sqrt(-gamma+gamma)
		#N_tau_mod =   i3 # (-gamma+gamma) * self.G_tau(tau)
		#print self.G_tau(-self.gamma)
		#print self.C_tau_mod(la, nu), self.N_tau(la, E, i2, i3),self.N_tau(nu, E, i2, i3), N_tau_mod
		#return self.C_tau_mod(la, nu) / sqrt(self.N_tau(la, E, i2, i3) * self.N_tau(nu, E, i2, i3) * -N_tau_mod)
		
	#def H(self, nu, la, E, i2, i3):
	
	def density_orbit_xy(self, x, y, E, i2, i3):
		nu, la = self.xy_to_conf(x, y)
		return self.density_orbit_conf(nu, la, E, i2, i3)
	
	def find_i3_max(self, E):
		def f(tau):
			#print tau,
			# make sure tau >= alpha
			r = abs(tau + self.alpha) 
			tau = r - self.alpha
			#print tau
			i3 = (E + self.G_tau(tau)) * (tau + self.gamma)
			return -i3
		tau_i3_max, neg_i3_max = scipy.optimize.fmin_powell(f, 1-self.alpha, full_output=True, disp=False)[:2]
		i3_max = -neg_i3_max
		r = abs(tau_i3_max + self.alpha) 
		tau_i3_max = r - self.alpha
		return i3_max, tau_i3_max
	
	def find_i2_max(self, E, i3):
		def f(tau):
			#print tau,
			# make sure tau >= alpha
			r = abs(tau + self.alpha) 
			tau = r - self.alpha
			#print tau
			i2 = (E + self.G_tau(tau) - i3 / (tau + self.gamma)) * (tau + self.alpha)
			return -i2
		tau_i2_max, neg_i2_max = scipy.optimize.fmin_powell(f, 1-self.alpha, full_output=True, disp=False)[:2]
		i2_max = -neg_i2_max
		r = abs(tau_i2_max + self.alpha) 
		tau_i2_max = r - self.alpha
		return i2_max, tau_i2_max
	
	def find_nu_max(self, E, i2, i3):
		def f(nu):
			#a = E - i2 / (tau+self.alpha) - i3 / (tau+self.gamma) + self.G_tau(tau)
			#i2 = (E + self.G_tau(tau) - i3 / (tau + self.gamma)) * (tau + self.alpha)
			return self.p_sq_tau(nu, E, i2, i3)
		eps = 1e-9
		nu_max, res = scipy.optimize.brentq(f, -self.gamma+eps, -self.alpha-eps, full_output=True, disp=False)
		print nu_max, self.p_sq_tau(nu_max, E, i2, i3)
		print "test", self.p_sq_tau(-self.gamma+1e-9, E, i2, i3)
		return nu_max
		#return a / (2*(tau+self.alpha))
		
	def find_la_range(self, E, i2, i3, la_max):
		def f(la):
			#a = E - i2 / (tau+self.alpha) - i3 / (tau+self.gamma) + self.G_tau(tau)
			#i2 = (E + self.G_tau(tau) - i3 / (tau + self.gamma)) * (tau + self.alpha)
			return self.p_sq_tau(la, E, i2, i3)
		eps = 1e-9
		la, res = scipy.optimize.brentq(f, -self.alpha+eps, la_max-eps, full_output=True, disp=False)
		print ">>", la_max, la, self.p_sq_tau(la_max, E, i2, i3), self.p_sq_tau(la, E, i2, i3)
		scale = 1.001
		if self.p_sq_tau(la/scale, E, i2, i3) < 0:
			print "peri", 
			la_max, res = scipy.optimize.brentq(f, la, 10**4, full_output=True, disp=False)
			la_min = la
			print "las", la_min, la_max
			assert self.p_sq_tau(la_min/scale, E, i2, i3) < 0
			assert self.p_sq_tau(la_min*scale, E, i2, i3) >= 0
			assert self.p_sq_tau(la_max/scale, E, i2, i3) >= 0
			assert self.p_sq_tau(la_max*scale, E, i2, i3) < 0
		else:
			print "apo"
			not implemented
		return la_min, la_max
		
		
		
		
class CommandFindIRange(object):
	def __init__(self, profile):
		self.profile = profile
		
	def run(self, args, opts, scope):
		re = 0.1
		E = self.profile.potential_xy(0.01, re)
		return i2/(tau+self.alpha) + i3 / (tau + self.gamma) - self.G_tau(tau)
		
class CoredPowerLawSpheroid(StackelProfileAxiOblate):
	def __init__(self, M, alpha, gamma, s=4):
		StackelProfileAxiOblate.__init__(self, alpha, gamma)
		self.M = M
		self.rho0 = 1
		self.c = sqrt(-self.gamma)
		self.s = s
		
	def psi_tau(self, tau):
		return self.rho0 * self.c ** self.s / tau**(self.s/2)
	
	def Psi_tau(self, tau):
		assert self.s != 2
		return 2. * self.rho0 * self.c**2 / (2.-self.s) * (self.c**(self.s-2.)/tau**((self.s-2.)/2.) - 1.)
	
class PerfectSpheroid(CoredPowerLawSpheroid):
	def __init__(self, M, alpha, gamma):
		CoredPowerLawSpheroid.__init__(self, M, alpha, gamma, 4)
		self.psiinf = self.Psi_tau(1e8)
		print "psiinf", self.psiinf
		self.fast_ = mab.gd.gdfast.PerfectSpheroidOblate(1., gamma, alpha, G);
		
	def fast(self):
		return self.fast_
		
	#def psi_tau(self, tau):
	#	return self.rho0 * self.c ** self.s / tau**(self.s/2)
	
	def rho_conf(self, nu, la):
		return self.rho0 * (self.alpha * self.gamma / (nu * la) ) ** 2
	
	def rho_conf(self, nu, la):
		return vectorize(self._rho_conf)(nu, la)
	def _rho_conf(self, nu, la):
		return self.fast_.density_conf(nu, la)
	def rho_xy(self, x, y):
		return vectorize(self._rho_xy)(x, y)
	def _rho_xy(self, x, y):
		#print "x,y",x, y
		return self.fast_.densityxz(x, y)
		msq = x**2/self.a1**2 + y**2/self.a2**2
		return self.rho0 * (1+msq)**-2
	
	def _G_tau(self, tau):
		#return 2 * pi * G * self.rho0 * -self.alpha * sqrt(-self.gamma/(tau+self.gamma)) * arctan(sqrt((tau+self.gamma)/-self.gamma))
		return self.fast_.G_tau(tau);
	
	def _G_tau_prime(self, tau):
		return self.fast_.G_tau_prime(tau);
		#x = (self.gamma+tau)/self.gamma
		#D = (-1./x)**(3./2) * (self.gamma*sqrt(-x) + tau * arctan(sqrt(-x)))/(2*self.gamma * tau)
		#return 2 * pi * G * self.rho0 * -self.alpha * D
		 
	def dphi_dla(self, nu, la):
		return self.fast_.dphi_dla(nu, la);
	
	def dphi_dnu(self, nu, la):
		return self.fast_.dphi_dnu(nu, la);
	
	def dphi_dx(self, x, y):
		return vectorize(self._dphi_dx)(x, y) 
	def dphi_dy(self, x, y):
		return vectorize(self._dphi_dy)(x, y) 
	def _dphi_dx(self, x, y):
		return self.fast_.dphidx(x, y);
	def _dphi_dy(self, x, y):
		return self.fast_.dphidz(x, y);
	
	def dphi_dx_num(self, x, y):
		return vectorize(self._dphi_dx_num)(x, y) 
	def dphi_dy_num(self, x, y):
		return vectorize(self._dphi_dy_num)(x, y) 
	def _dphi_dx_num(self, x, y):
		dx = 1e-6 
		return (self.potential_xy(x+dx,y) - self.potential_xy(x,y))/dx 
	def _dphi_dy_num(self, x, y):
		dy = 1e-6 
		return (self.potential_xy(x,y+dy) - self.potential_xy(x,y))/dy 
	
class TwoSlopeSpheroid(StackelProfileAxiOblate):
	def __init__(self, M, alpha, gamma, s1=-1, s2=-4):
		StackelProfileAxiOblate.__init__(self, alpha, gamma)
		self.M = M
		self.rho0 = 1
		self.c = sqrt(-self.gamma)
		self.s1 = s1
		self.s2 = s2
		#self.speed1 = speed1
		#self.speed2 = speed2
		self.rs1 = 0.01
		self.rs2 = 1.
		
	def __Psi_tau(self, tau):
		h2f1 = scipy.special.hyp2f1
		c1 = 1./(-2 + self.s1 - self.s2) * (self.rs2*(self.rs2 + self.gamma + tau)/(self.rs1 - self.rs2) + 0j)**(-self.s1/2)
		f1 = h2f1(-self.s1/2., 1/2.*(2 - self.s1 + self.s2), 1/2.*(4 - self.s1 + self.s2), -self.rs2/(self.rs1 - self.rs2) - 0.0000j)
		f2 = h2f1(-self.s1/2., 1/2.*(2 - self.s1 + self.s2), 1/2.*(4 - self.s1 + self.s2), -(self.rs2+self.gamma+tau)/(self.rs1 - self.rs2) - 0.0000j)
		a1 =  2*self.rs2**(1 + self.s2/2.)*(self.rs2 + self.gamma + tau)**(self.s1/2.)    * f1
		a2 = -2*self.rs2**(self.s1/2.)    *(self.rs2 + self.gamma + tau)**(1 + self.s2/2.)*f2
		#a2 = 1.
		#print "a1", a1
		#print "a2", a2
		#print "f1", f1
		#print "f2", f2
		#print "c1", c1
		#print "tau", tau
		#print -(self.s1/2.), 1/2.*(2 + self.s1 - self.s1), 1/2.*(4 + self.s1 - self.s2), -(self.rs2/(self.rs1 - self.rs2)) + 0j
		#print h2f1(-(self.s1/2.), 1/2.*(2 + self.s1 - self.s1), 1/2.*(4 + self.s1 - self.s2), -(self.rs2/(self.rs1 - self.rs2)) + 0j)
		#import pdb
		#pdb.set_trace()
		print imag(c1 * (a1 + a2))
		return real(c1 * (a1 + a2))
  
	
	def psi_tau(self, tau):
		assert all(tau >= self.gamma)
		ysq = (tau + self.gamma)
		y = sqrt(abs(ysq))
		#ysq[ysq > -1e-10] = 0
		#ysq = maxinum(ysq, 0)
		#y = sqrt(ysq)
		#assert not any(isnan(y)) 
		#print y
		#return self.rho0 * (1+y)**(self.s1) * (1.+y)**(self.s2-self.s1)
		#return self.rho0 * (self.c+y)**(self.s2)
		#return self.rho0 * (0.01+y**2)**((self.s1)/2)# * (1+y**2)**((self.s2-self.s1)/2)
		return self.rho0 * (0.01+y**2)**((self.s1)/2.) * (1+y**2)**((self.s2-self.s1)/2.)
	
class NFWLike(StackelProfileAxiOblate):
	def __init__(self, M, alpha, gamma):# s1=-1, s2=-3):
		StackelProfileAxiOblate.__init__(self, alpha, gamma)
		self.M = M
		self.rho0 = 1
		#self.c = sqrt(-self.gamma)
		self.rs1 = 0.01
		self.rs2 = 1.
		self.s1 = -1
		self.s2 = -3
		self.psiinf = self.Psi_tau(1e8)
		
	def Psi_tau(self, tau):
		return -2 * (arctan(sqrt(self.rs1/(self.rs2-self.rs1)))-arctan(sqrt(self.rs1+self.gamma + tau)))/(sqrt(self.rs2 - self.rs1))
  
	
	def psi_tau(self, tau):
		assert all(tau >= self.gamma)
		ysq = (tau + self.gamma)
		y = sqrt(abs(ysq))
		#ysq[ysq > -1e-10] = 0
		#ysq = maxinum(ysq, 0)
		#y = sqrt(ysq)
		#assert not any(isnan(y)) 
		#print y
		#return self.rho0 * (1+y)**(self.s1) * (1.+y)**(self.s2-self.s1)
		#return self.rho0 * (self.c+y)**(self.s2)
		#return self.rho0 * (0.01+y**2)**((self.s1)/2)# * (1+y**2)**((self.s2-self.s1)/2)
		return self.rho0 * (self.rs1+y**2)**((self.s1)/2.) * (self.rs2+y**2)**((self.s2-self.s1)/2.)
		
class EinastoLike(StackelProfileAxiOblate):
	def __init__(self, M, alpha, gamma, einasto_alpha, einasto_rs):# s1=-1, s2=-3):
		StackelProfileAxiOblate.__init__(self, alpha, gamma)
		self.M = M
		self.rho0 = 1
		#self.c = sqrt(-self.gamma)
		self.einasto_alpha = einasto_alpha
		self.einasto_rs = einasto_rs
		#self.rs1 = 0.01
		#self.rs2 = 1.
		#self.s1 = -1
		#self.s2 = -3
		self.psiinf = self.Psi_tau(1e8)
		
	#def Psi_tau(self, tau):
	#	return -2 * (arctan(sqrt(self.rs1/(self.rs2-self.rs1)))-arctan(sqrt(self.rs1+self.gamma + tau)))/(sqrt(self.rs2 - self.rs1))
  
	
	def psi_tau(self, tau):
		assert all(tau >= self.gamma)
		ysq = (tau + self.gamma)
		y = sqrt(abs(ysq))
		#ysq[ysq > -1e-10] = 0
		#ysq = maxinum(ysq, 0)
		#y = sqrt(ysq)
		#assert not any(isnan(y)) 
		#print y
		#return self.rho0 * (1+y)**(self.s1) * (1.+y)**(self.s2-self.s1)
		#return self.rho0 * (self.c+y)**(self.s2)
		#return self.rho0 * (0.01+y**2)**((self.s1)/2)# * (1+y**2)**((self.s2-self.s1)/2)
		x = y/self.einasto_rs
		return self.rho0 * exp(-2./self.einasto_alpha * (x**self.einasto_alpha-1))
		#return self.rho0 * (self.rs1+y**2)**((self.s1)/2.) * (self.rs2+y**2)**((self.s2-self.s1)/2.)		

class TwoSlopeSpheroidSpeed(StackelProfileAxiOblate):
	def __init__(self, M, alpha, gamma, s1=-1, s2=-4, speed1=1., speed2=1.):
		StackelProfileAxiOblate.__init__(self, alpha, gamma)
		self.M = M
		self.rho0 = 1
		self.c = sqrt(-self.gamma)
		self.s1 = s1
		self.s2 = s2
		self.speed1 = speed1
		self.speed2 = speed2
		
	def psi_tau(self, tau):
		assert all(tau >= self.gamma)
		ysq = (tau + self.gamma)
		y = sqrt(abs(ysq))
		#ysq[ysq > -1e-10] = 0
		#ysq = maxinum(ysq, 0)
		#y = sqrt(ysq)
		#assert not any(isnan(y)) 
		#print y
		#return self.rho0 * (1+y)**(self.s1) * (1.+y)**(self.s2-self.s1)
		#return self.rho0 * (self.c+y)**(self.s2)
		#return self.rho0 * (0.01+y**2)**((self.s1)/2)# * (1+y**2)**((self.s2-self.s1)/2)
		return self.rho0 * (0.01+y**self.speed1)**((self.s1)/self.speed1) * (1+y**self.speed2)**((self.s2-self.s1)/self.speed2)
	
