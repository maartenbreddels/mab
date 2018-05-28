# -*- coding: utf-8 -*-
#from go import const
import math
from mab.astrounits import *
import mab.constants
import scipy.integrate
import numpy


if 0:
	Omega_m = 0.3
	Omega_lambda = 0.7
	sigma_8 = 0.9
	h = 0.7 / hubbleH


if 0:
	Omega_m = 0.27
	Omega_lambda = 0.73
	sigma_8 = 0.9
	h = 0.7 / hubbleH


if 0:
	# condLF.dat
	Omega_m = 0.27 #       !  Omega_m
	Omega_lambda = 0.73 #       !  Omega_Lambda
	h = 0.7 / hubbleH  #      !  h (=H0/100.)
	sigma_8 = 0.9 #       !  Sigma_8
	n = 1.00  #     !  nspec
	#1          !  BBKS (1) or EFW (2) power spectrum
	#0.0        !  WDM filter m_ass (0.0 = CDM)
	#2          !  Peacock & Dodds (1) or Smith etal (2)
	#1          !  filter (1=top-hat, 2=Gaussian, 3=sharp-k)
	#2          !  method of computing S(M) (1=num, 2=fit)
	#'cmb'      !  model directory
	#2          !  Data sample (1=YMB, 2=ALL, 3=EARLY, 4=LATE)
	
#H0 = 100 * h * (KM/S) / (MPC)


class LCDMCosmology(object):
	def __init__(self, Omega_m, Omega_b, Omega_lambda, sigma8, n, h, X0=0.77, Y0=0.23):
		self.Omega_m = Omega_m
		self.Omega_b = Omega_b
		self.Omega_lambda = Omega_lambda
		self.sigma8 = sigma8
		self.n = n
		self.h = h
		self.H0 = 100 * hubbleH * (KM/S) / (MPC)
		#rho_crit_approx = 1e-26 # kg/m^3
		self.rho_crit = (3 / (8 * math.pi * mab.constants.G) * self.H0**2) #.as(MSOL/MPC**3)

		#self.X = 0.75
		self.X0 = X0
		self.Y0 = Y0
		
		self.D_norm = 1.
		self.D_norm = 1/self.Dz(0)

		#age_universe = 13.7
	def Hsqz(self, z):
		# H**2 as function of z
		return self.H0 ** 2 * (self.Omega_m*(1+z)**3 + (self.Omega_lambda))
	


	def Dz(self, z):
		C = 5 * self.Omega_m * self.H0.asNumber()**2 / 2. #* (Hsqz(z).asNumber()**(0.5))
		def f(z):
			return (1+z)/self.Hsqz(z).asNumber()**(3./2.)
		#i, ie = scipy.integrate.quad(f, z, scipy.integrate.inf)
		i, ie = scipy.integrate.quad(f, z, numpy.inf)
		#print i, ie
		# old D_norm = 1/0.760009686409
		return C * math.sqrt(self.Hsqz(z).asNumber()) * i * self.D_norm  # TODO/WARNING: where does this number come from???


#print D_norm, 1/D_norm
	def dtz(self,z1,z2):
		H0nr = (self.H0 * h).asNumber(1/S)
		def f(a):
			return (self.Omega_m/a+self.Omega_lambda*a**2)**-0.5
		#return (Omega_m/a+Omega_r/a**2+ Omega_lambda*a**2)**-0.5
	
		#_astart = 1/(1+1000000.)
		#_astart = 1e-18
		a1 = 1/(z1+1.)
		a2 = 1/(z2+1.)
		dtH, dtH_e = scipy.integrate.quad(f, a1, a2)
		#print dtH, dtH_e, H0nr
		return dtH/H0nr

	def dtdz(z):
		H0nr = (H0 * h).asNumber(1/S)
		return (1./numpy.sqrt(Omega_m*(1+z)**5+Omega_lambda*(1+z)**2)) / H0nr
	
	def tz(z=0):
		#print "--", scipy.integrate.inf, 1/scipy.integrate.inf
		return dtz(scipy.integrate.inf, z)
	
#def a_and_t(	
#	delta_a/(self.omega_m0/a+self.omega_r0/a**2+
#				self.omega_l0*a**2+(1 - self.omega_0))**0.5/H0

# WMAP 3 year, Spergel et al 2007

h = 0.732 / hubbleH
Omega_m = (0.1277 / h**2).asNumber()
Omega_b = (0.02229 / h**2).asNumber()
Omega_lambda = 1-Omega_m -0.001#(0.716)
sigma8 = 0.761
#sigma_8 = 0.8
#sigma_8 = 0.9
n = 1.0

wmap3y = LCDMCosmology(Omega_m, Omega_b, Omega_lambda, sigma8, n, h) 

#Omega_r = 4.6e-5 #Omega_0-Omega_m-Omega_lambda
#Omega_r = 4.6e-1 #Omega_0-Omega_m-Omega_lambda


if __name__ == "__main__":
	from go import const
	print "H0 = ",H0.asNumber(hubbleH * KM/S/MPC)
	print "t = ",tz(0)
	print "t_0 =", (tz(0) * S).__as(GYR)
	print "dt(z=inf,z=20) =", (dtz(scipy.integrate.inf,20) * S)._as(GYR)
	print "dt(z=inf,z=30) =", (dtz(scipy.integrate.inf,30) * S)._as(GYR)
	print "t(z=20) =", (dtz(20,0) * S)._as(GYR)
	print "t(z=15) =", (dtz(15,0) * S)._as(GYR)
	print "t(z=1000) =", (dtz(1000,0) * S)._as(GYR)
	print "dt(z=30,z=6) =", (dtz(30,6) * S)._as(GYR)
	print "dt(z=25,z=6) =", (dtz(25,6) * S)._as(GYR)
	print "dt(z=20,z=6) =", (dtz(20,6) * S)._as(GYR)
	print "dt(z=20,z=6.5) =", (dtz(20,6.5) * S)._as(GYR)
	print "dt(z=30,z=6.5) =", (dtz(30,6.5) * S)._as(GYR)
	print "dt(z=15,z=6.5) =", (dtz(15,6.5) * S)._as(GYR)
	print "dt(z=15,z=5) =", (dtz(15,5) * S)._as(GYR)
	print "dt(z=15,z=3) =", (dtz(15,3) * S)._as(GYR)
	print "dt(z=1000,z=200) =", (dtz(1000,200) * S)._as(GYR)
	print "dt(z=inf,z=200) =", (dtz(scipy.integrate.inf,200) * S)._as(GYR)
	print "dt(z=30,z=29) =", (dtz(30,29) * S)._as(GYR)
	print "dt(z=20,z=19.6) =", (dtz(20,19.6) * S)._as(GYR)
	print "dt(z=15,z=14.8) =", (dtz(15,14.8) * S)._as(GYR)
	print "t_H(z=6.5) =", (1/(Hsqz(6.5)**0.5*h))._as(GYR)
	#print "dt(z=20,z=6.5) =", (dtz(20,6.5) * S)._as(GYR)
	#print 1-Omega_m-Omega_lambda, H0nr
	
	print Omega_m, Omega_b, Omega_lambda
	print Omega_m + Omega_b + Omega_lambda
	#print 3 / (8 * math.pi * const.G) * H0**2
	print (3 / (8 * math.pi * const.G) * H0**2*h**2)._as(MSOL/MPC**3)
	print rho_crit
	print (rho_crit * Omega_b * h**2)._as(MH/CM**3)
	print (rho_crit * Omega_b * h**2)._as(MSOL/MPC**3), "%g" % (rho_crit * Omega_b * h**2).asNumber(MSOL/MPC**3)
	print (rho_crit * Omega_b * X*h**2)._as(MH/MPC**3)
	print "D(0)",Dz(0)
	print "D(10)",Dz(10)
	print "D(20)",Dz(20)
	print "D(30)",Dz(30)
	print "Omega_m", Omega_m
	print "Omega_L", Omega_lambda
