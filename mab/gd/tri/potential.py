# -*- coding: utf-8 -*-
from kaplot import *
import scipy

class NFWTri(object):
	def __init__(self, p, q, nfw_xaxis):
		self.p = p
		self.q = q
		self.nfw_xaxis = nfw_xaxis
		
	def run(self, args, opts, scope):
		global logrho_model
		N = 8
		#sigmas = zeros(N)  + arange(1,N)
		logsigmas = zeros(N) #(arange(N)-N/2.)*0.3
		logweights = zeros(N)
		#logweights = [0, 0]
		#logsigmas = [-2, 2]
		logxmin = -2.
		logxmax = 2.
		Nx = 100
		deltax = (logxmax-logxmin)/Nx
		logx = arange(logxmin, logxmax, deltax)
		x = 10**logx
		logrho_target = log(self.nfw_xaxis.densityr(x))
		logsigmas = arange(logxmin, logxmax, (logxmax-logxmin)/N)
		#logsigmas = log10(sigmas)
		#print self.nfw_xaxis.densityr(x)
		#dsa
		
		distance = 79/1000.
		parsec_km = 1.4959787068e8*648e3/pi
		conversion_factor = distance*1.0e6*tan(pi/(648e3))*parsec_km
		
		def f(params):
			global logrho_model
			sigmas = array([10**params[:N]]).T
			weights = array([10**params[N:]]).T
			#S = (-1/(2*sigmas**2)*x**2) + log(weights/(sigmas*sqrt(2*pi))**3)
			rho_model = zeros(Nx)
			for i in range(N):
				rho_model += exp((-x**2/(2*sigmas[i]**2)) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3))
			logrho_model = log(rho_model)
			logrho_model[isinf(logrho_model)] = -100
			#logrho_model = sum( S, axis=0 )
			#print S
			#print S.shape
			print "target  ", ("% 8.3f" * Nx) % tuple(logrho_target)
			print "model   ", ("% 8.3f" * Nx) % tuple(logrho_model)
			print "params  ", ("% 5.2f" * (N*2)) % tuple(params)
			print "sigmas  ", ("% 5.2f" * N) % tuple(sigmas.T[0])
			print "weights ", ("% 5.2f" * N) % tuple(weights.T[0])
			#print "weights ", weights.T
			value = log10(sum((logrho_model-logrho_target)**2))
			print "difference", value
			#raw_input()
			return value
			
		params = concatenate([logsigmas, logweights])
		print params
		bounds = None #[None] * N
		u = scipy.optimize.fmin_l_bfgs_b(f, params, None, bounds=bounds, approx_grad=True, iprint=1,factr=1e-10,maxfun=200000)[0]
		#u = params
		#u = scipy.optimize.fmin(f, params, maxiter=10000000, maxfun=1000000)
		print u
		sigmas = array(10.**u[:N])
		weights = array(10.**u[N:])
		logxmin = -3.
		logxmax = 3.
		Nx = 200
		deltax = (logxmax-logxmin)/Nx
		logx = arange(logxmin, logxmax, deltax)
		x = 10**logx
		logrho_target = log(self.nfw_xaxis.densityr(x))
		#logsigmas = arange(logxmin, logxmax, (logxmax-logxmin)/N)
		if 1:
			
			print sigmas
			rho_model = zeros(Nx)
			for i in range(N):
				rho_model += exp((-x**2/(2*sigmas[i]**2)) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3))
			logrho_model = log(rho_model)
			#logrho_model[isinf(logrho_model)] = -2
			#logrho_model[logrho_model<10] = -2
			print "target  ", ("% 8.3f" * Nx) % tuple(logrho_target)
			print "model   ", ("% 8.3f" * Nx) % tuple(logrho_model)
			print "params  ", ("% 5.2f" * (N*2)) % tuple(u)
			print "sigmas  ", ("% 5.2f" * N) % tuple(sigmas)
			print "weights ", ("% 5.2f" * N) % tuple(weights)
			#print "weights ", weights.T
			value = log10(sum((logrho_model-logrho_target)**2))
			print "difference", value
		logrho_model_total = logrho_model
		import pdb
		#pdb.set_trace()
		print logrho_model
		#box()
		mozaic(2,2,box)
		graph(logx, logrho_target)
		graph(logx, logrho_model, color="red", linestyle="dash")
		for i in range(N):
			logrho_model = (-1/(2*sigmas[i]**2)*x**2) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3)
			logrho_model[isinf(logrho_model)] = -10
			logrho_model[logrho_model<-10] = -10
			graph(logx, logrho_model, color="blue", linestyle="dot")
		#sigma, weight = 10**5, 1e20
		#logrho_model = (-1/(2*sigma**2)*x**2) + log(weight/(sigma*sqrt(2*pi))**3)
		#graph(logx, logrho_model, color="orange", linestyle="dot")
		ylim(min(logrho_target), max(logrho_target))
		select(1,0)
		graph(logx, 10**logrho_target)
		graph(logx, 10**logrho_model_total, color="red", linestyle="dash")
		for i in range(N):
			logrho_model = (-1/(2*sigmas[i]**2)*x**2) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3)
			logrho_model[isinf(logrho_model)] = -10
			logrho_model[logrho_model<-10] = -10
			graph(logx, 10**logrho_model, color="blue", linestyle="dot")
		#sigma, weight = 10**5, 1e20
		#logrho_model = (-1/(2*sigma**2)*x**2) + log(weight/(sigma*sqrt(2*pi))**3)
		#graph(logx, logrho_model, color="orange", linestyle="dot")
		#ylim(min(logrho_target), max(logrho_target))
		
		
		select(1,1)
		graph(10**logx, 10**logrho_target)
		graph(10**logx, 10**logrho_model_total, color="red", linestyle="dash")
		for i in range(N):
			logrho_model = (-1/(2*sigmas[i]**2)*x**2) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3)
			logrho_model[isinf(logrho_model)] = -10
			logrho_model[logrho_model<-10] = -10
			graph(10**logx, 10**logrho_model, color="blue", linestyle="dot")
		xlim(0, 10)
		
		select(0,1)
		graph(10**logx, logrho_target)
		graph(10**logx, logrho_model_total, color="red", linestyle="dash")
		for i in range(N):
			logrho_model = (-1/(2*sigmas[i]**2)*x**2) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3)
			logrho_model[isinf(logrho_model)] = -10
			logrho_model[logrho_model<-10] = -10
			graph(10**logx, logrho_model, color="blue", linestyle="dot")
		xlim(0, 10)
		ylim(min(logrho_target), max(logrho_target))
		
		print conversion_factor
		densities = []
		for i in range(N):
			# Msol/kpc**3
			
			mass = weights[i]
			#print mass
			#print "weight", w, "Msol/kpc^3"
			#mass = mass# /1000**3
			#print "weight", w, "Msol/pc^3"
			
			sigma = sigmas[i] # kpc
			density_central = mass / ( (sigma*1000*sqrt(2*pi))**3)
			#print "sigma", sigma, "kpc"
			#print "density", density_central, "Msol/kpc^3"
			#print "log", (-1/(2*sigmas[i]**2)*0**2) + log(weights[i]/(sigmas[i]*sqrt(2*pi))**3)
			sigma = sigma*1000*parsec_km # kpc to km
			sigma = sigma / conversion_factor # km to arcsec
			#print "sigma", sigma, "km"
			densities.append(density_central)

			# Msol/km^3
			#density = mass / ( (sigma*sqrt(2*pi))**2)
			#print "dens: ", mass / ( (sigma*sqrt(2*pi))**3), "Msol/km^3", 
			
			#print "sigma", sigma, "arcsec" # km to arcsec
			
			print density_central, sigma, "1.0 1.0"
			
		print "tot mass %e" % ( sum(weights))
		
		densities = array(densities) /parsec_km**3
		sigmas = sigmas*1000*parsec_km/conversion_factor * conversion_factor
		print densities, sigmas
		print "tot mass", 2*pi*(2*pi)**(1./2)*sum(densities*sigmas**3)
		
		print "m enc %e" % self.nfw_xaxis.enclosed_mass(10**2)
		print parsec_km, conversion_factor
		
		rmin, rmax = 10**-2, 10
		print "rlogmin, rlogmax:", log10(rmin*1000*parsec_km/conversion_factor), log10(rmax*1000*parsec_km/conversion_factor)
		
		draw()
			
			
		
		
		
		
		
		
		