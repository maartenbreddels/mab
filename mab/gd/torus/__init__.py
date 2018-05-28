import numpy
from numpy import *
import mab.gd.potential
import mab.gd.schw.galaxy
import mab.gd.torus
from scipy.optimize import fsolve, fmin, brentq, fmin_ncg, fsolve
import os

class TorusFit(object):
	def __init__(self, galaxy, toytorus, J1, J2, Nn=0, name="default"):
		self.galaxy = galaxy
		self.toytorus = toytorus
		self.J1 = J1
		self.J2 = J2
		self.Nn = Nn
		self.name = name
		self.Sn = zeros((self.Nn, self.Nn)) #numpy.reshape(x[2:], (Nn, Nn))
		
	def phasespace_toyangle(self, theta1_toy, theta2_toy):
		J1_toy, J2_toy = self.gettoyaction(theta1_toy, theta2_toy)
		return self.toytorus.phasespace(J1_toy, J2_toy, theta1_toy, theta2_toy)	

	def gettoyaction(self, theta1_toy, theta2_toy, Sn=None):
		if Sn is None:
			Sn = self.Sn
		#J1_toy = J1 + 2 * n1  * Sn sum(n*Sn*cos(n*theta_toy)
		J1_toy = self.J1
		J2_toy = self.J2
		for n1i in range(self.Nn):
			for n2i in range(self.Nn):
				#sprint n1, n2, Sn.shape
				n1 = n1i*2
				n2 = n2i*2
				J1_toy += n1 * Sn[n1i,n2i] * cos(n1 * theta1_toy + n2 * theta2_toy)
				J2_toy += n2 * Sn[n1i,n2i] * cos(n1 * theta1_toy + n2 * theta2_toy)

		return J1_toy, J2_toy
		
	def _enerychisquare(self, x):
		M, b = x[:2]
		Sn = reshape(x[2:], (self.Nn, self.Nn))
		theta1d = arange(0, 2*pi, 2*pi/9)
		self.toytorus.set_mass_and_scale(M, b)
		Es = []
		for theta1_toy in theta1d:
			for theta2_toy in theta1d:
				J1_toy, J2_toy = self.gettoyaction(theta1_toy, theta2_toy, Sn)
				r,phi,vr,vphi,Ekin,Epot,E = self.toytorus.phasespace(J1_toy, J2_toy, theta1_toy, theta2_toy)
				Epot = self.galaxy.potentialr(r)
				E = Epot + Ekin
				Es.append(E)
				# dEdb = dEpotdr * drdb # no 
				# 
				# dEdb = dEdr * drdb + 
				# dEdM = dEdr * drdM
		Es = array(Es)
		#chisq = sum(((Es-mean(Es))/mean(Es))**2)/len(Es)
		#print Es
		chisq = sum(((Es-mean(Es)))**2)/len(Es)
		print chisq,mean(Es)
		return chisq
		
		
	def save(self):
		filename = os.path.join("torus_" +self.name + ".npy")
		params = [self.toytorus.p.M, self.toytorus.p.scale] + list(self.Sn.flat)
		numpy.save(filename, params)
		
	def load(self):
		filename = os.path.join("torus_" +self.name + ".npy")
		x = numpy.load(filename)
		M, b = x[:2]
		self.Sn = reshape(x[2:], (self.Nn, self.Nn))
		print "loaded", M,b, self.Sn
		self.toytorus.set_mass_and_scale(M, b)
		
	def load_or_fit(self):
		filename = os.path.join("torus_" +self.name + ".npy")
		if os.path.exists(filename):
			self.load()
		else:
			self.fit()
			self.save()
		
	
	def fit(self):
		params = [self.toytorus.p.M, self.toytorus.p.scale] + list(self.Sn.flat)
		res = fmin(self._enerychisquare, params, full_output=True, disp=True, maxiter=500)
		x = res[0]
		M, scale = x[:2]
		print M, scale
		self.Sn = x[2:]
		self.Sn = reshape(self.Sn, (self.Nn, self.Nn))
		print self.Sn
		
		

class TorusModelFit(object):
	def __init__(self, galaxy, toytorus):
		self.galaxy = galaxy
		self.toytorus = toytorus
		
	def createtorus(self, J1, J2, name):
		return TorusFit(self.galaxy, self.toytorus, J1, J2, name=name)
		

class TorusModelIsochrone(object):
	def __init__(self, M, scale, distance, beta):
		self.p = mab.gd.potential.Isochrone(M, scale)
		self.galaxy = mab.gd.schw.galaxy.Galaxy1C_constant_anisotropy(self.p, distance, beta)
		self.k = self.p.M * self.p.G
		self.b = self.p.scale
		
	def set_mass_and_scale(self, M, scale):
		self.p = mab.gd.potential.Isochrone(M, scale)
		self.galaxy = mab.gd.schw.galaxy.Galaxy1C_constant_anisotropy(self.p, self.galaxy.distance, self.galaxy.beta)
		self.k = self.p.M * self.p.G
		self.b = self.p.scale
	
		
	def hamiltonian(self, J1, J2):
		L = J2
		return -2 * self.k**2 / (2*J1 + L + sqrt(4*self.b*self.k+L**2) )**2 

	def dphasespacedb(self, J1, J2, theta1, theta2):
	
		L = J2
		H = self.hamiltonian(J1, J2)
		print "energy", H #, "original", E
		# step 2
		a = -self.k / (2 * H) - self.b
		# step 3
		e = sqrt(1+L**2/(2*H*a**2))
		#wl = sqrt(self.k) / (2*(a+self.b)**(3./2.)) * (1+L/sqrt(4*self.b*self.k+L**2))
		#print "a,e,wl", a, e, wl, a*e/(a+b)
		print "a,e", a, e
		# step4
		def f_psi(psi):
			return psi - a*e/(a+self.b)*sin(psi) - theta1
		def f_psi_prime(psi):
			return 1 - a*e/(a+self.b)*cos(psi)
		res = fsolve(f_psi, pi/2, fprime=f_psi_prime, full_output=0)
		psi = res
		print "psi", psi
		r = a * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a))
		
		# a = -self.k / (2 * H) - self.b
		# H = -2 * self.k**2 / (2*J1 + L + sqrt(4*self.b*self.k+L**2) )**2
		# k = self.p.M * self.p.G
		dHdb = 2*2 * self.k**2 / (2*J1 + L + sqrt(4*self.b*self.k+L**2) )**3 * 0.5 / sqrt(4*self.b*self.k+L**2) * 4 * self.k
		dadb = self.k/(2*H**2) * dHdb -1
		dkdM = self.p.G
		dHdk = H * 2 * self.k
		dHdM = dHdk * dkdM
		dadM = -dkdM / (2* H)

		if 0:
			a0 = a
			b0 = self.b
			e0 = e

			epsilon = 1e-5
			de = epsilon
			e2 = e + de
			r2 = a * sqrt((1-e2*cos(psi))*(1-e2*cos(psi)+2*self.b/a))
			print "drde", (r2-r)/de
			del e2
			
			
			dpsi = epsilon
			psi2 = psi + dpsi
			r2 = a * sqrt((1-e*cos(psi2))*(1-e*cos(psi2)+2*self.b/a))
			print "drdpsi", (r2-r)/dpsi
			del psi2
			
			
			da = epsilon
			a2 = a + da
			r2 = a2 * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a2))
			print "drda", (r2-r)/da, r, r2
			del a2
			
			
			da = epsilon
			a = a + da
			psi2 = fsolve(f_psi, pi/2, fprime=f_psi_prime, full_output=0)
			a = a0
			print "dpsida", (psi2-psi)/da
			

			db = epsilon
			self.b = self.b + db
			H2 = self.hamiltonian(J1, J2)
			a = -self.k / (2 * H2) - self.b
			e = sqrt(1+L**2/(2*H2*a**2))
			psi2 = fsolve(f_psi, pi/2, fprime=f_psi_prime, full_output=0)
			print "dpsidb", (psi2-psi)/db
			print "dHdb", (H2-H)/db, dHdb
			print "dadb", (a-a0)/db, dadb
		
		
			dpsidb = (psi2-psi)/db
			self.b = b0
			e = e0
			a = a0
			db = epsilon
			b2 = self.b + db
			r2 = a * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*b2/a))
			print "drdb_partial", (r2-r)/db
			del b2
			
			a2 = a + da
			e2 = sqrt(1+L**2/(2*H*a2**2))
			print "deda", (e2-e)/da
				
			print "="* 40
		
		drde = 0.5 * a / sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) * ( -cos(psi) - cos(psi)  +2*e*cos(psi)**2 - cos(psi)*2*self.b/a)
		print "drde", drde
		
		drdpsi = 0.5 * a / sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) * ( +e*sin(psi) +e*sin(psi) - e**2*2*cos(psi)*sin(psi) + 2*e*self.b/a*sin(psi) )
		print "drdpsi", drdpsi
		#drda = sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) + a / sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) * (1-e*cos(psi))*(-2*self.b/a**2+2/a*dadb)
		
		drda = sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) + 0.5 * a / sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) * (1-e*cos(psi))*(-2*self.b/a**2)
		print "drda", drda

		drdb_partial = 0.5 * a / sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a)) * (1-e*cos(psi))*(2/a)
		print "drdb_partial", drdb_partial
		#
		#dpsida = 1./dpsida
		
		# deda
		deda_partial = -1/sqrt(1+L**2/(2*H*a**2)) * L**2/(2*H*a**3)
		dedH_partial = -0.5/sqrt(1+L**2/(2*H*a**2)) * L**2/(2*H**2*a**2)
		dedb_partial = dedH_partial * dHdb
		deda = deda_partial + dedb_partial / dadb
		print "deda", deda
		
		dedb = deda * dadb #* 0
		
		# dedM
		dedM = dedH_partial * dHdM
		
		
		# 
		u = a * e / (a+self.b) 
		duda_partial = ((a+self.b)*e - a*e)/(a+self.b)**2
		dudb_partial = -u / (a+self.b)
		dude_partial = a / (a+self.b)
		dudb = duda_partial * dadb + dude_partial*dedb + dudb_partial
		dpsidu = (psi* cos(psi)-sin(psi))/(u*cos(psi)-1) - (cos(psi)*(u*psi*cos(psi)-theta1-u*sin(psi)))/(u*cos(psi)-1)**2 
		#dpsida = dpsidu * duda_partial
		dpsidb = dpsidu * dudb
		#dpsidb = dpsida * dadb
		#print "dpsida", dpsida
		print "dpsidb", dpsidb, dpsidu*dudb, dpsidu, dudb
		dpsidM = dpsidu * dude_partial * dedM
		
		
		drdb = drda * dadb + drdpsi * dpsidb + drde * dedb + drdb_partial
		
		drdM = drda * dadM + drdpsi * dpsidM + drde * dedM
		print drda * dadM, drda, dadM
		print drdpsi * dpsidM, drdpsi, dpsidM
		print drde * dedM, drde, dedM
		#print drda * dadb, drda, dadb
		#print drdpsi * dpsidb, drdpsi, dpsidb
		#print drde * dedb, drde, dedb
		#print drdb_partial
		#b = self.b
		#q = sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a))
		
		#print q * dadb + a * ((1+2*b/a-e*cos(psi)) * e * sin(psi) * dpsidb + (1-cos(psi)*e)*(2/a-2*b*dadb/a**2+e*sin(psi)*dpsidb)) / (2*q)
		return drdb, drdM

	
	def phasespace(self, J1, J2, theta1, theta2):
		# step 1
		L = J2
		H = self.hamiltonian(J1, J2)
		#print "energy", H, "original", E
		# step 2
		a = -self.k / (2 * H) - self.b
		# step 3
		e = sqrt(1+L**2/(2*H*a**2))
		wl = sqrt(self.k) / (2*(a+self.b)**(3./2.)) * (1+L/sqrt(4*self.b*self.k+L**2))
		#print "a,e,wl", a, e, wl, a*e/(a+b)
		# step4
		def f_psi(psi):
			return psi - a*e/(a+self.b)*sin(psi) - theta1
		def f_psi_prime(psi):
			return 1 - a*e/(a+self.b)*cos(psi)
		res = fsolve(f_psi, pi/2, fprime=f_psi_prime, full_output=0)
		psi = res
		# step 5 (r, chi)
		r = a * sqrt((1-e*cos(psi))*(1-e*cos(psi)+2*self.b/a))
		#print "check", psi, psi - a*e/(a+b)*sin(psi) - theta1
		wratio = 0.5 * (1+J2/sqrt(4*self.k*self.b+J2**2))
		Lambda_psi = arctan(sqrt((1+e)/(1-e))*tan(psi/2)) +\
				L/sqrt(L**2+4*self.b*self.k)*\
				arctan( sqrt((a*(1+e)+2*self.b)/(a*(1-e)+2*self.b))*tan(0.5*psi) )
		# step 6: theta
		phi = theta2 - wratio*theta1 + Lambda_psi
		# step 7: pr,ptheta
		pr = sqrt(self.k/(a+self.b)) * a*e*sin(psi)/r
		vr = pr
		pphi = L#/cos(phi)
		vphi = pphi/r
		#print "psi, phi", psi, phi
		#print "p", pr, pphi
		#print "r,rorg", r, rorg, 1/(r/rorg)
		#print "r,phi", r, phi
		#print "vr,vphi", vr, vphi
		Epot = self.galaxy.potentialr(r)
		Ekin = 0.5 * (vphi**2 + vr**2)
		return r, phi, vr, vphi, Ekin, Epot, Ekin+Epot
		#print "E", Epot, Ekin, Ekin+Epot, H
		#print "J1, J2", J1, J2