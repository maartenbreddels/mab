import pyublas
from numpy import *
import numpy
import os
from mab.cvsfile import readcsv, writecsv
#from mab.gd.orbitreader import *
from scipy.integrate import dblquad

class SchwSet(object):
	def __init__(self, modelpath, setname):
		self.modelpath = modelpath
		self.setname = setname
		self.setdir = os.path.join(modelpath, "schw", setname)
		filename = os.path.join(self.setdir, "parameters.csv")
		self.parameters = readcsv(filename, noeval="name")
		
	def writeparameters(self, filename):
		filename = os.path.join(self.setdir, filename)
		writecsv(filename, self.parameters)
		
	def do(self, f, dirname=None, ignoreerror=False):
		if dirname:
			if dirname[-1] == "/":
				dirname = dirname[:-1]
		for i, name in enumerate(self.parameters.name):
			try:
				dirname_i = os.path.join(self.setdir, name)
				#print dirname_i, dirname 
				if (dirname is None) or (dirname == dirname_i):
					f(self.modelpath, dirname_i, name, self.parameters[i]) 
			except Exception, e:
				if ignoreerror:
					print "skipping", name, e
				else:
					raise
		
class NNLSoutput(object):
	def __init__(self, modelpath, setname, name, nnprefix="nn"):
		self.modelpath = modelpath
		self.name = name
		#self.setdir = os.path.join(modelpath, "schw", setname)
		self.filename = os.path.join(self.modelpath, "schw", setname, self.name, "datfil", nnprefix+"_nnls.out")
		lines = file(self.filename).readlines()
		self.redchi = float(lines[-2].split()[-1])
		self.chisq = float(lines[-2].split()[-3])
		self.chisqkin = float(lines[-1].split()[-3])
		#print "% 30s % 10.3f % 10.3f" %(name, chisq, redchi)
		
		
		
		
class IntrinsicMoments(object):
	def __init__(self, modelpath, dirname, nnprefix="nn", sphericalvelocities=True):
		self.modelpath = modelpath
		self.dirname = dirname
		self.filename = os.path.join(self.dirname, "datfil", nnprefix+"_intrinsic_moments.out")
		lines = file(self.filename).readlines()
		
		self.phiborders = array([float(k) for k in lines[3].split()])
		self.thetaborders = array([float(k) for k in lines[5].split()])
		
		self.phis = (self.phiborders[:-1] + self.phiborders[1:])/2
		self.thetas = (self.thetaborders[:-1] + self.thetaborders[1:])/2
		extralines = 0
		while not "phi" in lines[7+extralines+1]:
			extralines += 1
			#print lines[7+extralines]
		#print lines[7:7+extralines+1]
		self.rborders = array([float(k) for k in ("".join(lines[7:7+extralines+1])).split()])
		#self.rs = 10**( (log10(self.rborders[:-1]) + log10(self.rborders[1:]))/2 )
		self.rs = (self.rborders[:-1] + self.rborders[1:])/2
		#radii /= kpc_to_arcsec
		
		self.momentgrid = zeros((12, len(self.phis), len(self.thetas), len(self.rs))) # 3 velocities, 3 variances and 3 covariances, the order is xx, yy, zz, xy, yz, zx
		lines = lines[10+extralines+1:]
		for line in lines:
			values = line.split()
			phi, theta, r = [int(k) for k in values[0:3]]
			x, y, z , vx, vy, vz, xx, yy, zz, xy, yz, zx = [float(k) for k in values[6:9+9]]
			self.momentgrid[0, phi-1, theta-1, r-1] = vx
			self.momentgrid[1, phi-1, theta-1, r-1] = vy
			self.momentgrid[2, phi-1, theta-1, r-1] = vz
			self.momentgrid[3, phi-1, theta-1, r-1] = xx
			self.momentgrid[4, phi-1, theta-1, r-1] = yy
			self.momentgrid[5, phi-1, theta-1, r-1] = zz
			self.momentgrid[6, phi-1, theta-1, r-1] = xy
			self.momentgrid[7, phi-1, theta-1, r-1] = yz
			self.momentgrid[8, phi-1, theta-1, r-1] = zx
			
			self.momentgrid[9, phi-1, theta-1, r-1] = x
			self.momentgrid[10, phi-1, theta-1, r-1] = y
			self.momentgrid[11, phi-1, theta-1, r-1] = z
			
		self.varvphi_avg = self.momentgrid[0,0,0]*0.0
		self.varvtheta_avg = self.varvphi_avg * 0 
		self.varvr_avg = self.varvphi_avg * 0
		self.vr_avg = self.varvphi_avg * 0
		self.vphi_avg = self.varvphi_avg * 0
		self.vtheta_avg = self.varvphi_avg * 0
		self.varvrs = self.momentgrid[3] *0
		self.varvphis = self.momentgrid[3] *0
		self.varvthetas = self.momentgrid[3] *0
		
		N = 0
		i1 = 0
		i2 = 3
		if sphericalvelocities:
			for phii in range(len(self.phis)): #[i1:i1+1]:
				for thetai in range(len(self.thetas)): #[i2:i2+1]:#[2:-1]:
					vr = (self.momentgrid[0, phii, thetai])
					vphi = (self.momentgrid[1, phii, thetai])
					vtheta = (self.momentgrid[2, phii, thetai])
					varvr = (self.momentgrid[3, phii, thetai])
					varvphi = (self.momentgrid[4, phii, thetai])
					varvtheta = (self.momentgrid[5, phii, thetai])
					self.varvrs[phii, thetai] = varvr
					self.varvphis[phii, thetai] = varvphi
					self.varvthetas[phii, thetai] = varvtheta
					#print varvr
					#print varvphi
					#print varvtheta
					self.varvr_avg += varvr
					self.varvphi_avg += varvphi
					self.varvtheta_avg += varvtheta
					self.vr_avg += vr
					self.vphi_avg += vphi
					self.vtheta_avg += vtheta
					N += 1
			
			self.varvr_avg /= N #len(self.phis) * len(self.thetas)
			self.varvphi_avg /= N #len(self.phis) * len(self.thetas)
			self.varvtheta_avg /= N #len(self.phis) * len(self.thetas)
			self.vr_avg /= N
			self.vtheta_avg /= N
			self.vphi_avg /= N
		else:
			for phii in range(len(self.phis)): #[i1:i1+1]:
				for thetai in range(len(self.thetas)): #[i2:i2+1]:#[2:-1]:
					phi = self.phis[phii]
					theta = self.thetas[thetai]
					
					R1 = reshape(array([cos(phi) * sin(theta), -sin(phi), cos(phi)*cos(theta)]), (3,1))
					R2 = reshape(array([sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta)]), (3,1))
					R3 = reshape(array([cos(theta), 0.0, -sin(theta)]), (3,1))
					
					if 1:
						R1 = reshape(array([cos(phi) * sin(theta), sin(phi)*sin(theta), cos(theta)]), (3,1))
						R2 = reshape(array([-sin(phi), cos(phi), 0.0]), (3,1))
						R3 = reshape(array([cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta)]), (3,1))
						
					R1 = array([cos(phi) * sin(theta), -sin(phi), cos(phi)*cos(theta)])
					R2 = array([sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta)])
					R3 = array([cos(theta), 0.0, -sin(theta)])
						
					
					
					# vr, vphi, vtheta = vx * R1 + vy * R2 + vz * R3
					def T(x):
						#return reshape(x, (1, len(x)))
						return x
					vx = T(self.momentgrid[0, phii, thetai])
					vy = T(self.momentgrid[1, phii, thetai])
					vz = T(self.momentgrid[2, phii, thetai])
					varvx = T(self.momentgrid[3, phii, thetai])
					varvy = T(self.momentgrid[4, phii, thetai])
					varvz = T(self.momentgrid[5, phii, thetai])
					varvxvy = T(self.momentgrid[6, phii, thetai])
					varvyvz = T(self.momentgrid[7, phii, thetai])
					varvzvx = T(self.momentgrid[8, phii, thetai])
					if 0:
						varvr, varvphi, varvtheta = varvzvx* 0, varvzvx * 0, varvzvx * 0
						
						for i in range(len(self.rs)):
							#x = self.momentgrid[9, phii, thetai, i]
							#y = self.momentgrid[10, phii, thetai, i]
							#z = self.momentgrid[11, phii, thetai, i]
							#r = sqrt(x**2+y**2+z**2)
							#theta = arccos(z/r)
							#phi = arctan(y/x) 
							#R1 = array([cos(phi) * sin(theta), -sin(phi), cos(phi)*cos(theta)])
							#R2 = array([sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta)])
							#R3 = array([cos(theta), 0.0, -sin(theta)])
							A = transpose(array([R1, R2, R3]))
							X = array([[varvx[i], varvxvy[i], varvzvx[i]], [varvxvy[i], varvy[i], varvyvz[i]], [varvzvx[i], varvyvz[i], varvz[i]]])
							#S1 = tensordot(A, X, axes=([1], [0]))
							#S = tensordot(A, S1, axes=([0], [0]))
							S1 = tensordot(A, X, axes=([1], [0]))
							S = tensordot(S1, A, axes=([1], [1]))
							varvr[i], varvphi[i], varvtheta[i] = S[0,0], S[1,1], S[2,2]
					else:
						A = transpose(array([R1, R2, R3]))
						X = array([[varvx, varvxvy, varvzvx], [varvxvy, varvy, varvyvz], [varvzvx, varvyvz, varvz]])
						S1 = tensordot(A, X, axes=([1], [0]))
						S = tensordot(A, S1, axes=([1], [1]))
						varvr, varvphi, varvtheta = S[0,0], S[1,1], S[2,2]
						
					#A = transpose(array([R1, R2, R3]))
					#X = array([[varvx, varvxvy, varvzvx], [varvxvy, varvy, varvyvz], [varvzvx, varvyvz, varvz]])
					
					#S1 = tensordot(A, X, axes=([1], [0]))
					#S = tensordot(A, S1, axes=([0], [0]))
					#import pdb; pdb.set_trace()
					
					
					#import pdb; pdb.set_trace()
					#varvr, varvphi, varvtheta = varvx * R1**2 + R2**2*varvy + R3**2*varvz + 2*R1*R2*varvxvy + 2*R2*R3*varvyvz + 2*R3*R1*varvzvx
					#varvr, varvphi, varvtheta = S[0,0], S[1,1], S[2,2]
					#print S[:,:,6]
					#vr, vphi, vtheta = R1*vx + R2*vy + R3*vz
					
					vr, vphi, vtheta = dot(A, [vx, vy, vz])
					
					self.varvrs[phii, thetai] = varvr
					self.varvphis[phii, thetai] = varvphi
					self.varvthetas[phii, thetai] = varvtheta
					#print varvr
					#print varvphi
					#print varvtheta
					self.varvr_avg += varvr
					self.varvphi_avg += varvphi
					self.varvtheta_avg += varvtheta
					self.vr_avg += vr
					self.vphi_avg += vphi
					self.vtheta_avg += vtheta
					N += 1
			
			self.varvr_avg /= N #len(self.phis) * len(self.thetas)
			self.varvphi_avg /= N #len(self.phis) * len(self.thetas)
			self.varvtheta_avg /= N #len(self.phis) * len(self.thetas)
			self.vr_avg /= N
			self.vtheta_avg /= N
			self.vphi_avg /= N
		
		
class OrbitWeights(object):
	def __init__(self, dirname, schwmodel, nnprefix="nn"):
		#self.modelpath = modelpath
		#self.name = name
		#self.setdir = os.path.join(modelpath, "schw", setname)
		self.filename = os.path.join(dirname, "datfil", nnprefix+"_orb.out")
		#orbfilename = os.path.join(self.dirname, "datfil", "nn_orb.out")
		lines = file(self.filename).readlines()[1:]
		self.orbitweights = ones((schwmodel.nI1, schwmodel.nI2)) * 0.0
		
		for line in lines:
			values = [float(v) for v in line.split()]
			#print values
			nr, E, I2, I3, _, _, weight = values
			#if E > 0 and I2 > 0 and I3 > 0:
			self.orbitweights[abs(E)-1, abs(I2)-1] += weight
			#grid2IM[E-1, I2-1] = weight
		#return grid2IM

class Configuration(object):
	def __init__(self, modelpath, dirname):
		self.modelpath = modelpath
		self.dirname = dirname
		self.orblibfilename = os.path.join(self.dirname, "datfil", "orblib.dat")
		self.orbitreader = Orblib(self.orblibfilename)
		self.orbitreader.read()
		self.noI1 = self.orbitreader.noI1
		self.noI2 = self.orbitreader.noI2
		self.noI3 = self.orbitreader.noI3
		self.noI = self.noI1*self.noI2*self.noI3
		self.noConstraints = self.orbitreader.noConstraints
		self.noMaxVelocityHistograms = self.orbitreader.noMaxVelocityHistograms

	def nnls(self, A, y):
		import scipy.optimize
		x, err = scipy.optimize.nnls(A, y)
		return x, err
	
		
	def readorbitmatrix(self, aperture=None):
		orbitmatrix = zeros((self.noI1*self.noI2*self.noI3, self.noConstraints, self.noMaxVelocityHistograms))
		#print orbitmatrix.shape
		self.orbitreader.fillorbit(orbitmatrix)
		#areas = zeros(self.noConstraints)
		if aperture:
			for c in range(1, self.noConstraints+1):
				area = aperture.area(c)
				#orbitmatrix[:,c-1,:] /= area
				#areas[c-1] = area
		A = reshape(transpose(orbitmatrix), (self.noMaxVelocityHistograms*self.noConstraints, self.noI))
		return A
	
	def createtargetvector(self, aperture):
		m = self.createtargetmatrix(aperture)
		y = reshape(transpose(m), (self.noMaxVelocityHistograms*self.noConstraints))
		return y
	
	def readtrueorbitweightsvector(self, weightname=""):
		return ravel(self.readtrueorbitweights(weightname))
	
	def readtrueorbitweightskinvector(self, weightname=""):
		return ravel(self.readtrueorbitweightskin(weightname))
	
	def readtrueorbitweights(self, weightname=""):
		filename = os.path.join(self.modelpath, "orbitweights" +weightname +".npy")
		orbitweights = numpy.load(filename)
		return orbitweights
	
	def readtrueorbitweightskin(self, weightname=""):
		filename = os.path.join(self.modelpath, "orbitweights" +weightname +"_kin.npy")
		orbitweights = numpy.load(filename)
		return orbitweights
	
	def loglikelihoodbycount(self, pvector, counts):
		mask = pvector > 0
		return sum(log10(pvector[mask])*counts[mask])
		
		
	def createtargetmatrix(self, aperture):
		targetmatrix = zeros((self.noConstraints, self.noMaxVelocityHistograms))
		stars2d = readcsv(os.path.join(self.modelpath, "vlos.csv"))
		for star in stars2d:
			eta = star.eta * 60 * 60
			xi = star.xi * 60 * 60
			contraint = aperture.contraint(eta, xi)
			vlos = star.vlos
			vindex = int(vlos / self.orbitreader.dvhist+int(self.noMaxVelocityHistograms/2))
			#print contraint, vlos, vindex, orbitreader.dvhist
			if (vindex > 0) and (vindex < self.noMaxVelocityHistograms):
				targetmatrix[contraint-1, vindex] += 1
				#print vindex
				#assert vindex <
			else:
				print vindex, contraint
		print "targetmatrix", sum(targetmatrix)
		return targetmatrix
			
	def readfittedorbitweightsvector(self):
		return ravel(self.readfittedorbitweights())
	
	def readfittedorbitweights(self):
		orbfilename = os.path.join(self.dirname, "datfil", "nn_orb.out")
		lines = file(orbfilename).readlines()[1:]
		grid2IM = ones((self.noI1,self.noI2)) * 0.0
		
		for line in lines:
			values = [float(v) for v in line.split()]
			#print values
			nr, E, I2, I3, _, _, weight = values
			#if E > 0 and I2 > 0 and I3 > 0:
			grid2IM[abs(E)-1, abs(I2)-1] += weight
			#grid2IM[E-1, I2-1] = weight
		return grid2IM
		
	def solutionvector2moments(self, y, moments=4):
		m = reshape(y, (self.noMaxVelocityHistograms, self.noConstraints))
		#momentlist = []
		grid = zeros((moments, self.noConstraints))
		v = (arange(self.noMaxVelocityHistograms)-(self.noMaxVelocityHistograms-1)/2) * self.orbitreader.dvhist
		#print v
		#int(vlos / self.orbitreader.dvhist+int(self.noMaxVelocityHistograms/2))
		for j in range(self.noConstraints):
			for i in range(moments):
				grid[i,j] = sum(v**(i+1)*m[:,j])/sum(m[:,j])
		return grid
		
		
		
		

class Aperture(object):
	def __init__(self, dirname, masses=False, schwpostfix=""):
		self.dirname = dirname
		aperturefilename = os.path.join(self.dirname, "aperture" +schwpostfix+ ".dat")
		#print dirname, aperturefilename
		if os.path.exists(aperturefilename):
			binsfilename = os.path.join(self.dirname, "bins" +schwpostfix+ ".dat")
		else:
			aperturefilename = os.path.join(self.dirname, "infil", "aperture" +schwpostfix+ ".dat")
			binsfilename = os.path.join(self.dirname, "infil", "bins" +schwpostfix+ ".dat")
		lines = file(binsfilename).readlines()
		noBins = int(lines[1])
		#self.bins = zeros(noBins, dtype=int)
		self.bins = array([int(k) for k in (" ".join(lines[2:])).split()])
		self.min = min(self.bins)
		self.max = max(self.bins)
		
		lines = file(aperturefilename).readlines()
		#print lines
		self.x0, self.y0 = [float(k) for k in lines[1].strip().split()]
		self.width, self.height = [float(k) for k in lines[2].strip().split()]
		self.gridx, self.gridy = [int(k) for k in lines[4].strip().split()]
		self.bins.resize(self.gridx, self.gridy)
		
		if masses:
			massfilename = os.path.join(self.dirname, "datfil", "mass_aper.dat")
			lines = file(massfilename).readlines()
			self.masses = zeros(int(lines[0].strip()))
			for line in lines[1:]:
				index, mass = line.strip().split()
				index, mass = int(index), float(mass)
				self.masses[index-1] = mass
			

	def constraint(self, x, y):
		xgrid = int((x+self.x0)/self.width * self.gridx)
		ygrid = int((y+self.y0)/self.height * self.gridy)
		return self.bins[xgrid, ygrid]
	
	def mass(self, galaxy):
		#print self.width/2+self.x0
		assert abs(self.width/2  +self.x0) < 1e-5
		assert abs(self.height/2 +self.y0) < 1e-5
		centerx = ((arange(self.gridx) + 0.5) - self.gridx/2.)/(self.gridx) * self.width
		centery = ((arange(self.gridy) + 0.5) - self.gridy/2.)/(self.gridy) * self.height
		x, y = meshgrid(centerx, centery)
		r = sqrt(x**2+y**2)
		masses = []
		dA = galaxy.arcsec_to_kpc(self.width/self.gridx) * galaxy.arcsec_to_kpc(self.height/self.gridy)
		for i in range(self.max):
			mask = self.bins == (i+1)
			rs = ravel(r[mask])
			masses.append(dA * sum(galaxy.stellar_profile.densityR(galaxy.arcsec_to_kpc(rs))))
		return masses
	
	def mass_precise(self, galaxy):
		#print self.width/2+self.x0
		assert abs(self.width/2  +self.x0) < 1e-5
		assert abs(self.height/2 +self.y0) < 1e-5
		#centerx = ((arange(self.gridx) + 0.5) - self.gridx/2.)/(self.gridx) * self.width
		#centery = ((arange(self.gridy) + 0.5) - self.gridy/2.)/(self.gridy) * self.height
		#x, y = meshgrid(centerx, centery)
		#r = sqrt(x**2+y**2)
		#masses = []
		#dA = galaxy.arcsec_to_kpc(self.width/self.gridx) * galaxy.arcsec_to_kpc(self.height/self.gridy)
		J = galaxy.arcsec_to_kpc(1) * galaxy.arcsec_to_kpc(1)
		#for i in range(self.max):
		#	mask = self.bins == (i+1)
		#	rs = ravel(r[mask])
		#	masses.append(dA * sum(galaxy.stellar_profile.densityR(galaxy.arcsec_to_kpc(rs))))
		masses = zeros(self.max)
		for xi in range(self.gridx):
			#print xi, "out of", self.gridx
			for yi in range(self.gridy):
				constraint = self.bins[xi, yi]
				if constraint > 0: 
					x0 = self.x0 + xi * self.width/self.gridx
					y0 = self.y0 + yi * self.height/self.gridy
					x1 = self.x0 + (xi+1) * self.width/self.gridx
					y1 = self.y0 + (yi+1) * self.height/self.gridy
					def f(x, y):
						r = sqrt(x*x+y*y)
						return galaxy.stellar_profile.densityR(galaxy.arcsec_to_kpc(r))
					mass, error = dblquad(f, x0, x1, lambda x: y0, lambda x: y1)
					masses[constraint-1] += mass/J
		return masses
	

	def area(self, contraint):
		return sum(self.bins == (contraint))
	
	def areas(self):
		return array([self.area(c+1) for c in range(max(self.bins))])

