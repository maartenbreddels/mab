from numpy import *
import numpy
import os

class SchwHelper(object):
	def __init__(self):
		pass
	
	@classmethod
	def getrs(cls, galaxy, dither=1, dE=False, physical=False):
		logr1, logr2 = cls.logE1, cls.logE2
		kpc_to_arcsec = galaxy.kpc_to_arcsec(1.)
		if physical:
			logrmin, logrmax = logr1-log10(kpc_to_arcsec), logr2-log10(kpc_to_arcsec)
		else: 
			logrmin, logrmax = logr1, logr2 
		nE = cls.nI1
		logrs = arange(nE*dither, dtype=float) / (nE*dither-1) * (logrmax-logrmin) + logrmin
		rs = 10**logrs
		return rs
	
	@classmethod
	def getrborders(cls, galaxy, dither=1, dE=False, physical=False):
		logr1, logr2 = cls.logE1, cls.logE2
		kpc_to_arcsec = galaxy.kpc_to_arcsec(1.)
		if physical:
			logrmin, logrmax = logr1-log10(kpc_to_arcsec), logr2-log10(kpc_to_arcsec)
		else: 
			logrmin, logrmax = logr1, logr2 
		nE = cls.nI1
		logrs = (arange(nE*dither+1, dtype=float) -0.5) / (nE*dither-1) * (logrmax-logrmin) + logrmin
		rs = 10**logrs
		return rs
	
	@classmethod
	def getEs(cls, galaxy, dither=1, dE=False):
		logr1, logr2 = cls.logE1, cls.logE2
		kpc_to_arcsec = galaxy.kpc_to_arcsec(1.)
		
		logrmin, logrmax = logr1-log10(kpc_to_arcsec), logr2-log10(kpc_to_arcsec) 
		nE = cls.nI1
		logrs = arange(nE*dither, dtype=float) / (nE*dither-1) * (logrmax-logrmin) + logrmin
		rs = 10**logrs
		Es = galaxy.potentialr(rs)
		if dE:
			dEs = concatenate( ([Es[1] - Es[0]], (Es[2:] - Es[0:-2])/2, [Es[-1] - Es[-2]]) )
			return Es, dEs
		else: 
			return Es
	
	
	@staticmethod
	def index_to_orbitnr(i1, i2, i3):
		return i1*nI2*nI3 + i2*nI3 + i3

class SchwSolution(object):
	def __init__(self, dirname, n_moments, n_constraints, modelpath, weightname="", fitted=True, addLz=0):
		self.dirname = dirname
		if not fitted:
			filename = os.path.join(modelpath, "orbitweights" +weightname +".npy")
			#orbitweights = ravel(numpy.load(filename))
			orbitweights = array(numpy.load(filename).flat)
			if addLz:
				allorbitweights = zeros((len(orbitweights), addLz))
				for i in range(addLz/2):
					allorbitweights[:,i] = orbitweights 
					allorbitweights[:,i+addLz/2+1] = orbitweights 
				allorbitweights[:,addLz/2] = 1*orbitweights
				orbitweights = ravel(allorbitweights)/(addLz)
		else:
			filename = os.path.join(dirname, "orbitweights" +weightname +".npy")
			#orbitweights = ravel(numpy.load(filename))
			orbitweights = array(numpy.load(filename).flat)
		
		filename = os.path.join(dirname, "projectedmoments.npy")
		projectedmoments = array(memmap(filename, dtype='float64', mode='readonly', shape=(len(orbitweights), n_moments, n_constraints)))
		self.orblibmoments = projectedmoments 
		#projectedmoments = load()
		mask = projectedmoments[:,0,:] > 0
		#print projectedmoments.shape 
		#print mask.shape 
		#mask = mask[:,newaxis,:]
		#print mask.shape
		f = (10000*5*5*25)
		projectedmoments /= f
		#for i in range(1, projectedmoments.shape[1]):  
		#	projectedmoments[:,i,:][mask] /= projectedmoments[:,0,:][mask]
		self.projectedmoments = tensordot(projectedmoments, orbitweights, axes=([0], [0]))
		
		densities = array(projectedmoments[:,0,:])
		self.rho2d = sum(orbitweights * transpose(densities), axis=1)

		