import mab.cvsfile
import mab.gd.logging as logging
import numpy

import os
import math
from numpy import *
import random as pyrandom
import mab.random
import mab.parallelize

logger = logging.getLogger("gd.galaxysampling")

class GalaxySampling_SphericalNonRotating(object):
	def __init__(self, **kwargs):
		pass
	
def export_csv(filename, records, seperator=","):
	#import pdb
	#pdb.set_trace()
	#N = len(records.dtype)
	f = open(filename, "w")
	#for name in records.dtype.names:
	print >>f, seperator.join(records.dtype.names)
	for record in records:
		values = [repr(value) for value in record]
		print >>f, seperator.join(values)
	
class DataSamplingMock(object):
	def __init__(self, modelpath, light_model, input, output_npy, N, velocity_error=0, seed=0, rmax=None, output_csv=None, multi_sets=None, multi_sets_seed=None):
		self.modelpath = modelpath
		self.light_model = light_model
		self.input = input
		self.output_npy = output_npy
		self.output_csv = output_csv
		self.N = N
		self.velocity_error = velocity_error
		self.seed = seed
		self.rmax = rmax
		self.multi_sets_seed = multi_sets_seed
		self.multi_sets = multi_sets
		
		self.dirname = os.path.join(modelpath, "data")
		
		names = "eta,xi,rc,Rdeg,vlos_true,vlos,e_vlos".split(",")
		names = [k.strip() for k in names]
		types = ["f8"] * len(names)
		# concatenate the values to the intrinsic values
		self.dtype = self.input.dtype + zip(names, types)
		
	def run(self, args, opts, scope):
		self.init()
		self.sample(scope)
		self.save()
		
	def init(self):
		self.input.load()
		
	def sample(self, scope):
		if self.multi_sets:
			self.stars_sets = mab.utils.numpy.mmapzeros((self.multi_sets, self.N), dtype=self.dtype).view(recarray)
		else:
			self.stars_sets = zeros((1, self.N), dtype=self.dtype).view(recarray)
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def wrapper(index):
			#print "set", index
			seed = self.multi_sets_seed(self.seed, index)
			random.seed(seed)
			pyrandom.seed(seed)
			numpy.random.seed(seed)
			self.sample_set(self.stars_sets[index])
		if self.multi_sets:
			indices = range(self.multi_sets)
			wrapper(indices)
		else:
			random.seed(self.seed)
			pyrandom.seed(self.seed)
			numpy.random.seed(self.seed)
			self.sample_set(self.stars_sets[0])
			
	def save(self):
		
		filename = os.path.join(self.dirname, self.output_npy)
		logger.info("output line of sight velocities(npy): %s" % filename)
		save(filename, self.stars_sets)
		if self.output_csv:
			if self.multi_sets:
				logger.warning("multisets and csv not supported, output to mutiple files")
				for i in range(self.multi_sets):
					filename = os.path.join(self.dirname, self.output_csv+"_%03i.csv" % i )
					logger.info("output line of sight velocities(csv): %s" % filename)
					export_csv(filename, self.stars_sets[i])
			else:
				filename = os.path.join(self.dirname, self.output_csv)
				logger.info("output line of sight velocities(csv): %s" % filename)
				export_csv(filename, self.stars_sets[0])
			
			
		#mab.cvsfile.writecsv(outputfilename, stars)
		
		if 0:
			outputfilename = os.path.join(modelpath, opts.obsname+".asc")
			print "output line of sight velocities (asc):", outputfilename
			indices = numpy.argsort(stars.Rdeg)
			f = file(outputfilename, "w")
			for i in indices:
				i = int(i)
				print >>f, stars[i].vlos, stars[i].e_vlos
			f.close()
			#writecsv(outputfilename, stars)
		
		if 0:
			stars = allstars[:opts.Nphoto]
			
			for star in stars:
				star.eta = degrees(atan(star.z/self.galaxy.distance))
				star.xi = degrees(atan(star.y/self.galaxy.distance))
				star.R = numpy.sqrt(star.eta**2 + star.xi**2)
				star.attributes = star.attributes + ["eta", "xi", "R"]   	
			outputfilename = os.path.join(modelpath, "photo.csv")
			print "output photometry:", outputfilename
			writecsv(outputfilename, stars)
			
	def sample_set(self, stars):
		#print (self.multi_sets, self.N)
		#stars = zeros((self.multi_sets, self.N), dtype=self.dtype).view(recarray)
		
		#stars = mab.cvsfile.readcsv(os.path.join(self.modelpath, self.input))
		#stars = self.input.df_samples
		
		samples = self.input.df_samples
		if self.rmax:
			r2d = numpy.sqrt(samples.y**2+samples.z**2)
			samples = samples[r2d<self.rmax]
			
		rotation_matrix = mab.random.randomSO3()
		
		samples.x, samples.y, samples.z = tensordot(rotation_matrix, [samples.x, samples.y, samples.z], axes=[(1,), (0,)])
		samples.vx, samples.vy, samples.vz = tensordot(rotation_matrix, [samples.vx, samples.vy, samples.vz], axes=[(1,), (0,)])
		#print x.shape
		#print z/self.light_model.distance
		#print numpy.degrees(numpy.arctan(z/self.light_model.distance)) * 60 **2
		#dsadsa
		
		samples.eta = numpy.degrees(numpy.arctan(samples.z/self.light_model.distance)) * 60 **2
		samples.xi = numpy.degrees(numpy.arctan(samples.y/self.light_model.distance)) * 60 **2
		
		#allstars = stars.clone()
		#print len(stars)
		#stars = stars[:opts.Nvel]
		#random.seed(self.seed)
		samples = self._picksample(samples)
		logger.info("%d samples" % len(samples))
		
		stars.eta = numpy.degrees(numpy.arctan(samples.z/self.light_model.distance)) * 60 **2
		stars.xi = numpy.degrees(numpy.arctan(samples.y/self.light_model.distance)) * 60 **2
		stars.rc = numpy.sqrt(stars.eta**2 + stars.xi**2)
		stars.Rdeg = numpy.sqrt(stars.eta**2 + stars.xi**2) / 60**2
		
		for i in range(len(samples)):
			e_vlos = self.velocity_error # km/s
			#print self.input.dtype
			for name in samples.dtype.names:
				stars[name][i] = samples[name][i]
			stars.vlos_true[i] = samples.vx[i]
			if e_vlos > 0.0:
				stars.vlos[i] = stars.vlos_true[i] + numpy.random.normal(0, e_vlos)
			else:
				stars.vlos[i] = stars.vlos_true[i]
			stars.e_vlos[i] = e_vlos
			#star.attributes = star.attributes + ["eta", "xi", "rc", "Rdeg", "vlos_true", "vlos", "e_vlos"]
		
			
	
	def _picksample(self, stars):
		#stars = stars.filter(lambda star: sqrt(star.y**2+star.z**2)  < self.rmax)
		mask = sqrt(stars.y**2+stars.z**2)  < self.rmax
		total = sum(mask)
		starsindices = pyrandom.sample(range(total), self.N)
		return stars[mask][starsindices]

			
class DataSamplingMockBiased(DataSamplingMock):
	def __init__(self, modelpath, light_model, input, output_npy, sampling, N, velocity_error=0, seed=0, rmax=None, output_csv=None, multi_sets=None, multi_sets_seed=None):
		DataSamplingMock.__init__(self, modelpath, light_model, input, output_npy, N, velocity_error=velocity_error, seed=seed, rmax=rmax, output_csv=output_csv, multi_sets=multi_sets, multi_sets_seed=multi_sets_seed)
		self.sampling = sampling
		logger.info("biased sampling")
		
	def _picksample(self, stars):
		source_stars = mab.cvsfile.readcsv(os.path.join(self.modelpath, "data", self.sampling))
		indices = []
		if len(source_stars) < self.N:
			logger.info("source population is smaller than target, using nearest sample without duplicates")
			for j in range(self.N):
				i = random.randint(0, len(source_stars)-1)
				xi = source_stars.xi[i]
				eta = source_stars.eta[i]
				#R = source_stars.Rdeg[i]
				#nearest_index = numpy.argsort(abs(stars.Rdeg-R))
				distances = numpy.sqrt((xi-stars.xi)**2 + (eta-stars.eta)**2)
				nearest_index = numpy.argsort(distances)
				for k in nearest_index:
					if k not in indices:
						indices.append(k)
						break
		else:
			source_indices = pyrandom.sample(range(len(source_stars)), self.N)
			for i in source_indices:
				#print i, len(source_indices)
				xi = source_stars.xi[i]
				eta = source_stars.eta[i]
				distances = numpy.sqrt((xi-stars.xi)**2 + (eta-stars.eta)**2)
				nearest_index = numpy.argmin(distances)
				#R = source_stars.Rdeg[i]
				#nearest_index = argmin(abs(stars.Rdeg-R))
				indices.append(nearest_index)
		#print len(stars)
		#print indices
		return stars[indices]
			
		
