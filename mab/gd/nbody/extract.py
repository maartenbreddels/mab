# -*- coding: utf-8 -*-
import glob
import mab.gd.logging as logging
from numpy import *
import os

logger = logging.getLogger("gd.nbody.gadget")

class TimeSeries(object):
	def __init__(self, filepattern, dt, snapshotclass, extractions=[]):
		self.filepattern = filepattern
		self.dt = dt
		self.snapshotclass = snapshotclass
		self.extractions = extractions
		
	def run(self, args, opts, scope):
		logger.info("filepattern: %s" % self.filepattern)
		filenames = glob.glob(self.filepattern)
		filenames.sort()
		t = 0
		for filename in filenames:
			snapshot = self.snapshotclass(filename=filename)
			snapshot.load()
			snapshot.center()
			for extraction in self.extractions:
				extraction.extract(t, snapshot)
			t += self.dt
			print filename
		for extraction in self.extractions:
			extraction.save()

class RadialPath(object):
	def __init__(self, dirname, name, componentname, index=0):
		self.filename = os.path.join(dirname, name+"_radial.npy")
		self.filename_t = os.path.join(dirname, name+"_radial_t.npy")
		self.componentname = componentname
		self.index = index
		self.times = []
		self.rs = []

	def extract(self, t, snapshot):
		component = snapshot.componentmap[self.componentname]
		x, y, z = q = component.q[:,self.index]
		r = sqrt(x**2+y**2+z**2)
		self.times.append(t)
		self.rs.append(r)
		#print q
		
	def save(self):
		self.times = array(self.times)
		self.rs = array(self.rs)
		save(self.filename, self.rs)
		save(self.filename_t, self.times)
		logger.info("saving to %s and %s" % (self.filename, self.filename_t))
		logger.debug("radii: %r" % (list(self.rs)))
		
	def load(self):
		self.rs = load(self.filename)
		self.times = load(self.filename_t)
		
		
	
