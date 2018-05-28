# -*- coding: utf-8 -*-
import mab.asciifile

class HaloList(object):
	def __init__(self, filename):
		self.filename = filename
		
	def load(self):
		self.halos = mab.asciifile.readsimple(self.filename)
		#print self.stars
		return self.halos
		
	def run(self, args, opts, scope):
		self.load()
		for halo in self.halos:
			print halo.halo, halo.treepos, halo.rho0, halo.rs
		
		