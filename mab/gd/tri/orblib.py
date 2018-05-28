# -*- coding: utf-8 -*-
import mab
from numpy import *
class Orblib(object):
	def __init__(self, filename):
		self.filename = filename
		
		
	def load(self):
		self.orblib = mab.gd.gdfast.Orblib(self.filename)
		self.orblib.read()
		self.noOrbits = self.orblib.noI1 * self.orblib.noI2 * self.orblib.noI3
		self.grid = zeros((self.noOrbits, self.orblib.noConstraints, self.orblib.noMaxVelocityHistograms))
		self.orblib.fillorbit(self.grid)
		
		