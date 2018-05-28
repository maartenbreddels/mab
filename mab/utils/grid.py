import numpy
import scipy
from numpy import *


class Gridder2d(object):
	def __init__(self, xmin, xmax, ymin, ymax, nx, ny):
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.nx = nx
		self.ny = ny
		self.resize = (xmin, ymin), (xmax, ymax)
		
	def meshgrid(self):
		x = (numpy.arange(self.nx) + 0.5) / (self.nx)
		y = (numpy.arange(self.ny) + 0.5) / (self.ny)
		x, y = numpy.meshgrid(x, y)
		return x, y
		
		
	def __call__(self, x, y, **kwargs):
		return self.histogram(x, y, **kwargs)
		
	def histogram(self, x, y, **kwargs):
		data, _, _ = scipy.histogram2d(x, y, bins=[self.nx, self.ny], range=[[self.xmin, self.xmax], [self.ymin, self.ymax]], **kwargs)
		return data
