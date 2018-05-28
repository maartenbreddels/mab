from __future__ import absolute_import
from numpy import *

def medianboxfilter2d(x, y, values, scale):
	assert len(x) == len(y) == len(values)
	values_filtered = zeros_like(values)
	for i in range(len(x)):
		xi = x[i]
		yi = y[i]
		mask = (x > xi - scale/2) & (x < xi + scale/2) & \
				(y > yi - scale/2) & (y < yi + scale/2)
		values_filtered[i] = median(values[mask])
	return values_filtered

def mediancirclefilter2d(x, y, values, scale):
	assert len(x) == len(y) == len(values)
	values_filtered = zeros_like(values)
	for i in range(len(x)):
		xi = x[i]
		yi = y[i]
		r = sqrt((xi-x)**2+(yi-y)**2)
		mask = r < scale
		values_filtered[i] = median(values[mask])
	return values_filtered
	
	
class MedianBoxFilter2d(object):
	def __init__(self, scale):
		self.scale = scale
		
	def filter(self, x, y, values):
		return medianboxfilter2d(x, y, values, self.scale)
	
class MedianCircleFilter2d(object):
	def __init__(self, scale):
		self.scale = scale
		
	def filter(self, x, y, values):
		return mediancirclefilter2d(x, y, values, self.scale)
