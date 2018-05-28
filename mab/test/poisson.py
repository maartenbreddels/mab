# -*- coding: utf-8 -*-
from kaplot import *

import scipy


class Light(object):
	def __init__(self, light_distribution, N=100000, M=50, Nr=100, rmax=1.5):
		self.light_distribution = light_distribution
		self.rmax = rmax
		self.Nr = Nr
		self.N = N
		self.M = M
		
		
	def run(self, args, opts, scope):
		print scipy.integrate.quad(lambda r: self.light_distribution.densityr(r, M=1) * r**2 * 4 * pi, 0, 1.5)
		print scipy.integrate.quad(lambda R: self.light_distribution.densityR(R, M=1) * R**1 * 2 * pi, 0, 1.5)
		dsa
		counts = []
		for i in range(self.M):
			rs = self.light_distribution.sample_r(N=self.N, rmax=1.5)
			count, bins = numpy.histogram(rs, self.Nr, [0, self.rmax])
			counts.append(count)
		counts = array(counts)
		count = mean(counts, axis=0)
		box()
		histogramline(bins, count)
		ylim(0, max(count)*1.1)
		centers = (bins[0:-1] + bins[1:])/2
		print centers
		errorbars(centers, count, yerr=sqrt(count))
		errorbars(centers, count, yerr=var(counts, axis=0)**0.5, color="red")
		
		
		
		draw()
		