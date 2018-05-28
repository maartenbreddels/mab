# -*- coding: utf-8 -*-
from kaplot import *


class Plot4(object):
	def __init__(self, f):
		self.f = f
		
	def run(self, args, opts, scope):
		box()
		
		dp = 0.001
		p = arange(dp/2, 1+dp/2, dp)
		x = self.f.R(p)
		y = self.f(p)
		print x
		print y
		graph(x, y)
		#ylim(0, 1)
		draw()
		
		
class PlotEdgeworth(object):
	def __init__(self, f, Nsigma=4, N=100, ymax=None):
		self.f = f
		self.Nsigma = Nsigma
		self.N = N
		self.ymax = ymax
		
	def run(self, args, opts, scope):
		x = arange(self.N) / (self.N-1.) - 0.5
		x *= self.f.sigma * self.Nsigma * 2
		y = gaussian(x, 0, self.f.sigma)
		graph(x, y, color="red")
		print x
		
		y = self.f(x)
		print y
		graph(x, y)
		if self.ymax:
			ylim(0, self.ymax)
		
		
		draw()
		
		


class PlotEdgeworths(object):
	def __init__(self, fs, Nsigma=4, N=100000, ymax=None):
		self.fs = fs
		self.Nsigma = Nsigma
		self.N = N
		self.ymax = ymax
		
	def run(self, args, opts, scope):
		i = 0
		document()
		page(fontsize="15pt")
		mozaic(2,2,box)
		for f, color in zip(self.fs, nicecolors):
			x = arange(self.N) / (self.N-1.) - 0.5
			x *= f.sigma * self.Nsigma * 2
			dx = x[1] - x[0]
			#print x
			
			y = f(x)
			print "m0", sum(y*dx), f.moment(0) if hasattr(f, "moment") else 0.
			#print y
			select(0, 0)
			graph(x, y, color=color)
			select(1, 0)
			mask = y > 0
			graph(x[mask], log10(y[mask]), color=color, addlegend=False)
			ylim(-4, 0)
			if sum(~mask) > 0:
				graph(x[~mask], log10(-y[~mask]), color=color, addlegend=False, linestyle="dot")
			i += 1
			y = f(x)*x**2
			print "m2", sum(y*dx), f.moment(2) if hasattr(f, "moment") else 0.
			#import pdb;
			#pdb.set_trace()
			#print y
			select(0, 1)
			graph(x, y, color=color, addlegend=False, )
			y = f(x)*x**4
			print "m4", sum(y*dx), f.moment(4) if hasattr(f, "moment") else 0.
			#print y
			select(1, 1)
			graph(x, y, color=color, addlegend=False, )
			
		#y = gaussian(x, 0, f.sigma)
		#graph(x, y, color=nicecolors[i])
		select(0,0)
		alabels = [f.label for f in self.fs]
		#labels.append("normal") 
		autolegend(*alabels)
		select(1,0)
		labels("x", "log prob. dens")
		select(0,0)
		labels("x", "prob. dens")
		select(0,1)
		labels("x", "f(x)x<sup>2</sup>")
		select(1,1)
		labels("x", "f(x)x<sup>4</sup>")
		if self.ymax:
			ylim(0, self.ymax)
		
		draw()
		