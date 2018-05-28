# -*- coding: utf-8 -*-
from kaplot import *
from mab.utils.numpy import range_and_borders

#class FunctionAdapeter(object):
#	def __init__(self, object, methodname,  


class Series(object):
	def __init__(self, series, xattr, yattr, colors=[], labels=[], xlim=None, ylim=None):
		self.series = series
		self.xattr = xattr
		self.yattr = yattr
		self.colors = colors
		self.labels = labels
		self.xlim = xlim
		self.ylim = ylim
		
		
	def run(self, args, opts, scope):
		box()
		index = 0
		for series in self.series:
			if hasattr(series, "load"):
				series.load()
			x = getattr(series, self.xattr)
			y = getattr(series, self.yattr)
			print len(x)
			kwargs = {}
			if index < len(self.colors):
				kwargs["color"] = self.colors[index]
			graph(x, y, **kwargs)
			index += 1
		if self.labels:
			autolegend(*self.labels)
		if self.xlim:
			xlim(*self.xlim)
		if self.ylim:
			ylim(*self.ylim)
		draw()

class Functions(object):
	def __init__(self, xmin, xmax, N=100, functions=[], logx=False, logy=False, colors=nicecolors, ymin=None, ymax=None):
		self.functions = functions
		self.xmax = xmax
		self.xmin = xmin
		self.N = N
		self.functions = functions
		self.colors = colors
		self.logx = logx
		self.logy = logy
		self.ymin = ymin
		self.ymax = ymax
		
	def run(self, args, opts, scope):
		self.x, self.x_borders =  range_and_borders(self.N, self.xmin, self.xmax)
		box()
		self.plot(scope)
		draw()
		
	def plot(self, scope):
		u = 10**self.x if self.logx else self.x
		for function, color in zip(self.functions, self.colors):
			y = function(u)
			if self.logy:
				y = log10(y)
			print y
			mask = ~(isnan(y) | isinf(y))
			graph(self.x[mask], y[mask], color=color)

		ylim(self.ymin, self.ymax)
		
		
		
class Bars(object):
	def __init__(self, values, labels, groupnames, **kwargs):
		self.values = values
		self.kwargs = kwargs
		self.labels = labels
		self.groupnames = groupnames
		
		
	def run(self, args, opts, scope):
		#mozaic(2,2,box)
		document()
		page(fontsize="20pt", symbolsize="16pt")
		box()
		barchart(self.values, **self.kwargs)
		ax = current.container.xaxis
		def label(x):
			index = int(x)
			print x, index
			return self.labels[index]
		ax.label = label
		ax.start = 0
		ax.interval = 1
		ax.subticks = 0
		#ax.
		ylabel("relative evidence")
		if 0:
			import pdb
			pdb.set_trace()
		contexts = self.kwargs["contexts"]
		count = len(contexts)
		objects = [kaplot.objects.PlotObject(current.container, **context) for context in contexts]
		legend(["squaresolid"] * count, self.groupnames, objects)
		grow(top=1.4)
		draw()
		
		