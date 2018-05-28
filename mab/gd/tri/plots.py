# -*- coding: utf-8 -*-
from kaplot import *
import mab

class PlotOrblib(object): # abstract
	def __init__(self, orblib):
		self.orblib = orblib
		
	def draw(self, axis=None, aperture=None):
		nx = self.orblib.orblib.noI1
		if axis is None:
			ny = self.orblib.orblib.noI2*self.orblib.orblib.noI3
		else:
			ny = self.orblib.orblib.noI2#*self.orblib.orblib.noI3
		grid = self.orblib.grid.reshape((self.orblib.orblib.noI1, self.orblib.orblib.noI2, self.orblib.orblib.noI3, self.orblib.orblib.noConstraints, self.orblib.orblib.noMaxVelocityHistograms))
		mozaic(nx, ny, container)
		colormap = "blackwhite"
		if aperture is not None:
			area = [sum(aperture.aperture == k+1) for k in range(grid.shape[3])]
			print "area", area
			area = array(area)
		else:
			area = None
		
		for i in range(nx):
			for j in range(ny):
				select(i, j)
				if axis is None:
					indexedimage(grid[i,j/self.orblib.orblib.noI2, j%self.orblib.orblib.noI2].T, colormap=colormap)
				else:
					border()
					spacer()
					for k in range(self.orblib.orblib.noI3):
						mass = grid[i,j, k].sum(axis=axis)
						if area is not None:
							mass = mass/area
						graph(mass, color=nicecolors[k])
		
		#draw()

class PlotDensity(PlotOrblib):
	def __init__(self, orblib, aperture):
		PlotOrblib.__init__(self, orblib)
		self.aperture = aperture

	def run(self, args, opts, scope):
		self.orblib.load()
		self.aperture.load()
		self.draw(axis=1, aperture=self.aperture)
		draw()

class PlotVelocity(PlotOrblib):
	def __init__(self, orblib):
		PlotOrblib.__init__(self, orblib)

	def run(self, args, opts, scope):
		self.orblib.load()
		self.draw(axis=0)
		draw()

class PlotOrbits(PlotOrblib):
	def __init__(self, orblib):
		PlotOrblib.__init__(self, orblib)

	def run(self, args, opts, scope):
		self.orblib.load()
		self.draw(axis=None)
		draw()
		
class PlotDensity2D(PlotOrblib):
	def __init__(self, orblib, aperture):
		PlotOrblib.__init__(self, orblib)
		self.aperture = aperture

	def run(self, args, opts, scope):
		self.orblib.load()
		self.aperture.load()
		nx = self.orblib.orblib.noI1
		ny = self.orblib.orblib.noI2*self.orblib.orblib.noI3 #/2
		grid = self.orblib.grid.reshape((self.orblib.orblib.noI1, self.orblib.orblib.noI2, self.orblib.orblib.noI3, self.orblib.orblib.noConstraints, self.orblib.orblib.noMaxVelocityHistograms))
		document(size="15cm,25cm")
		mozaic(nx, ny, container)
		colormap = "blackwhite"
		if 0:
			if aperture is not None:
				area = [sum(aperture.aperture == k+1) for k in range(grid.shape[3])]
				print "area", area
				area = array(area)
			else:
				area = None
		print "%" * 70
		for i in range(nx):
			for j in range(ny):
				select(i, j)
				if isinstance(self.aperture, mab.gd.tri.aperture.ApertureNative):
					density = grid[i,j/self.orblib.orblib.noI2, j%self.orblib.orblib.noI2].sum(axis=1)
					density2d = density.reshape((self.aperture.Nx, self.aperture.Ny))
				else:
					density = grid[i,j/self.orblib.orblib.noI2, j%self.orblib.orblib.noI2].sum(axis=1)
					density2d = zeros_like(self.aperture.aperture)
					N = self.aperture.aperture.max()
					print "max", N
					for k in range(N):
						#import pdb
						#pdb.set_trace()
						density2d[self.aperture.aperture==k+1] = density[k]
						if (i==3) and (j==3):
							print self.aperture.aperture==k+1
							print density[k]
				indexedimage(density2d**0.25, colormap=colormap)
		draw()

		