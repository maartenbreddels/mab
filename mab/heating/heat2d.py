# -*- coding: utf-8 -*-
import PIL
import sys
from numpy import *
from kaplot import *

class Material(object):
	def __init__(self, color, conductivity, heat_capacity, density):
		self.color = color
		self.conductivity = conductivity
		self.heat_capacity = heat_capacity
		self.density = density
		"""
			conductivity=0.14 # W/(m*K)
			heat_capacity=1.7 # J/(g*K)
			density=800 # kg/m3
		
		"""
		#self.alpha = self.conductivity * 
		
		
class Grid(object):
	def __init__(self, image, materials, initial_T, gridsize=0.001): # gridsize in meters
		self.image = image
		self.materials = materials
		self.material_map = {}
		for material in self.materials:
			self.material_map[material.color] = material
		self.initial_T = initial_T
		self.gridsize = gridsize

	def run(self, args, opts, scope):
		self.image.load()
		#print "size", self.image.image.size
		width, height = self.image.image.size
		missing = False
		for color in self.image.colors:
			if color not in self.material_map:
				print "missing material for color", color
				missing = True
		if missing:
			sys.exit(-1)
		T = zeros((width, height))
		dT2 = zeros((width, height))
		dE = zeros((width, height)) # energy transport, J/m
		#material_index = zeros((width, height))
		alpha = zeros((width, height))
		beta = zeros((width, height))
		
		for x in range(width):
			for y in range(height):
				color = self.image.image.getpixel((x, y))[:3]
				material = self.material_map[color]
				#material.conductivity = 0.04 # * 5
				Rvalue = self.gridsize / material.conductivity # m2 * K / W
				area = self.gridsize # area is actually length
				alpha[x, y] = area / Rvalue  # W/K/m
				volume = self.gridsize ** 2 # volume is actually surface area, or see it as 'per meter'
				beta[x, y] = material.heat_capacity * material.density*1000. * volume # J/K/m
				T[x,y] = self.initial_T.get(material, 0)
				#print alpha[x, y], beta[x, y], color
				#import pd
				#print material
		T0 = T * 1.
		dt = 1.
		N = 150000*2
		
		for i in range(N):
		#if 1:
			dE = dE * 0
			if 1:
				dTrl = T[:-1,:] - T[1:,:]
				alphas = (alpha[:-1,:] + alpha[1:,:])/2
				dE[1:] += alphas * dTrl * dt
				dE[:-1] -= alphas * dTrl * dt
				#T[1:] += alphas * dTrl * dt / beta[1:]
				#T[:-1] += alphas * dTrl * dt / beta[:-1] 
				
				dTtb = T[:,:-1] - T[:,1:]
				alphas = (alpha[:,:-1] + alpha[:,1:])/2
				dE[:,1:] += alphas * dTtb * dt
				dE[:,:-1] -= alphas * dTtb * dt
				#T[:,1:] += alphas * dTtb * dt / beta[:,1:]
				#T[:,:-1] += alphas * dTtb * dt / beta[:,:-1] 
			
			if 0:
				for x in range(width):
					for y in range(height):
						if x > 0: # do left
							dT = T[x-1, y] - T[x, y]
							a = (alpha[x-1, y] + alpha[x, y])/2
							dE[x,y] += a * dT * dt
							#if dT > 0:
								
								#import pdb
								#pdb.set_trace()
						if x < width-1: # do right
							dT = T[x+1, y] - T[x, y]
							a = (alpha[x+1, y] + alpha[x, y])/2
							dE[x,y] += a * dT * dt
						if y > 0: # do bottom
							dT = T[x, y-1] - T[x, y]
							a = (alpha[x, y-1] + alpha[x, y])/2
							dE[x,y] += a * dT * dt
						if y < height-1: # do top
							dT = T[x, y+1] - T[x, y]
							a = (alpha[x, y+1] + alpha[x, y])/2
							dE[x,y] += a * dT * dt
			print "total energy", sum(dE)
			print sum(T)
			for x in range(width):
				for y in range(height):
					dT2[x, y] = dE[x, y] / beta[x,y]
					T[x, y] += dE[x, y] / beta[x,y]
					
				
		#box()
		xs = [27, 34, 50]
		mozaic(2,2,box)
		indexedimage(T.T, colormap="blackwhite")
		contour(T.T, 10., color="red")
		for color, x in zip(nicecolors, xs):
			vline(x, color=color)
		color = "red"
		hline(3.5, color=color)
		hline(7.5, color=color)
		hline(17.5, color=color)
		select(1,0)
		for color, x in zip(nicecolors, xs):
			graph(T[x,:], color=color)
			scatter(range(len(T[x,:])), T[x,:])
		vline(3.5)
		vline(7.5)
		vline(17.5)
		hline(10.0, color="red")
		hline(12.4, color="red")
		select(0, 1)
		indexedimage(T0.T)
		contour(T.T, 10.)
		contour(T.T, [10., 12.4], color="blue")
		hline(3.5, color=color)
		hline(7.5, color=color)
		hline(17.5, color=color)
		
		if 0:
			select(0,1)
			indexedimage(dE.T)
			for x in xs:
				vline(x)
			select(1,1)
			for x in xs:
				graph(dE[x,:])
			select(0,2)
			indexedimage(dT2.T)
			for x in xs:
				vline(x)
			select(1,2)
			for x in xs:
				graph(dT2[x,:])
		draw()
		
		
	
class Image(object):
	def __init__(self, filename):
		self.filename = filename
		
	def run(self, args, opts, scope):
		self.load()
		
	def load(self):
		self.image = PIL.Image.open(self.filename)
		self.colors = []
		for data in self.image.getcolors():
			self.colors.append(data[1][:3])
		for i, color in enumerate(self.colors):
			print "color", i, color
		
		