# -*- coding: utf-8 -*-
from kaplot import *
#from scipy.interpolate import griddata
#from scipy.interpolate import griddata




class TestInter(object):
	def __init__(self):
		pass
	
	
	def run(self, args, opts, scope):
		def func(x, y):
			return x*(1-x)*cos(4*pi*x) * sin(4*pi*y**2)**2
		grid_x, grid_y = mgrid[0:1:100j, 0:1:200j]
		points = random.rand(1000, 2)
		values = func(points[:,0], points[:,1])
		grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
		indexedimage(grid_z2)
		draw()

class TestGaussian(object):
	def __init__(self, tree):
		self.tree = tree
		
	def run(self, args, opts, scope):
		print self.tree.dimension
		self.tree.init()
		self.tree.read()
		colormap = "whiterainbow"
		if 1:
			sigmax, sigmay = 0.1, 0.02
			
			def f(*point):
				x, y = point
				xp = x - 0.5
				yp = y - 0.5
				angle = radians(40)
				xr = xp * cos(angle) + yp * sin(angle)
				yr = xp * sin(angle) - yp * cos(angle)
				z = gaussian(xr, 0.0, sigmax) *  gaussian(yr, 0.0, sigmay)
				return log(z)
				
			self.tree.update(f)
			self.tree.write_solved()
			g = exp(self.tree.grid(100))
			dx = 1./100
			print g
			
			print g.shape
			mozaic(2,2,box)
			
			indexedimage(g.T, colormap=colormap, resize=[(0, 0), (1,1)])
			select(1,0)
			x = arange(100)/(100.)
			y = g.sum(axis=0) * dx
			graph(x, y / y.max())
			x = arange(0, 1, 0.001)
			y = gaussian(x, 0.5, sigmax)
			print "TEST", sum(y*0.001)
			y /= y.max()
			graph(x, y, color="red")
			
			grid_x, grid_y = mgrid[0:1:100j, 0:1:200j]
			x = arange(100)/(100.)
			y = exp(f(grid_x, grid_y)).sum(axis=1) * dx
			y /= y.max()
			graph(x, y, color="blue")
			
			select(0, 0)
			points = array(self.tree.points)
			x, y = points.T
			scatter(x, y)
			
			
			select(1,1)
			def f(x, y):
				return \
					1 * (1-x) * (1-y) +\
					0 * (x) * (1-y) + \
					1 * (x) * (y) + \
					0 * (1-x) * (y)
			#x, y = mgrid[0:1:0.01, 0:1:0.01]x
			dx = 0.01
			x = arange(0, 1+dx/2, dx)
			y = arange(0, 1+dx/2, dx)
			x, y = numpy.meshgrid(x, y)
			z = f(x, y)
			print z.shape
			indexedimage(z.T, colormap=colormap)
			draw()
	
class TestGrid(object):
	def __init__(self, tree):
		self.tree = tree
		
	def run(self, args, opts, scope):
		print self.tree.dimension
		self.tree.init()
		self.tree.read()
		colormap = "whiterainbow"
		if 1:
			def f(*point):
				print point
				if point == (0.5, 0.5):
					print "yep"
					return 1.
				elif point == (1., 1.):
					print "yep"
					return 0.
				else:
					return 0.
			self.tree.update(f)
			self.tree.write_solved()
			g = 10**self.tree.grid(100)
			print g
			
			print g.shape
			mozaic(2,2,box)
			
			indexedimage(g.T, colormap=colormap)
			select(1,0)
			graph(g.sum(axis=0))
			
			
			select(1,1)
			def f(x, y):
				return \
					1 * (1-x) * (1-y) +\
					0 * (x) * (1-y) + \
					1 * (x) * (y) + \
					0 * (1-x) * (y)
			#x, y = mgrid[0:1:0.01, 0:1:0.01]x
			dx = 0.01
			x = arange(0, 1+dx/2, dx)
			y = arange(0, 1+dx/2, dx)
			x, y = numpy.meshgrid(x, y)
			z = f(x, y)
			print z.shape
			indexedimage(z.T, colormap=colormap)
			
			
		else:
			def f(*point):
				print point
				if point == (0.5, 0.5, 0.5):
					print "yep"
					return 1.
				elif point == (0.75, 0.75, 0.75):
					print "yep"
					return 1.
				else:
					return 0.
			self.tree.update(f)
			self.tree.write_solved()
			g = 10**self.tree.grid(50)
			#print g
			
			print g.shape
			mozaic(3,3,box)
			
			select(0,0)
			indexedimage(g.sum(0).T, colormap=colormap)
			select(0,1)
			indexedimage(g.sum(1).T, colormap=colormap)
			select(0,2)
			indexedimage(g.sum(2).T, colormap=colormap)

			select(1,0)
			graph(g.sum(axis=0).sum(axis=0))
			select(1,1)
			graph(g.sum(axis=0).sum(axis=1))
			select(1,2)
			graph(g.sum(axis=1).sum(axis=1))
			
		draw()