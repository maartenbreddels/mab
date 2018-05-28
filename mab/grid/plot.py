from kaplot import *

class PlotFunction(object):
	def __init__(self, f, grid):
		self.f = f
		self.grid = grid
		
		
	def run(self, args, opts, scope):
		n = 100
		x = arange(0, n) / (n-1.)
		y = arange(0, n) / (n-1.)
		xg, yg = numpy.meshgrid(x, y)
		I = self.f(xg, yg).T
		
		box()
		#indexedimage(I)
		probimage2d(I, 0, 1, x, y, drawcontourlines=True, fill=False, colormap=None)
		#draw()
		def f(x):
			y = log(self.f(x))
			if y < -300:
				return -300
			#if y 
			#3print y
			return y
		self.grid.rootNode.split(2)
		self.grid.evaluate(f)
		if 1:
			for i in [50,50,50,50,50,50,50,50]:
				self.grid.optimize(i, 8)
				self.grid.evaluate(f)
		
		#for x, y in zip(xg, yg):
		g = zeros_like(xg)
		if 1:
			for i in range(n):
				for j in range(n):
					g[i,j] = exp(self.grid.interpolate((yg[i,j], xg[i,j])))
		
				
		#print g.mean()
		probimage2d(g, 0, 1, x, y, drawcontourlines=True, fill=False, color="green", colormap=None)
		points = array(self.grid.points).T
		print "points", len(points[0])
		#print g.max()
		self.grid.calculate_volumes()
		print "integral:", self.grid.integrate()
		dy = dx = x[1] - x[0]
		print " ", sum(I*dx*dy)
		#print points[:10]
		scatter(points[0], points[1], symbolsize="10pt")
		
		draw()
		