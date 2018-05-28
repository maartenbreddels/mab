from numpy import *
"""
Grids are build up from 1d grids.
Each 1d Grid has the following members:
	* umin, umax: the original space' min and max (e.g. logrmin, logrmax)
	* xmin, xmax: the target space' min and max  (eg 10**logrmin, 10**logrmax)
	* x, x_borders: array with values between xmin and xmax (excluding them), and their border (including them)
	* u, u_borders:                           umin and umax
	* N: number of points between xmin and xmax (same length as x and u array)
"""
class Grid1dAperture(object):
	def __init__(self, aperture, name="logr", label="nolabel"):
		aperture.load()
		self.xmin = aperture.aperture_rcenters_kpc[0]
		self.xmax = aperture.aperture_rcenters_kpc[-1]
		self.aperture = aperture
		self.label = label
		self.name = name
		
		self.centers = self.aperture.aperture_rcenters_kpc
		self.u_centers = self.centers
		self.borders = self.aperture.aperture_rborders_kpc
		self.u_borders = self.borders
		
		self.N = len(self.centers)
			
		
class Grid1dLog(object):
	def __init__(self, logxmin, logxmax, N, name="logr", label="nolabel"):
		self.umin = logxmin
		self.umax = logxmax
		self.logxmin = logxmin
		self.logxmax = logxmax
		self.xmax = 10**logxmin
		self.xmax = 10**logxmax
		self.label = label
		self.name = name
		self.N = N
		
		
		uniform = (arange(N)+0.5) / N
		self.u = self.logx = uniform * (logxmax - logxmin) + logxmin
		self.x = 10**self.logx
		
		uniform = (arange(N+1.)) / (N)
		logx = uniform * (logxmax - logxmin) + logxmin
		self.x_borders = 10**logx
		self.u_borders = logx
		
	def subdivide(self, M):
		return Grid1dLog(self.logxmin, self.logxmax, self.N * M, name=self.name, label=self.label)
	
class Grid1dLinear(object):
	def __init__(self, xmin=0., xmax=1., N=10, name="r", label="nolabel"):
		self.xmin = xmin
		self.xmax = xmax
		self.label = label
		
		u = (arange(N)+0.5) / N
		self.u_centers = self.centers = u * (xmax - xmin) + xmin
		
		u = (arange(N+1.)) / (N)
		self.u_borders = self.borders = u * (xmax - xmin) + xmin
		
		self.name = name
		self.N = N
		
	def subdivide(self, M):
		return Grid1dLinear(self.xmin, self.xmax, self.N * M, name=self.name, label=self.label)
		
class Grid2d(object):
	def __init__(self, gridx, gridy):
		self.gridx = gridx
		self.gridy = gridy
		self.shape = (self.gridx.N, self.gridy.N)
		
	def subdivide(self, Mx, My):
		return Grid2d(self.gridx.subdivide(Mx), self.gridy.subdivide(My))
	
class Grid2dDithered(Grid2d):
	def __init__(self, gridx, gridy, ditherx=2, dithery=2):
		Grid2d.__init__(self, gridx, gridy)
		self.subgrid = self.subdivide(ditherx, dithery)
		
		
class Grid3d(object):
	def __init__(self, gridx, gridy, gridz):
		self.gridx = gridx
		self.gridy = gridy
		self.gridz = gridz
		self.shape = (self.gridx.N, self.gridy.N, self.gridz.N)
		
	def subdivide(self, Mx, My, Mz):
		return Grid3d(self.gridx.subdivide(Mx), self.gridy.subdivide(My), self.gridz.subdivide(Mz))
	
class Grid3dDithered(Grid3d):
	def __init__(self, gridx, gridy, gridz, ditherx=2, dithery=2, ditherz=2):
		Grid3d.__init__(self, gridx, gridy, gridz)
		self.subgrid = self.subdivide(ditherx, dithery, ditherz)
