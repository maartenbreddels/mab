from mab.constants import *
from mab.astrounits  import *
from numpy import *
from kaplot import *
import mab.parallelize
import mab.gd.logging as logging
import scipy.ndimage
import mab.cvsfile
import scipy.optimize.nnls

kpc_to_km = (1*KPC).asNumber(KM)
s_to_gyr = (S/GYR).asNumber()

logger = logging.getLogger("gd.schw.uv")


class UVGridder(object):
	def __init__(self, xmin, xmax, ymin, ymax, nx, ny):
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.nx = nx
		self.ny = ny
		self.resize = (xmin, ymin), (xmax, ymax)
		
	def __call__(self, x, y):
		return self.histogram(x, y)
		
	def histogram(self, x, y):
		data, _, _ = histogram2d(x, y, bins=[self.nx, self.ny], range=[[self.xmin, self.xmax], [self.ymin, self.ymax]])
		return data

class UVData(object):
	def __init__(self, filename, filters=[]):
		self.filename = filename
		self.filters = filters
		
	#def run(self, args, opts, scope):
	def load(self):
		stars = mab.cvsfile.readcsv(self.filename)
		for filter in self.filters:
			stars = filter(stars)
		print len(stars)
		return stars

class UVPlane(object):
	def __init__(self, filename, uv_data, uv_gridder, sigma=1.):
		self.filename = filename
		self.uv_data = uv_data
		self.uv_gridder = uv_gridder
		self.sigma = sigma
		
	def load(self):
		self.plane = load(self.filename)
	
	def run(self, args, opts, scope):
		stars = self.uv_data.load()
		print stars.Vr.min(), stars.Vr.mean(), stars.Vr.max()
		print stars.Vphi.min(), stars.Vphi.mean(), stars.Vphi.max()
		data = self.uv_gridder(stars.Vr, stars.Vphi)
		logger.info("writing UV plane to: %s" % self.filename)
		numpy.save(self.filename, data)
		box()
		I = data.T
		I = scipy.ndimage.gaussian_filter(I, self.sigma)
		indexedimage(I, resize=self.uv_gridder.resize, colormap="whiteblack")
		I /= I.max()
		print 
		for value in [0.9, 0.8, 0.7, 0.6, 0.5, 0.2, 0.1, 0.05]:
			level = value#10**(value) 
			contour(I, [level], resize=self.uv_gridder.resize, color="red")
		flipx()
		flipy()
		draw()
		
		
class UVFit(object):
	def __init__(self, uv_plane, uv_model, sigma=1.):
		self.uv_plane = uv_plane
		self.uv_model = uv_model
		self.sigma = sigma
		
	def run(self, args, opts, scope):
		self.fit()
		
	def fit(self):
		self.uv_model.load()
		self.uv_plane.load()
		uv_planes = self.uv_model.uv_planes
		uv_gridder = self.uv_model.uv_gridder
		dfgrid = self.uv_model.dfgrid
		
		uv_planes = uv_planes.reshape((dfgrid.nx * dfgrid.ny, uv_gridder.nx * uv_gridder.ny))
		target_uv_plane = self.uv_plane.plane.reshape((uv_gridder.nx * uv_gridder.ny))
		target_uv_plane = scipy.ndimage.gaussian_filter(target_uv_plane, self.sigma)
		for i in range(len(uv_planes)):
			uv_planes[i] =  scipy.ndimage.gaussian_filter(uv_planes[i], self.sigma)
		print uv_planes.shape, target_uv_plane.shape
		solution, error = scipy.optimize.nnls(uv_planes.T, target_uv_plane)
		print error
		print solution
		print solution.shape
		uv = dot(solution, uv_planes)
		solution = solution.T.reshape(dfgrid.nx, dfgrid.ny)
		#box()
		mozaic(2,2,box)
		indexedimage(solution.T, resize=dfgrid.resize, colormap="whiteblackred")
		
		uv = uv.reshape(uv_gridder.nx, uv_gridder.ny)
		labels("vr", "vphi")
		select(1,0)
		indexedimage(uv.T, resize=uv_gridder.resize, colormap="whiteblackred")
		labels("vr", "vphi")
		select(1,1)
		indexedimage(scipy.ndimage.gaussian_filter(self.uv_plane.plane.T, self.sigma), resize=uv_gridder.resize, colormap="whiteblackred")
		labels("vr", "vphi")
		
		
		draw()
		
		

class UVPlot(object):
	def __init__(self, uv_model, sigma=0):
		self.uv_model = uv_model
		self.sigma = sigma
		
	def run(self, args, opts, scope):
		self.uv_model.load()
		
		dfgrid = self.uv_model.dfgrid
		document(size="25cm,25cm")
		mozaic(dfgrid.nx, dfgrid.ny)
		for i in range(dfgrid.nx):
			for j in range(dfgrid.ny):
				select(i, j)
				border()
				spacer("1mm")
				#indexedimage(self.uv_model.density_planes[i,j])
				I = self.uv_model.uv_planes[i,j].T
				I = scipy.ndimage.gaussian_filter(I, self.sigma)
				vr = self.uv_model.dfgrid.x[i]
				vphi = self.uv_model.dfgrid.y[j]
				I /= I.max()
				I = log10(I)
				I[I<-4] = -4
				
				indexedimage(I, resize=self.uv_model.uv_gridder.resize, colormap="whiteblackred")
				symbol(vr, vphi, color="blue", symbolsize="30pt")
		draw()

class DensityPlot(object):
	def __init__(self, uv_model, sigma=0):
		self.uv_model = uv_model
		self.sigma = sigma
		
	def run(self, args, opts, scope):
		self.uv_model.load()
		
		dfgrid = self.uv_model.dfgrid
		document(size="25cm,25cm")
		mozaic(dfgrid.nx, dfgrid.ny)
		for i in range(dfgrid.nx):
			for j in range(dfgrid.ny):
				select(i, j)
				border()
				spacer("1mm")
				#indexedimage(self.uv_model.density_planes[i,j])
				I = self.uv_model.density_planes[i,j].T
				I = scipy.ndimage.gaussian_filter(I, self.sigma)
				indexedimage(I, colormap="whiteblackred")
		draw()

class UVModel(object):
	def __init__(self, dfgrid, orbit_integrator, r0, dr, phi0, dphi, filename, uv_gridder, density_gridder):
		self.dfgrid = dfgrid
		self.orbit_integrator = orbit_integrator
		self.r0 = r0
		self.dr = dr
		self.phi0 = phi0
		self.dphi = dphi
		self.filename = filename
		self.uv_gridder = uv_gridder
		self.density_gridder = density_gridder
		
	def run(self, args, opts, scope):
		
		self.uv_planes = mab.utils.numpy.mmapzeros((self.dfgrid.nx, self.dfgrid.ny, self.uv_gridder.nx, self.uv_gridder.ny))
		self.density_planes = mab.utils.numpy.mmapzeros((self.dfgrid.nx, self.dfgrid.ny, self.density_gridder.nx, self.density_gridder.ny))
		
		@mab.parallelize.parallelize(cores=opts.cores, info=opts.progressbar)
		def do(i, j):
			dither1 = self.dfgrid.dither
			dither2 = self.dfgrid.dither
			dither = dither1*dither2
			
			Ninput = dither1*dither2 #self.dfgrid.subgrid.n
			#N = self.orbitintegrator.orbital_points
			q0 = zeros((Ninput, 2))
			p0 = zeros((Ninput, 2))
			x0, y0 = self.r0 * cos(self.phi0), self.r0 * sin(self.phi0)
			#print s_to_gyr
			dt = ones(Ninput) * 1e14#/100
			# / s_to_gyr
			for d1 in range(dither1):
				for d2 in range(dither2):
					#d = d1*dither2+d2
					q0[d1*dither2+d2,:] = (x0, y0)
					#print i, self.dfgrid.subgrid.x[i]
					vr = self.dfgrid.subgrid.x[i*dither1+d1]
					vphi = self.dfgrid.subgrid.y[j*dither2+d2]
					vx = (vr * x0  - vphi * y0)/self.r0
					vy = (vphi * x0 + vr * y0)/self.r0
					#vr*r = (x*vx + y*vy)
					#vt*r = (x*vy - y*vx)
					
					p0[d1*dither2+d2,:] = (vx, vy)
					#print p0[d1*dither2+d2]
					
			q, p = self.orbit_integrator.integrate(dt, q0, p0) #, xout, yout, vxout, vyout)
			#print p
			#print q
			#print p.shape, q.shape
			density = zeros((100, 100))
			uvdensity = zeros((100, 100))
			#box()
			plot = False
			if plot:
				mozaic(2,1,box)
			vphi0 = 220
			for d1 in range(dither1):
				for d2 in range(dither2):
					d = d1*dither2 + d2
					x = q[0,d]
					y = q[1,d]
					vx = p[0,d]
					vy = p[1,d]
					
					r = numpy.sqrt(x**2+y**2)
					phi = (arctan2(y, x) + 2 * pi) % (2*pi)
					vr = (x*vx + y*vy)/r
					vt = (x*vy - y*vx)/r
					
					#density_d, _, _ = histogram2d(x, y, bins=100, range=[[-10,10], [-10, 10]])
					density_d = self.density_gridder(x, y)
					density += density_d
					self.density_planes[i,j] += density_d 
					
					#print vr, vt
					mask = (r >= self.r0 - self.dr/2) & (r < self.r0 + self.dr/2) &\
							(phi >= self.phi0 - self.dphi/2) & (phi < self.phi0 + self.dphi/2)
					scatter(vr[mask], vt[mask])
					s = 75
					#uvdensity_d, _, _ = histogram2d(vr[mask], vphi[mask], bins=100, range=[[-s, s], [vphi0-s,vphi0+s]])
					uvdensity_d = self.uv_gridder(vr[mask], vt[mask])
					uvdensity += uvdensity_d
					self.uv_planes[i,j] += uvdensity_d
					if plot:
						select(0, 0)
						scatter(x, y)
					#print len(x[mask])
					if plot:
						select(1,0)
						scatter(vr[mask], vt[mask])
						select(0,0)
						scatter(x[mask], y[mask])
						squarelim(-10,10)
			
			if plot:
				select(1,0)
				indexedimage(uvdensity.T, resize=self.uv_gridder.resize, colormap="whiteblackred")
				vx0, vy0 = self.dfgrid.x[i], self.dfgrid.y[j]
				symbol(vx0, vy0, color="blue", symbolsize="20pt")
				#xlim(-150, 150)
				#y#lim(220-150, 220+150)
				select(0, 0)
				#indexedimage(density.T, resize=((-10, -10), (10,10)))
				select(1,0)
				draw()
		#i = range(self.dfgrid.nx)
		#do(i)
		#do(1,15)
		#return
		i1s = [i1 for i1 in range(self.dfgrid.nx) for i2 in range(self.dfgrid.ny)]
		i2s = [i2 for i1 in range(self.dfgrid.nx) for i2 in range(self.dfgrid.ny)]
		do(i1s, i2s)
		
		filename = self.filename + "_uv_plane.npy"
		logger.info("saving uv planes to: %s" % filename)
		numpy.save(filename, self.uv_planes)
		
		filename = self.filename + "_density_plane.npy"
		logger.info("saving density planes to: %s" % filename)
		numpy.save(filename, self.density_planes)
		
	def load(self):
		filename = self.filename + "_uv_plane.npy"
		logger.info("loading uv planes from: %s" % filename)
		self.uv_planes = numpy.load(filename)
		
		filename = self.filename + "_density_plane.npy"
		logger.info("loading density planes from: %s" % filename)
		self.density_planes = numpy.load(filename)
		
		