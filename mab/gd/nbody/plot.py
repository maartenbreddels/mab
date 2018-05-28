# -*- coding: utf-8 -*-
from kaplot import *
import mab
from mab.binningtools import bingrid, binrange

class Movie(object):
	def __init__(self, gadgetdirname, snapshotclass, plot, name="frame"):
		self.snapshotclass = snapshotclass
		self.gadgetdirname = gadgetdirname
		self.plot = plot
		self.name = name
	
	def run(self, args, opts, scope):
		snapshots = []
		for i in range(1000):
			filename = os.path.join(self.gadgetdirname, "output", "snapshot_%03d.hdf5" % i)
			if os.path.exists(filename):
				print filename
				snapshot = self.snapshotclass(filename)
				snapshot.nr = i
				snapshots.append(snapshot)
		snapshots[0].load()
		for snapshot in snapshots:
			snapshot.total = 500#len(snapshots)
		def do(i, dohardcopy=True):
			snapshot = snapshots[i]
			filename = os.path.join(self.gadgetdirname, self.name +"_%03d.gif" % snapshot.nr)
			if os.path.exists(filename):
				print "skipping", filename
			else:
				snapshot.load()
				self.plot.snapshot = snapshot
				self.plot.snapshot0 = snapshots[0]
				clear()
				self.plot.plot()
				if dohardcopy:
					#filename = os.path.join(self.gadgetdirname, self.name +"_%03d.gif" % snapshot.nr)
					hardcopy(filename)
				#draw()
		#time.sleep(0.5)
		if opts.iteration:
			do(opts.iteration, False)
			draw()
		else:
			do = mab.parallelize.parallelize(cores=opts.cores)(do)
			do(range(len(snapshots)))
class PlotXYZ(object):
	def __init__(self, snapshot, orbit=None, galaxy=None, dt=None):
		self.snapshot = snapshot
		self.galaxy = galaxy
		self.snapshot0 = None
		self.orbit = orbit
		self.dt = dt
		
	def run(self, args, opts, scope):
		self.snapshot.load()
		self.plot()
		draw()
		
	def plot(self):
		s = 200
		document(size="35cm,25cm")
		hsplit(container, 2, [0.97, 0.03])
		info_container = select(1)
		select(0)
		##mozaic(1,2,box)
		mozaic(6,4,box)
		#print current.container_array
		#print current.container_array[1]
		current.container_array = current.container_array[0][1]
		#print current.container_array
		i = 0
		colormap = {"disk":"red", "bulge":"blue", "halo":"black"}
		for name in "halo bulge disk".split():
			if name in self.snapshot.componentmap:
				component = self.snapshot.componentmap[name]
				color = colormap[name]
				if len(component.q.flat) > 0:
					#print "q", q, len(q), q.shape
					self.plot_component(s, component, 0, color=color)
		if self.orbit:
			self.orbit.Tintegrate = 10.
			p, q = self.orbit.integrate()
			x_orbit, y_orbit = p
			select(2,0)
			graph(x_orbit, y_orbit, color="lightgreen", linewidth="2px")
			
		s = 1.5
		self.snapshot.q_center, self.snapshot.p_center = self.snapshot.center()
		x, y, z = self.snapshot.q_center
		rc = sqrt(x**2+y**2+z**2)
		if rc == 0:
			self.snapshot.q_center_normal = array([1., 0., 0.])
		else:
			self.snapshot.q_center_normal = self.snapshot.q_center/rc
		if self.orbit:
			fraction = (rc-self.orbit.rp)/(self.orbit.ra-self.orbit.rp)#*0.5 + 0.5
		else:
			fraction = 1.
		fraction = self.snapshot.nr / float(self.snapshot.total)
		#print fraction
		rectangle(0, 0, fraction, 1, solid=True, color="red", container=info_container)
		if hasattr(self.snapshot, "nr"):
			if self.dt:
				text("snapshot T=%05.2f Gyr" % (self.snapshot.nr*self.dt), 0.5, 0.5, container=info_container)
			else:
				text("snapshot %03d" % self.snapshot.nr, 0.5, 0.5, container=info_container)
		xlim(0, 1, container=info_container)
		ylim(0, 1, container=info_container)
		if self.snapshot0:
			self.snapshot0.center()
		#for name in "halo bulge disk".split():
		for name in "bulge disk".split():
			if name in self.snapshot.componentmap:
				component = self.snapshot.componentmap[name]
				color = colormap[name]
				if len(component.q.flat) > 0:
					#print "q", q, len(q), q.shape
					self.plot_component(s, component, 2, color=color)
					
					
		def do(f, names="halo bulge disk".split()):
			for name in names:
				if name in self.snapshot.componentmap:
					component = self.snapshot.componentmap[name]
					f(component, color=colormap[name])

		select(3, 1)
		do(self.plotv)
		xlim(-2, 2)
		ylim(0, 30)
		labels("logr", "v")
		
		select(3, 2)
		do(self.plotsigmas)
		xlim(-2, 2)
		ylim(0, 30)
		labels("logr", "sigma v")
		
		select(3, 3)
		do(self.plotgamma2r)
		xlim(-2, 2)
		ylim(1, 6)
		hline(3)
		labels("logr", "gamma2")
		
		select(3, 0)
		do(self.plotani)
		xlim(-2, 2)
		ylim(-3, 1)
		hline(0)
		labels("logr", "anisotropy")
		
		select(4,0)
		do(self.plotencmass)
		xlim(-2, 2)
		ylim(-2, 0.1)
		vline(0)
		
		#ylim(-3, 1)
		labels("logr", "Menc fraction")
		
		select(4,1)
		do(self.plotencmasslin)
		xlim(-2, 2)
		ylim(0., 1.)
		vline(0)
		
		#ylim(-3, 1)
		labels("logr", "Menc fraction")
		
		
		select(5,0)
		do(self.plot_Rho, names=["disk"])
		select(5,1)
		do(self.plot_Rho, names=["bulge"])
		select(5,2)
		do(self.plot_vlos, names=["disk", "bulge"])
		xlim(0, 2.5)
		ylim(0, 30)
		
		select(5,3)
		do(self.plot_gamma2, names=["disk", "bulge"])
		xlim(0, 2.5)
		ylim(1, 6)
		hline(3)
		#xlim(-2, 2)
		#ylim(0., 1.)
		#vline(0)
		
		#ylim(-3, 1)
		#labels("logr", "Menc fraction")
		
	def plot_vlos(self, component, **kwargs):
		x, y, z = q = component.q
		r = sqrt(x**2+y**2+z**2)
		mask = r < 15
		x, y, z = q = component.q[:,mask] + self.snapshot.q_center.reshape(3,1)
		r = sqrt(x**2+y**2+z**2)
		
		up = array([0, 0, 1])
		right = numpy.cross(self.snapshot.q_center_normal, up)
		up = numpy.cross(self.snapshot.q_center_normal, right)
		X = dot(right, q)
		Y = dot(up, q)
		R = sqrt(X**2+Y**2)
		p = component.p[:,mask] + self.snapshot.p_center.reshape(3,1)
		vlos = dot(self.snapshot.q_center_normal, p)
		Rs = []
		vvarlos = []
		for n, Rbin, vlosbin in binrange(2500, R, vlos):
			Rs.append(mean(Rbin))
			vvarlos.append(var(vlosbin))
		vvarlos = array(vvarlos)
		Rs = array(Rs)
		graph(Rs, vvarlos**0.5, **kwargs)
	
	
	def plot_gamma2(self, component, **kwargs):
		x, y, z = q = component.q
		r = sqrt(x**2+y**2+z**2)
		mask = r < 15
		x, y, z = q = component.q[:,mask] + self.snapshot.q_center.reshape(3,1)
		r = sqrt(x**2+y**2+z**2)
		
		up = array([0, 0, 1])
		right = numpy.cross(self.snapshot.q_center_normal, up)
		up = numpy.cross(self.snapshot.q_center_normal, right)
		X = dot(right, q)
		Y = dot(up, q)
		R = sqrt(X**2+Y**2)
		p = component.p[:,mask] + self.snapshot.p_center.reshape(3,1)
		vlos = dot(self.snapshot.q_center_normal, p)
		Rs = []
		vvarlos = []
		m4los = []
		def m4(x):
			xc = x - mean(x)
			return sum(xc**4)/len(xc)
		for n, Rbin, vlosbin in binrange(2500, R, vlos):
			Rs.append(mean(Rbin))
			vvarlos.append(var(vlosbin))
			m4los.append(m4(vlosbin))
		vvarlos = array(vvarlos)
		m4los = array(m4los)
		Rs = array(Rs)
		gamma2 = m4los/vvarlos**2
		graph(Rs, gamma2, **kwargs)
	
	def plot_Rho(self, component, **kwargs):
		x, y, z = q = component.q
		r = sqrt(x**2+y**2+z**2)
		mask = r < 15
		x, y, z = q = component.q[:,mask] + self.snapshot.q_center.reshape(3,1)
		r = sqrt(x**2+y**2+z**2)
		
		e_los = q/r
		#v = array([vx[i], vy[i], vz[i]])
		#v_los = dot(e_los, v)
		# 'distances'
		#d = dot(q, self.center_normal)
		print "center", self.snapshot.q_center_normal
		up = array([0, 0, 1])
		print "up", up
		right = numpy.cross(self.snapshot.q_center_normal, up)
		up = numpy.cross(self.snapshot.q_center_normal, right)
		print "new up", up
		print "right", right
		X = dot(right, q)
		Y = dot(up, q)
		R = sqrt(X**2+Y**2)
		#import pdb
		#pdb.set_trace()
		#scatter(X, Y, **kwargs)
		contourlevels=[0.95, 0.5, 0.1, 0.05, 0.01, 0.001, 0.0001]#, 0.001]#, 0.0001]#, 0.00001]
		bincount=150
		s = 2.5
		density2d(X, Y, bincount=bincount, xmin=-s, xmax=s, ymin=-s, ymax=s, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, **kwargs)
		
		#select(4,3)
		
		
		#p_projected = q - d * self.center_normal
		#r = sqrt(p_projected[0]**2 + p_projected[1]**2 + p_projected[2]**2)
		
	def plotencmass(self, component, **kwargs):
		def do(component, **kwargs):
			x, y, z = component.q
			r = sqrt(x**2 + y**2 + z**2)
			mask = r > 0
			r = r[mask]
			logr = log10(r)
			def f(n):
				return log10(maximum(1e-6, n))
			cumhistogram(logr, datamin=-3., datamax=3., bincount=6*50, normalize=True, function=f, **kwargs)
		if self.snapshot0:
			do(self.snapshot0.componentmap[component.name], linestyle="dot", **kwargs)
		do(component, **kwargs)
	
	def plotencmasslin(self, component, **kwargs):
		def do(component, **kwargs):
			x, y, z = component.q
			r = sqrt(x**2 + y**2 + z**2)
			mask = r > 0
			r = r[mask]
			logr = log10(r)
			def f(n):
				return log10(maximum(1e-6, n))
			cumhistogram(logr, datamin=-3., datamax=3., bincount=6*50, normalize=True, **kwargs)
		if self.snapshot0:
			do(self.snapshot0.componentmap[component.name], linestyle="dot", **kwargs)
		do(component, **kwargs)
	
	
	def plotani(self, component, **kwargs):
		logrs, mr, mphi, mtheta = component.moments_spherical(5000, 2)
		graph(logrs, 1-(mphi[2]+mtheta[2])/(2*mr[2]), **kwargs)
		graph(logrs, 1-(mphi[2])/mr[2], linestyle="dash", **kwargs)
		graph(logrs, 1-(mtheta[2])/mr[2], linestyle="dot", **kwargs)
	
	def plotsigmas(self, component, **kwargs):
		logrs, mr, mphi, mtheta = component.moments_spherical(5000, 2)
		graph(logrs, mr[2]**0.5, **kwargs)
		graph(logrs, mphi[2]**0.5, linestyle="dash", **kwargs)
		graph(logrs, mtheta[2]**0.5, linestyle="dot", **kwargs)
		
	def plotgamma2r(self, component, **kwargs):
		logrs, mr, mphi, mtheta = component.moments_spherical(5000, 4)
		graph(logrs, mr[4]/mr[2]**2, **kwargs)
		graph(logrs, mphi[4]/mphi[2]**2, linestyle="dash", **kwargs)
		graph(logrs, mtheta[4]/mtheta[2]**2, linestyle="dot", **kwargs)
		
	def plotv(self, component, **kwargs):
		logrs, mr, mphi, mtheta = component.moments_spherical(5000, 1)
		graph(logrs, mr[1], **kwargs)
		graph(logrs, mphi[1], linestyle="dash", **kwargs)
		graph(logrs, mtheta[1], linestyle="dot", **kwargs)
		
	
	def plot_component(self, s, component, offset, **kwargs):
		q, p = component.q, component.p
		x, y, z = q
		
		#box()
		def do(i, x1, x2):
			select(i, 0+offset)
			scatter(x1, x2, **kwargs)
			squarelim(-s, s)
		
			select(i, 1+offset)
			contourlevels=[0.95, 0.5, 0.1, 0.05, 0.01, 0.001, 0.0001]#, 0.001]#, 0.0001]#, 0.00001]
			bincount=150
			density2d(x1, x2, bincount=bincount, xmin=-s, xmax=s, ymin=-s, ymax=s, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, **kwargs)
		
		do(2, x, y)
		do(1, x, z)
		do(0, y, z)
		
class PlotConfiguration(object):
	def __init__(self, snapshot, galaxy=None):
		self.snapshot = snapshot
		self.galaxy = galaxy
		
	def run(self, args, opts, scope):
		self.snapshot.load()
		#print self.snapshot.q.shape
		#x, y, z = self.snapshot.q #[1:4,:]
		#s = max(abs(x).max(), abs(y).max(), abs(z).max())
		s = 100
		mozaic(4,3,box)
		i = 0
		colormap = {"disk":"red", "bulge":"blue", "halo":"black"}
		#for component, color in zip(self.snapshot.components[::-1], "black blue red".split()):
		for name in "halo bulge disk".split():
			if name in self.snapshot.componentmap:
				component = self.snapshot.componentmap[name]
				color = colormap[name]
				if len(component.q.flat) > 0:
					#print "q", q, len(q), q.shape
					self.plot(s, component, color=color)
					i += 0
				#self.plot(s, self.snapshot.q[:,:10000])
		if self.galaxy:
			logr = arange(-1, 4, 0.2)
			r = 10**logr
			jeans = self.galaxy.jeans()
			sigmar = jeans.sigmar(r)
			select(2,1)
			graph(logr, sigmar, color="orange")
			select(3,2)
			#r = 
			import mab
			vcirc = (self.galaxy.dphidr(r)*r/mab.gd.potential.G)**0.5
			import pdb;
			#pdb.set_trace()
			#print vcirc
			graph(logr, vcirc, color="orange")
		draw()
	
	def plot(self, s, component, **kwargs):
		q, p = component.q, component.p
		x, y, z = q
		
		#box()
		select(0, 0)
		scatter(x, y, **kwargs)
		squarelim(-s, s)
		
		select(1,0)
		scatter(z, y, **kwargs)
		squarelim(-s, s)
		
		select(0, 1)
		scatter(x, z, **kwargs)
		squarelim(-s, s)
		
		r = sqrt(x**2+y**2+z**2)
		logr = log10(r)
		R = sqrt(x**2+y**2)
		logR = log10(R)
		select(1,1)
		obj = histogram(logr, bincount=50., datamin=-3., datamax=3., **kwargs)
		data = obj.data
		bins = obj.bins
		print "bins", bins
		centers = (bins[1:-1] + bins[0:-2])/2
		print "centers", centers
		print data, data.shape
		print bins, bins.shape
		r = 10**centers
		print r.shape
		dens = data/r**3
		select(2,0)
		print dens
		mask = dens > 0
		y = log10(dens[mask])
		print y, centers[mask]
		graph(centers[mask], y, **kwargs)
		
		select(2,1)
		if 0:
			logr, vr, vrvar = component.moments_spherical(N=500)
			print "logr", logr
			print "vdisp", vrvar**0.5
			graph(logr, vrvar**0.5, **kwargs)
	
		if 0:
			logR, vR, vRvar, vphi, vphivar, vz, vzvar = component.moments_polar(N=500)
			print "logr", logr
			print "vdisp", vrvar**0.5
			select(3,0)
			graph(10**logR, log10(vRvar**0.5), **kwargs)
			graph(10**logR, log10(vphivar**0.5), linestyle="dash", **kwargs)
			graph(10**logR, log10(vzvar**0.5), linestyle="dot", **kwargs)
			xlim(0, 6)
			#ylim(-2, -0.5)
			select(3,1)
			graph(10**logR, log10(vphivar**0.5), **kwargs)
			xlim(0, 6)
			select(3,2)
			graph(logR, vphi, **kwargs)
			vline(log10(8.0))
		#hline(220)		
		
		select(0,2)
		#draw()
		obj = cumhistogram(logR, bincount=50., datamin=-2., datamax=2., **kwargs)
		data = obj.data
		bins = obj.bins
		print "bins", bins
		centers = (bins[1:-1] + bins[0:-2])/2
		print "centers", centers
		print data, data.shape
		print bins, bins.shape
		R = 10.**centers
		print r.shape
		select(1,2)
		graph(R, sqrt(data/R), **kwargs)
		xlim(0, 10*3.5)#)
		
		
		select(2, 2)
		x, y, z = q
		R = sqrt(x**2+y**2)
		#logR = log10(R)
		#select(1,1)
		data, bins = numpy.histogram(R, bins=50, range=(0., 7.), normed=True) #bincount=50., datamin=0., datamax=7., **kwargs)
		#data = obj.data
		#bins = obj.bins
		print "bins", bins
		centers = (bins[1:] + bins[0:-1])/2
		print "centers", centers
		print data, data.shape
		print bins, bins.shape
		R = centers #10**centers
		print r.shape
		dens = data/R**1
		#select(2,0)
		print dens
		mask = dens > 0
		y = log10(dens[mask])
		print y, centers[mask]
		if len(y)>0:
			graph(R[mask], y, **kwargs)
		
		
class PlotConfigurationFromFile(object):
	def __init__(self, snapshot_class, plot_snapshot, scale=None):
		self.snapshot_class = snapshot_class
		self.plot_snapshot = plot_snapshot
		self.scale = scale
		
	def run(self, args, opts, scope):
		self.snapshot = self.snapshot_class(filename=args[1])
	
		self.plot_snapshot.snapshot = self.snapshot
		self.plot_snapshot.run(args, opts, scope)		
