from kaplot import *
try:
	import h5py
except:
	pass
import mab.gd.logging as logging
import mab.parallelize
logger = logging.getLogger("gd.gadget.plot")
from mab.binningtools import bingrid, binrange

def find_center(x, y, vx, vy):
	xc, yc = mean(x), mean(y)
	distances = sqrt((x-xc)**2+(y-yc)**2)
	N = len(x)
	indices = argsort(distances)[:N/4]
	xsub = x[indices]
	ysub = y[indices]
	vxsub = vx[indices]
	vysub = vy[indices]
	xc, yc = mean(xsub), mean(ysub)
	vxc, vyc = mean(vxsub), mean(vysub)
	if 0:
		for i in range(4):
			distances = sqrt((xsub-xc)**2+(ysub-yc)**2)
			
			N = len(xsub)
			indices = argsort(distances)[:N/2]
			xsub = xsub[indices]
			ysub = ysub[indices]
			vxsub = vxsub[indices]
			vysub = vysub[indices]
			xc, yc = mean(xsub), mean(ysub)
			vxc, vyc = mean(vxsub), mean(vysub)
		
	return xc, yc, vxc, vyc


class Gadget1(object):
	def __init__(self, galaxy, galaxy_host, orbit, working_directory, cores):
		self.galaxy = galaxy
		self.galaxy_host = galaxy_host
		self.orbit = orbit
		self.working_directory = working_directory
		self.cores = cores
		
	def run(self, args, kwargs, scope):
		oldcwd = os.getcwd()
		try:
			if self.orbit:
				self.orbit.find_orbit()
			#t = arange(N)/(N+1.) * T
				self.orbit.Tintegrate = 10.
				p, q = self.orbit.integrate()
				#print "p", p
				#print q
				x_orbit, y_orbit = p
			
			logger.info("changing working directory to: %s" % self.working_directory)
			os.chdir(self.working_directory)
			if 0:
				filename = os.path.join("output", "snapshot_%03d.hdf5" % 0)
				if not os.path.exists(filename):
					return
				
				print filename
				f = h5py.File(filename)
				if self.orbit:
					w = self.orbit.ra+25
				else:
					w = 10
				#rectangle(-w, -w, w, w, solid=True, color="lightgrey")
				pos = f["/PartType1/Coordinates"]
				vel = f["/PartType1/Velocities"]
				pot_dm = f["/PartType1/Potential"]
				x_dm, y_dm, z_dm = x, y, z = transpose(pos)
				r = (x_dm**2+y_dm**2)**0.5
				rmax = max(r)
				print "rmax", rmax
			
			#for i in range(100)[:]:
			@mab.parallelize.parallelize(cores=self.cores, info=False)
			def do(i):
				clear()
				document(size="42cm,30cm")
				mozaic(4,3,box)
				#box()
				select(0,0)
				filename = os.path.join("output", "snapshot_%03d.hdf5" % i)
				if not os.path.exists(filename):
					return
				
				print filename
				try:
					f = h5py.File(filename)
					if self.orbit:
						w = self.orbit.ra+25
					else:
						w = 10
						#w = 3.
					#rectangle(-w, -w, w, w, solid=True, color="lightgrey")
					if 1:
						ids = array(f["/PartType3/ParticleIDs"])
						mask = (ids > 15000) & (ids < 20000)
						mask = ids > 0
						pos = array(f["/PartType3/Coordinates"])[mask,:]
						vel = f["/PartType3/Velocities"][mask,:]
						pot_dm = f["/PartType3/Potential"][mask]
						x_dm, y_dm, z_dm = x, y, z = transpose(pos)
						vx_dm, vy_dm, vz_dm = vx, vy, vz = transpose(vel)
						
					ids = array(f["/PartType2/ParticleIDs"])
					print ids.min(), ids.max()
					mask = ids > 5000
					mask = ids > 0
					print ids
					print mask
					pos = f["/PartType2/Coordinates"][mask,:]
					vel = f["/PartType2/Velocities"][mask,:]
					pot = f["/PartType2/Potential"][mask]
					x, y, z = transpose(pos)
					vx, vy, vz = transpose(vel)
					if 0:
						x_dm = x
						y_dm = y
						z_dm = z
						vx_dm = vx
						vy_dm = vy
						vz_dm = vz
						pot_dm = pot
				except:
					print "error for", filename
					raise
					return
				
				if i == 0:
					logger.info("# particles = %d" % len(x))
				
				#indices = arg
				r = sqrt(x**2+y**2+z**2)
				center_index = argmin(pot+self.galaxy_host.profile_model.potentialr(r))
				center_index = argmin(pot_dm)
				x0 = x_dm[center_index]
				y0 = y_dm[center_index]
				z0 = z_dm[center_index]
				#x0, y0, z0 = 0.1, 0.1, 0
				print "zero", x0, y0, z0
				#vx0, vy0 = 0, 0
				#x0, y0, vx0, vy0 = find_center(x, y, vx, vy)
				#x0, y0, vx0, vy0 = 0, 0, 0, 0
				xc = x-x0
				yc = y-y0
				zc = z-z0
				xc_dm = x_dm - x0
				yc_dm = y_dm - y0
				zc_dm = z_dm - z0
				
				r = (xc**2+yc**2+zc**2)**0.5
				mask = r < 0.3
				vx0 = mean(vx[mask])
				vy0 = mean(vy[mask])
				vz0 = mean(vz[mask])
					
				STEP = 1
				w = 1.70
				scatter(xc_dm[::STEP], zc_dm[::STEP], color="black")
				scatter(xc[::STEP], zc[::STEP], color="red")
				xlim(-w, w)
				ylim(-w, w)
				labels("x (kpc)", "z (kpc)")
				select(0,1)
				scatter(yc_dm[::STEP], zc_dm[::STEP], color="black")
				scatter(yc[::STEP], zc[::STEP], color="red")
				xlim(-w, w)
				ylim(-w, w)
				labels("y (kpc)", "z (kpc)")
				select(0,2)
				scatter(xc_dm[::STEP], yc_dm[::STEP], color="black")
				scatter(xc[::STEP], yc[::STEP], color="red")
				xlim(-w, w)
				ylim(-w, w)
				labels("x (kpc)", "y (kpc)")
				select(3,0)
				labels("x (kpc)", "y (kpc)")
				scatter(x_dm[::STEP], y_dm[::STEP], color="black")
				scatter(x[::STEP], y[::STEP], color="red")
				title("snapshot %03d" % i)
				if self.orbit:
					graph(x_orbit, y_orbit, color="green")
				xlim(-200, 200)
				ylim(-200, 200)
				select(1,0)
				labels("x (kpc)", "z (kpc)")
				bincount = 650
				#density2d(x, y, xmin=-w, xmax=w, ymin=-w, ymax=w, normalize=True, colormap="whiterainbow", f=lambda x: maximum(log10(x), -6))
				w = 200
				bincount = 200
				w = w2 = 1.70
				contourlevels=[0.95, 0.5, 0.1, 0.05, 0.01, 0.001]#, 0.0001]#, 0.001]#, 0.0001]#, 0.00001]
				density2d(xc_dm, zc_dm, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="black")
				density2d(xc, zc, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="red")
				
				#density2d(x_dm, y_dm, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=[0.1, 0.01, 0.001], smooth=True)
				
				
				#x0, y0, vx0, vy0 = find_center(x, y, vx, vy)
				#symbol(x0, y0, color="orange")
				#x0_dm, y0_dm, vx0_dm, vy0_dm = find_center(x_dm, y_dm, vx_dm, vy_dm)
				#symbol(x0_dm, y0_dm, color="blue")
				#x0 = 0
				#y0 = 0
				#x0, y0 = 0., 0.
				#vx0, vy0 = 0., 0.
				symbol(x0, y0, color="purple", symbolName="squaresolid", symbolsize="10pt")
				symbol(0, 0, color="green", symbolName="triangle", symbolsize="10pt")
				xlim(x0-w2, x0+w2)
				ylim(y0-w2, y0+w2)
				xlim(-w2, +w2)
				ylim(-w2, +w2)
				vline(0)
				hline(0)
				
				if 1:
					select(1,1)
					labels("y (kpc)", "z (kpc)")
					#print z.min(), z.max()
					density2d(yc_dm, zc_dm, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="black")
					density2d(yc, zc, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="red")

					xlim(-w2, +w2)
					ylim(-w2, +w2)
					vline(0)
					hline(0)				
					select(1,2)
					labels("x (kpc)", "y (kpc)")
					#print z.min(), z.max()
					density2d(xc_dm, yc_dm, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="black")
					density2d(xc, yc, bincount=bincount, xmin=-w, xmax=w, ymin=-w, ymax=w, drawimage=False, drawcontour=True, contourlevels=contourlevels, smooth=True, color="red")

					xlim(-w2, +w2)
					ylim(-w2, +w2)
					vline(0)
					hline(0)				
				N = len(x)
				print x0, y0, vx0, vy0
				#x0, y0 = 0., 0.
				#x0, y0, vx0, vy0 = find_center(x, y, vx, vy)
				def do(xp, yp, zp, vx, vy, **kwargs):
					#x0, y0, vx0, vy0 = find_center(x, y, vx, vy)
					#r = ((xp-x0)**2 + (yp-y0)**2)**0.5
					r = ((xp-x0)**2 + (yp-y0)**2+(zp-z0)**2)**0.5
					print "r min", r.min()
					cumhistogram(log10(r), function=lambda x: log10(x/N), **kwargs)
					xlim(-3, 3)
				if 0:
					select(2,0)
					do(x, y, z, vx, vy, color="red")
					do(x_dm, y_dm, z_dm, vx_dm, vy_dm, color="black")
				
					rmax = 100*0.5
					logrs = arange(-3, log10(rmax), 0.1)
					if 0:
						m = self.galaxy.profile_model.dm_profile.enclosed_mass(10**logrs)
						m /= max(m)
						#y *= len(x)
						graph(logrs, log10(m), color="blue")
					
					logrs = arange(-3, log10(rmax), 0.1)
					p = mab.gd.potential.NFW(1., 10.)
					m = p.enclosed_mass(10**logrs)
					m /= max(m)
					#y *= len(x)
					graph(logrs, log10(m), color="blue", linestyle="dot")
				
				if 0:
					logrs = arange(-3, 3, 0.1)
					m = self.galaxy.light_model.light_profile.enclosedmass(10**logrs)
					m /= max(m)
					#y *= len(x)
					graph(logrs, log10(m), color="red")
				
					ylim(-4, 0) #log10(len(x)))
				if 0:
					select(3,0)
					histogram(zc, datamin=-4, datamax=4., bincount=50)
				
				def kin(xp, yp, zp, vx, vy, vz, first, **kwargs):
					Nperbin = 3000
					r3d = sqrt(xp**2+yp**2+zp**2)
					print "r3d min", r3d.min()
					#scale = float(i%5)/4.
					#color = Color(scale, 0, 0)
					#color = 
					r = sqrt(xp*xp+yp*yp+zp*zp).astype(float64)
					rhosq = (xp*xp+yp*yp);
					rho = sqrt(rhosq);
					vr = (vx*xp + vy*yp + vz*zp)/r;
					vphi = (vy*xp - vx*yp)/rho;
					vtheta = (vx*zp*xp/rho + vy*zp*yp/rho - vz*rho)/r;
					vr[isnan(vr)] = 0
					vphi[isnan(vphi)] = 0
					vtheta[isnan(vtheta)] = 0
					vr[isinf(vr)] = 0
					vphi[isinf(vphi)] = 0
					vtheta[isinf(vtheta)] = 0
					xs = []
					betas = []
					betasphi = []
					betastheta = []
					varvrs = []
					varvphis = []
					varvthetas = []
					vphis = []
					vthetas = []
					vrs = []
					rhos = []
					varvlos = []
					Rs = []
					pseudo_f_densities = []
					for n, r, vr, vphi, vtheta, x, y, vz in binrange(Nperbin, r3d, vr, vphi, vtheta, xp, yp, vz):
						#print n, mean(vr), std(vr), std(vphi)
						vphis.append(mean(vphi))
						vthetas.append(mean(vtheta))
						vrs.append(mean(vr))
						varvphi = var(vphi)# - mean(vphi)**2
						varvtheta = var(vtheta)# - mean(vtheta)**2
						varvr = var(vr)# - mean(vr)**2
						beta = 1 - (varvphi+varvtheta)/(2*varvr)
						betatheta = 1 - (varvtheta)/(varvr)
						betaphi = 1 - (varvphi)/(varvr)
						#print mean(r), beta
						dr = max(r) - min(r)
						volume = max(r) ** 3 - min(r)**3
						r = mean(r)
						logr = log10(r)
						xs.append(logr)
						rho = n/volume#(r**2*dr)
						assert not isnan(varvr), "error: %d %f %f %r" % (n, mean(vr), var(vr), vr) 
						assert n > 0
						assert rho > 0
						assert all(varvr > 0)
						assert all(varvtheta > 0)
						assert all(varvphi > 0)
						rhos.append(rho)
						pseudo_f_densities.append(rho/varvr**(3*0.5))
						betas.append(beta)
						betasphi.append(betaphi)
						betastheta.append(betatheta)
						varvrs.append(varvr)
						varvphis.append(varvphi)
						varvthetas.append(varvtheta)
						varvlos.append(var(vz))
						Rs.append(mean(sqrt(x**2+y**2)))
						
					rhos = array(rhos)
					pseudo_f_densities = array(pseudo_f_densities)
					varvlos = array(varvlos)
					
					if 0:
						select(0,2)
						graph(xs, sqrt(varvrs), color=color, alpha=alpha)
						graph(xs, sqrt(varvphis), color=color, linestyle="dash", alpha=alpha)
						graph(xs, sqrt(varvthetas), color=color, linestyle="dot", alpha=alpha)
					logrs = arange(min(xs), max(xs), 0.01)
					#if isinstance(galaxy, Galaxy_constant_anisotropy)
						
					
					#ylim(0, 20)
					xs = array(xs)
					rs = 10**xs
					drhodr = (rhos[1:] - rhos[:-1])/(rs[1:]-rs[:-1])
					logslope = drhodr * (rs[1:]+rs[:-1])/(rhos[1:]+rhos[:-1])
					print "slopes", logslope
					if 1:
						select(3,2)
						if first:
							labels("log r (kpc)", "logslope")
						graph((xs[1:]+xs[:-1])/2, logslope, **kwargs)
						xlim(-2, 2)
						ylim(-5, 1)
					#select(1,0)
					#print len(xs), len(betas)
					if 1:
						select(2, 2)
						graph(xs, betas, **kwargs)
						graph(xs, betasphi, linestyle="dash", **kwargs)
						graph(xs, betastheta, linestyle="dot", **kwargs)
						if first:
							labels("log r (kpc)", "anisotropy")
						xlim(-2, 2)
						ylim(-3, 1)
						hline(0)
					
					if 0:
						select(1, 2)
						graph(10**xs, betas, **kwargs)
						ylim(-3, 1)
						xlim(0, 3.0)
						hline(0)
				
					select(2, 0)
					#print ">>>", xs, varvrs
					graph(xs, vrs, linestyle="normal", **kwargs)
					graph(xs, vthetas, linestyle="dot", **kwargs)
					graph(xs, vphis, linestyle="dash", **kwargs)
					if first:
						labels("log r (kpc)", "vr,vtheta,vphi (km/s)")
					vmax = 15
					ylim(-vmax, vmax)
					xlim(-2, 2.0)
					hline(0)
					if 0:
						select(1, 2)
						graph(10**xs, vrs, linestyle="normal", **kwargs)
						graph(10**xs, vthetas, linestyle="dot", **kwargs)
						graph(10**xs, vphis, linestyle="dash", **kwargs)
						ylim(-vmax, vmax)
						xlim(0, 5.0)
						hline(0)
					
					varvrs = array(varvrs)
					varvthetas = array(varvthetas)
					varvphis = array(varvphis)
					select(2, 1)
					graph(xs, varvrs**0.5, linestyle="normal", **kwargs)
					graph(xs, varvthetas**0.5, linestyle="dot", **kwargs)
					graph(xs, varvphis**0.5, linestyle="dash", **kwargs)
					if first:
						labels("log r (kpc)", "vel disp (km/s)")
					vsigmamax = 15
					ylim(0, vsigmamax)
					xlim(-2, 2.0)
					if 0:
						select(2, 2)
						graph(10**xs, varvrs**0.5, linestyle="normal", **kwargs)
						graph(10**xs, varvthetas**0.5, linestyle="dot", **kwargs)
						graph(10**xs, varvphis**0.5, linestyle="dash", **kwargs)
						ylim(0, vsigmamax)
						xlim(0, 3.0)
				
					if 1:
						select(3, 1)
						if first:
							labels("log r (kpc)", "log rho")
						graph(xs, log10(rhos), **kwargs)
						xlim(-2, 2.0)
						ylim(-5, 10)
					
					
					select(3, 1)
					if 0:
						mask = pseudo_f_densities > 0
						graph(xs[mask], log10(pseudo_f_densities[mask]), **kwargs)
						xlim(-2, 2.0)
						ylim(-5, 5)
					else:
						if 0:
							graph(log10(Rs), varvlos**0.5, **kwargs)
							ylim(0, vsigmamax)
							xlim(-3, 3)
						
				#x0, y0 = 0, 0
				#vx0, vy0 = 0, 0
				kin(x-x0, y-y0, z-z0, vx-vx0, vy-vy0, vz-vz0, True, color="red")
				##kin(x-x0, y-y0, z, vx-vx0, vy-vy0, vz, color="red")
				##kin(x_dm-x0, y_dm-y0, z_dm, vx_dm-vx0, vy_dm-vy0, vz_dm, color="black")
				kin(x_dm-x0, y_dm-y0, z_dm-z0, vx_dm-vx0, vy_dm-vy0, vz_dm-vz0, False, color="black")
				
				#draw(
				hardcopy("frame_%03d.png" % i)
				f.close()
			#for i in range(100):
			#	do(i)
			do(range(100))
			#do([99])
		finally:
			os.chdir(oldcwd)		
		