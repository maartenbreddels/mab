# -*- coding: utf-8 -*-
from kaplot import *

class PlotEffectivePotential(object):
	def __init__(self, profile):
		self.profile = profile
		
	def run(self, argv, opts, scope):
		mozaic(2,2,box)
		re = float(argv[1])
		E = self.profile.potential_xy(0.01, re)
		print "E=", E
		#i2 = float(argv[2]) * abs(E)
		#i3 = float(argv[3]) * abs(E)
		i3_max, t = self.profile.find_i3_max(E)
		i3 = i3_max * float(argv[2])
		i2_max, t = self.profile.find_i2_max(E, i3)
		i2 = i2_max * float(argv[3])
		
		print i3_max, t, t-1
		print "I2,I3", i2, i3
		#dsa
		
		#re = float(argv[1])
		#E = self.profile.potential_xy(0.2, re)
		#print "E=", E
		#i2 = float(argv[2])
		#i3 = float(argv[3])
		tau = arange(1e-4, 1., 0.01) - self.profile.gamma
		pot = self.profile.V_eff_tau(tau, i2, i3)
		U = -self.profile.G_tau(tau)
		#potinf = self.profile.G_tau(1e8)
		#print potinf
		#dsa
		#print tau
		#print pot
		pot[isnan(pot)] = 0
		pot[isinf(pot)] = 0
		graph(tau, pot)
		graph(tau, U, color="blue")
		graph(tau, i2/(tau+self.profile.alpha), color="orange")
		graph(tau, i3 / (tau + self.profile.gamma), color="purple")
		#print "U", U
		#dsa
		hline(E, color="green")
		hline(0)
		#hline(-E, color="green")
		
		vline(-self.profile.alpha, color="red")
		vline(-self.profile.gamma, color="red")
		xlim(-0.1-self.profile.gamma, tau.max())
		s = 0.0001/3
		ylim(-s, s)
		select(1, 0)
		
		n_mu = 150
		n_lambda = 250
		lambda_max = 3.
		u = (arange(n_mu)+0.5) / (n_mu)
		#u = (arange(n_mu)) / (n_mu-1.)
		nus = self.profile.uniform_to_nu(u)
		u = (arange(n_lambda)+0.5) / (n_lambda)
		#u = (arange(n_lambda)) / (n_lambda-1.)
		lambdas = self.profile.uniform_to_lambda(u, lambda_max) + 0.001
		pot = self.profile.V_eff_tau(nus, i2, i3)
		#print "nus", nus
		#print "pot", pot, E 
		#print pot>E
		
		nus = nus[pot>E] 
		
		pot = self.profile.V_eff_tau(lambdas, i2, i3)
		lambdas = lambdas[pot<E]
		#print pot, E 
		#print pot<E
		#print lambdas
		
		
		for nu in nus:
			x, y = self.profile.conf_to_xy(nu, lambdas)
			if len(x) > 0:
				pass
				graph(x, y)
				#scatter(x, y, color="red")
				#print x, y
		for l in lambdas:
			x, y = self.profile.conf_to_xy(nus, l)
			if len(x) > 0:
				pass
				graph(x, y)
				#scatter(x, y, color="red")
		if 1:
			draw()
			sys.exit(0)
		elif 1:
			nu, la = nus.max(), lambdas.max()
			nu, la = nus.mean(), lambdas.mean()
			x, y = self.profile.conf_to_xy(nu, la)
			symbol(x, y, color="red")
			print "nu, lambda", nu, la
			dla = dnu = 1e-10
			#dla = 1e-8
			
			x0, y0 = x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu, la+dla)
			print "dx/dl", dla/(x2-x1)
			
			dx = 1e-8
			#nu1, la1 = self.profile.xy_to_conf(x1, y1)
			#nu2, la2 = self.profile.xy_to_conf(x1+dx, y1)
			
			#print dxdl - dxdn * dnu/   
			
	
			#draw()
			x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu, la+dla)
			a = (x2-x1)/dla
			x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu, la+dla)
			b = (y2-y1)/dla
			x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu+dnu, la)
			c = (x2-x1)/dnu
			x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu+dnu, la)
			d = (y2-y1)/dnu
			#m = numpy.matrix([[a,b],[c,d]])
			m = numpy.matrix([[a,c],[b,d]])
			print "m", "d(x,y)/d(la,nu)\n", m
			print "i", m.I
			
			#print [a,b],[c,d]
			print "inv", 1/(d/(a*d-b*c))
			
			#print "dl/dx", dla/(x2-x1)
			#print "dl/dy", dla/(y2-y1)
			dy = dx = 1e-8
			#print "x,y", x1, y1
			nu1, la1 = self.profile.xy_to_conf(x0, y0)
			nu2, la2 = self.profile.xy_to_conf(x0, y0+dy)
			#a = (
			
			nu1, la1 = self.profile.xy_to_conf(x0, y0)
			nu2, la2 = self.profile.xy_to_conf(x0, y0+dy)
			
			nu1, la1 = self.profile.xy_to_conf(x0, y0)
			nu2, la2 = self.profile.xy_to_conf(x0+dx, y0)
			
			nu1, la1 = self.profile.xy_to_conf(x0, y0)
			nu2, la2 = self.profile.xy_to_conf(x0+dx, y0)
			
			
			print x0, y0, x1, y1
			print nu, la, nu1, la1
			#print "nu2, la2", nu2,la2
			#print nu, la, nu1, la1
			
			print "> dl/dx", (la2-la1)/(dx)
			dy = dx
			nu1, la1 = self.profile.xy_to_conf(x1, y1)
			nu2, la2 = self.profile.xy_to_conf(x1, y1+dy)
			print "> dl/dy", (la2-la1)/(dy)
			print 2 * x * (self.profile.alpha-self.profile.gamma) / (nu+self.profile.alpha)
			
			alpha = self.profile.alpha
			gamma = self.profile.gamma
			Psq = (la - nu) / (4*(la+alpha) * (la + gamma))
			Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma))
			print "m", x / (2*(la+alpha)*Psq)
			print "m", x / (2*(nu+alpha)*Qsq)
			print "m", y / (2*(la+gamma)*Psq)
			print "m", y / (2*(nu+gamma)*Qsq)
			
			
			p1 = self.profile.G_tau(la)
			p2 = self.profile.G_tau(la+dla)
			print "Gprime ", (p2-p1)/dla
			print "Gprime ", self.profile.G_tau_prime(la)
			p1 = self.profile.G_tau(nu)
			p2 = self.profile.G_tau(nu+dnu)
			print "Gprime ", (p2-p1)/dnu
			print "Gprime ", self.profile.G_tau_prime(nu)
			
			p1 = self.profile.potential_conf(nu, la)
			p2 = self.profile.potential_conf(nu, la+dla)
			print "dphidla ", (p2-p1)/dla
			print "dphidla", self.profile.dphi_dla(nu, la)
			print "dphidla", self.profile.dphi_dla(nu, la+dla)
			print
			
			p1 = self.profile.potential_conf(nu, la)
			p2 = self.profile.potential_conf(nu+dnu, la)
			print "dphidnu ", (p2-p1)/dnu
			print "dphidnu", self.profile.dphi_dnu(nu, la)
			
			p1 = self.profile.potential_xy(x, y)
			p2 = self.profile.potential_xy(x+dx, y)
			print "dphidx ", (p2-p1)/dx
			print "dphidx", self.profile.dphi_dx(x, y)
			
			p1 = self.profile.potential_xy(x, y)
			p2 = self.profile.potential_xy(x, y+dy)
			print "dphidy ", (p2-p1)/dy
			print "dphidy", self.profile.dphi_dy(x, y)
			
			P = sqrt(Psq)
			Q = sqrt(Qsq)
 			
			dxdt = 1.
			dydt = 1.4
			dt = 1e-6
			x1 = x
			x2 = x1+dxdt*dt
			nu1, la1 = self.profile.xy_to_conf(x1, y)
			nu2, la2 = self.profile.xy_to_conf(x2, y+dydt*dt)
			print "dladt", (la2-la1)/dt
			print "dnudt", (nu2-nu1)/dt
			v1 = t = dxdt * a + dydt * b
			print t, t/Psq
			v2 = t = dxdt * c + dydt * d
			print t, t/Qsq
			print "s", (dxdt**2+dydt**2)**0.5
			#print ((v1)**2+(v2)**2)**0.5
			print ((v1/P)**2+(v2/Q)**2)**0.5
			#print ((v1/Psq)**2+(v2/Qsq)**2)**0.5
			
			
			#la1 = la
			dladt = 1.
			dnudt = 1.4
			#la2 = la+dladt*dt
			x1, y1 = self.profile.conf_to_xy(nu, la)
			x2, y2 = self.profile.conf_to_xy(nu+dnudt*dt, la+dladt*dt)
			print "vx", (x2-x1)/dt
			#print P*dladt
			print "vy", (y2-y1)/dt
			
			a = x / (2*(la+alpha)*P)
			b = x / (2*(nu+alpha)*Q)
			c = y / (2*(la+gamma)*P)
			d = y / (2*(nu+gamma)*Q)
			v1 = dladt * a * P + dnudt * b * Q
			#print dladt, a,b,c,d,Qsq, Psq # * c * Qsq
			v2 = dladt * P * c + dnudt * d * Q
			print a,b,c,d,Qsq, Psq, P, Q
			print "vx", v1
			print "vx", v2
			print "s", (v1**2+v2**2)**0.5
			print ((dladt*P)**2+(dnudt*Q)**2)**0.5
			#import pdb
			#pdb.set_trace()
			sys.exit(0)
			
			r = 1.5
			theta = 0.1
			x = r * cos(theta)
			y = r * sin(theta)
			a = x/r
			b = y/r
			c = -y/r**2
			d = x/r**2
			m = numpy.matrix([[a,c],[b,d]])
			print "m", m
			print "i", m.I
			
			
			a = cos(theta)
			b = -y
			c = sin(theta)
			d = x
			m = numpy.matrix([[a,c],[b,d]])
			print "m", m
			print "i", m.I
			
			
			dsa
 
		xlim(0, lambda_max**0.5)
		ylim(0, lambda_max**0.5)
		xlim(0, 2)
		ylim(0, 2)
		draw()

class PlotOrbitFromIC(object):
	def __init__(self, profile, orbit_integrator, gridder, initial_conditions, orbit_library=None):
		self.profile = profile
		self.orbit_integrator = orbit_integrator
		self.gridder = gridder
		self.initial_conditions = initial_conditions
		self.orbit_library = orbit_library
		
	def run(self, argv, opts, scope):
		self.initial_conditions.load()
		if self.orbit_library:
			self.orbit_library.load()
		document(size="35cm,25cm")
		mozaic(3,2,box)
		i,j,k = int(argv[1]), int(argv[2]), int(argv[3])
		print i,j,k
		ic = self.initial_conditions.initial_conditions
		x, y, vx, vy, E, i2, i3, i2_max, i3_max, Torbit , la_max = ic[i,j,k]
		print "x,y", x, y
		print "v", vx, vy
		print "integrals", E, i2, i3
		print "i2/i3 max", i2_max, i3_max
		print "Torbit", Torbit
		
		eps = 1e-9
		nu_min = -self.profile.gamma + eps
		nu_max = self.profile.find_nu_max(E, i2, i3)
		la_min, la_max = self.profile.find_la_range(E, i2, i3, la_max)
		scale = 1.00001
		la_mid = (la_max + la_min) / 2
		nu_mid = (nu_max + nu_min) / 2
		print "nu", nu_min, nu_mid, nu_max
		print "la", la_min, la_mid, la_max
		if 1:
			la_min *= scale
			la_max /= scale
			nu_min *= scale
			nu_max /= scale
		print "> mid", self.profile.density_orbit_conf(nu_mid, la_mid, E, i2, i3)
		for nu in [nu_min, nu_max]:
			for la in [la_min, la_max]:
				print nu, la, self.profile.density_orbit_conf(nu, la, E, i2, i3)
				
		import scipy.integrate
		eps = 0.001
		if 0:
			N = 100
			x = arange(N) / (N-1.)
			x = x * (la_max - la_min) + la_min
			y = self.profile.density_orbit_conf(nu_min, x, E, i2, i3)
			print y
			graph(x, y)
			#print scipy.integrate.quad(lambda la: self.profile.density_orbit_conf(nu_min, la, E, i2, i3), la_min, la_max)
			#print scipy.integrate.dblquad(lambda nu, la: self.profile.density_orbit_conf(nu, la, E, i2, i3)/1e8, nu_min, nu_max, lambda _: la_min, lambda _: la_max,epsrel=1e-3)
			draw()
			dsa
		
		
		#dsa
		if self.orbit_library:
			select(1,1)
			moments = self.orbit_library.momentgrid[i,j,k]
			mass = moments[0]
			mass[isnan(mass)] = 0
			mass /= mass.max()
			I = log10(mass)
			I[isnan(I)] = 0
			I[I<-4] = -4
			print "imax", I.max(), I.min()
			
			indexedimage(I.T, resize=self.gridder.resize)
			contour(I.T, levels=40, resize=self.gridder.resize, color="red")
			if 0:
				las = arange(-self.profile.alpha, 4, 0.01)
				all = array([self.profile.conf_to_xy(nu_max, la) for la in las])
				x, y = all.T
				graph(x, y, color="red")
				nus = arange(-self.profile.gamma, -self.profile.alpha, 0.01)
				all = array([self.profile.conf_to_xy(nu, la_max) for nu in nus])
				x, y = all.T
				graph(x, y, color="red")
				nus = arange(-self.profile.gamma, -self.profile.alpha, 0.01)
				all = array([self.profile.conf_to_xy(nu, la_min) for nu in nus])
				x, y = all.T
				graph(x, y, color="red")
				draw()
			
		 
		
		Lz = sqrt(2*i2)
		Lzs = array([Lz])
		#p0 = array([[0., 0.]])
		import pyublas
		p0 = array([[vx, vy]])
		q0 = array([[x, y]])
		Norbits = 2000*5
		self.orbit_integrator.orbital_points = 100000
		dt = array([Torbit/self.orbit_integrator.orbital_points*Norbits])
		print dt, Lzs, q0, p0
		print ">>>", dt.dtype, Lzs.dtype, q0.dtype, p0.dtype
		print ">>>", dt.shape, Lzs.shape, q0.shape, p0.shape
		q1, p1 = self.orbit_integrator.integrate(dt, Lzs, q0, p0)
		x = ravel(q1[0,:])
		y = ravel(q1[1,:])
		vx = ravel(p1[0,:])
		vy = ravel(p1[1,:])
		select(0, 0)
		graph(x, abs(y))
		xo, yo = x, y
		s = 1
		xlim(0, s)
		ylim(0, s)
		
		
		x = arange(0, 1, 1./50) + 1/50./2
		y = arange(0, 1, 1./50) + 1/50./2
		
		x, y = numpy.meshgrid(x, y)
		import time
		t1 = time.time()
		rho = self.profile.density_orbit_xy(x, y, E, i2, i3)
		print time.time() - t1
		dsa
		rho[isnan(rho)] = 0
		select(1, 0)
		mass = rho# * x
		mass /= mass.max()
		I = log10(mass)
		I[I<-4] = -4
		print I.max()
		indexedimage(I, resize=self.gridder.resize)
		#contour(I, levels=20, resize=self.gridder.resize)
		select(2,0)
		Ia = I
		indexedimage(I, resize=self.gridder.resize)
		graph(xo, yo, alpha=0.3)
		
		select(0,1)
		I = self.gridder(xo, abs(yo))
		density_orbit = I
		I /= I.max()
		I = log10(I)
		I[I<-4] = -4
		indexedimage(I.T, resize=self.gridder.resize)
		contour(I.T, levels=30, resize=self.gridder.resize)
		contour(Ia, levels=40, resize=self.gridder.resize, color="red")
		
		
		select(2,1)
		vsq = (vx**2 + vy**2)#**0.5
		densityvsq = self.gridder(xo, abs(yo), weights=vsq, normed=True)
		I = densityvsq/density_orbit
		I[isnan(I)] = 0
		I /= I.max()
		#I = log10(I)
		#I[I<-4] = -4
		indexedimage(I.T, resize=self.gridder.resize)
		contour(I.T, levels=10, resize=self.gridder.resize)
		
		vx, vy = self.profile.v_xy(x, y, E, i2, i3)
		vsq = (vx**2 + vy**2)**0.5
		I = vsq
		I[isnan(I)]= 0
		I /= I.max()
		#I = log10(I)
		#I[I<-4] = -4
		#indexedimage(I, resize=self.gridder.resize)
		contour(I, levels=10, resize=self.gridder.resize, color="red")
		#graph(xo, yo, alpha=0.3)
		
		
		
		draw()
			

class PlotOrbitLibrary(object):
	def __init__(self, orbit_library):
		self.orbit_library = orbit_library
		
	def run(self, args, opts, scope):
		self.orbit_library.load()
		grid = self.orbit_library.initial_conditions.grid
		document(size="25cm,25cm")
		if len(args) != 3:
			print "usage: <E index> <moment index>"
		i = int(args[1])
		moment = int(args[2])
		mozaic(grid.gridy.N, grid.gridz.N, container)
		for j in range(grid.gridy.N):
			for k in range(grid.gridz.N):
				select(j, k)
				border()
				spacer("1mm")
				data = self.orbit_library.momentgrid[i,j,k,moment]
				if moment > 0:
					data /= self.orbit_library.momentgrid[i,j,k,0]
				data[isnan(data)] = 0
				print data.max(), data.min()
				indexedimage(data.T, colormap="whiterainbow")
		draw()
		
		
		
class PlotOrbitAxi(object):
	def __init__(self, profile, orbit_integrator, gridder):
		self.profile = profile
		self.orbit_integrator = orbit_integrator
		self.gridder = gridder
		
	def run(self, argv, opts, scope):
		document(size="35cm,25cm")
		mozaic(3,2,box)
		if 0:
			re = float(argv[1])
			E = self.profile.potential_xy(0.2, re)
			print "E=", E
			i2 = float(argv[2]) * abs(E)
			i3 = float(argv[3]) * abs(E)
			print "I2,I3", i2, i3
		
		re = float(argv[1])
		E = self.profile.potential_xy(0.01, re)
		print "E=", E
		#i2 = float(argv[2]) * abs(E)
		#i3 = float(argv[3]) * abs(E)
		i3_max, t = self.profile.find_i3_max(E)
		i3 = i3_max * float(argv[2])
		i2_max, t = self.profile.find_i2_max(E, i3)
		i2 = i2_max * float(argv[3])
		
		
		n_y = 155
		n_x = 155
		if 1:
			y = (arange(n_y) + 0.5) / (n_y) * 2.
			x = (arange(n_x) + 0.5) / (n_x) * 2.
			x, y = numpy.meshgrid(x, y)
			xs, ys = x, y
			rho = self.profile.density_orbit_xy(x, y, E, i2, i3)
			print rho
			if 0:
				rho[isnan(rho)] = -10
				rho = log10(rho/rho.max())
				rho[rho<-3] = 0
				rho[isnan(rho)] = -3
			indexedimage(rho, colormap="whiteblackred", resize=[(0, 0), (2,2)])
			
		if 1:
			#Lz = i2*1100
			#Lz = i2#*350
			Lz = sqrt(2*i2)
			timescale = 1e20*20#/10*1*0.725
			#t1 = self.profile.dphi_dx(x, y)
			#t2 = self.profile.dphi_dx_num(x, y)
			#print "error", (abs(t1-t2)/t1).max()
			#indexedimage(t1, resize=[(0, 0), (2,2)])
			pot = self.profile.potential_xy(x, y) + Lz*Lz/(2*x*x) - E
			#import pdb
			#pdb.set_trace()
			#indexedimage(pot, resize=[(0, 0), (2,2)])
			#contour(pot, levels=10, resize=[(0, 0), (2,2)])
			contour(pot, levels=[0], resize=[(0, 0), (2,2)], color="red", linewidth="1px")
			nu, la = self.profile.xy_to_conf(x, y)
			vnu = self.profile.V_eff_tau(nu, i2, i3)
			vla = self.profile.V_eff_tau(la, i2, i3)
			contour(vnu-E, levels=[0], resize=[(0, 0), (2,2)], color="orange", linewidth="1px")
			contour(vla-E, levels=[0], resize=[(0, 0), (2,2)], color="orange", linewidth="1px")
			
			
			print "pot", pot.min(), pot.max()
			dt = array([timescale/self.orbit_integrator.orbital_points*10.1])
			#x0, z0 = 0.81, 1.64
			#x0, z0 = 0.726754, 1.12584
			x0, z0 = 0.239, 0.73
			x0, z0 = 0.885, 0.97
			x0, z0 = 1.13367, 0.379
			#x0, z0 = 0.501, 1.090
			x0, z0 = 0.501, 1.090
			x0, z0 = 1.0, 0.8
			x0, z0 = 0.4, 0.2
			x0, z0 = 0.351902,  0.0903412
			if 0:
				zs = arange(0, 2, 0.001)
				pot = self.profile.potential_xy(x0, zs) + Lz*Lz/(2*x0*x0) - E
				#mask = pot<E
				indices = argsort(pot)
				best = None
				for i in indices:
					if pot[i] < 0:
						best = i
						break
				#import pdb
				#pdb.set_trace()
				z0 = zs[best]

			#nu0, la0 = self.profile.xy_to_conf(x0, z0)
			vx0, vz0 = self.profile.v_xy(x0, z0, E, i2, i3)
			vx0, vz0 = -vx0, -vz0
			print vx0, vz0
			#dsa
			#x0, z0 = 0.63, 1.06
			q0 = array([[x0, z0]])
			#eps = 1e-16
			#v = vz
			#angle = math.radians(float(args[3]))
			#p0 = array([[v*sin(angle), v*cos(angle)]])
			Lzs = array([Lz])
			#p0 = array([[0., 0.]])
			p0 = array([[vx0, vz0]])
			print ">>>", dt, Lzs, q0, p0
			print ">>>", dt.dtype, Lzs.dtype, q0.dtype, p0.dtype
			print ">>>", dt.shape, Lzs.shape, q0.shape, p0.shape
			q1, p1 = self.orbit_integrator.integrate(dt, Lzs, q0, p0)
			x = ravel(q1[0,:])
			z = ravel(q1[1,:])
			vx = ravel(p1[0,:])
			vz = ravel(p1[1,:])
			#vx = abs(vx)
			#vz = abs(vz)
			mask = ~(isnan(x) & isnan(z))
			step = 1150
			graph(x[mask], z[mask], color="blue", alpha=0.7)
			scatter(x[::step], z[::step], color="black", symbolsize="20pt")
			
			
			#zmax, xmax
			
			select(2,0)
			v = sqrt(vx**2+vz**2)
			densityv = self.gridder(x, z, weights=abs(v), normed=True)
			densityr = self.gridder(x, z, normed=False)
			density = densityv / densityr
			density[isnan(density)] = 0
			indexedimage(density.T, colormap="whiterainbow", resize=self.gridder.resize)
			d = density.T
			#d /= d.max()
			#d = log10(d)
			#d[d<-2] = -2
			contour(d, levels=20, resize=self.gridder.resize)
			select(2,1)
			draw()
			
			#d = self.profile.density_orbit_xy(xs, ys, E, i2, i3)
			dx, dy = self.profile.v_xy(xs, ys, E, i2, i3)
			d = sqrt(dx**2+dy**2)
			print d.shape
			print d.max()
			d[isnan(d)] = 0
			#d /= d.max()
			#d = log10(d)
			#d[d<-2] = -2
			indexedimage(d, colormap="whiterainbow", resize=self.gridder.resize)
			contour(d, levels=20, resize=self.gridder.resize)
			select(2,0)
			contour(d, levels=20, resize=self.gridder.resize, linestyle="dash", color="red")
			
			if 0:
				nu, la = self.profile.xy_to_conf(abs(x), abs(z))
				select(0,1)
				graph(nu)
				graph(la)
				labels("t", "tau")
				select(1,1)
				graph(vx, color="red")
				graph(vz, color="red", linestyle="dot")
				scatter(arange(len(vx))[::step], vx[::step], symbolsize="20pt")
				for t in arange(len(vx))[::step]:
					vline(t)
				v = sqrt(vx**2+vz**2)
				print "v", v
				graph(v, color="orange")
				#graph(vz, color="red")
				
				alpha = self.profile.alpha
				gamma = self.profile.gamma
				x0, y0 = x, z
				@vectorize
				def conv(vx, vz, nu, la, x, y):
					#print "...",vx, vz, nu, la
					#dsadsa
					Psq = (la - nu) / (4*(la+alpha) * (la + gamma))
					Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma))
					#matrix = numpy.matrix([[a,b], [c,d]])
					a = x / (2*(la+alpha)*Psq)
					b = x / (2*(nu+alpha)*Qsq)
					c = y / (2*(la+gamma)*Psq)
					d = y / (2*(nu+gamma)*Qsq)
					#print "1", nu, la
					#print "2", a, b, vx, vz
					return (vx * a + vz * b) #* sqrt(Psq)
				#vla = vx * a + vz * b
				print vx.shape, vz.shape, nu.shape, la.shape, x.shape, z.shape
				beta = alpha
				mu =- alpha
				row1a = sqrt((la+beta) * (la+gamma) * (nu + alpha) / ((alpha-gamma) * (la-mu) * (la-nu)))
				row1b = -(mu+gamma) * (la+alpha) * (nu + alpha) / ((alpha-gamma) * (mu-nu) * (mu-la))
				row1c = -sqrt((nu+beta) * (nu+gamma) * (la+alpha) / ((alpha-gamma) * (nu-la) * (nu-mu)))
				
				print row1c
				#dsa
				#da
				vla = conv(vx, vz, nu, la, x, z)
				#vla = array([conv(vx_, vz_, nu_, la_) for vx_, vz_, nu_, la_ in zip(vx, vz, nu, la)])
				print vla.shape 
				select(1,0)
				labels("t", "v_tau")
				#graph(abs(vla), color="red")
				Psq = (la - nu) / (4*(la+alpha) * (la + gamma))
				Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma))
				P = sqrt(Psq)
				Q = sqrt(Qsq)
				#s2, s3 = 0.90, 0.9
				#s2, s3 = 1.10, 1.1
				s2, s3 = 1.0, 1.
				vla = self.profile.p_sq_tau(la, E, i2*s2, i3*s3)**0.5 / P
				mask = ~isnan(vla)
				graph(arange(len(vla))[mask], vla[mask], color="blue")
				vnu = self.profile.p_sq_tau(nu, E, i2*s2, i3*s3)**0.5 / Q
				#print vnu.max()
				mask = ~isnan(vnu)
				graph(arange(len(vnu))[mask], vnu[mask], color="blue", linestyle="dot")
				for t in arange(len(vx))[::step]:
					vline(t)
				
				
				select(1,1)
				labels("t", "vx, vy")
				vx1 = vla * row1a
				vx2 = vnu * row1c
				x, y = x0, y0
				a = x / (2*(la+alpha)*P)
				b = x / (2*(nu+alpha)*Q)
				c = y / (2*(la+gamma)*P)
				d = y / (2*(nu+gamma)*Q)
				#vx = vla * a / P + vnu * b / Q#* Qsq
				vx = vla * a + vnu * b
				vy = vla * c + vnu * d
				vx = abs(vx)
				vy = abs(vy)
				
				D = abs(a*d-c*b)
				print "dets", D.mean(), D.max(), D.min()
				
				
				#3print dladt, a,b,c,d,Qsq, Psq # * c * Qsq
				#vy dladt * Psq * c
				
				#graph(vx, linestyle="dot") 
				#graph(vx1, linestyle="dot") 
				#graph(vx2, linestyle="dot") 
				#vx = vx1+vx2
				#draw()
				mask = ~isnan(vx)
				graph(arange(len(vx))[mask], vx[mask])#, linestyle="dash")
				mask = ~isnan(vy)
				graph(arange(len(vy))[mask], vy[mask], linestyle="dot")
				v = sqrt(vla**2+vnu**2)
				v2 = sqrt(vx**2+vy**2)
				print (v-v2)
				#dsa
				mask = ~isnan(v)
				print "v", v
				graph(arange(len(v))[mask], v[mask], color="purple", linestyle="dash")
				mask = ~isnan(v2)
				graph(arange(len(v2))[mask], v2[mask], color="yellow", linestyle="dot")
				#vx = vla * row1a + vnu * row1c
	
				
				#vx = self.profile.vx(nu, la)
				#graph(vx, color="blue")
			
			
		draw()
			
			
		
		#for i in range(0, 10):
		#	x = i/3.
		#	print x, self.profile.density_orbit_xy(x, re, E, i2, i3)
		#print self.profile.density_orbit_xy(0.1, re, E, i2, i3)


class PlotCoordinatesAxi(object):
	def __init__(self, profile):
		self.profile = profile
		
	def run(self, argv, opts, scope):
		
		n_nu = 40
		n_lambda = 100
		lambda_max = 2.
		#u = (arange(n_nu)+0.5) / (n_nu)
		u = (arange(n_nu)) / (n_nu-1.)
		nus = self.profile.uniform_to_nu(u)
		u = (arange(n_lambda)) / (n_lambda-1.)
		lambdas = self.profile.uniform_to_lambda(u, lambda_max)
		print lambdas
		
		box()
		for nu in nus:
			x, y = self.profile.conf_to_xy(nu, lambdas)
			graph(x, y)
			scatter(x, y, color="red")
			
		n_nu = 100
		n_lambda = 40
		lambda_max = 2.
		#u = (arange(n_nu)+0.5) / (n_nu)
		u = (arange(n_nu)) / (n_nu-1.)
		nus = self.profile.uniform_to_nu(u)
		u = (arange(n_lambda)) / (n_lambda-1.)
		lambdas = self.profile.uniform_to_lambda(u, lambda_max)
			
		for la in lambdas:
			x, y = self.profile.conf_to_xy(nus, la)
			graph(x, y)
			scatter(x, y, color="blue")
		grow(1.1)
		draw()
		
		
class PlotDensity(object):
	def __init__(self, profile, fit=None):
		self.profile = profile
		self.fit = fit
		
	def run(self, argv, opts, scope):
		if self.fit:
			self.fit.load()
		document(size="30cm,25cm")
		page(fontsize="18pt")
		mozaic(2,2,box)
		n_y = 100
		n_x = 100
		#y = (arange(n_y) + 0.5) / (n_y) * 10
		
		logy = (arange(n_y)) / (n_y-1.) * 3.5 - 2.0
		y = 10**logy
		x = 0.
		nus, las = self.profile.xy_to_conf(x, y)
		rho = self.profile.rho_conf(nus, las)
		print "x", x
		print "y", y
		print "nu", nus
		print "la", las
		print "rho", rho
		#dsa
		labels("log r", "log rho") 
		
		graph(log10(y), log10(rho))
		if self.fit:
			rho_fit = self.fit.rho_xy(x, y)
			print "x,y",x, y
			print rho_fit
			#graph(log10(y), log10(rho_fit), linestyle="dot", color="green")
			graph(log10(y), log10(rho_fit), linestyle="dot")
			
			rho_fit = self.fit.density_dir(y, pi/20)
			graph(log10(y), log10(rho_fit), linestyle="dotdash", color="blue")
			rho_fit = self.fit.density_dir(y, pi/2-pi/20)
			graph(log10(y), log10(rho_fit), linestyle="dotdash", color="blue")
			#rho_fit = self.profile.density_dir(y, pi/2-pi/20)
			#graph(log10(y), log10(rho_fit), linestyle="dotdash", color="red")
			vline(-0.75)
			vline(0)
		
		
		
		tau = y**2 - self.profile.gamma
		#graph(log10(y), log10(self.profile.psi_tau(tau)), color="red", linestyle="dash")
		print tau
		
		
		dlogr = logy[1]-logy[0]
		logslope = (log10(rho[1:]) - log10(rho[0:-1])) / dlogr
		select(0, 1)
		graph(log10(y), logslope)
		labels("log r", "dlog rho / dlogr") 
		ylim(-5, 1)
		
		
		if 0:
			rhoinf = self.profile.potential_xy(0.1, 100)
			select(2,0)
			logy = (arange(n_y)) / (n_y-1.) * 2.5 - 1.5
			y = 10**logy
			x = 0.01
			#nus, las = self.profile.xy_to_conf(x, y)
			rho = self.profile.potential_xy(x, y)
			print "x", x
			print "y", y
			print "nu", nus
			print "la", las
			#rho = -rho
			rho = rho+rhoinf
			print "rho", rho
			rho = -rho
			print "rho", rho
			#dsa
			
			graph(log10(y), log10(rho))
			
			logx = (arange(n_x)) / (n_x-1.) * 2.5 - 1.5
			x = 10**logx
			y = 0.01
			rho = self.profile.potential_xy(x, y) + rhoinf
			rho = -rho
			print "x", x
			print "y", y
			print "nu", nus
			print "la", las
			print "rho", rho
			graph(log10(x), log10(rho), color="green")
			
			
		
		select(0, 0)
		
		if 1:
			#x = (arange(n_x) + 0.5) / (n_x) * 10
			logx = (arange(n_x)) / (n_x-1.) * 2.5 - 1.5
			x = 10**logx
			y = 0.
			nus, las = self.profile.xy_to_conf(x, y)
			rho = self.profile.rho_conf(nus, las)
			print "x", x
			print "y", y
			print "nu", nus
			print "la", las
			print "rho", rho
			graph(log10(x), log10(rho), color="green")
			if self.fit:
				rho_fit = self.fit.rho_xy(x, y)
				print "x,y",x, y
				print rho_fit
				graph(log10(x), log10(rho_fit), linestyle="dot", color="green")
			
			
			dlogr = logx[1]-logx[0]
			logslope = (log10(rho[1:]) - log10(rho[0:-1])) / dlogr
			select(0, 1)
			graph(log10(x), logslope, color="green")
			graph(log10(x), -2 * x / (x + 1) -1, color="blue")
			#x = r/self.r_2
			alpha = 0.27
			logslope = -2 * x**alpha
			graph(log10(x), logslope, color="green", linestyle="dot")
			
			
			hline(0)
			hline(-1)
			hline(-3)
			
		n_y = 30
		n_x = 30
		if 0:
			y = (arange(n_y) + 0.5) / (n_y) * 1
			x = (arange(n_x) + 0.5) / (n_x) * 1
			x, y = numpy.meshgrid(x, y)
			#nus, las = self.profile.xy_to_conf(x, y)
			rho = self.profile.potential_xy(x, y)
			rho = rho - rho.min() + 0.1
			print "pot", rho
			print rho.shape
			select(1,1)
			rho = rho/rho.max()
			#indexedimage(log10(rho), resize=[(0, 0), (1.0, 1.0)])
			n = 10
			#levels = (arange(n) + 0.5) / n * (rho.max() - rho.min()
			#contour(log10(rho.T), levels=levels, resize=[(0, 0), (10, 10)])
			contour(log10(rho), levels=40, resize=[(0, 0), (1.0, 1.0)])
				
			#, colormap 
		if 1:
			nx = n_y = 100
			y = (arange(n_y) + 0.5) / (n_y) * 2 - 1
			x = (arange(n_x) + 0.5) / (n_x) * 2 - 1
			x, y = numpy.meshgrid(x, y)
			nus, las = self.profile.xy_to_conf(x, y)
			rho = self.profile.rho_xy(x, y)
			#rho = x
			print rho.shape
			select(1,0)
			rho = rho#/rho.max()
			#indexedimage(log10(rho), resize=[(0, 0), (1.0, 1.0)])
			n = 10
			#levels = (arange(n) + 0.5) / n * (rho.max() - rho.min()
			#contour(log10(rho.T), levels=levels, resize=[(0, 0), (10, 10)])
			levels = 12
			c = contour(log10(rho), levels=levels, resize=[(-1, -1), (1.0, 1.0)])
			print log10(rho).min(), log10(rho).max()
			print "levels", c.levels
			levels = c.levels
			if self.fit:
				rho = self.fit.rho_xy(x, y)
				print log10(rho).min(), log10(rho).max()
				contour(log10(rho), levels=levels, resize=[(-1, -1), (1.0, 1.0)], color="blue", linestyle="dot")
				#contour(log10(rho), levels=[-0.11848672], resize=[(-1, -1), (1.0, 1.0)], color="blue", linestyle="dot")
			circle(0,0, 10**-0.75, color="black", linewidth="3px")
			circle(0,0, 10**0.0, color="black", linewidth="3px")
			labels("x", "y")
			#, colormap 
			
			
		draw()