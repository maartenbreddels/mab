# -*- coding: utf-8 -*-
from mab.gd.schw.cmdutils import *
import sys
import mab.gd.potential

import mab.gd.logging as logging
logger = logging.getLogger("gd.utils.mass_concentration")

def c_maccio(Mvir, h=0.72):
	return 10**(1.020 - 0.109 * (log10(Mvir*h) - 12))

class Maccio(object):
	def __init__(self, mass_range, scale_range, dirname):
		self.mass_range = mass_range
		self.scale_range = scale_range
		self.dirname = dirname
		self.filename_concentration = os.path.join(self.dirname, "concentration.npy")
		self.filename_M200 = os.path.join(self.dirname, "M200.npy")
		self.filename_pdf = os.path.join(self.dirname, "mass-concentration-pdf.npy")
		self.h = 0.72# hubble param

	def run(self, args, opts, scope):
		self.calculate()
		self.save()
	
	def save(self):
		logger.info("saving concentration to: %s" % self.filename_concentration)
		save(self.filename_concentration, self.c_grid)
		logger.info("saving M200 to: %s" % self.filename_M200)
		save(self.filename_M200, self.M200_grid)
		logger.info("saving pdf to: %s" % self.filename_pdf)
		save(self.filename_pdf, self.pdf)
	
	def load(self):
		logger.info("loading concentration from: %s" % self.filename_concentration)
		self.c_grid = load(self.filename_concentration)
		logger.info("loading M200 from: %s" % self.filename_M200)
		self.M200_grid = load(self.filename_M200)
		logger.info("loading pdf from: %s" % self.filename_pdf)
		self.pdf = load(self.filename_pdf)
		
	def calculate(self):
		logM1kpc_min, logM1kpc_max =  self.mass_range.min, self.mass_range.max
		logrs_min, logrs_max =  self.scale_range.min, self.scale_range.max
		n_M1kpc = 100
		n_rs = 100
		logM1kpcs = arange(n_M1kpc) / (n_M1kpc+1.0) * (logM1kpc_max - logM1kpc_min) + logM1kpc_min
		M1kpcs = 10**logM1kpcs
		logrss = arange(n_rs) / (n_rs+1.0) * (logrs_max - logrs_min) + logrs_min
		rss = 10**logrss
		#if not os.path.exists(self.filename_concentration):
		if 1:
			c_grid = zeros((n_M1kpc, n_rs))
			M200_grid = zeros((n_M1kpc, n_rs))
			for i in range(n_M1kpc):
				for j in range(n_rs):
					M1kpc = M1kpcs[i]
					rs = rss[j]
					nfw = mab.gd.potential.NFW1kpc(M1kpc, rs)
					M200_grid[i,j] = nfw.M200
					c_grid[i,j] = nfw.c
			#print "saving", filename
		#else:
		#	print "loading", filename
		#	c_grid = load(self.filename_concentration)
		self.c_grid = c_grid
		self.M200_grid = M200_grid
		print "cgrid", c_grid.min(), c_grid.max()
		#lalala
		#7.6, 8.3)),  ("nfw.rs", "rs_%06e", lambda x: 10**x, (-1, 1)
		
		#Mvir = 3e10
		#print "Mvir: %e" % Mvir
		h = 0.72
		
		def c_maccio(Mvir):
			return 10**(1.020 - 0.109 * (log10(Mvir*self.h) - 12))
		def c_maccio_up(Mvir):
			return exp(log(c_maccio(Mvir)) + 0.4)
		def c_maccio_down(Mvir):
			return exp(log(c_maccio(Mvir)) - 0.4)
		
		def p_c_maccio(c, Mvir):
			return exp( -0.5 * (log(c) - log(c_maccio(Mvir)))**2/0.33**2)
		
		p_c_grid = zeros((n_M1kpc, n_rs))
		#filename = os.path.join(extra.modelpath, "data", "p_c_grid.npy")
		#if not os.path.exists(self.filename_pdf):
		if 1:
			for i in range(n_M1kpc):
				for j in range(n_rs):
					#p_c_grid[i,j] = p_c_maccio(c_grid[i,j], mab.gd.potential.NFW1kpc(M1kpcs[i], rss[j]).M200)
					p_c_grid[i,j] = p_c_maccio(c_grid[i,j], M200_grid[i,j])
		#	save(filename, p_c_grid)
		#else:
		#	p_c_grid = load(self.filename_pdf)

		def c_mc2011(Mvir):
			#Mvir = mab.gd.potential.NFW1kpc(M1kpc, 0.5).M200
			z = 0
			m = 0.097
			w = 0.029
			gamma = 16.885
			alpha = -110.001
			beta = 2469.720
			bz = alpha / (z+gamma) + beta / (z+gamma)**2
			az = w*z -m
			print bz, az 
			c = 10**(az * log10(Mvir*h) + bz)
			#c = e**(az * log(Mvir*h) + bz)
			return c
			#print "c=",c, 10**c
			#print "az", az
			#sys.exit(0)
		print "c =", c_mc2011(1e8)
		self.pdf = p_c_grid
	def c(self, Mvir):
		return 10**(1.020 - 0.109 * (log10(Mvir*self.h) - 12))
		
if 0:
	parametersetcollect, extra = program("parametersetcollect", kaplot=True, iteration=0, schw=True)
	parametersetcollect.load(int(extra.opts.iteration))

	galaxy = None
	try:
		galaxy = mab.gd.configure.createobj_fromini(extra.ini, "galaxy", scope=extra.scope)
	except:
		print "no galaxy model known"
		
	light_model = mab.gd.configure.createobj_fromini(extra.ini, "light_model", scope=extra.scope)
	#galaxy = scope["galaxy"]

	#box()'

	def test(grid, resize):
		(xmin, ymin), (xmax, ymax) = resize
		nx, ny = grid.shape
		
		contours = zeros((3, 2, nx))
		medians = zeros((nx))
		for i in range(nx):
			ps = grid[i]
			cumps = cumsum(ps)
			cumps /= cumps.max()
			#if i == nx/2:
			#	print cumps
			#	import pdb
			#	pdb.set_trace() 
			#print cumps[-1]
			x = xmin + float(i)/nx * (xmax-xmin)
			def do(value, **kwargs):
				median_index = argmin(abs(cumps-value))
				if cumps[median_index] < value:
					i1 = median_index
				else:
					i1 = max(0, median_index-1)
				
				i2 = where(cumps>min(1, cumps[i1]+1e-16))[0][0]
				deltay = (cumps[i2] - cumps[i1])
				if median_index == 0: 
					deltax = (cumps[0]-cumps[i1])/deltay * (i2-i1)
				else:
					deltax = (value-cumps[i1])/deltay * (i2-i1)
				#print deltax, (i2-i1), cumps[i1], cumps[i2], i1, i2
				if deltax > 1e30:
					import pdb
					pdb.set_trace()
				assert deltax >= 0
				if i1 == i2 == 0:
					y = ymin
				else:
					#print float(median_index)/nx
					y = ymin + float(median_index+deltax)/ny * (ymax-ymin)
				#print y
				return y
				#symbol(x, y, **kwargs)
				
			
			sigmas = [0.682689492137 , 0.954499736104, 0.997300203937]
			#color = 
			for j, sigma in enumerate(sigmas):
				contours[j, 0,i] = do(0.5-sigma/2, color="yellow")
				contours[j, 1,i] = do(0.5+sigma/2, color="yellow")
			medians[i] = do(0.5, color="blue")
		x = arange(nx)/(float(nx)-1) * (xmax-xmin) + xmin
		alphas = [0.9, 0.6, 0.3]
		for i in range(3)[::-1]:
			color = Color.blue.blend(Color.white, alphas[i])
			obj = fillrange(x, contours[i,0], contours[i,1], color=color)
			if i == 0:
				ci_fillrange = obj  
		#for i in range(3):
		#	graph(x, contours[i,0], color="lightgrey")
		#	graph(x, contours[i,1], color="lightgrey")
		#print contours[0,1]
		#print medians
		mediangraph = graph(x, medians, color="green")
		return mediangraph, ci_fillrange 
		
			
		
		
	if parametersetcollect.dimension == 2:
		document(size="15cm,20cm")
		mozaic(2,3,box)
		
		def makeresize(i, j):
			return (parametersetcollect.probability_range[i][0], parametersetcollect.probability_range[j][0]), (parametersetcollect.probability_range[i][1], parametersetcollect.probability_range[j][1])
		
		select(0, 0)
		
		if galaxy is None:
			filename = os.path.join(extra.modelpath, "data", "c_grid.npy")
			M1kpc_min, M1kpc_max =  parametersetcollect.probability_range[0]
			rs_min, rs_max =  parametersetcollect.probability_range[1]
			n_M1kpc = 100
			n_rs = 100
			logM1kpcs = arange(n_M1kpc) / (n_M1kpc+1.0) * (M1kpc_max - M1kpc_min) + M1kpc_min
			M1kpcs = 10**logM1kpcs
			logrss = arange(n_rs) / (n_rs+1.0) * (rs_max - rs_min) + rs_min
			rss = 10**logrss
			if not os.path.exists(filename):
				c_grid = zeros((n_M1kpc, n_rs))
				for i in range(n_M1kpc):
					for j in range(n_rs):
						M1kpc = M1kpcs[i]
						rs = rss[j]
						c_grid[i,j] = mab.gd.potential.NFW1kpc(M1kpc, rs).c
				print "saving", filename
				save(filename, c_grid)
			else:
				print "loading", filename
				c_grid = load(filename)
			print c_grid.min(), c_grid.max()
			#lalala
			#7.6, 8.3)),  ("nfw.rs", "rs_%06e", lambda x: 10**x, (-1, 1)
			
			#Mvir = 3e10
			#print "Mvir: %e" % Mvir
			
			h = 0.72
			def c_maccio(Mvir):
				return 10**(1.020 - 0.109 * (log10(Mvir*h) - 12))
			def c_maccio_up(Mvir):
				return exp(log(c_maccio(Mvir)) + 0.4)
			def c_maccio_down(Mvir):
				return exp(log(c_maccio(Mvir)) - 0.4)
			
			def p_c_maccio(c, Mvir):
				return exp( -0.5 * (log(c) - log(c_maccio(Mvir)))**2/0.33**2)
			
			p_c_grid = zeros((n_M1kpc, n_rs))
			filename = os.path.join(extra.modelpath, "data", "p_c_grid.npy")
			if not os.path.exists(filename):
				for i in range(n_M1kpc):
					for j in range(n_rs):
						p_c_grid[i,j] = p_c_maccio(c_grid[i,j], mab.gd.potential.NFW1kpc(M1kpcs[i], rss[j]).M200)
				save(filename, p_c_grid)
			else:
				p_c_grid = load(filename)
			
				
			def c_mc2011(Mvir):
				#Mvir = mab.gd.potential.NFW1kpc(M1kpc, 0.5).M200
				z = 0
				m = 0.097
				w = 0.029
				gamma = 16.885
				alpha = -110.001
				beta = 2469.720
				bz = alpha / (z+gamma) + beta / (z+gamma)**2
				az = w*z -m
				print bz, az 
				c = 10**(az * log10(Mvir*h) + bz)
				#c = e**(az * log(Mvir*h) + bz)
				return c
				#print "c=",c, 10**c
				#print "az", az
				#sys.exit(0)
			print "c =", c_mc2011(1e8)
			
			
		resize = makeresize(0,1)
		#probimage2d(parametersetcollect.probability_grid, 0, 1, resize=resize, colormap="whiteblue", drawcontourlines=True)
		fill = False
		probimage2d(parametersetcollect.probability_grid, 0, 1, resize=resize, colormap=None, color="blue", drawcontourlines=True, premultiply=True, fill=fill)
		#probimage2d(p_c_grid, 0, 1, resize=resize, color="green", drawcontourlines=True, fill=False, colormap=None, linewidth="2pt", linestyle="dot")
		if galaxy is None:
			probimage2d(parametersetcollect.probability_grid * p_c_grid, 0, 1, resize=resize, color="black", drawcontourlines=True, fill=False, colormap=None, linewidth="1pt", alpha=1.0)
		#scatter(parametersetcollect.parameter_points[0], parametersetcollect.parameter_points[1], color="black", symbolsize="12pt")
		#levels = [5, 10, 15, 20, 35, 50, 100, 150]
		levels = arange(10, 60, 10)
		#linestyle = ["normal", "dot",
		if not galaxy: 
			clearautolegend()
			for level, linestyle, color in zip(levels, alllinestyles*3, (goodcolors*3)[1:]):
				contour(c_grid.transpose(), levels=[level], resize=resize, color=color, linestyle=linestyle, addlegend=level==10 or level==50)
			
			n_M200 = 100
			M200_min, M200_max = 5, 12
			#10**arange(
			logM200s = arange(n_M200) / (n_M200+1.0) * (M200_max - M200_min) + M200_min
			logM1kpcs = []
			rss = []
			logM1kpcs_up = []
			rss_up = []
			logM1kpcs_down = []
			rss_down = []
			print "c(1e10) =", c_mc2011(1e10)
			print "c(1e12) =", c_mc2011(1e12)
			print log10(mab.gd.potential.NFW1kpc(10**7.9, 1.).M200)
			#fsjklsa
			for i in range(n_M200):
				M200 = 10**logM200s[i]
				c = c_mc2011(M200)
				c = c_maccio(M200)
				pot = mab.gd.potential.NFW(M200, 1.)
				rs = pot.r200/c
				rss.append(rs)
				pot = mab.gd.potential.NFW(M200, rs)
				M1kpc = pot.enclosed_mass(1.0)
				logM1kpcs.append(log10(M1kpc))
				
				c = c_maccio_up(M200)
				pot = mab.gd.potential.NFW(M200, 1.)
				rs = pot.r200/c
				rss_up.append(rs)
				pot = mab.gd.potential.NFW(M200, rs)
				M1kpc = pot.enclosed_mass(1.0)
				logM1kpcs_up.append(log10(M1kpc))
				
				c = c_maccio_down(M200)
				pot = mab.gd.potential.NFW(M200, 1.)
				rs = pot.r200/c
				rss_down.append(rs)
				pot = mab.gd.potential.NFW(M200, rs)
				M1kpc = pot.enclosed_mass(1.0)
				logM1kpcs_down.append(log10(M1kpc))

				
			rss = array(rss)
			logM1kpcs = array(logM1kpcs)
			graph(logM1kpcs, log10(rss), linewidth="3pt", linestyle="dot")
			#graph(logM1kpcs_up, log10(rss_up), linewidth="3pt", linestyle="dot")
			#graph(logM1kpcs_down, log10(rss_down), linewidth="3pt", linestyle="dot")
			xlim(resize[0][0], resize[1][0])
			ylim(resize[0][1], resize[1][1])
			pot = mab.gd.potential.NFW1kpc(10**8.0, 10**1.00327)
			print "c = %e %e" % (pot.c, pot.M200)
			#sys.exit(0)
			
			autolegend(*["c=10", "c=50", "Eq. 27", ], location="left,top", edgespacing="3mm")
		#if galaxy:
			#if hasattr(galaxy.profile_model.dm_profile, "M1kpc"):
			#	symbol(log10(galaxy.profile_model.dm_profile.M1kpc), log10(galaxy.profile_model.dm_profile.rs), color="red") #symbolName="plus")
			#	print log10(galaxy.profile_model.dm_profile.rs), galaxy.profile_model.dm_profile.rs
		labels("M<sub>1kpc</sub> - solar mass" , "log10(<i>r</i><sub>s</sub>) - kpc")
		
		
		print parametersetcollect.probability_grid.max()
		if 1:
			select(0, 1)
			probgraph(parametersetcollect.probability_grid, 0, resize=parametersetcollect.probability_range[0], name="M<sub>1kpc</sub>")
			labels("M<sub>1kpc</sub> - solar mass" , "p")
			if galaxy:	
				if hasattr(galaxy.profile_model.dm_profile, "M1kpc"):
					vline(log10(galaxy.profile_model.dm_profile.M1kpc), color="red", linestyle="dash")
				
			select(0, 2)
			probgraph(parametersetcollect.probability_grid, 1, resize=parametersetcollect.probability_range[1], name="log10(<i>r</i><sub>s</sub>)")
			labels("log10(<i>r</i><sub>s</sub>) - kpc" , "p")	
			if galaxy:	
				if hasattr(galaxy.profile_model.dm_profile, "rs"):
					vline(log10(galaxy.profile_model.dm_profile.rs), color="red", linestyle="dash")
			
		
		if 1:
			select(1,0)
			resize = parametersetcollect.anisotropy_range
			indexedimage(transpose(parametersetcollect.anisotropy_grid), resize=resize)
			#test(parametersetcollect.anisotropy_grid, resize=resize)
			xlim(resize[0][0], resize[1][0])
			ylim(resize[0][1], resize[1][1])
			if galaxy:	
				hline(galaxy.beta, color="red", linestyle="dash")
			
			labels("r - kpc", "anisotropy")
				
			select(1,1)
			
			resize = parametersetcollect.mass_enclosed_range
			x = transpose((parametersetcollect.mass_enclosed_grid))
			
			if 1: #galaxy is not None:
				logrmin = resize[0][0]
				logrmax = resize[1][0]
				n = 100
				logr = arange(n)/(float(n)-1) * (logrmax-logrmin) + logrmin
				r = 10**logr
				rs = r
			
			
			#print x.max()
			#limit = -5
			#x[x-x.max() < limit] = limit 
			#indexedimage(x, resize=resize)
			test(parametersetcollect.mass_enclosed_grid, resize=resize)
			if galaxy is not None:
				graph(logr, log10(galaxy.profile_model.dm_profile.enclosed_mass(r)), color="red", linestyle="dash")
			y = array([log10(light_model.cumdensityr(0, r_)) for r_ in rs])
			clearautolegend()
			graph(logr, y, color="black", linestyle="dot")
			ylim(3, 10)
			xlim(-3,1)
			autolegend("Stellar mass", location="left,top", edgespacing="3mm")
			labels("log <i>r</i> - kpc", "log10 M(&lt;r) - solar mass")
			
			select(1,2)
			
			resize = parametersetcollect.logslope_range
			#indexedimage(transpose(parametersetcollect.logslope_grid), resize=resize)
			obj1, obj2 = test(parametersetcollect.logslope_grid, resize=resize)
			if galaxy is not None:
				obj1, obj2 = test(parametersetcollect.logslope_grid, resize=resize)
				obj0 = graph(logr, galaxy.profile_model.dm_profile.logslope(r), color="red", linestyle="dash")
				legend(["line", "line", "squaresolid"], ["true value", "median", "confidence intervals"], [obj0, obj1, obj2], location="left, top")
			else:
				legend(["line", "squaresolid"], ["median", "confidence intervals"], [obj1, obj2], location="left, top")
			xlim(resize[0][0], resize[1][0])
			ylim(resize[0][1], resize[1][1])
			labels("log <i>r</i> - kpc", "dlog &rho; / dlog <i>r</i>") 
			#symbol(8, 0.5, color="red") #
		
	if parametersetcollect.dimension == 3:
		document(size="20cm,20cm")
		mozaic(3,3,box)
		
		def makeresize(i, j):
			return (parametersetcollect.probability_range[i][0], parametersetcollect.probability_range[j][0]), (parametersetcollect.probability_range[i][1], parametersetcollect.probability_range[j][1])
		
		select(0, 1)
		resize = makeresize(0,1)
		probimage2d(parametersetcollect.probability_grid, 0, 1, resize=resize, colormap="whiteblue", fill=False, color="blue", drawcontourlines=True, premultiply=True)
		scatter(parametersetcollect.parameter_points[0], parametersetcollect.parameter_points[1], color="black", symbolsize="5pt")
		if galaxy:
			#symbol(log10(galaxy.profile_model.dm_profile.M1kpc), log10(galaxy.profile_model.dm_profile.rs), color="red") #symbolName="plus")
			#print log10(galaxy.profile_model.dm_profile.rs), galaxy.profile_model.dm_profile.rs
			print
		labels("M<sub>1kpc</sub> - solar mass" , "log10(<i>r</i><sub>s</sub>) - kpc")
		
		select(0, 0)
		resize = makeresize(0,2)
		probimage2d(parametersetcollect.probability_grid, 0, 2, resize=resize, colormap="whiteblue", fill=False, color="blue", drawcontourlines=True, premultiply=True)
		scatter(parametersetcollect.parameter_points[0], parametersetcollect.parameter_points[2], color="black", symbolsize="5pt")
		labels("M<sub>1kpc</sub> - solar mass" , "alpha")
		
		select(1, 0)
		resize = makeresize(1,2)
		probimage2d(parametersetcollect.probability_grid, 1, 2, resize=resize, colormap="whiteblue", fill=False, color="blue", drawcontourlines=True, premultiply=True)
		scatter(parametersetcollect.parameter_points[1], parametersetcollect.parameter_points[2], color="black", symbolsize="5pt")
		labels("rs" , "alpha")
		
		print parametersetcollect.probability_grid.max()
		if 1:
			select(0, 2)
			probgraph(parametersetcollect.probability_grid, 0, resize=parametersetcollect.probability_range[0], name="M<sub>1kpc</sub>")
			labels("M<sub>1kpc</sub> - solar mass" , "p")
			#if galaxy:	
			#	vline(log10(galaxy.profile_model.dm_profile.M1kpc), color="red", linestyle="dash")
				
			select(1, 1)
			probgraph(parametersetcollect.probability_grid, 1, resize=parametersetcollect.probability_range[1], name="log10(<i>r</i><sub>s</sub>)")
			labels("log10(<i>r</i><sub>s</sub>) - kpc" , "p")	
			#if galaxy:	
			#	vline(log10(galaxy.profile_model.dm_profile.rs), color="red", linestyle="dash")
			select(2, 0)
			probgraph(parametersetcollect.probability_grid, 2, resize=parametersetcollect.probability_range[2], name="alpha")
			labels("alpha" , "p")	
			#if galaxy:	
			#	vline(log10(galaxy.profile_model.dm_profile.rs), color="red", linestyle="dash")
			
		
		if 1:
			select(2,1)
			resize = parametersetcollect.anisotropy_range
			#indexedimage(transpose(parametersetcollect.anisotropy_grid), resize=resize)
			test(parametersetcollect.anisotropy_grid, resize=resize)
			xlim(resize[0][0], resize[1][0])
			ylim(resize[0][1], resize[1][1])
			if galaxy:	
				hline(galaxy.beta, color="red", linestyle="dash")
			
			labels("r - kpc", "anisotropy")
				
			select(2,2)
			
			resize = parametersetcollect.mass_enclosed_range
			x = transpose((parametersetcollect.mass_enclosed_grid))
			
			if 1: #galaxy is not None:
				logrmin = resize[0][0]
				logrmax = resize[1][0]
				n = 100
				logr = arange(n)/(float(n)-1) * (logrmax-logrmin) + logrmin
				r = 10**logr
				rs = r
			
			
			#print x.max()
			#limit = -5
			#x[x-x.max() < limit] = limit 
			#indexedimage(x, resize=resize)
			test(parametersetcollect.mass_enclosed_grid, resize=resize)
			if galaxy is not None:
				graph(logr, log10(galaxy.profile_model.dm_profile.enclosed_mass(r)), color="red", linestyle="dash")
			y = array([log10(light_model.cumdensityr(0, r_)) for r_ in rs])
			clearautolegend()
			graph(logr, y, color="black", linestyle="dot")
			ylim(3, 10)
			xlim(-3,1)
			autolegend("Stellar mass", location="left,top", edgespacing="3mm")
			labels("log <i>r</i> - kpc", "log10 M(&lt;r) - solar mass")
			
			select(1,2)
			
			resize = parametersetcollect.logslope_range
			#indexedimage(transpose(parametersetcollect.logslope_grid), resize=resize)
			obj1, obj2 = test(parametersetcollect.logslope_grid, resize=resize)
			if galaxy is not None:
				obj1, obj2 = test(parametersetcollect.logslope_grid, resize=resize)
				obj0 = graph(logr, galaxy.profile_model.dm_profile.logslope(r), color="red", linestyle="dash")
				legend(["line", "line", "squaresolid"], ["true value", "median", "confidence intervals"], [obj0, obj1, obj2], location="left, top")
			else:
				legend(["line", "squaresolid"], ["median", "confidence intervals"], [obj1, obj2], location="left, top")
			xlim(resize[0][0], resize[1][0])
			ylim(resize[0][1], resize[1][1])
			labels("log <i>r</i> - kpc", "dlog &rho; / dlog <i>r</i>") 
			#symbol(8, 0.5, color="red") #
		



	draw()

