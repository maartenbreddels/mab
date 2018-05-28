# -*- coding: utf-8 -*-
from kaplot import *
from mab.binningtools import bingrid, binrange
import mab.gd.logging as logging
import scipy.stats
logger = logging.getLogger("gd.schw.plot")

class OrbitLibraryLosvdPlot(object):
	def __init__(self, storage_2d):
		self.storage_2d = storage_2d
		
	def run(self, args, opts, scope):
		self.storage_2d.load()
		#OrbitLibraryLosvdPlot
		nE, nL = self.storage_2d.dfgrid.n_I1, self.storage_2d.dfgrid.n_I2
		print self.storage_2d.losvds.shape
		#mozaic(
		
		
class ProjectionMatrix(object):
	def __init__(self, projection_matrix):
		self.projection_matrix = projection_matrix
		
	def run(self, args, opts, scope):
		self.projection_matrix.load()
		I = self.projection_matrix.matrices[0]
		I = log10(I)
		I[I<-5] = -5
		bins = I.shape[0]
		print bins
		resize = (0, scope["schw_logr_min"]), (bins-1, scope["schw_logr_max"])
		box() 
		indexedimage(I.T, colormap="whiteblack", resize=resize)
		draw()
		#document(size="40cm,25cm")
	
class OrbitMoments(object):
	def __init__(self, storage_3d):
		self.storage_3d = storage_3d
		
	def run(self, args, opts, scope):
		self.storage_3d.load()
		document(size="40cm,25cm")
		moments = self.storage_3d.moments3d
		moment_nr = self.storage_3d.v02index
		moment_indices = [self.storage_3d.v00index, self.storage_3d.v20index, self.storage_3d.v02index]
		
		dfgrid = self.storage_3d.dfgrid
		nE, nL = dfgrid.n_I1, dfgrid.n_I2
		#mozaic(nE, nL, container)
		mozaic(nL, len(moment_indices), box)
		if 1:
			for j in range(nL):
				moment_indices = [self.storage_3d.v00index, self.storage_3d.v20index, self.storage_3d.v02index]
				for m, moment_index in enumerate(moment_indices):
					select(j, m)
					for i in range(nE)[::4]:
						#for color, m in zip(goodcolors, ):
						y = moments[i+j*nE][m]
						if any(y>0):
							y /= y.max()
						color = "red"
						x = self.storage_3d.x
						graph(x, y, color=color)
		else:
			for i in range(nE):
				for j in range(nL):
					select(i, j)
					border()
					spacer()
					for color, m in zip(goodcolors, ):
						x = moments[i+j*nE][m]
						if any(x>0):
							x /= x.max()
						graph(x, color=color)
		draw()
		print moments.shape

class Density2dCheck(object):
	def __init__(self, storage2d, weights, photometry):
		self.storage2d = storage2d
		self.weights = weights
		self.photometry = photometry
		
	def run(self, args, opts, scope):
		self.storage2d.load()
		self.weights.load()
		
		mozaic(2,2, box)
		select(0, 0)
		print dir(self.weights)
		w = self.weights.orbitweights
		print w.shape
		masses = self.storage2d.masses
		print masses.shape
		shape = masses.shape
		newshape = (shape[0] * shape[1], ) + shape[2:]
		print "shape:", shape, newshape
		masses = masses.reshape(newshape) * 1.
		
		masses = dot(masses.T, w).T
		
		max_order = -3
		
		
		if 0:
			if 0:
				masses = sum(masses, axis=0)
			else:
				masses = masses[2]
		print masses.shape
		resize = self.storage2d.masses_resize
		mass_schw = masses.reshape((self.storage2d.NR, self.storage2d.NR))
		logmass_schw = log10(mass_schw)
		logmass_schw -= logmass_schw.max()
		logmass_schw[logmass_schw<max_order] = max_order
		
		indexedimage(logmass_schw.T, resize=resize)
		contour(logmass_schw.T, levels=10, resize=resize)
		if 0:
			x,y = mgrid[-self.storage2d.Rmax:self.storage2d.Rmax: self.storage2d.NR*1j,-self.storage2d.Rmax:self.storage2d.Rmax: self.storage2d.NR*1j]
			#print x
			#print x.shape
			R = sqrt(x**2+y**2)
			mass = self.light_model.densityR(R)
		else:
			self.photometry.load()
			mass = self.photometry.grid2d
		select(1,0)
		logmass = log10(mass)
		logmass -= logmass.max()
		logmass[logmass<max_order] = max_order
		indexedimage(logmass.T, resize=resize)
		contour(logmass.T, resize=resize, levels=10)
		select(0,1)
		contour(logmass.T, resize=resize, levels=10)
		contour(logmass_schw.T, resize=resize, levels=10, color="red", linestyle="dash")
		print sum(mass_schw), sum(mass)
		draw()
		

class PlotSchwSolutionLosvdData(object):
	def __init__(self, solution, schwmodel, storage_2d, storage_losvd, binned_data_m2, binned_data_m4, df_samples=None, galaxy=None, ratios=None, foreground_model=None, mean_vsys=None, simplemodel=None, member_filter=None):
		self.solution = solution
		self.schwmodel = schwmodel
		self.galaxy = galaxy
		self.storage_2d = storage_2d
		self.storage_losvd = storage_losvd
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.df_samples = df_samples
		self.member_ratios = ratios
		self.foreground_model = foreground_model
		self.mean_vsys = mean_vsys
		self.simplemodel = simplemodel
		self.member_filter = member_filter
		
	def run(self, args, opts, scope):
		print self.solution
		self.solution.load()
		self.member_ratios.load()
		#self.storage_2d.aperture.load()
		#self.schwmodel.storage_2d.aperture.load()
		dirname = os.path.join(self.binned_data_m2.modelpath, "schw", "aperture")
		cached_filename= os.path.join(dirname, "observation_group" +self.binned_data_m2.postfix + "_cache" +".npy")
		self.binned_data_m2.observation.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		self.binned_data_m2.aperture.load()
		self.binned_data_m4.aperture.load()
		self.storage_losvd.load()
		if self.df_samples:
			self.df_samples.load()
			
		if 0:
			numpy.random.seed(2)
			N = 100.
			x = numpy.random.normal(0., 1, N)
			#print scipy.stats.kstest(x, 'norm')
			D, p = scipy.stats.kstest(x, scipy.stats.norm(0., 1.2).cdf)
			print D, p
			print scipy.stats.ksone.cdf(D, N)
			print (1-scipy.stats.ksone.cdf(D, N))
			
			print "D =", D
			print (1-scipy.stats.ksone.cdf(D, N))*2
			print (1-scipy.stats.kstwobign.cdf(D*sqrt(N)))
			sys.exit(0)
			
			box()
			x = arange(0.0001, 1., 0.01)
			y = scipy.stats.ksone.pdf(x, N)
			print y
			graph(x, y)
			draw()
			
			import pdb;
			pdb.set_trace()
			
		stars_group = []
		no_groups = self.binned_data_m2.moments.shape[1]
		stars = self.binned_data_m2.observation.stars
		stars, indices = self.binned_data_m2.aperture.stars_to_indices(stars)
		for group in range(no_groups):
			stars_group.append([])
		clean_stars = []
		for star, index in zip(stars, indices):
			print ">", index
			if(index != -1) and (star.vlos > -self.storage_losvd.vmax) and (star.vlos < self.storage_losvd.vmax) :
				stars_group[index].append(star)
				clean_stars.append(star)
		import mab
		for group in range(no_groups):
			stars_group[group] = mab.cvsfile.CsvObject(stars_group[group])
		clean_stars = mab.cvsfile.CsvObject(clean_stars)
		if 0:
			if not os.path.exists(cached_filename):
				for index in range(no_groups):
					print index
					stars = self.binned_data_m2.get_stars_from_aperture(index)
					print max([star.vlos for star in stars])
					print len(stars)
					stars_group.append(stars)
				
				logger.info("writing cache: " + cached_filename)
				numpy.save(cached_filename, stars_group)
			else:
				logger.info("reading cached: " + cached_filename)
				stars_group = numpy.load(cached_filename).tolist()
				print stars_group
		#import pickle
		#pickle.
		#import pdb
		#pdb.set_trace()
		#da
		
		
		Rborders = arange(self.storage_losvd.NR+1) / (0.0+self.storage_losvd.NR) * (self.storage_losvd.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		#print Rborders
		rho2d_target = array([self.schwmodel.light_model.cumdensityR(R1, R2) for R1, R2 in zip(R1s, R2s)])
		
		
		losvds = self.storage_losvd.losvds
		print self.solution.orbitweights.shape, losvds.shape
		
		losvd = dot(losvds.T, self.solution.orbitweights).T
		mozaic(3,3,box)
		
		print self.storage_losvd.vborders
		print "*" * 70
		print self.storage_losvd.vcenters
		cumulative = True
		cumulative = False
		
		mean_sigma = 2.
		v1 = self.mean_vsys-self.storage_losvd.vmax
		v2 = self.mean_vsys+self.storage_losvd.vmax
		print self.simplemodel
		
		f_m_in = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v1, v2)[0]
		f_n_in = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v1, v2)[0]
		f_m_out = scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.galaxy_velocity_model(v, mean_sigma), v2, inf)[0]
		f_n_out = scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), -inf, v1)[0] + \
			scipy.integrate.quad(lambda v: self.simplemodel.foreground_velocity_model(v, mean_sigma), v2, inf)[0]
		weight = f_m_in/f_n_in # additional weight, since we are interested in the vmin,vmax region
		
		
		for j in range(3):
			for i in range(3):
				index = i + j * 3
				if index < no_groups:
					select(i, j)
					title("bin = %d" % index)
					stars = stars_group[index]
					#scatter(stars.rc, stars.vlos)
					sigma_v = 2.
					vs = self.storage_losvd.vcenters
					evs = vs * 0 + sigma_v
					dv = delta_v = self.storage_losvd.delta_v
					#losvd_fg = array([exp(self.foreground_model.logL(v+self.mean_vsys, e_v)) for v, e_v in zip(vs, evs)])
					losvd_fg = array([self.simplemodel.foreground_velocity_model(v+self.mean_vsys, e_v) for v, e_v in zip(vs, evs)])
					losvd_fg /= sum(losvd_fg*dv)
					if self.member_ratios:
						y = losvd[:,0] * 0.
						y_fg = y * 0.
						N = 0
						for star in stars:
							Rkpc = self.storage_losvd.light_model.arcsec_to_kpc(star.rc)
							Rindex = Rkpc * (0.0+self.storage_losvd.NR) /self.storage_losvd.Rmax
							Rindex = int(Rindex)
							#print Rindex, Rkpc
							#ratio = self.member_ratios.ratios1(Rkpc)
							#self.ratios.append(ratio)
							Rarcsec = star.rc
							ratio = 0
							for i in range(self.simplemodel.catalogues):
								if int(star.catalogue_mask) & (1<<i):
									f_member = self.simplemodel.memberRe(Rarcsec, weight=weight, catalogue=i) # i)
									f_non_member = self.simplemodel.non_memberRe(Rarcsec, weight=weight, catalogue=i)#, i)
									#graph(Rarcsec, log10(f_member), color=color)
									#graph(Rarcsec, log10(f_non_member), color=color, linestyle="dash")
									ratio = max(f_member/f_non_member, ratio)
									
							
							losvd_at_R = losvd[:,Rindex]
							losvd_at_R = scipy.ndimage.gaussian_filter(losvd_at_R, [sigma_v/delta_v], mode='constant')
							losvd_at_R = losvd_at_R / sum(losvd_at_R * delta_v)
							p_member = ratio / (1+ratio)
							#print ratio, sum(losvd_fg*dv), sum(losvd_at_R*dv)
							#print Rkpc, p_member
							y += losvd_at_R * p_member + losvd_fg * (1-p_member)
							y_fg += losvd_fg * (1-p_member)
							N += 1
						y_fg /= N
						y /= N
					else:
						Rcenter = self.binned_data_m2.aperture.aperture_rcenters_kpc[index]
						
						Rindex = Rcenter * (0.0+self.storage_losvd.NR) /self.storage_losvd.Rmax
						Rindex = int(Rindex)
						losvd_at_R = losvd[:,Rindex]
						delta_v = self.storage_losvd.delta_v
						losvd_at_R = scipy.ndimage.gaussian_filter(losvd_at_R, [sigma_v/delta_v], mode='constant')
						y = losvd_at_R/losvd_at_R.max()
						y /= sum(y)
						dv = self.storage_losvd.vcenters[1] - self.storage_losvd.vcenters[0]
						y /= dv
					
					if cumulative:
						ycum = cumsum(y)
						ycum /= max(ycum)
						graph(self.storage_losvd.vcenters, ycum)
						print self.storage_losvd.vcenters.shape, ycum.shape
						cdf = scipy.interpolate.interp1d(self.storage_losvd.vcenters, ycum)
						x = self.storage_losvd.vcenters
						graph(x, cdf(x), color="red", linestyle="dash")
					else:
						graph(self.storage_losvd.vcenters, y)
						graph(self.storage_losvd.vcenters, y_fg, color="blue")
						v = self.storage_losvd.vcenters
						print "SIGMA", (sum(v**2*y)/sum(y))**0.5
					if cumulative:
						cumhistogram(stars.vlos, datamin=-50, datamax=50, binwidth=0.5, normalize=True)
						if 0:
							D, p = scipy.stats.kstest(stars.vlos, cdf)
							print index, "-" * 70
							N = len(stars.vlos)
							print "D, p", D, p
							print (1-scipy.stats.ksone.cdf(D, N))*2
							print (1-scipy.stats.kstwobign.cdf(D*sqrt(N)))
							title("D=%f p=%f" % (D, p))
							print "-" * 70
					else:
						histogram(stars.vlos, datamin=-50, datamax=50, binwidth=2., normalize=True)
					#print index, Rcenter, Rindex
					x = self.storage_losvd.vcenters
					print x.shape, y.shape
					m2 = sum(x**2*y*dv)
					m4 = sum(x**4*y*dv)
					sigma = m2**0.5
					gamma = m4/m2**2
					graph(x, gaussian(x, 0, sigma), color="red")
					print "sigma", sigma, "gamma", gamma
					
					
		select(1,2)
		losvd = losvd * 1
		for i in range(losvd.shape[1]):
			losvd[:,i] /= sum(losvd[:,i])
			
		indexedimage(losvd)
		
		
		select(2,2)
		if 0:
			allstars = []
			for group in stars_group:
				allstars.extend(group)
			print len(allstars)
			_, Rindices = self.storage_losvd.stars_to_apertureindices(allstars)
			#print max(indices)
			Vindices = [self.storage_losvd.velocity_to_index(star.vlos) for star in allstars]
			N = len(Vindices)
			logL = 0
			for i in range(N):
				logL += log10(losvd[Vindices[i], Rindices[i]])
				
			print logL
			logLobs = logL
			
			if 0:
				M = 1000
				#grid = losvd * 0.
				logLs = []
				for j in range(M):
					logL = 0
					for i in range(N):
						Rindex = Rindices[i]
						cdf = cumsum(losvd[:,Rindex])
						cdf /= cdf.max()
						u = random.uniform()
						index = 0
						while u >= cdf[index]:
							index +=1
							#print index, u >= cdf[index], u, cdf[index]
						logL += log10(losvd[index, Rindex])
						if losvd[index, Rindex] == 0:
							import pdb;
							pdb.set_trace()
						#grid[index, i] += 1
					logLs.append(logL)
					
				logLs = array(logLs)
				print logLs.mean(), logLs.std()
				histogram(logLs, datamin=logL-60, datamax=logL+60, bincount=100)
				vline(logLobs, color="green")
				#indexedimage(grid)
			
			
			
		else:
			N = 100000
			grid = losvd * 0.
			for i in range(losvd.shape[1]):
				cdf = cumsum(losvd[:,i])
				cdf /= cdf.max()
				#print cdf
				for j in range(N):
					u = random.uniform()
					index = 0
					while u >= cdf[index]:
						index +=1
						#print index, u >= cdf[index], u, cdf[index]
					grid[index, i] += 1
				
				#print cdf
			indexedimage(grid)
			
		
		
		
		draw()
		return
		
		select(0, 0)
		indexedimage(losvd)
		select(0, 2)
		for i in range(losvd.shape[1])[:4]:
			y = losvd[:,i]
			if y.max() > 0:
				graph(y/y.max())
				
		select(0, 1)
		indexedimage(losvd/rho2d_target)

		if self.df_samples:
			df_samples = self.df_samples.df_samples
		else:
			df_samples = None

		if df_samples:
			R = sqrt(df_samples.x**2 + df_samples.y**2)
			vlos = df_samples.vz
			#numpy.histogram2d(
			select(1,0)
			data = numpy.histogram2d(R, vlos, bins=[30, 30], range=[(0, 1.5), (-40, 40)])
			data = data[0].T
			indexedimage(data)
			#data = density2d(R, vlos, xmin=0, xmax=1.5, ymin=-40., ymax=40., bincountx=30, bincounty=30)
			#data = data.data2d
		
		select(1,1)
		if df_samples:
			contour(data, levels=10)
		contour(losvd, levels=10, color="red")
		
		print losvds.shape
		#losvd = 
		draw()
	
class PlotSchwSolutionLosvd(object):
	def __init__(self, solution, schwmodel, storage_2d, storage_losvd, binned_data_m2, binned_data_m4, df_samples=None, galaxy=None, kinematics=None):
		self.solution = solution
		self.schwmodel = schwmodel
		self.galaxy = galaxy
		self.storage_2d = storage_2d
		self.storage_losvd = storage_losvd
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.df_samples = df_samples
		self.kinematics = kinematics
		
	def run(self, args, opts, scope):
		print self.solution
		self.solution.load()
		self.storage_2d.aperture.load()
		#self.schwmodel.storage_2d.aperture.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		self.binned_data_m2.aperture.load()
		self.binned_data_m4.aperture.load()
		self.storage_losvd.load()
		if self.df_samples:
			self.df_samples.load()
		
		
		Rborders = arange(self.storage_losvd.NR+1) / (0.0+self.storage_losvd.NR) * (self.storage_losvd.Rmax)
		R1s = Rborders[0:-1]
		R2s = Rborders[1:]
		
		dRs = R2s - R1s
		delta_R = R2s[0] - R1s[0]
		assert all(abs(dRs - delta_R) < 1e-10), "no constant dR" 
		
		#print Rborders
		rho2d_target = array([self.schwmodel.light_model.cumdensityR(R1, R2) for R1, R2 in zip(R1s, R2s)])
		
		
		losvds = self.storage_losvd.losvds
		print self.solution.orbitweights.shape, losvds.shape
		
		losvd = dot(losvds.T, self.solution.orbitweights).T
		mozaic(4,4,box)
		select(0, 0)
		indexedimage(losvd)
		select(0, 3)
		for i in range(losvd.shape[1])[:4]:
			y = losvd[:,i]
			if y.max() > 0:
				graph(y/y.max())
				
		select(0, 2)
		indexedimage(self.solution.orbitweights.reshape(8, 20))
		select(0, 1)
		indexedimage(losvd/rho2d_target)
		
		if self.kinematics:
			self.kinematics.load()
			select(1,0)
			indexedimage(self.kinematics.true_losvd)
			for index in range(self.kinematics.true_losvd.shape[1]):
				i1 = index%4
				i2 = index/4
				xi = i2 / 4
				yi = i2 % 4
				print i1, xi, yi
				select(xi+2, yi)
				y = self.kinematics.true_losvd[:,index]
				graph(y/y.sum())
				y = losvd[:,index][self.kinematics.index_min:self.kinematics.index_max]
				graph(y/y.sum(), color="red", linestyle="dot")
				
		if self.df_samples:
			df_samples = self.df_samples.df_samples
		else:
			df_samples = None

		if df_samples:
			R = sqrt(df_samples.x**2 + df_samples.y**2)
			vlos = df_samples.vz
			#numpy.histogram2d(
			select(1,0)
			data = numpy.histogram2d(R, vlos, bins=[30, 30], range=[(0, 1.5), (-40, 40)])
			data = data[0].T
			indexedimage(data)
			#data = density2d(R, vlos, xmin=0, xmax=1.5, ymin=-40., ymax=40., bincountx=30, bincounty=30)
			#data = data.data2d
		
		select(1,1)
		if df_samples:
			contour(data, levels=10)
		contour(losvd, levels=10, color="red")
		
		print losvds.shape
		#losvd = 
		draw()

class PlotSchwSolution(object):
	def __init__(self, solution, schwmodel, storage_2d, binned_data_m2, binned_data_m4, galaxy=None, observation=None, storage_losvd=None):
		self.solution = solution
		self.schwmodel = schwmodel
		self.galaxy = galaxy
		self.storage_2d = storage_2d
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.observation = observation
		self.storage_losvd = storage_losvd
		
		
	def run(self, args, opts, scope):
		self.solution.load()
		#import pdb
		#pdb.set_trace()
		if 0:
			c = self.solution.orbitweights
			c = c.reshape((8, 20))
			c[4:,:]  = 0
			c = c.reshape(8*20)
			self.solution.orbitweights = c
		if self.observation:
			self.observation.load()
		print self.solution

		if self.storage_losvd:
			self.storage_losvd.load()
			losvds = self.storage_losvd.losvds
			self.losvd = dot(losvds.T, self.solution.orbitweights).T
			#print self.losvd.shape
			#import pdb
			#pdb.set_trace()
			#dsa
		
		self.storage_2d.aperture.load()
		#self.schwmodel.storage_2d.aperture.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		self.binned_data_m2.aperture.load()
		self.binned_data_m4.aperture.load()
		self.moments_observation_m2 = self.binned_data_m2.moments
		self.moments_observation_m4 = self.binned_data_m4.moments
		self.storage_3d_log = scope["storage_3d_log"]
		self.storage_3d_log.load()
		print self.storage_3d_log.moments3d
		
		self.moments_3d_log = self.solution.calculate_solution_moments(self.storage_3d_log.moments3d)
		self.light_profile = scope["light_profile"]

		
		self.markers = [10**-1.5, 0.05, 0.1, 0.25]
		self.markers = []
		self.markerslog10 = log10(self.markers)
		self.storage_2d_m0 = scope["storage_2d_m0"]
		self.storage_2d_m2 = scope["storage_2d_m2"]
		self.moments_m0 = scope["storage_2d_m0"]
		self.moments_m0.load()
		self.moments_m0.aperture.load()
		print dir(self.moments_m0.aperture)
		self.xborders_2d = self.moments_m0.aperture.aperture_rborders_kpc
		self.momentsprojected_m0 = self.solution.calculate_solution_moments(self.moments_m0.projectedmoments)
		self.moments_m0 = self.momentsprojected_m0[0]
		#self.momentsprojected_m0 = self.solution.calculate_solution_moments(self.moments_m0.projectedmoments)
		#self.masses2d = self.solution.calculate_solution_moments(self.storage_2d.masses)
		#print "!2d density sum:", sum(self.moments_m0), sum(self.solution.moments_solution_m0[0]), self.solution.moments_solution_m0.shape
		#assert self.solution.binned_data.aperture == self.schwmodel.storage_2d.aperture
		
		self.moments_m2 = moments_m2 = self.solution.moments_solution_m2
		self.moments_m4 = moments_m4 = self.solution.moments_solution_m4
		self.moments3d = moments3d = self.solution.moments3d_solution
		
		self.varvr = moments3d[1] 
		self.varvphi = moments3d[2]/2
		#self.varvtheta = moments3d[6]
		self.anisotropy  = 1 - (self.varvphi)/(self.varvr)
		self.x_2d_m2 = self.solution.binned_data_m2.aperture.aperture_rcenters_kpc
		self.x_2d_m4 = self.solution.binned_data_m4.aperture.aperture_rcenters_kpc
		self.x_3d = self.schwmodel.storage_3d.x
		self.xborders_3d = self.schwmodel.storage_3d.xborders
		
		
		#beta2 = 1 - varvtheta/varvr
		#beta3 = 1 - varvphi/varvr
		galaxy = None
		try:
			galaxy = scope["galaxy"]
		except:
			pass
		
		#document("25cm,15cm")
		#page(fontsize="12pt")
		#document(size="30cm,15cm")
		try:
			short = scope["short"]
		except:
			short = False
		try:
			vdisp = scope["vdisp"]
		except:
			vdisp = False
			
		try:
			df = scope["df"]
		except:
			df = False
		try:
			full = scope["full"]
		except:
			full = False
			

		if df:
			document("10cm,10cm")
			page(fontsize="15pt")
			mozaic(1,1,box)
			select(0,0)
			self.plot_orbitweights()
		elif vdisp:
			document("10cm,10cm")
			page(fontsize="15pt")
			mozaic(1,1,box)
			select(0,0)
			#self.plot_anisotropy(galaxy)
			#select(0, 1)
			#self.plot_vdisp_hr(galaxy)
			self.plot_pvdisp()
		elif short:
			document("10cm,15cm")
			page(fontsize="15pt")
			mozaic(1,2,box)
			select(0,0)
			self.plot_anisotropy(galaxy)
			select(0, 1)
			self.plot_vdisp(galaxy)
			
		elif full:
			document("30cm,20cm")
			page(fontsize="15pt")
			mozaic(3,3,box)
			select(1,0)
			self.plot_anisotropy(galaxy)
			select(1,1)
			self.plot_mass_3d()
			#select(1,2)
			#self.plot_mass_3d_log()
			select(0,2)
			#self.plot_orbitweights_real()
			select(0,1)
			self.plot_orbitweights()
			select(0, 0)
			self.plot_vdisp(galaxy)
			select(2, 0)
			self.plot_pvdisp()
			select(2, 1)
			self.plot_kurtosis()
			#self.plot_mass_2d()
			#select(0,2)
			select(2,2)
			self.plot_mass_3d_log()
			select(1,2)
			self.plot_logslope_3d()
			select(0,2)
			#self.plot_jeanscheck(galaxy)
			
			self.schwmodel.storage_2d.load()
			if 0:
				moments = self.schwmodel.storage_2d.projectedmoments
				n = moments.shape[0]
				grid = zeros((n, moments.shape[2]))
				for i in range(n):
					grid[i] = moments[i,2,:] * self.solution.orbitweights[i]
				grid = log10(grid)
				grid -= grid.max()
				grid[grid<-4] = -4
				indexedimage(grid.T, colormap="whiteblack")
			
				select(1,2)
				
				self.schwmodel.storage_2d.load()
				self.schwmodel.storage_2d.storage_3d.load()
				storage_3d = self.schwmodel.storage_2d.storage_3d
				print storage_3d.moments3d.shape, self.solution.orbitweights.shape
				moments3d= tensordot(self.solution.orbitweights, storage_3d.moments3d, axes=[(0,), (0,)])
				print moments3d.shape
				mask = moments3d[0]>0
				moments3d[1:,mask] /= moments3d[0,mask]
				x = storage_3d.x
				graph(x, moments3d[storage_3d.v20index]**0.5, color="red")
				graph(x, (moments3d[storage_3d.v02index]/2)**0.5, color="green")
				
				select(2,2)
				rho = moments3d[storage_3d.v00index]/(10**x)**3
				mask = rho > 0
				graph(x[mask], log10(rho[mask]))
				
				select(3,2)
				N = len(rho)
				print "N=",N
				print rho, x
				M = 5
				rho = sum(rho.reshape(N/M,M), axis=1)/M
				x = sum(x.reshape(N/M,M), axis=1)/M
				print rho, x
				drho = log10(rho[1:]) - log10(rho[:-1])
				dx   = x[1:] - x[:-1]
				mrho = (rho[1:] + rho[:-1])/2
				mx   = (x[1:] + x[:-1])/2
				mask = rho > 0
				slope = drho/dx
				#mask
				#mask = (mx > 0) 
				mask = ~(isnan(slope) | isinf(slope))
				graph(mx[mask], slope[mask])
				
				print slope[mask]
				ylim(-5.5,1.5)
		else:
			document("20cm,15cm")
			page(fontsize="15pt")
			mozaic(2,2,box)
			select(1,0)
			self.plot_anisotropy(galaxy)
			select(0,0)
			self.plot_mass_3d()
			select(0,1)
			self.plot_orbitweights()
			select(1, 1)
			self.plot_vdisp(galaxy)
			if 0:
				select(2, 0)
				self.plot_pvdisp()
				select(2, 1)
				self.plot_kurtosis()
				select(0,2)
				
				self.schwmodel.storage_2d.load()
				moments = self.schwmodel.storage_2d.projectedmoments
				n = moments.shape[0]
				grid = zeros((n, moments.shape[2]))
				for i in range(n):
					grid[i] = moments[i,2,:] * self.solution.orbitweights[i]
				grid = log10(grid)
				grid -= grid.max()
				grid[grid<-4] = -4
				indexedimage(grid.T, colormap="whiteblack")
				
				select(1,2)
				
				self.schwmodel.storage_2d.load()
				self.schwmodel.storage_2d.storage_3d.load()
				storage_3d = self.schwmodel.storage_2d.storage_3d
				print storage_3d.moments3d.shape, self.solution.orbitweights.shape
				moments3d= tensordot(self.solution.orbitweights, storage_3d.moments3d, axes=[(0,), (0,)])
				print moments3d.shape
				mask = moments3d[0]>0
				moments3d[1:,mask] /= moments3d[0,mask]
				x = storage_3d.x
				graph(x, moments3d[storage_3d.v20index]**0.5, color="red")
				graph(x, (moments3d[storage_3d.v02index]/2)**0.5, color="green")
				
				select(2,2)
				rho = moments3d[storage_3d.v00index]/(10**x)**3
				mask = rho > 0
				graph(x[mask], log10(rho[mask]))
				
				select(3,2)
				N = len(rho)
				print "N=",N
				print rho, x
				M = 5
				rho = sum(rho.reshape(N/M,M), axis=1)/M
				x = sum(x.reshape(N/M,M), axis=1)/M
				print rho, x
				drho = log10(rho[1:]) - log10(rho[:-1])
				dx   = x[1:] - x[:-1]
				mrho = (rho[1:] + rho[:-1])/2
				mx   = (x[1:] + x[:-1])/2
				mask = rho > 0
				slope = drho/dx
				#mask
				#mask = (mx > 0) 
				mask = ~(isnan(slope) | isinf(slope))
				graph(mx[mask], slope[mask])
				
				print slope[mask]
				ylim(-5.5,1.5)
		
		draw()
		
	def plot_jeanscheck(self, galaxy):
		logrs = self.storage_3d_log.x
		rs = 10**logrs
		logr_borders = self.storage_3d_log.xborders
		r_borders = 10**logr_borders

		#return
		print type(galaxy.profile_model)
		print dir(galaxy.profile_model.dm_profile)
		#Gm_over_r = array([galaxy.profile_model.dm_profile.G * galaxy.profile_model.enclosed_mass(0, r)/r for r in rs])
		
		enclosed_mass_stars = cumsum(self.moments_3d_log[0]) * 1e6
		#print sum(m0)
		#print m0
		#dsa
		#print cumsum(m0).shape
		#Gm_over_r = array([(galaxy.profile_model.dm_profile.G * (galaxy.profile_model.dm_profile.enclosed_mass(rs[i]) + enclosed_mass_stars[i]))/rs[i] for i in range(len(rs))])
		#Gm_over_r = array([((galaxy.profile_model.dm_profile.enclosed_mass(rs[i]) + enclosed_mass_stars[i])) for i in range(len(rs))])
		Gm_over_r = array([((galaxy.profile_model.dm_profile.enclosed_mass(rs[i]) + enclosed_mass_stars[i]*0)) for i in range(len(rs))])
		print Gm_over_r.shape
		print enclosed_mass_stars[-1]
		#dsa
		#sdsa
		graph(logrs, log10(Gm_over_r), linestyle="dot")
		
		
		sigmarsq = self.moments_3d_log[1]/self.moments_3d_log[0]
		sigmatsq = self.moments_3d_log[2]/self.moments_3d_log[0]
		beta = 1 - sigmatsq/2/sigmarsq
		beta_borders= (beta[:-1] + beta[1:]) / 2
		Gm_over_r_borders = (Gm_over_r[:-1] + Gm_over_r[1:])/2
		sigmarsq_borders = (sigmarsq[:-1] + sigmarsq[1:])/2
		#alpha = dlog sigma_r^2/dlog r
		dlogsigmarsq = log10(sigmarsq[1:]) - log10(sigmarsq[:-1])
		dlogr = logrs[1:] - logrs[:-1]
		alpha = dlogsigmarsq/dlogr
		#gamma = -galaxy.profile_model.light_profile.logslope(r_borders[1:-1])
		#print self.moments_3d_log.shape
		#dsa
		dlogr = logrs[1] - logrs[0]
		#dr = dlogr * rs
		rho = self.moments_3d_log[0]/(4*pi*rs**3*dlogr)
		#dlogmass = log(self.moments_3d_log[0,1:]) - log(self.moments_3d_log[0,:-1])
		dlogrho = log10(rho[1:]) - log10(rho[:-1])
		gamma = -(dlogrho)/dlogr
		mask = ~(isnan(gamma) | isinf(gamma))
		#graph(logr_borders[1:-1][mask], gamma[mask], color="red", linestyle="dash")
		print gamma
		#dsa
		#draw()
		
		print sigmarsq_borders.shape, gamma.shape, 2*beta_borders.shape, alpha.shape
		y = log10(sigmarsq_borders*(gamma-2*beta_borders-alpha) * rs[:-1] / galaxy.profile_model.dm_profile.G)
		mask = ~isnan(y)
		print y
		graph(logr_borders[1:-1][mask], y[mask], color="red", linestyle="dash")
		
		mask = ~isnan(beta)
		#print beta
		#graph(logrs[mask], beta[mask])
		if 0:
			y = log10(self.moments_3d_log[0])
			y[y<-10] = -10
			dy = y[1:] - y[:-1]
			logr = self.storage_3d_log.x
			dlogr = logr[1:] - logr[:-1]
			y = dy/dlogr-3
			#mask = ~isinf(y)
			graph(self.storage_3d_log.xborders[1:-1], y)
			
			
			print log10(Gm_over_r)
		for x in self.markerslog10:
			vline(x, color="blue")
		#dsa
		labels("log r/kpc", "GM/r")
		
	def plot_anisotropy(self, galaxy):
		vfill(0, 0.1, color="lightgrey")
		#vline(0.05, linestyle="dash")
		graph(self.x_3d, self.anisotropy, color="red", linestyle="dash")
		if self.galaxy and hasattr(self.galaxy, "beta"): 
			hline(galaxy.beta, color="black")
		ylim(-2, 1)
		labels("3d radius (kpc)", "anisotropy") 
		if galaxy:
			import mab
			if isinstance(galaxy, mab.gd.schw.galaxy2.Galaxy_constant_anisotropy):
				hline(galaxy.beta)
		#else:
		if 0: #self.observation:
			#import pdb; pdb.set_trace()
			
			xs = []
			betas = []
			varvrs = []
			varvphis = []
			varvthetas = []
			Nperbin = 100
			samples = self.observation.stars
			for n, r, vr, vphi, vtheta in binrange(Nperbin, samples.r3d, samples.vr, samples.vphi, samples.vtheta):
				#print n, mean(vr), std(vr), std(vphi)
				varvphi = var(vphi)
				varvtheta = var(vtheta)
				varvr = var(vr)
				beta = 1 - (varvphi+varvtheta)/(2*varvr)
				#print mean(r), beta
				xs.append(log10(mean(r)))
				betas.append(beta)
				varvrs.append(varvr)
				varvphis.append(varvphi)
				varvthetas.append(varvtheta)
			
			xs = array(xs)
			graph(10**xs, betas, color="green")
			
			
		autolegend("output", "true value")
		clearautolegend()
		for x in self.markers:
			vline(x, color="blue")
		
		
	def plot_logslope_3d(self):
		print self.storage_3d_log.x
		#box()
		y = log10(self.moments_3d_log[0])
		y[y<-10] = -10
		dy = y[1:] - y[:-1]
		logr = self.storage_3d_log.x
		dlogr = logr[1:] - logr[:-1]
		y = dy/dlogr-3
		#mask = ~isinf(y)
		graph(self.storage_3d_log.xborders[1:-1], y)
		ylim(-5, 2)
		
		y = self.light_profile.logslope(10**logr)
		graph(logr, y, color="red", linestyle="dash")
		
		for x in self.markerslog10:
			vline(x, color="blue")
		#borders = 10**self.storage_3d_log.xborders
		#mass = [self.light_profile.cumdensityr(r1, r2, M=1.) for r1, r2 in zip(borders[:-1], borders[1:])]
		#y = log10(mass)
		#y[y < -10] = -10
		#histogramline(self.storage_3d_log.xborders, y, color="red", linestyle="dash")
		#kaplot.grid(xinterval=1, yinterval=1)
		labels("log r/kpc", "dlog &rho; / dlog r")
		
	def plot_mass_3d_log(self):
		print self.storage_3d_log.x
		#box()
		y = log10(self.moments_3d_log[0])
		y[y < -10] = -10
		
		borders = 10**self.storage_3d_log.xborders
		mass = [self.light_profile.cumdensityr(r1, r2, M=1.) for r1, r2 in zip(borders[:-1], borders[1:])]
		if 1:
			histogramline(self.storage_3d_log.xborders, y)
			y = log10(mass)
			y[y < -10] = -10
			histogramline(self.storage_3d_log.xborders, y, color="red", linestyle="dash")
			labels("log r/kpc", "d M(r)/dlogr")
			kaplot.grid(xinterval=1, yinterval=1)
		else:
			histogramline(self.storage_3d_log.xborders, y-log10(mass), color="red")
			ylim(-0.1, 0.1)
		for x in self.markerslog10:
			vline(x, color="blue")
		
		
	def plot_mass_3d(self):
		kpc_to_arcsec = 1/self.schwmodel.light_model.arcsec_to_kpc(1.)
		#print kpc_to_arcsec 
		#histogramline(schwmodel.storage_3d.xborders - log10(kpc_to_arcsec), moments3d[0], binned=True, fill=False, drawverticals=False)
		#print moments3d[0]
		y = self.moments3d[0]
		normalize = True
		if normalize:
			y = y /sum(y)
		plotlog = False
		print y
		#dsa

		if plotlog:
			histogramline(self.xborders_3d, log10(y), binned=True, fill=False, drawverticals=False, color="red", linestyle="dash")
		else:
			histogramline(self.xborders_3d, y, binned=True, fill=False, drawverticals=False, color="red", linestyle="dash")
		print "3d density sum:", sum(self.moments3d[0])
		
		r1 = self.xborders_3d[:-1]
		r2 = self.xborders_3d[1:]
		rho3d_true = []
		for r1, r2 in zip(r1, r2):
			#r1 = schwmodel.galaxy.arcsec_to_kpc(r1)
			#r2 = schwmodel.galaxy.arcsec_to_kpc(r2)
			#print r1, r2
			rho3d_true.append(self.schwmodel.light_model.light_profile.cumdensityr(r1, r2, M=1.))
		rho3d_true = array(rho3d_true)
		if normalize:
			rho3d_true = rho3d_true/sum(rho3d_true)
		if plotlog:
			histogramline(self.schwmodel.storage_3d.xborders, log10(rho3d_true), binned=True, fill=False, drawverticals=False, alpha=0.7)
		else:
			histogramline(self.schwmodel.storage_3d.xborders, (rho3d_true), binned=True, fill=False, drawverticals=False, alpha=0.7)
			
		autolegend("output", "true value")
		
		labels("3d radius (kpc)", "relative mass (M<sub>&solar;</sub>)")
		clearautolegend()
		for x in self.markers:
			vline(x, color="blue")
		
	def plot_mass_2d(self):
		kpc_to_arcsec = 1/self.schwmodel.light_model.arcsec_to_kpc(1.)
		#print kpc_to_arcsec 
		#histogramline(schwmodel.storage_3d.xborders - log10(kpc_to_arcsec), moments3d[0], binned=True, fill=False, drawverticals=False)
		#print moments3d[0]
		y = self.moments_m0#@[0]
		print "2d density sum:", sum(self.moments_m0)
		normalize = True
		if normalize:
			y = y /sum(y)
		
		histogramline(self.xborders_3d, y, binned=True, fill=False, drawverticals=False, color="red", linestyle="dash")
		
		R1 = self.xborders_2d[:-1]
		R2 = self.xborders_2d[1:]
		rho3d_true = []
		for R1, R2 in zip(R1, R2):
			#r1 = schwmodel.galaxy.arcsec_to_kpc(r1)
			#r2 = schwmodel.galaxy.arcsec_to_kpc(r2)
			#print r1, r2
			rho3d_true.append(self.schwmodel.light_model.light_profile.cumdensityR(R1, R2, M=1.))
		rho3d_true = array(rho3d_true)
		if normalize:
			rho3d_true = rho3d_true/sum(rho3d_true)
		histogramline(self.schwmodel.storage_3d.xborders, rho3d_true, binned=True, fill=False, drawverticals=False, alpha=0.7)
		autolegend("output", "true value")
		
		
		#histogramline(self.solution.Rborders, self.masses2d, binned=True, fill=False, drawverticals=False, alpha=0.7)
		
		
		labels("2d radius (kpc)", "relative mass (M<sub>&solar;</sub>)")
		clearautolegend()		
		
	def plot_orbitweights_real(self):
		stars = self.observation.stars
		e = log10(-stars.E)
		l = stars.L/stars.Lmax * 8
		#print e, 
		#scatter(e, l)
		histogram(l, binwidth=1./5)
		
	def plot_orbitweights(self):
		dfgrid = self.schwmodel.dfgrid
		n = 200
		x = (arange(n)+0.5)/(n)
		logrs = dfgrid.logrmin + x * (dfgrid.logrmax - dfgrid.logrmin)
		ls = 1 * x
		c  = self.solution.orbitweights
		df = array([[dfgrid(c, logr, l) for logr in logrs] for l in ls])
	
		resize = ((0, 0), (dfgrid.n_I1, dfgrid.n_I2))
		#resize = ((0, 0), (dfgrid.n_I1, dfgrid.n_I2))
		if 0:
			df = log10(df)
			df -= df.max()
			df[df<-6] = -6
		im = indexedimage((df), colormap="whiterainbow", resize=resize)
		#im = indexedimage((df), colormap="whiteblack", resize=resize)
		c = c.reshape(8, 20)
		print "DF SHAPE", c.shape
		y = sum(c, axis=0)
		y = y / max(y) * c.shape[0]/2
		x = arange(0, c.shape[1]) + 0.5
		graph(x, y)
		y = sum(c, axis=1)
		y = y / max(y) * c.shape[1]/2
		x = arange(0, c.shape[0]) + 0.5
		y = c.shape[1]-y
		graph(y, x)
		labels("energy index", "angular mom. index")
		levels = [-1,1e6]#0, im.datamax/2, im.datamax]
		innercolorbar("weight", im, location="left, bottom", edgespacing="6mm", size="0.5cm,3cm", colormap=im.colormap, levels=levels)
		print im.datamin, im.datamax
		clearautolegend()
		
		
		
	def plot_vdisp_hr(self, galaxy):
		#sigma = (self.momentsprojected_m0[2]/self.momentsprojected_m0[0])**0.5
		self.storage_2d_m2.load()
		self.momentsprojected_m0 = self.solution.calculate_solution_moments(self.storage_2d_m2.projectedmoments)
		#print self.storage_2d_m2
		#print self.storage_2d_m0
		Rmax = self.storage_2d_m2.aperture.rmax
		N = len(sigma)
		dR = Rmax*1./N
		R = arange(0, N)/(1.*N) * Rmax + dR/2
		#graph(
		#import pdb;
		#pdb.set_trace()
		graph(R, sigma)
		xlim(0, 1.5)
		ylim(0, 15)
		
		#print self.moments_m0.shape
		#dsa
	def plot_vdisp(self, galaxy):
		print "vr", sqrt(self.varvr)
		graph(self.x_3d, sqrt(self.varvr), color="red", linestyle="dash")
		#graph(schwmodel.storage_3d.x, sqrt(varvphi), color="green")
		#graph(schwmodel.storage_3d.x, sqrt(varvtheta), color="blue")
		#graph(schwmodel.storage_3d.x, sqrt(), color="blue")
		#graph(self.x_3d, sqrt((self.varvphi+self.varvtheta)/2), color="green", linestyle="dash")
		graph(self.x_3d, sqrt(self.varvphi), color="green", linestyle="dash")
		

		if galaxy:
			import mab
			if isinstance(galaxy, mab.gd.schw.galaxy2.Galaxy_constant_anisotropy):
				jeans = galaxy.jeans()
				sigmar = jeans.sigmar(self.x_3d)
				sigmat = (1-galaxy.beta)**0.5*sigmar
				graph(self.x_3d, sigmar)
				graph(self.x_3d, sigmat)
		
		if self.galaxy:
			jeans = galaxy.jeans()
			r3d = schwmodel.storage_3d.x
			sigmar = jeans.sigmar(r3d)
			graph(schwmodel.storage_3d.x, sigmar, color="black")
			sigmat = sqrt((1-galaxy.beta))*sigmar
			graph(schwmodel.storage_3d.x, sigmat, color="black")
			#graph(schwmodel.storage_3d.x, sigmar, color="blue", linestyle="dash")
			#r = arange(
			
			
		autolegend("output &sigma;<sub>r</sub>", "output &sigma;<sub>t</sub>", "true value", location="left,bottom")
		labels("3d radius (kpc)", "&sigma; - km/s")
		ylim(0,15)
		
		#draw()
		
	def plot_pvdisp(self):
		if self.storage_losvd:
			#losvd = dot(losvds.T, self.solution.orbitweights).T
			dv = self.storage_losvd.delta_v
			p = self.losvd# / sum(self.losvd * dv)
			vs = self.storage_losvd.vcenters
			vm0 = sum(      p.T * dv, axis=1)
			vm  = sum(vs   *p.T * dv, axis=1)/vm0
			vm2 = sum(vs**2*p.T * dv, axis=1)/vm0
			print vm
			print vs
			print "sigmas", vm2**0.5
			graph(self.storage_losvd.Rcenters, vm2**0.5, color="blue", linestyle="dash")
			#print sum(vs * p * dv, axis=1)
			#print sum(vs * p * dv, axis=0)
			#print self.losvd
			#dsa
			if 0:
				for i in range(p.shape[1]):
					graph(p[:,i])
					print (sum(vs**2*p[:,i]) / sum(p[:,i]))**0.5
			#return
		
		scatter(self.x_2d_m2, self.moments_observation_m2[2]**0.5, symbolsize="28pt")
		graph(self.x_2d_m2, sqrt(self.moments_m2[2]), color="red", linestyle="dash")
		print "m2", self.x_2d_m2, len(self.x_2d_m2)
		print "moments", sqrt(self.moments_m2[2]), len(sqrt(self.moments_m2[2]))
		
		ylim(0, 15)
		xlim(0, 1.0)
		labels("R (kpc)", "&sigma;<sub>los</sub> (km/s)")
		#xlim(0, 1.5)
	
	def plot_kurtosis(self):
		if 0:
			scatter(self.x_2d_m4, self.moments_observation_m4[4]**(1./4), symbolName="squaresolid")
			graph(self.x_2d_m4, self.moments_m4[4]**(1./4), color="red", linestyle="dash")
			if self.galaxy and hasattr(self.galaxy, "jeans"):
				filename = os.path.join(self.dirname, "sigma_los.npy")
				n = 100
				R = arange(n)/(n-1.) * self.Rmax
				sigma_los = load(filename) 
				filename = os.path.join(self.dirname, "m4_los.npy")
				m4_los = load(filename) 
				graph(R, m4_los**(1./4), color="blue")
			ylim(0,20)
			xlim(0, 1.5)
		else:
			if self.storage_losvd:
				#losvd = dot(losvds.T, self.solution.orbitweights).T
				dv = self.storage_losvd.delta_v
				p = self.losvd# / sum(self.losvd * dv)
				vs = self.storage_losvd.vcenters
				vm0 = sum(      p.T * dv, axis=1)
				vm2  = sum(vs**2   *p.T * dv, axis=1)/vm0
				vm4 = sum(vs**4*p.T * dv, axis=1)/vm0
				graph(self.storage_losvd.Rcenters, vm4/vm2**2, color="blue", linestyle="dash")
			scatter(self.x_2d_m4, self.moments_observation_m4[4]/self.moments_observation_m4[2]**2, symbolName="squaresolid")
			graph(self.x_2d_m4, self.moments_m4[4]/self.moments_m4[2]**2, color="red", linestyle="dash")
			#print self.x_2d, len(self.x_2d)
			#print sqrt(self.moments[2]), len(sqrt(self.moments[2]))
			#ylim(0, 15)
			#xlim(0, 1.5)
			if self.galaxy and hasattr(self.galaxy, "jeans"):
				filename = os.path.join(self.dirname, "sigma_los.npy")
				n = 100
				R = arange(n)/(n-1.) * self.Rmax
				sigma_los = load(filename) 
				filename = os.path.join(self.dirname, "m4_los.npy")
				m4_los = load(filename) 
				graph(R, m4_los/sigma_los**4, color="blue")
			
			ylim(1,5)
			xlim(0, 1.5)
class PlotSchwSolution2(PlotSchwSolution):
	def __init__(self, solution, schwmodel, galaxy=None):
		self.solution = solution
		self.schwmodel = schwmodel
		self.galaxy = galaxy
		
		
	def run(self, args, opts, scope):
		self.solution.load()
		document("20cm,15cm")
		page(fontsize="15pt")
		#document(size="30cm,15cm")
		mozaic(2,2,box)
		self.plot()
		draw()
		
	def plot(self):
		pass
		
