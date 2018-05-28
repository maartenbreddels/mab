# -*- coding: utf-8 -*-
import cvxopt.solvers
import os
import numpy
from numpy import *
import mab.gd.logging as logging
import mab.utils.numpy
import sys
import scipy.optimize
from mab.gd import gdfast_schw
import numpy.linalg
import pyublas

logger = logging.getLogger("gd.schw.solution2")
class OrbitWeightCalculatorSeperable(object):
	def __init__(self, modelpath, dfgrid, galaxy, df, name="orbitweights", postfix="", plot=True):
		self.modelpath = modelpath
		self.dfgrid = dfgrid
		self.plot = plot
		self.galaxy = galaxy
		self.name = name
		self.df = df
		self.postfix = postfix
		
	def run(self, args, opts, scope):
		self.calculate(scope)
		self.save()
		
	def calculate(self, scope):
		if 0:
			filename = os.path.join(self.modelpath, "df", "fE" +self.postfix +".npy")
			logger.info("loading f1(E) as: %s" % filename)
			dff = load(filename)
			
			filename = os.path.join(self.modelpath, "df", "fE-E.npy")
			logger.info("loading corresponding energy as: %s" % filename)
			dfEs = load(filename)
		else:
			
			self.df.load()
			dff = self.df.fE
			dfEs = -self.df.es
			#dff_interpol = scipy.interpolate.interp1d(log10(-dfEs[::-1]), log10(dff[::-1]))
		#orbitweights = mab.utils.numpy.mmapzeros((self.dfgrid.n_I1, self.dfgrid.n_I2))
		orbitweights_inner_products = mab.utils.numpy.mmapzeros((self.dfgrid.dof))
		df_inner_products = mab.utils.numpy.mmapzeros((self.dfgrid.dof))
		@mab.parallelize.parallelize(cores=scope["cores"], info=scope["progressbar"])
		def orbitweight(I1, I2):
			dither1 = self.dfgrid.dither
			dither2 = self.dfgrid.dither
			for i in range(dither1):
				E = self.dfgrid.subgrid.Es[I1*dither1 + i]
				dE = (self.dfgrid.subgrid.E_borders[I1*dither1 + i + 1] - self.dfgrid.subgrid.E_borders[I1*dither1 + i])
				Lmax = self.dfgrid.profile_model.Lmax_at_E(E)
				#for j in range(len(dfgrid.subgrid.ls)):
				for j in range(dither2):
					L = (self.dfgrid.subgrid.ls[I2*dither2 + j])*Lmax
					dL = (self.dfgrid.subgrid.l_borders[I2*dither2 + j + 1] - self.dfgrid.subgrid.l_borders[I2*dither2 + j])*Lmax
					
					#if opts.true:
					#	g = galaxy.gEL(E, L)
					#	mass = g * galaxy.fEL(E, L) * dE * dL
					#if opts.truenumerical:
					if 1:
						Ei = argmin(abs(dfEs-E))
						fEL = dff[Ei] * self.galaxy.fL(L) #* L ** (-2*self.anisotropy_beta)
						#fEL = 10**dff_interpol(log10(-E)) * self.galaxy.fL(L) #* L ** (-2*self.anisotropy_beta)
						#print E, dfEs[Ei]
						g = self.dfgrid.profile_model.gEL(E, L)
						mass = g * fEL * dE * dL
						#countG = g * dE *  dL
					u = float(i + 0.5) / dither1
					v = float(j + 0.5) / dither2
					#print u, v
					for k in range(self.dfgrid.dof_per_cell):
						orbitweights_inner_products[self.dfgrid.dof_index(I1, I2, k)] += mass * self.dfgrid.basis(k, u, v)
						df_inner_products[self.dfgrid.dof_index(I1, I2, k)] += fEL * dE * dL * self.dfgrid.basis(k, u, v)
				
		nE = self.dfgrid.n_I1
		nL = self.dfgrid.n_I2
		I1s = [i for i in range(nE) for j in range(nL)] #[:1]
		I2s = [j for i in range(nE) for j in range(nL)] #[:1]
		#print len(I1s), len(dfgrid.subgrid.E_borders), len(dfgrid.subgrid.Es)
		orbitweight(I1s, I2s)
		#print orbitweights_inner_products
		self.orbitweights = self.dfgrid.solve_coordinates(orbitweights_inner_products)
		self.df = self.dfgrid.solve_coordinates(df_inner_products)
		#print self.orbitweights
		if 0:
			for I1 in I1s:
				for I2 in I2s:
					orbitweight(I1, I2)
					
	def load(self):
		filename_base = os.path.join(self.modelpath, "df", self.name+self.postfix) 
		filename_npy = filename_base + ".npy"
		self.orbitweights = numpy.load(filename_npy)
		#filename = os.path.join(self.modelpath, "df", "df.npy")
		#self.df = numpy.load(filename)
		return self.orbitweights 
		
	def save(self):
		filename = os.path.join(self.modelpath, "df", "df" +self.postfix+".npy")
		logger.info("storing distribution function in file: %s" % filename)
		numpy.save(filename, self.df)
		filename_base = os.path.join(self.modelpath, "df", self.name+self.postfix)
		filename_png = filename_base + ".png" 
		filename_npy = filename_base + ".npy"
		logger.info("storing orbitweights in file: %s" % filename_npy) 
		numpy.save(filename_npy, self.orbitweights)
		import kaplot
		kaplot.box()
		n = 200
		x = (arange(n)+0.5)/(n)
		us = self.dfgrid.umin + x * (self.dfgrid.umax - self.dfgrid.umin)
		ls = 1 * x
		#logrs = self.dfgrid.subgrid.logrs
		#logrs = arange(len(logrs)) / (len(logrs)-1.)
		#ls = self.dfgrid.subgrid.ls
		#print ls, logrs
		#print self.orbitweights
		df = array([[self.dfgrid(self.orbitweights, u, l) for u in us] for l in ls])
		print ls
		#print df
		#print df.shape
		kaplot.indexedimage((df), colormap="whiterainbow")
		print self.orbitweights.shape
		#kaplot.indexedimage((self.orbitweights.reshape(10, 30)))
		kaplot.hardcopy(filename_png)

class RegularizationTester(object):
	def __init__(self, solution):
		self.solution = solution
		
	def run(self, args, opts, scope):
		self.solution.init()
		from kaplot import *
		box()
		regs = arange(-3, 3, 0.1)
		chisqs = []
		for l in regs:
			self.solution.regularization.regularization_delta = 10**l
			self.solution.solve()
			chisqs.append(self.solution.chisqs["m2"]/self.solution.dof)
		ylim(0, 3)
		graph(regs, chisqs)
		hline(1)
		draw()
			
class FakeMoments(object):
	def __init__(self, solution, binned_data_m2, binned_data_m4):
		self.solution = solution
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		
	def run(self, args, opts, scope):
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		self.solution.load()
		
		
		self.binned_data_m2.moments[2] = self.solution.moments_solution_m2[2]
		self.binned_data_m2.moments[4] = self.solution.moments_solution_m2[4]
		
		self.binned_data_m4.moments[2] = self.solution.moments_solution_m4[2]
		self.binned_data_m4.moments[4] = self.solution.moments_solution_m4[4]
		self.binned_data_m2.save()
		self.binned_data_m4.save()
		
class SolutionMoments(object):
	
	def calculate_solution_moments(self, moments):
		return tensordot(self.orbitweights, moments, axes=[(0,), (0,)])
		
	def load(self):
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		self.orbitweights = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		self.moments_solution_m2 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		self.moments_solution_m4 = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		self.moments3d_solution = numpy.load(filename)
		logger.debug("loading %s" % filename)
		
	def solve(self):
		x = self.orbitweights = orbitweights = self.findsolution()
		
		moments_m0 = tensordot(orbitweights, self.storage_2d_m0.projectedmoments, axes=[(0,), (0,)])
		moments_m2 = tensordot(orbitweights, self.storage_2d_m2.projectedmoments, axes=[(0,), (0,)])
		moments_m4 = tensordot(orbitweights, self.storage_2d_m4.projectedmoments, axes=[(0,), (0,)])
		moments3d = tensordot(orbitweights, self.storage_3d.moments3d, axes=[(0,), (0,)])
		if 0:
			moments_m0[1:] /= moments_m0[0]
			moments_m2[1:] /= moments_m2[0]
			moments_m4[1:] /= moments_m4[0]
		else:
			moments_m0[1:] /= moments_m0[0]
			moments_m2[1:] /= self.binned_data_m2.moments[0]
			moments_m4[1:] /= self.binned_data_m4.moments[0]
		moments3d[1:] /= moments3d[0]
		
		import pdb
		#pdb.set_trace()
		
		#self.moments2d_solution = moments_m2
		self.moments2d_solution_m2 = moments_m2
		self.moments2d_solution_m4 = moments_m4
		self.moments3d_solution = moments3d
		
		dirname = filename = os.path.join(self.dirname, "results")
		if not os.path.exists(dirname):
			os.makedirs(dirname)
		filename = os.path.join(self.dirname, "results/orbitweights" +self.postfix +".npy")
		numpy.save(filename, orbitweights)
		logger.debug("saving %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m0.npy")
		numpy.save(filename, moments_m0)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m2.npy")
		numpy.save(filename, moments_m2)
		logger.debug("saving %s" % filename)
		filename = os.path.join(self.dirname, "results/solution_projectedmoments" +self.postfix +"_m4.npy")
		numpy.save(filename, moments_m4)
		logger.debug("saving %s" % filename)
		
		filename = os.path.join(self.dirname, "results/solution_moments3d" +self.postfix +".npy")
		numpy.save(filename, moments3d)
		logger.debug("saving %s" % filename)
		
		
		self.logger.info("sum of orbitweight %f" % sum(orbitweights))
		self.logger.info("total light (2d) m0 %f" % sum(moments_m0[0]))
		self.logger.info("total light (2d) m2 %f" % sum(moments_m2[0]))
		self.logger.info("total light (2d) m4 %f" % sum(moments_m4[0]))
		
		
		f = file(os.path.join(self.dirname, "results/statistics_qp" +self.postfix +".txt"), "w")
		#print "lenghts", len(kinematics.m2), len(varvlos)
		dof_m2 = len(self.binned_data_m2.moments[2])
		dof_m4 = len(self.binned_data_m4.moments[2])
		dof = self.dof = dof_m2 + dof_m4
		chisq_m2_check = sum( ((self.binned_data_m2.moments[2]-moments_m2[2])/self.binned_data_m2.e_moments[2])**2 )
		chisq_m4_check = sum( ((self.binned_data_m4.moments[4]-moments_m4[4])/self.binned_data_m4.e_moments[4])**2 )
		#print self.binned_data_m2.moments[2]**0.5, moments_m2[2]**0.5
		self.logger.info("(check) chisq 2nd moment: %f (reduced: %f, dof: %d)" % (chisq_m2_check, chisq_m2_check/dof_m2, dof_m2))
		self.logger.info("(check) chisq 4th moment: %f (reduced: %f, dof: %d)" % (chisq_m4_check, chisq_m4_check/dof_m4, dof_m4))
		#print "chisq2m2", ((self.binned_data_m2.moments[2]-moments_m2[2])/self.binned_data_m2.e_moments[2])**2
		#print "chisq2m4", ((self.binned_data_m4.moments[4]-moments_m4[4])/self.binned_data_m4.e_moments[4])**2
		#chisq = ((self.binned_data.moments[2]-moments[2])/self.binned_data.e_moments[2])**2
		#print chisq
		#print argsort(chisq)
		#print chisq[argsort(chisq)]
		#print chis**2
		#print kinematic_m2
		#print moments[2]
		#print e_kinematic_m2
		chisq_tot, chisqs = self.calc_chisq(x, print_info=True, returnall=True)
		self.chisqs = dict(chisqs)
		print >>f, "chisq, redchsiq, dof"
		#chisq = sum(chis**2)
		print >>f, chisq_tot, ",", chisq_tot/dof,",", dof
				
		#print "chisq:",chisq_tot, chisq
		
		chisq = chisq_m2_check + chisq_m4_check
		f = file(os.path.join(self.dirname, "results/probability" +self.postfix +".txt"), "w")
		print >>f, "%e" % exp(-0.5*chisq)
		f = file(os.path.join(self.dirname, "results/logprobability" +self.postfix +".txt"), "w")
		print >>f, "%e" % (-0.5*chisq)
		
		logps = [(name, -0.5*chisq) for name, chisq in chisqs] 
		f = file(os.path.join(self.dirname, "results/logprobability_seperate" +self.postfix +".txt"), "w")
		print >>f, repr(logps)
		
		return orbitweights
		
	def calc_chisq(self, x, print_info=False, returnall=False):
		dof_m2 = len(self.binned_data_m2.moments[2])
		dof_m4 = len(self.binned_data_m4.moments[2])
		dof_reg = 1. # TODO, what is a reasonable value
		chisq_tot = 0
		chisqs = []
		
		chisq_m2 = dot(dot(x, self.P_m2),x)
		chisq_tot += chisq_m2
		chisqs.append(("m2", chisq_m2))
		if print_info:
			self.logger.info("chisq 2nd moment: %f (reduced: %f, dof: %d)" % (chisq_m2, chisq_m2/dof_m2, dof_m2))
		 
		chisq_m4 = dot(dot(x, self.P_m4),x)
		chisq_tot += chisq_m4
		chisqs.append(("m4", chisq_m4))
		if print_info: 
			self.logger.info("chisq 4th moment: %f (reduced: %f, dof: %d)" % (chisq_m4, chisq_m4/dof_m4, dof_m4))
		if self.regularization: 
			chisq_reg = dot(dot(x, self.P_reg),x) 
			chisqs.append(("reg", chisq_reg))
			chisq_tot += chisq_reg
			if print_info: 
				self.logger.info("chisq regular.  : %f (reduced: %f)" % (chisq_reg, chisq_reg/dof_reg))
			#self.logger.info("chisq 2nd moment: %f (reduced: %f)" % (chisq_m22, chisq_m22/dof))
		else:
			chisq_reg = 0
			
		
		if self.fitdensity3d:
			dof_light = len(self.masses3d_target)
			chisq_light = sum(((dot(x, self.schw_masses3d) - self.masses3d_target)/self.e_masses3d)**2)
			chisq_light += chisq_light
			self.logger.info("chisq light.  : %f (reduced: %f) dof=%s" % (chisq_light, chisq_light/dof_light, dof_light))
		else:
			chisq_light = 0
			
		if print_info:
			self.logger.info("chisq total     : %f (reduced: %f)" % (chisq_tot, chisq_tot/dof_reg))
		#filename = os.path.join(self.dirname, ...
		if self.use_fourth_moment:
			chisqs.append(("kin", chisq_m2+chisq_m4+chisq_light))
		else:
			chisqs.append(("kin", chisq_light))
		if returnall:
			return chisq_tot, chisqs
		else:
			return chisq_tot
		
	def init(self):
		pass
	
	def run(self, args, opts, scope):
		self.init()
		self.solve()


class KnownWeights2(SolutionMoments):
	def __init__(self, modelpath):
		self.modelpath = modelpath
	
	def load(self):
		pass
		#self.storage_2d.load()
		#self.storage_3d.load()
		#self.binned_data.load()
	
	def solve(self):
		filename = os.path.join(self.modelpath, "df", "orbitweights.npy")
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		return c
	

class QP(SolutionMoments):
	def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d, storage_3d, fitdensity2d, fitdensity3d, binned_data, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix=""):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.binned_data = binned_data
		self.dfgrid = dfgrid
		self.regularization = regularization
		#self.regularization_delta = regularization_delta
		self.require_stability = require_stability
		self.use_jeans = use_jeans
		self.use_fourth_moment = use_fourth_moment
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.qp")
		self.postfix = postfix

		
	def init(self):
		self.storage_2d.load()
		self.storage_3d.load()
		self.binned_data.load()
		
	
	def findsolution(self):
		logger.info("using quadratic programming to find a solution")
		n_nodes = self.storage_2d.projectedmoments.shape[0]
		nI1nodes = self.dfgrid.n_I1
		nI2nodes = self.dfgrid.n_I2
		
		constraints = len(self.binned_data.moments[0])
		dof_3d = self.storage_3d.moments3d.shape[2]
		logger.info("dof 3d: %d" % dof_3d)
		logger.info("unknown weights: %d" % n_nodes)
		logger.info("constraints: %d" % constraints)
		logger.debug("shape of projectedmoments: %r" % (self.storage_2d.projectedmoments.shape,))
		rescale = 1.
		
		Gpos = cvxopt.spmatrix(-1, range(n_nodes), range(n_nodes))
		#print "Gpos", array(Gpos).shape, Gpos
		hpos = cvxopt.matrix(0, size=(n_nodes, 1))
		assert not (self.fitdensity2d and self.fitdensity3d), "cannot fit both 2d and 3d density"
		
		
		indices1 = []
		indices2 = []
		values = []
		h_stability = []
		#self.require_stability = True
		if self.require_stability: #self.dfgrid.order == 0:
			filename = os.path.join(self.modelpath, "data", "DoS.npy")
			DoS = load(filename)
			
			assert self.dfgrid.order == 0
			n = 0
			for i2 in range(nI2nodes):
				for i1 in range(nI1nodes-1):
						lowE = self.dfgrid.index_to_orbitnr(i1,i2, 0)
						highE = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
						# -lowE + highE <= 0
						indices1.append(lowE)
						indices2.append(n)
						values.append(-1./DoS[lowE])
						
						indices1.append(highE)
						indices2.append(n)
						values.append(1./DoS[highE])
						h_stability.append(0.0)
						n += 1
			G_stability_matrix = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		
		
		
		densitiesprojected = self.storage_2d.projectedmoments[:,0,:]
		densityprojected_target = self.binned_data.moments[0]
		
		moments3d = self.storage_3d.moments3d
		
		rborders = self.storage_3d.xborders
		r = self.storage_3d.x
		
		#light_model
		#density_target = self.storage_3d.moments3d[:,0,:]
		#print density_target.shape
		#density_target
		
		kinematic_m2 = self.binned_data.moments[2]
		schw_m2 = self.storage_2d.projectedmoments[:,2,:]
		e_kinematic_m2 = self.binned_data.e_moments[2]
		kinematic_m4 = self.binned_data.moments[4]
		schw_m4 = self.storage_2d.projectedmoments[:,4,:]
		e_kinematic_m4 = self.binned_data.e_moments[4]
		logger.debug("shape of densitiesprojected: %r" % (densitiesprojected.shape,))
		logger.debug("shape of densityprojected_target: %r" % (densityprojected_target.shape,))
		
		light_model = self.light_model
		Erborders = self.dfgrid.r_borders
		#density_Etarget = array([galaxy.stellar_profile.cumdensityr(galaxy.arcsec_to_kpc(Erborders[i]), galaxy.arcsec_to_kpc(Erborders[i+1])) for i in range(self.dfgrid.n_I1)])
		density_Etarget = array([light_model.light_profile.cumdensityr((Erborders[i]), (Erborders[i+1]), M=1.) for i in range(self.dfgrid.n_I1)])
		if 0:
			from kaplot import *
			box()
			#graph(arange(len(density_Etarget)), density_Etarget)
			l = 0.8
			E = self.dfgrid.E_borders[len(self.dfgrid.E_borders)/2]
			Lmax = self.profile_model.Lmax_at_E(E)
			#Ls = [l * self.profile_model.Lmax_at_E(E) for E in self.dfgrid.Es]
			v = []
			#for i, (L, E) in enumerate(zip(Ls, self.dfgrid.Es)):
			#	deltaE = self.dfgrid.E_borders[i+1] - self.dfgrid.E_borders[i]
			#	print deltaE
			#	v.append(self.profile_model.gEL(E,L) / deltaE)
			for l in arange(0, 1, 0.01):
				L = l * Lmax
				v.append(self.profile_model.gEL(E,L))
				print v
			graph(arange(len(v)), v)
			draw()
	
		# this part constrainst the density to within a certain range
		if self.fitdensity3d: # this part is for testing, it uses the 3d density
			diff = (density_target*0.01).max()
			Gmin = cvxopt.matrix(-transpose(densities))
			hmin = cvxopt.matrix(minimum(-(density_target-diff/2), 0), size=(Nr3d, 1))
			Gmax = cvxopt.matrix(transpose(densities))
			hmax = cvxopt.matrix(density_target+diff/2, size=(Nr3d, 1))
		if self.fitdensity2d:
			diff = max(densityprojected_target*0.01)
			Gmin = cvxopt.matrix(-1/rescale*transpose(densitiesprojected))
			hmin = cvxopt.matrix(minimum(-(densityprojected_target-diff/2), 0), size=(constraints, 1))
			Gmax = cvxopt.matrix(1/rescale*transpose(densitiesprojected))
			hmax = cvxopt.matrix(densityprojected_target+diff/2, size=(constraints, 1))
			
		if self.use_jeans:
			dr = rborders[1:] - rborders[:-1]
			rc = (r[1:] + r[:-1])/2
			drc = r[1:] - r[:-1]
			nu = moments3d[:,0,:] / (4*pi*r**2*dr)
			mask = moments3d[:,0,:] > 0
			
			varvr = moments3d[:,4,:] * 0
			varvr[mask] = moments3d[:,4,:][mask]/moments3d[:,0,:][mask]
			
			varvphi = moments3d[:,5,:] * 0
			varvphi[mask] = moments3d[:,5,:][mask]/moments3d[:,0,:][mask]
			
			varvtheta = moments3d[:,6,:] * 0
			varvtheta[mask] = moments3d[:,6,:][mask]/moments3d[:,0,:][mask]
			
			#varvphi = moments3d[:,5,:]/moments3d[:,0,:]
			#varvtheta = moments3d[:,6,:]/moments3d[:,0,:]
			assert all(~isnan(varvr))
			assert all(~isnan(varvphi))
			assert all(~isnan(varvtheta))
			#print varvr.shape, nu.shape
			rho_varvr = varvr * nu
			#print rho_varvr.shape, drc.shape
			jeans_term1 = (rho_varvr[:,1:] - rho_varvr[:,:-1]) / drc
			
			varvrc = (varvr[:,:-1] + varvr[:,1:])/2
			varvphic = (varvphi[:,:-1] + varvphi[:,1:])/2
			varvthetac = (varvtheta[:,:-1] + varvtheta[:,1:])/2
			
			nuc = (nu[:,:-1] + nu[:,1:])/2
			#jeans_term2 = 2 * -0.5 * nuc * varvrc/rc
			jeans_term2 = nuc * (2*varvrc - varvphic - varvthetac) / rc
			#print sum(moments3d[0,0])
			
			#dphidr = rc * 0.
			#for i in range(len(rc)):
			dphidr = self.profile_model.dphidr(rc)
			#print jeans_term1.shape, jeans_term2.shape
			jeans =  jeans_term1 + jeans_term2
			jeans_term_rhs = -self.light_model.densityr(rc)/self.light_model.light_profile.M * self.profile_model.dphidr(rc)
			#print nu.shape, jeans.shape, jeans_term_rhs.shape
			#print jeans,jeans_term_rhs
			J = (jeans/jeans_term_rhs)
			#print J.max(), J.min()
			
			
			JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs))-1)
			Jhmin = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			
			#JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs)))
			#Jhmin = cvxopt.matrix(1.25, size=(constraints-1, 1))
			
			JGmax = cvxopt.matrix(-(transpose(jeans/jeans_term_rhs)-1))
			Jhmax = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			#JGmax = cvxopt.matrix((transpose(jeans/jeans_term_rhs)-1))
			#Jhmax = cvxopt.matrix(1.51, size=(constraints-1, 1))
			#print Jhmax
			#print jeans
			
			#print array(JGmin).shape
			#print array(JGmax).shape
			#print " "
			#print array(Jhmin).shape
			#print array(Jhmax).shape
			
			
			#target = graph(xc, log10(nuc * dphidr), color="orange", linestyle="dash")
			import sys
			#sys.exit(0)
			
		if self.regularization:
			Greg, hreg = self.regularization.assemble_system()
			delta = self.regularization.regularization_delta
		if 0: #regularisation:
			G = cvxopt.matrix([Gpos, Gmin, Gmax, Greg])
			h = cvxopt.matrix([hpos, hmin, hmax, hreg])
		else:
			Glist = []
			hlist = []
			if self.use_jeans:
				#G = cvxopt.matrix([Gpos, Gmin, Gmax])
				#h = cvxopt.matrix([hpos, hmin, hmax])
				G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin, JGmax])
				h = cvxopt.matrix([hpos, hmin, hmax, Jhmin, Jhmax])
				
				print array(Gpos).shape
				print array(Gmin).shape
				print array(Gmax).shape
				print array(JGmin).shape
				print array(JGmax).shape
				print array(G).shape
				print " "
				print array(hpos).shape
				print array(hmin).shape
				print array(hmax).shape
				print array(Jhmin).shape
				print array(Jhmax).shape
				print array(h).shape
				
				#G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin])
				#h = cvxopt.matrix([hpos, hmin, hmax, Jhmin])
			else:
				Glist.extend([Gpos, Gmin, Gmax])
				hlist.extend([hpos, hmin, hmax])
				if self.require_stability:
					Glist.append(cvxopt.matrix(transpose(G_stability_matrix)))
					hlist.append(cvxopt.matrix(h_stability))
					
				G = cvxopt.matrix(Glist)
				h = cvxopt.matrix(hlist)
		
		P = cvxopt.spmatrix(0., range(n_nodes), range(n_nodes))
		q = cvxopt.matrix(0., size=(n_nodes, 1))
		if 0:
			P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
			x = array(sigmar**2)
			mask = density_target > 0
			x[mask] *= density_target[mask] 
			q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
		else:
			if 0:
				Pt = hstack((varvloss/kinematics.em2, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
				x *= hstack((densityprojected_target2, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			elif 0:
				print varvloss.shape, n, n_nodes, Greg.shape
				Pt = hstack((Greg/delta, varvloss/kinematics.em2))#, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((array(hreg)/delta, kinematics.m2/kinematics.em2))#, kinematics.m4/kinematics.em4))
				x *= hstack((ones(n), densityprojected_target2))#, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			else:
				#print varvloss.shape, n, orbits, Greg.shape
				#print Greg.shape
				#print kinematic_m2.shape
				#print schw_m2.shape
				Pt_list = []
				if self.regularization:
					Pt_list.append(Greg/delta)
				Pt_list.append((kinematic_m2- schw_m2/densityprojected_target)/(e_kinematic_m2))
				#Pt_list.append((kinematic_m2*densityprojected_target- schw_m2)/(e_kinematic_m2))
				if self.use_fourth_moment:
					Pt_list.append((kinematic_m4- schw_m4/densityprojected_target)/(e_kinematic_m4))
					#Pt = 1/rescale*hstack((Greg/delta, kinematic_m2/e_kinematic_m2*densitiesprojected- schw_m2/e_kinematic_m2))
					#Pt = 1/rescale*hstack((Greg/delta,\
					#(kinematic_m2*densitiesprojected- schw_m2)/(densityprojected_target*e_kinematic_m2),\
					#(kinematic_m4*densitiesprojected - schw_m4)/(densityprojected_target*e_kinematic_m4) ))
					#Pt = 1/rescale*hstack((Greg/delta,\
					#(kinematic_m2- schw_m2/densityprojected_target)/(e_kinematic_m2),\
					#(kinematic_m4- schw_m4/densityprojected_target)/(e_kinematic_m4) ))
					
					#, kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
				#else:
					#Pt = hstack(((kinematic_m2/e_kinematic_m2*densityprojected_target- schw_m2/e_kinematic_m2),\
					#	(kinematic_m4/e_kinematic_m4*densityprojected_target- schw_m4/e_kinematic_m4)))
				#	Pt = hstack((	(kinematic_m2 - schw_m2/densityprojected_target)/e_kinematic_m2,\
				#					(kinematic_m4 - schw_m4/densityprojected_target)/e_kinematic_m4 ))
							# kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
				Pt = hstack(tuple(Pt_list))
				test = tensordot(Pt, Pt, axes=([1],[1]))
				P = cvxopt.matrix(test)
				x_list = []
				if self.regularization:
					x_list.append(array(hreg)/delta)
					#x = hstack((array(hreg)/delta, zeros(constraints), zeros(constraints)))#, zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
					#x *= hstack((ones(n), zeros(constraints), zeros(constraints)))
				#else:
					#x = hstack((zeros(constraints), zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
					#x *= hstack((zeros(constraints), zeros(constraints)))
				x_list.append(zeros(constraints))
				if self.use_fourth_moment:
					x_list.append(zeros(constraints))
				x = hstack(tuple(x_list))
				#print Pt.shape, x.shape
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
				#Pt = ((kinematics.m2/kinematics.em2*densitiesprojected - varvloss/kinematics.em2))
				#P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				#x = ((kinematics.m2/kinematics.em2))
				#x *= ((densityprojected_target2))
			
			
		A = cvxopt.matrix(1/rescale, size=(1, n_nodes))
		b = cvxopt.matrix([1.0])
		
		orbitweights, chisq_solution = self.solve_problem(P, q, G, h, A, b)
		#cvxopt.solvers.options["maxiters"] = 500
		#solution = cvxopt.solvers.qp(P, q, G, h, A, b)
		#chisq_solution = solution['primal objective']
		#import pdb
		#pdb.set_trace()
		
		#print solution['x']
		#orbitweights = array(solution['x']).flatten()
		if self.regularization:
			self.P_reg = tensordot(Greg/delta, Greg/delta, axes=[(1,), (1,)])
		Ph = (kinematic_m2 - schw_m2/densityprojected_target)/e_kinematic_m2
		self.P_m2 = tensordot(Ph, Ph, axes=[(1,), (1,)])
		Ph = (kinematic_m4 - schw_m4/densityprojected_target)/e_kinematic_m4
		self.P_m4 = tensordot(Ph, Ph, axes=[(1,), (1,)])
		P = tensordot(Pt, Pt, axes=[(1,), (1,)])
		chisq_tot = tensordot(orbitweights, tensordot(orbitweights, P, axes=[(0,), (0,)]), axes=[(0,),(0,)]) +\
				tensordot(orbitweights, q, axes=[(0,), (0,)])
		logger.info("chisq_tot: %s %s" % (chisq_tot, chisq_solution*2))
		self.hessian = P
		
		
		return orbitweights
		
	def solve_problem(self, P, q, G, h, A, b):
		cvxopt.solvers.options["maxiters"] = 500
		cvxopt.solvers.options['show_progress'] = False
		solution = cvxopt.solvers.qp(P, q, G, h, A, b)
		chisq_solution = solution['primal objective']
		orbitweights = array(solution['x']).flatten()
		return orbitweights, chisq_solution

class KnownWeights(QP):
	"""def __init__(self, modelpath, schwsetname, schwmodelname, dfgrid, storage_2d, storage_3d, binned_data, postfix=""):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dfgrid = dfgrid
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.binned_data = binned_data
		self.light_model = dfgrid.light_model
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.postfix = postfix
		self.logger = logging.getLogger("gd.schw.solution.known")"""
		
		
	"""def load(self):
		self.binned_data.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_3d.load()"""
	
	def _load(self):
		self.storage_2d.load()
		self.storage_3d.load()
		self.binned_data.load()
	
	def solve_problem(self, P, q, G, h, A, b):
		filename = os.path.join(self.modelpath, "df", "orbitweights" +self.postfix +".npy")
		logger.info("using known solution from: %s" % filename)
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		return c, "?"
	
	

class QPGeneric(SolutionMoments):
	#def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d, storage_3d, fitdensity2d, fitdensity3d, binned_data, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix=""):
	def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d_m0, storage_2d_m2, storage_2d_m4, storage_3d, storage_3d_log, fitdensity2d, fitdensity3d, aperture_light, binned_data_m2, binned_data_m4, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix="", dfname=None):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d_m0 = storage_2d_m0
		self.storage_2d_m2 = storage_2d_m2
		self.storage_2d_m4 = storage_2d_m4
		self.storage_3d = storage_3d
		self.storage_3d_log = storage_3d_log
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.aperture_light = aperture_light
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.dfgrid = dfgrid
		self.regularization = regularization
		#self.regularization_delta = regularization_delta
		self.require_stability = require_stability
		self.use_jeans = use_jeans
		self.use_fourth_moment = use_fourth_moment
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.qp")
		self.postfix = postfix
		self.dfname = dfname

	def init(self):
		self.storage_2d_m0.load()
		self.storage_2d_m2.load()
		self.storage_2d_m4.load()
		self.storage_3d.load()
		self.storage_3d_log.load()
		self.aperture_light.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		
	
	def findsolution(self):
		logger.info("using quadratic programming to find a solution")
		n_nodes = self.storage_2d_m0.projectedmoments.shape[0]
		nI1nodes = self.dfgrid.n_I1
		nI2nodes = self.dfgrid.n_I2
		
		constraints_m0 = len(self.aperture_light.radial_surface_densities)
		constraints_m2 = len(self.binned_data_m2.moments[0])
		constraints_m4 = len(self.binned_data_m4.moments[0])
		dof_3d = self.storage_3d.moments3d.shape[2]
		logger.info("dof 3d: %d" % dof_3d)
		logger.info("unknown weights: %d" % n_nodes)
		logger.info("constraints m0: %d" % constraints_m0)
		logger.info("constraints m2: %d" % constraints_m2)
		logger.info("constraints m4: %d" % constraints_m4)
		logger.debug("shape of projectedmoments m0: %r" % (self.storage_2d_m0.projectedmoments.shape,))
		logger.debug("shape of projectedmoments m2: %r" % (self.storage_2d_m2.projectedmoments.shape,))
		logger.debug("shape of projectedmoments m4: %r" % (self.storage_2d_m4.projectedmoments.shape,))
		rescale = 1.
		
		Gpos = cvxopt.spmatrix(-1, range(n_nodes), range(n_nodes))
		#print "Gpos", array(Gpos).shape, Gpos
		hpos = cvxopt.matrix(0, size=(n_nodes, 1))
		#assert not (self.fitdensity2d and self.fitdensity3d), "cannot fit both 2d and 3d density"
		
		
		indices1 = []
		indices2 = []
		values = []
		h_stability = []
		#self.require_stability = True
		if self.require_stability: #self.dfgrid.order == 0:
			filename = os.path.join(self.modelpath, "data", "DoS.npy")
			DoS = load(filename)
			
			assert self.dfgrid.order == 0
			n = 0
			for i2 in range(nI2nodes):
				for i1 in range(nI1nodes-1):
						lowE = self.dfgrid.index_to_orbitnr(i1,i2, 0)
						highE = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
						# -lowE + highE <= 0
						indices1.append(lowE)
						indices2.append(n)
						values.append(-1./DoS[lowE])
						
						indices1.append(highE)
						indices2.append(n)
						values.append(1./DoS[highE])
						h_stability.append(0.0)
						n += 1
			G_stability_matrix = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		
		
		
		densitiesprojected_m0 = self.storage_2d_m0.projectedmoments[:,0,:]
		densitiesprojected_m2 = self.storage_2d_m2.projectedmoments[:,0,:]
		densitiesprojected_m4 = self.storage_2d_m4.projectedmoments[:,0,:]
		densityprojected_target_m0 = self.aperture_light.radial_surface_densities
		densityprojected_target_m2 = self.binned_data_m2.moments[0]
		densityprojected_target_m4 = self.binned_data_m4.moments[0]
		
		moments3d = self.storage_3d.moments3d
		
		rborders = self.storage_3d.xborders
		r = self.storage_3d.x
		
		#light_model
		#density_target = self.storage_3d.moments3d[:,0,:]
		#print self.storage_3d
		#print rborders
		#dsa
		density_target  = 1
		#print density_target.shape
		#density_target
		
		kinematic_m2 = self.binned_data_m2.moments[2]
		schw_m2 = self.storage_2d_m2.projectedmoments[:,2,:]
		e_kinematic_m2 = self.binned_data_m2.e_moments[2] #* 10
		
		#print kinematic_m2
		#print e_kinematic_m2
		
		kinematic_m4 = self.binned_data_m4.moments[4]
		schw_m4 = self.storage_2d_m4.projectedmoments[:,4,:]
		e_kinematic_m4 = self.binned_data_m4.e_moments[4] #* 10
		#print kinematic_m4
		#print e_kinematic_m4
		
		schw_m0 = self.storage_2d_m0.projectedmoments[:,0,:]
		m0 = densityprojected_target_m0[:]
		e_m0 = 1e-6 * m0**0.5 # * 0 + 0.001
		#print m0, m0.shape
		
		
		logger.debug("shape of densitiesprojected m0: %r" % (densitiesprojected_m0.shape,))
		logger.debug("shape of densitiesprojected m2: %r" % (densitiesprojected_m2.shape,))
		logger.debug("shape of densitiesprojected m4: %r" % (densitiesprojected_m4.shape,))
		logger.debug("shape of densityprojected_target m0: %r" % (densityprojected_target_m0.shape,))
		logger.debug("shape of densityprojected_target m2: %r" % (densityprojected_target_m2.shape,))
		logger.debug("shape of densityprojected_target m4: %r" % (densityprojected_target_m4.shape,))
		
		light_model = self.light_model
		Erborders = self.dfgrid.r_borders
		#density_Etarget = array([galaxy.stellar_profile.cumdensityr(galaxy.arcsec_to_kpc(Erborders[i]), galaxy.arcsec_to_kpc(Erborders[i+1])) for i in range(self.dfgrid.n_I1)])
		density_Etarget = array([light_model.light_profile.cumdensityr((Erborders[i]), (Erborders[i+1]), M=1.) for i in range(self.dfgrid.n_I1)])
		if 0:
			from kaplot import *
			box()
			#graph(arange(len(density_Etarget)), density_Etarget)
			l = 0.8
			E = self.dfgrid.E_borders[len(self.dfgrid.E_borders)/2]
			Lmax = self.profile_model.Lmax_at_E(E)
			#Ls = [l * self.profile_model.Lmax_at_E(E) for E in self.dfgrid.Es]
			v = []
			#for i, (L, E) in enumerate(zip(Ls, self.dfgrid.Es)):
			#	deltaE = self.dfgrid.E_borders[i+1] - self.dfgrid.E_borders[i]
			#	print deltaE
			#	v.append(self.profile_model.gEL(E,L) / deltaE)
			for l in arange(0, 1, 0.01):
				L = l * Lmax
				v.append(self.profile_model.gEL(E,L))
				print v
			graph(arange(len(v)), v)
			draw()
	
		# this part constrainst the density to within a certain range
		if self.fitdensity3d: # this part is for testing, it uses the 3d density
			binmin, binmax = 0, 249#250./6 * 1.1, 250./6*(6-1.75)
			binmin, binmax = 250./6 * (-2.8+3), 250./6*(1.2+3.)
			#binmin = 0
			#print self.storage_3d_log.x.shape
			#print self.storage_3d_log.x[binmin], self.storage_3d_log.x[binmax]
			#print self.storage_3d_log.moments3d.shape
			schw_masses3d = self.storage_3d_log.moments3d[:,0,binmin:binmax]
			logborders = self.storage_3d_log.xborders[binmin:binmax+1]
			borders = 10**logborders
			masses3d_target = array([self.light_model.cumdensityr(r1, r2, M=1.) for r1, r2 in zip(borders[:-1], borders[1:])])
			N_light = 20000
			e_masses3d = (masses3d_target/N_light)**0.5 #**1.5 
			self.schw_masses3d  = schw_masses3d
			self.masses3d_target = masses3d_target
			self.e_masses3d = e_masses3d
			
			
			#print schw_masses3d.shape, masses3d_target.shape
			#dsa
			if 0:
				diff = (density_target*0.01).max()
				Gmin = cvxopt.matrix(-transpose(densities))
				hmin = cvxopt.matrix(minimum(-(density_target-diff/2), 0), size=(Nr3d, 1))
				Gmax = cvxopt.matrix(transpose(densities))
				hmax = cvxopt.matrix(density_target+diff/2, size=(Nr3d, 1))
		skipbins = 0
		if self.fitdensity2d:
			densityprojected_target_m0 = densityprojected_target_m0[skipbins:]
			constraints_m0 = constraints_m0 - skipbins
			#print densitiesprojected_m0.shape
			densitiesprojected_m0 = densitiesprojected_m0[:,skipbins:]
			#print densitiesprojected_m0.shape
			#dsa
			#densityprojected_target_m0[0] *= 0.85
			#densityprojected_target_m0[1] *= 0.90
			diff = max(densityprojected_target_m0*0.01)
			Gmin = cvxopt.matrix(-1/rescale*transpose(densitiesprojected_m0))
			hmin = cvxopt.matrix(minimum(-(densityprojected_target_m0-diff/2), 0), size=(constraints_m0, 1))
			
			#Gmax = cvxopt.matrix(transpose(densitiesprojected_m0)*1.)
			#hmax = cvxopt.matrix(densityprojected_target_m0*1.05, size=(constraints_m0, 1))
			Gmax = cvxopt.matrix(1/rescale*transpose(densitiesprojected_m0))
			hmax = cvxopt.matrix(densityprojected_target_m0+diff/2, size=(constraints_m0, 1))
			
		if self.use_jeans:
			dr = rborders[1:] - rborders[:-1]
			rc = (r[1:] + r[:-1])/2
			drc = r[1:] - r[:-1]
			nu = moments3d[:,0,:] / (4*pi*r**2*dr)
			mask = moments3d[:,0,:] > 0
			
			varvr = moments3d[:,1,:] * 0
			varvr[mask] = moments3d[:,1,:][mask]/moments3d[:,0,:][mask]
			
			varvphi = moments3d[:,2,:] * 0
			varvphi[mask] = moments3d[:,2,:][mask]/moments3d[:,0,:][mask]/2
			
			varvtheta = moments3d[:,2,:] * 0
			varvtheta[mask] = moments3d[:,2,:][mask]/moments3d[:,0,:][mask]/2
			
			#varvphi = moments3d[:,5,:]/moments3d[:,0,:]
			#varvtheta = moments3d[:,6,:]/moments3d[:,0,:]
			assert all(~isnan(varvr))
			assert all(~isnan(varvphi))
			assert all(~isnan(varvtheta))
			#print varvr.shape, nu.shape
			rho_varvr = varvr * nu
			#print rho_varvr.shape, drc.shape
			jeans_term1 = (rho_varvr[:,1:] - rho_varvr[:,:-1]) / drc
			
			varvrc = (varvr[:,:-1] + varvr[:,1:])/2
			varvphic = (varvphi[:,:-1] + varvphi[:,1:])/2
			varvthetac = (varvtheta[:,:-1] + varvtheta[:,1:])/2
			
			nuc = (nu[:,:-1] + nu[:,1:])/2
			#jeans_term2 = 2 * -0.5 * nuc * varvrc/rc
			jeans_term2 = nuc * (2*varvrc - varvphic - varvthetac) / rc
			#print sum(moments3d[0,0])
			
			#dphidr = rc * 0.
			#for i in range(len(rc)):
			dphidr = self.profile_model.dphidr(rc)
			#print jeans_term1.shape, jeans_term2.shape
			jeans =  jeans_term1 + jeans_term2
			jeans_term_rhs = -self.light_model.densityr(rc)/self.light_model.light_profile.M * self.profile_model.dphidr(rc)
			#print nu.shape, jeans.shape, jeans_term_rhs.shape
			#print jeans,jeans_term_rhs
			J = (jeans/jeans_term_rhs)
			#print J.max(), J.min()
			
			
			JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs))-1)
			Jhmin = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			
			#JGmin = cvxopt.matrix((transpose(jeans/jeans_term_rhs)))
			#Jhmin = cvxopt.matrix(1.25, size=(constraints-1, 1))
			
			JGmax = cvxopt.matrix(-(transpose(jeans/jeans_term_rhs)-1))
			Jhmax = cvxopt.matrix(self.jeans_fraction, size=(dof_3d-1, 1))
			#JGmax = cvxopt.matrix((transpose(jeans/jeans_term_rhs)-1))
			#Jhmax = cvxopt.matrix(1.51, size=(constraints-1, 1))
			#print Jhmax
			#print jeans
			
			#print array(JGmin).shape
			#print array(JGmax).shape
			#print " "
			#print array(Jhmin).shape
			#print array(Jhmax).shape
			
			
			#target = graph(xc, log10(nuc * dphidr), color="orange", linestyle="dash")
			import sys
			#sys.exit(0)
			
		if self.regularization:
			Greg, hreg = self.regularization.assemble_system()
			delta = self.regularization.regularization_delta
		if 0: #regularisation:
			G = cvxopt.matrix([Gpos, Gmin, Gmax, Greg])
			h = cvxopt.matrix([hpos, hmin, hmax, hreg])
		else:
			Glist = []
			hlist = []
			if self.use_jeans:
				#G = cvxopt.matrix([Gpos, Gmin, Gmax])
				#h = cvxopt.matrix([hpos, hmin, hmax])
				#G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin, JGmax])
				#h = cvxopt.matrix([hpos, hmin, hmax, Jhmin, Jhmax])
				G = cvxopt.matrix([Gpos, Gmax, JGmin, JGmax])
				h = cvxopt.matrix([hpos, hmax, Jhmin, Jhmax])
				G = cvxopt.matrix([Gpos, JGmin, JGmax])
				h = cvxopt.matrix([hpos, Jhmin, Jhmax])
				
				print array(Gpos).shape
				print array(Gmin).shape
				print array(Gmax).shape
				print array(JGmin).shape
				print array(JGmax).shape
				print array(G).shape
				print " "
				print array(hpos).shape
				print array(hmin).shape
				print array(hmax).shape
				print array(Jhmin).shape
				print array(Jhmax).shape
				print array(h).shape
				
				#G = cvxopt.matrix([Gpos, Gmin, Gmax, JGmin])
				#h = cvxopt.matrix([hpos, hmin, hmax, Jhmin])
			else:
				if 0:
					Glist.extend([Gpos*1., Gmin, Gmax, ])
					hlist.extend([hpos*1., hmin, hmax])
				else:
					Glist.extend([Gpos*1., Gmax, Gmin])
					hlist.extend([hpos*1., hmax, hmin])
					#Glist.extend([Gpos*1.])
					#hlist.extend([hpos*1.])
					pass
				#Glist.extend([Gpos])
				#hlist.extend([hpos])
				if self.require_stability:
					Glist.append(cvxopt.matrix(transpose(G_stability_matrix)))
					hlist.append(cvxopt.matrix(h_stability))
					
				G = cvxopt.matrix(Glist)
				h = cvxopt.matrix(hlist)
				#G = None
		
		P = cvxopt.spmatrix(0., range(n_nodes), range(n_nodes))
		q = cvxopt.matrix(0., size=(n_nodes, 1))
		if 0:
			P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
			x = array(sigmar**2)
			mask = density_target > 0
			x[mask] *= density_target[mask] 
			q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
		else:
			if 0:
				Pt = hstack((varvloss/kinematics.em2, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
				x *= hstack((densityprojected_target2, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			elif 0:
				print varvloss.shape, n, n_nodes, Greg.shape
				Pt = hstack((Greg/delta, varvloss/kinematics.em2))#, m4loss/kinematics.em4))
				P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				x = hstack((array(hreg)/delta, kinematics.m2/kinematics.em2))#, kinematics.m4/kinematics.em4))
				x *= hstack((ones(n), densityprojected_target2))#, densityprojected_target2))
				q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
			else:
				#print varvloss.shape, n, orbits, Greg.shape
				#print Greg.shape
				#print kinematic_m2.shape
				#print schw_m2.shape
				if 0:
					Pt_list = []
					rho = 0
					rho = 0.867819669542 # numerical result for a gaussian
					if self.regularization:
						Pt_list.append(Greg/delta)
					Pt_list.append((kinematic_m2- schw_m2/densityprojected_target_m2)/(e_kinematic_m2))
					if self.use_fourth_moment:
						Pt_list.append((kinematic_m4- schw_m4/densityprojected_target_m4)/(e_kinematic_m4)/sqrt(1-rho**2))
					Ptl_list = list(Pt_list)
					Ptr_list = list(Pt_list)
					Ptl_list.append((kinematic_m2- schw_m2/densityprojected_target_m2)/(e_kinematic_m2)/sqrt(1-rho**2) * sqrt(2*rho))
					Ptr_list.append((kinematic_m4- schw_m4/densityprojected_target_m4)/(e_kinematic_m4)/sqrt(1-rho**2) * -sqrt(2*rho))
					Ptl = hstack(tuple(Ptl_list))
					Ptr = hstack(tuple(Ptl_list))
					test = tensordot(Ptl, Ptr, axes=([1],[1]))
					P = cvxopt.matrix(test)
					x_list = []
					if self.regularization:
						x_list.append(array(hreg)/delta)
					x_list.append(zeros(constraints_m2))
					if self.use_fourth_moment:
						x_list.append(zeros(constraints_m4))
					x_list.append(zeros(constraints_m4))
					x = hstack(tuple(x_list))
					#print Pt.shape, x.shape
					#import pdb
					#pdb.set_trace()
					#print constraints_m2, constraints_m4
					assert (constraints_m2) == (constraints_m4)
					#print n_nodes, x.shape, Ptl.shape
					q = cvxopt.matrix(-tensordot(x*-0, Ptl, axes=([0], [1])), size=(n_nodes, 1))
				else:
					Pt_list = []
					#self.regularization = True
					if self.regularization:
						Pt_list.append(Greg/delta)
					Pt_list.append((kinematic_m2- schw_m2/densityprojected_target_m2)/(e_kinematic_m2))
					#Pt_list.append(schw_m2/densityprojected_target_m2/(e_kinematic_m2))
					if self.fitdensity3d:
						#Pt_list.append((masses3d_target - schw_masses3d)/(e_masses3d))
						Pt_list.append(schw_masses3d/e_masses3d)
					#print ">", sum(masses3d_target)
					#dsa
					#if self.fitdensity2d:
					#	Pt_list.append((m0 - schw_m0)/(e_m0))
					#Pt_list.append((schw_m0)/(e_m0))
					#print kinematic_m2.shape
					#print m0.shape, e_m0.shape, ((m0 - schw_m0)/(e_m0)).shape
					#print "dsa"
					#print schw_m0.shape
					#print schw_m2.shape
					#print Greg.shape
					#dsa
					#Pt_list.append((kinematic_m2*densityprojected_target- schw_m2)/(e_kinematic_m2))
					if self.use_fourth_moment:
						Pt_list.append((kinematic_m4- schw_m4/densityprojected_target_m4)/(e_kinematic_m4))
						#Pt_list.append(schw_m4/densityprojected_target_m4/e_kinematic_m4)
						
						#Pt = 1/rescale*hstack((Greg/delta, kinematic_m2/e_kinematic_m2*densitiesprojected- schw_m2/e_kinematic_m2))
						#Pt = 1/rescale*hstack((Greg/delta,\
						#(kinematic_m2*densitiesprojected- schw_m2)/(densityprojected_target*e_kinematic_m2),\
						#(kinematic_m4*densitiesprojected - schw_m4)/(densityprojected_target*e_kinematic_m4) ))
						#Pt = 1/rescale*hstack((Greg/delta,\
						#(kinematic_m2- schw_m2/densityprojected_target)/(e_kinematic_m2),\
						#(kinematic_m4- schw_m4/densityprojected_target)/(e_kinematic_m4) ))
						
						#, kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
					#else:
						#Pt = hstack(((kinematic_m2/e_kinematic_m2*densityprojected_target- schw_m2/e_kinematic_m2),\
						#	(kinematic_m4/e_kinematic_m4*densityprojected_target- schw_m4/e_kinematic_m4)))
					#	Pt = hstack((	(kinematic_m2 - schw_m2/densityprojected_target)/e_kinematic_m2,\
					#					(kinematic_m4 - schw_m4/densityprojected_target)/e_kinematic_m4 ))
								# kinematics.m4/kinematics.em4*densitiesprojected - m4loss/kinematics.em4))
					Pt = hstack(tuple(Pt_list))
					Ptl, Ptr = Pt, Pt
					test = tensordot(Pt, Pt, axes=([1],[1]))
					P = cvxopt.matrix(test)
					x_list = []
					if self.regularization:
						x_list.append(array(hreg)/delta)
						#x = hstack((array(hreg)/delta, zeros(constraints), zeros(constraints)))#, zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
						#x *= hstack((ones(n), zeros(constraints), zeros(constraints)))
					#else:
						#x = hstack((zeros(constraints), zeros(constraints))) #kinematics.m2/kinematics.em2, kinematics.m4/kinematics.em4))
						#x *= hstack((zeros(constraints), zeros(constraints)))
					x_list.append(kinematic_m2/e_kinematic_m2*0)
					if self.fitdensity3d:
						#x_list.append(zeros(len(masses3d_target)))
						x_list.append(masses3d_target/e_masses3d)
					#if self.fitdensity2d:
					#	x_list.append(zeros(len(m0)))
					#x_list.append(-m0/e_m0)
					
					if self.use_fourth_moment:
						x_list.append(kinematic_m4/e_kinematic_m4*0)
					x = hstack(tuple(x_list))
					#print Pt.shape, x.shape
					q = cvxopt.matrix(-tensordot(x, Pt, axes=([0], [1])), size=(n_nodes, 1))
				#Pt = ((kinematics.m2/kinematics.em2*densitiesprojected - varvloss/kinematics.em2))
				#P = cvxopt.matrix(tensordot(Pt, Pt, axes=([1],[1])) )
				#x = ((kinematics.m2/kinematics.em2))
				#x *= ((densityprojected_target2))
			
			
		A = cvxopt.matrix(1/rescale, size=(1, n_nodes))
		#print sum(m0)
		#dsa
		#b = cvxopt.matrix([sum(m0)])
		b = cvxopt.matrix([1.0])
		#A = None
		#b = None
		
		#print "P", array(P).shape
		#print "q", array(q).shape
		#print "G", array(G).shape
		#print "h", array(h).shape
		#G = None
		#h = None

		orbitweights, chisq_solution = self.solve_problem(P, q, G, h, A, b)
		
		#print masses3d_target, sum(masses3d_target)
		#x = tensordot(orbitweights, schw_masses3d, axes=[(0,), (0,)])
		#print x, sum(x)
		
		
		#cvxopt.solvers.options["maxiters"] = 500
		#solution = cvxopt.solvers.qp(P, q, G, h, A, b)
		#chisq_solution = solution['primal objective']
		#import pdb
		#pdb.set_trace()
		
		#print solution['x']
		#orbitweights = array(solution['x']).flatten()
		if 1:
			if self.regularization:
				self.P_reg = tensordot(Greg/delta, Greg/delta, axes=[(1,), (1,)])
			Ph = (kinematic_m2 - schw_m2/densityprojected_target_m2)/e_kinematic_m2
			self.P_m2 = tensordot(Ph, Ph, axes=[(1,), (1,)])
			Ph = (kinematic_m4 - schw_m4/densityprojected_target_m4)/e_kinematic_m4
			self.P_m4 = tensordot(Ph, Ph, axes=[(1,), (1,)])
			P = tensordot(Ptl, Ptr, axes=[(1,), (1,)])
		#P = tensordot(Ptl, Ptr, axes=[(1,), (1,)])
		P = array(P)
		self.hessian = P
		#print P.shape
		chisq_tot = 0.5*tensordot(orbitweights, tensordot(orbitweights, P, axes=[(0,), (0,)]), axes=[(0,),(0,)]) +\
				tensordot(orbitweights, q, axes=[(0,), (0,)])
		logger.info("chisq_tot: %s %s" % (chisq_tot, chisq_solution))
		#import pdb
		#pdb.set_trace()
		if 0:
			from kaplot import *
			#box()
			mozaic(2,2,box)
			x0 = orbitweights * 1.
			chisq_tot0 = tensordot(x0, tensordot(x0, P, axes=[(0,), (0,)]), axes=[(0,),(0,)]) +\
				tensordot(x0, q, axes=[(0,), (0,)])
			s = zeros((20, 8))
			dx = 1e-6
			for x in range(20):
				for y in range(8):
					index = x + y*20
					x1 = orbitweights * 1.
					x1[index] += dx
					x1 /= sum(x1)
					chisq_tot1 = tensordot(x1, tensordot(x1, P, axes=[(0,), (0,)]), axes=[(0,),(0,)]) +\
						tensordot(x1, q, axes=[(0,), (0,)])
					d = (chisq_tot1-chisq_tot0)/dx
					s[x, y] = d
			select(0, 0)
			i = log10(s**2)
			i -= i.max()
			i[i < -5] = -5
			indexedimage(i.T, colormap="whiteblack")
			select(0, 1)
			spos = s *1
			spos[spos < 0] = 1e-2
			spos = log10(spos)
			indexedimage(spos.T, colormap="whiteblack")
			select(1, 1)
			sneg = s * -1
			sneg[sneg < 0] = 1e-2
			sneg = log10(sneg)
			indexedimage(sneg.T, colormap="whiteblack")
			print s
			print s.min()
			draw()
		if 0:
			from kaplot import *
			#box()
			mozaic(2,2,box)
			m = P.max()
			P /= m
			i = log10(P)
			i[P < 1e-7] = -7
			print i
			indexedimage(i)
			s = zeros((20, 8))
			for x in range(20):
				for y in range(8):
					index = x + y*20
					s[x, y] = i[index,index]
			select(0, 1)
			indexedimage(s.T, colormap="whiteblack")
			
			draw()
		
		
		return orbitweights
		
	def solve_problem(self, P, q, G, h, A, b):
		cvxopt.solvers.options["maxiters"] = 500
		cvxopt.solvers.options['show_progress'] = False
		#initvals = zeros(160) + 1
		if 0:
			initvals = zeros((20, 8)) * 1.
			x = arange(8)
			for i in range(20):
				initvals[i,:] = x**3
			initvals /= sum(initvals)
			#initvals = initvals.reshape(160)
		solution = cvxopt.solvers.qp(P, q, G, h, A, b) #, initvals=initvals)
		chisq_solution = solution['primal objective']
		orbitweights = array(solution['x']).flatten()
		return orbitweights, chisq_solution
	

class KnownWeightsGeneric(QPGeneric):
	"""def __init__(self, modelpath, schwsetname, schwmodelname, dfgrid, storage_2d, storage_3d, binned_data, postfix=""):
		self.modelpath = modelpath
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.dfgrid = dfgrid
		self.storage_2d = storage_2d
		self.storage_3d = storage_3d
		self.binned_data = binned_data
		self.light_model = dfgrid.light_model
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.postfix = postfix
		self.logger = logging.getLogger("gd.schw.solution.known")"""
		
		
	"""def load(self):
		self.binned_data.load()
		self.storage_2d.init()
		self.storage_2d.load()
		self.storage_3d.load()"""
	
	def _load(self):
		self.storage_2d.load()
		self.storage_3d.load()
		self.binned_data.load()
	
	def solve_problem(self, P, q, G, h, A, b):
		filename = os.path.join(self.modelpath, "df", "orbitweights_" +self.dfname +".npy")
		logger.info("using known solution from: %s" % filename)
		#filename = os.path.join(self.modelpath, "df", "orbitweights_tang.npy")
		orbitweights = load(filename)
		c = orbitweights.flatten() 
		c /= sum(c)
		#print sum(c)
		#dsa
		return c, "?"
	



class OptGeneric(SolutionMoments):
	#def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d, storage_3d, fitdensity2d, fitdensity3d, binned_data, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix=""):
	def __init__(self, modelpath, light_model, profile_model, schwsetname, schwmodelname, storage_2d_m0, storage_2d_m2, storage_2d_m4, storage_3d, storage_3d_log, fitdensity2d, fitdensity3d, aperture_light, binned_data_m2, binned_data_m4, dfgrid,  regularization=None, require_stability=False, use_jeans=False, use_fourth_moment=False, jeans_fraction = 0.05, postfix="", dfname=None):
		self.modelpath = modelpath
		self.light_model = light_model
		self.profile_model = profile_model
		self.schwsetname = schwsetname
		self.schwmodelname = schwmodelname
		self.storage_2d_m0 = storage_2d_m0
		self.storage_2d_m2 = storage_2d_m2
		self.storage_2d_m4 = storage_2d_m4
		self.storage_3d = storage_3d
		self.storage_3d_log = storage_3d_log
		self.fitdensity2d = fitdensity2d
		self.fitdensity3d = fitdensity3d
		self.aperture_light = aperture_light
		self.binned_data_m2 = binned_data_m2
		self.binned_data_m4 = binned_data_m4
		self.dfgrid = dfgrid
		self.regularization = regularization
		#self.regularization_delta = regularization_delta
		self.require_stability = require_stability
		self.use_jeans = use_jeans
		self.use_fourth_moment = use_fourth_moment
		self.jeans_fraction = jeans_fraction 
		self.dirname = os.path.join(self.modelpath, "schw", self.schwsetname, self.schwmodelname)
		self.logger = logging.getLogger("gd.schw.solution.qp")
		self.postfix = postfix
		self.dfname = dfname

	def init(self):
		self.storage_2d_m0.load()
		self.storage_2d_m2.load()
		self.storage_2d_m4.load()
		self.storage_3d.load()
		self.storage_3d_log.load()
		self.aperture_light.load()
		self.binned_data_m2.load()
		self.binned_data_m4.load()
		
	def make_opts(self):
		densitiesprojected_m0 = self.storage_2d_m0.projectedmoments[:,0,:]
		densitiesprojected_m2 = self.storage_2d_m2.projectedmoments[:,0,:]
		densitiesprojected_m4 = self.storage_2d_m4.projectedmoments[:,0,:]
		densityprojected_target_m0 = self.aperture_light.radial_surface_densities
		densityprojected_target_m2 = self.binned_data_m2.moments[0]
		densityprojected_target_m4 = self.binned_data_m4.moments[0]


		kinematic_m2 = self.binned_data_m2.moments[2]
		schw_m2 = self.storage_2d_m2.projectedmoments[:,2,:]
		e_kinematic_m2 = self.binned_data_m2.e_moments[2] #* 10
		
		kinematic_m4 = self.binned_data_m4.moments[4]
		schw_m4 = self.storage_2d_m4.projectedmoments[:,4,:]
		e_kinematic_m4 = self.binned_data_m4.e_moments[4] #* 10
		
		schw_m0 = self.storage_2d_m0.projectedmoments[:,0,:]
		m0 = densityprojected_target_m0[:]
		e_m0 = 1e-6 * m0**0.5 # * 0 + 0.001

		
		opts = []
		As = []
		bs = []
		self.model_m2 = schw_m2/densityprojected_target_m2/e_kinematic_m2
		As.append(schw_m2/densityprojected_target_m2/e_kinematic_m2)
		bs.append(kinematic_m2/e_kinematic_m2)
		self.opt_m2 = mab.gd.gdfast_schw.OptimizationMatrixChiSquare(self.model_m2, kinematic_m2/e_kinematic_m2)
		self.opt_m4 = mab.gd.gdfast_schw.OptimizationMatrixChiSquare(schw_m4/densityprojected_target_m4/e_kinematic_m4, kinematic_m4/e_kinematic_m4)
		As.append(schw_m4/densityprojected_target_m4/e_kinematic_m4)
		bs.append(kinematic_m4/e_kinematic_m4)
		opts.append(self.opt_m2)
		opts.append(self.opt_m4)
		self.opts = opts
		if self.fitdensity3d: # this part is for testing, it uses the 3d density
			binmin, binmax = 0, 249#250./6 * 1.1, 250./6*(6-1.75)
			binmin, binmax = 250./6 * (-2.80+3), 250./6*(1.2+3.)
			#binmin = 0
			print self.storage_3d_log.x.shape
			print self.storage_3d_log.x[binmin], self.storage_3d_log.x[binmax]
			print self.storage_3d_log.moments3d.shape
			schw_masses3d = self.storage_3d_log.moments3d[:,0,binmin:binmax]
			logborders = self.storage_3d_log.xborders[binmin:binmax+1]
			borders = 10**logborders
			masses3d_target = array([self.light_model.cumdensityr(r1, r2, M=1.) for r1, r2 in zip(borders[:-1], borders[1:])])
			e_masses3d = 0.1*1.1 * masses3d_target**0.5 #**1.5
			
			self.opt_m0 = mab.gd.gdfast_schw.OptimizationMatrixChiSquare(schw_masses3d/e_masses3d, masses3d_target/e_masses3d)
			opts.append(self.opt_m0)
			As.append(schw_masses3d/e_masses3d)
			bs.append(masses3d_target/e_masses3d)
		if self.regularization:
			Greg, hreg = self.regularization.assemble_system()
			delta = self.regularization.regularization_delta
			self.opt_reg = mab.gd.gdfast_schw.OptimizationMatrixChiSquare(Greg/delta, array(hreg)/delta)
			opts.append(self.opt_reg)
			As.append(Greg/delta)
			bs.append(array(hreg)/delta)
		self.As = As
		self.bs = bs
		
	def findsolution(self):
		self.make_opts()
		opts = self.opts
		As = numpy.hstack(self.As)
		bs = numpy.hstack(self.bs)
		print "As", As.shape
		print "bs", bs.shape
		x, bla = scipy.optimize.nnls(transpose(As), bs)
		print x
		print "x", x.shape
		#dsa
		return x
		dsa
		
		class Dummy:
			pass
		counter = Dummy()
		counter.count = 0
		
		def f_and_g(x):
			#print x
			grad = x * 0
			gradn = x * 0
			logp = sum([opt.logp(x) for opt in opts]) #opt_m2.logp(x) * 1
			for opt in opts:
				opt.dlogpdx(x, grad)
			#self.photometry.optimizer.dlogpdx(x, grad)
			#opt_norm.dlogpdx(x, grad)
			counter.count += 1
			if (counter.count % 100) == 0: #debug:
				#print "%10f %10f %10f %10f %10f" % (logp,  sum(x), dot(totalmass_matrix, x/sum(x)), dot(ptotalmass_matrix, x/sum(x)), dot(losvds.T, x).sum() * delta_R * delta_v / dot(ptotalmass_matrix, x/sum(x)))
				print "%5.15f %10f" % (logp,  sum(x))
				pass
			
			return -logp, -grad
		
		Nparams = self.model_m2.shape[0]
		x = ones(Nparams)
		x /= sum(x)
		logger.info("%d free parameters" % Nparams)
		bounds = [(1e-10, 1) for i in range(Nparams)]
		#print schw_m2.shape
		
		x = scipy.optimize.fmin_l_bfgs_b(f_and_g, x, None, bounds=bounds, approx_grad=False, iprint=-1,factr=1e-2,maxfun=200000)[0]
		#x /= sum(x)
		self.logL = 1
		u = log(x)
		print u
		#dsa

		#print "bla"
		#ddsa
		return x
		
import emcee
class OptGenericExplore(object):
	def __init__(self, solution):
		self.solution = solution
		
	def run(self, args, opts, scope):
		print "bla"
		from kaplot import *
		
		self.solution.init()
		self.solution.load()
		self.solution.make_opts()
		mozaic(3,3,box)
		select(0, 0)
		#indexedimage(self.solution.orbitweights.reshape)
		
		I = self.solution.orbitweights * 1.
		dfgrid = self.solution.dfgrid
		I = I.reshape((dfgrid.n_I2, dfgrid.n_I1))
		#indexedimage(I)
		
		
		c  = self.solution.orbitweights
		Nparam = len(c)
		#for i in range(N):
		grad = c * 0.
		for opt in self.solution.opts:
			opt.dlogpdx(c, grad)
		I = grad
		I = I.reshape((dfgrid.n_I2, dfgrid.n_I1))
		print grad.min(), grad.max()
		I /= I.max()
		print I.max(), I.min()
		I[I<10**-22] = 10**-22
		I = log10(I)
		print I.max(), I.min()
		#I -= I.max()
		I[I<-4] = -4
		indexedimage(I)#, colormap="whiteblack")# mask=I>0)
		draw()
		dsa
			
			
		
		if 0:
		
			n = 200
			x = (arange(n)+0.5)/(n)
			logrs = dfgrid.logrmin + x * (dfgrid.logrmax - dfgrid.logrmin)
			ls = 1 * x
			df = array([[dfgrid(c, logr, l) for logr in logrs] for l in ls])
		
			resize = ((0, 0), (dfgrid.n_I1, dfgrid.n_I2))
			if 0:
				df = log10(df)
				df -= df.max()
				df[df<-6] = -6
			im = indexedimage((df), colormap="whiterainbow", resize=resize)
			
		#for i in range(Nparam):
		#	for j in 
		i, j = 3, 2
		x0 = c * 1.
		x1 = c * 1.
		dx = 1e-8
		x1[i] += dx
		grad0 = x0 * 0
		grad1 = x0 * 0
		
		self.solution.opt_m2.dlogpdx(x0, grad0)
		self.solution.opt_m2.dlogpdx(x1, grad1)
		print "grad", (grad1[j] - grad0[j])/dx
		#print -self.solution.model_m2[j,i]**2
		M = numpy.matrix(self.solution.model_m2)
		hessian = M * M.T
		print -hessian[i,j]
		print hessian
		u, s, v = numpy.linalg.svd(hessian)
		print s
		print u.shape
		print v.shape
		#hessian = hessian.I
		var = c * 0
		for i in range(Nparam):
			var[i] = hessian[i,i]
		var = var.reshape((dfgrid.n_I2, dfgrid.n_I1))
		print var.min(), var.max(), var.mean()
		select(1,0)
		#indexedimage(var)
		#dsa
		if 0:
			for j in range(16):# + [-1, -2, -3, -4, -5]:
				select(j/4,j%4)
				var = c * 0
				for i in range(Nparam):
					var[i] = v[j,i]
				var = var.reshape((dfgrid.n_I2, dfgrid.n_I1))
				print var.min(), var.max(), var.mean()
				#select(1,1)
				indexedimage(var)
			
		print "test"
		for i in [0, 1, -1]:
			var = c * 0
			for k in range(Nparam):
				var[k] = v[i,k]
			x1 = x0*1.
			x1 += var * dx
			dx = 1e-7
			print (self.solution.opt_m2.logp(x1)-self.solution.opt_m2.logp(x0))/dx
			
		nwalkers = len(s)*2
		print nwalkers
		print v.shape, x0.shape
		print array(v[0]).reshape(-1).shape
		#dsa
		#import pdb
		#pdb.set_trace()
		v = array(v)
		#p0 = array([(log10(x0) + x0*(v[:,walker/2].reshape(-1) * 1/s[walker/2] * random.normal(0,1) * 1e-18)).reshape(-1) for walker in range(nwalkers)])
		#print v[:,walker/2].reshape(-1)
		p0 = array([(log10(x0) + x0*(random.normal(0,maximum(1e-10, v[:,walker/2].reshape(-1)**2) * 1/s[walker/2]) * 1e-10)).reshape(-1) for walker in range(nwalkers)])
		#dsa
		print p0.shape, x0.shape
		ndim = len(s)
		p0 = p0.reshape(ndim*2, ndim)
		print p0.shape, x0.shape
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,args=[self])
		pos, prob, state = sampler.run_mcmc(p0, 15)
		print prob
		print prob.shape
		#sampler.reset()
		#pos, prob, state = sampler.run_mcmc(pos, 20)
		
		
		#histogram(prob[len(prob)/2:])
		histogram(sampler.lnprobability[:,:].reshape(-1), color="red", bincount=100)
		#histogram(sampler.lnprobability[:,100:].reshape(-1), bincount=100)
		print sampler.acceptance_fraction
		select(2,0)
		for k in range(nwalkers):
			graph(sampler.lnprobability[k])
		#draw()
		
		print sampler.flatchain.shape
		
		storage_3d = scope["storage_3d"]
		moments = storage_3d.moments3d
		print moments.shape
		print sampler.chain.shape
		
		select(0,1)
		for w in range(nwalkers):
			for i in range(sampler.chain.shape[1]):
				m = tensordot(10**sampler.chain[walker,i,:], moments, axes=[(0,), (0,)])
				sigmarsq = m[1]/m[0]
				sigmatsq = m[2]/m[0]
				beta = 1 - sigmatsq/2/sigmarsq
				#print beta.shape
				graph(beta)
			#self.solution.calculate_solution_moments
		#dsa
		ylim(-3,1)
		vline(0)
		
		
		
		draw()
		
	def logL(self, u):
		logL = 0
		x = 10**u
		for opt in self.solution.opts:
			logL += opt.logp(x)
		return logL
		
def lnprob(x, fitter):
	return fitter.logL(x)