# -*- coding: utf-8 -*-
from numpy import *
import mab.gd.logging as logging
import cvxopt
logger = logging.getLogger("gd.schw.regularization")
			
class RegularizationQP_toiso(object):
	def __init__(self, modelpath, regularization_delta, dfgrid, light_model, orbitweight_calculator, scale=8):
		self.modelpath = modelpath
		self.regularization_delta = regularization_delta
		self.dfgrid = dfgrid
		self.light_model = light_model
		self.orbitweight_calculator = orbitweight_calculator
		self.scale = scale
		
	def assemble_system(self):
		weights = self.orbitweight_calculator.load().T
		#print weights.shape
		#Erborders = self.dfgrid.r_borders
		#density_Etarget = array([self.light_model.light_profile.cumdensityr((Erborders[i]), (Erborders[i+1]), M=1.) for i in range(self.dfgrid.n_I1)])
		
		n_nodes = self.dfgrid.get_dof()
		nI1nodes = self.dfgrid.n_I1
		nI2nodes = self.dfgrid.n_I2
		
		weights = weights.reshape((nI2nodes, nI1nodes))
		weights = weights.T
		if 0:
			from kaplot import *
			box()
			indexedimage(weights)
			draw()
		#print weights.shape
		
		indices1 = []
		indices2 = []
		values = []
		hreg = []
		#delta = 10000.1
		logger.info("regularization delta: %.5f" % self.regularization_delta)
		delta =  self.regularization_delta #.regdelta #sqrt(delta / (2*orbits))
		#delta = 1.0
		#delta = sqrt(delta / (2*orbits))
		logger.info("internal regularization delta: %.5f" % delta)
		n = 0
		s = 1.
		scale_i1 = 1 * s
		scale_i2 = self.scale * s
		rscale = 1e-4
		
		#s = 1
		#scale_i1 = 1 * s
		#scale_i2 = 8 * s
		#rscale = 1e-4
		
		#s = 0.01
		#scale_i1 = 1 * s
		#scale_i2 = 1 * s
		#rscale = 1e-4
		if 0:
			for i in range(4):
				print "=" * 10
				for vi in range(7):
					for ui in range(7):
						u = ui/6.
						v = vi/6.
						print "% 8.6e" % self.dfgrid.basis(i, u, v),
					print
				print "=" * 10
		#sys.exit(0)
		if self.dfgrid.order == 0:
			for i1 in range(1,nI1nodes-1):
				for i2 in range(nI2nodes):
					extras = 1 #/ ((i2 + 0.5) / (nI2nodes))
					o1 = self.dfgrid.index_to_orbitnr(i1-1,i2, 0) 
					o2 = self.dfgrid.index_to_orbitnr(i1+0,i2, 0)
					o3 = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
					# minimize -x_{i-1,j} + 2 x_{i,j} -x_{i+1,j} 
					r1 = 1 + (random.random() - 0.5) * rscale
					r2 = 1 + (random.random() - 0.5) * rscale
					r3 = 1 + (random.random() - 0.5) * rscale
					values.append(-1.*r1/weights[i1-1,i2]/scale_i1*extras)
					assert o1 < n_nodes
					assert o2 < n_nodes
					assert o3 < n_nodes, "o3 = %d, n_nodes = %d" % (o3, n_nodes)
					indices1.append(o1)
					indices2.append(n)
					
					values.append(2.*r2/weights[i1,i2]/scale_i1*extras)
					indices1.append(o2)
					indices2.append(n)
					
					values.append(-1.*r3/weights[i1+1,i2]/scale_i1*extras)
					indices1.append(o3)
					indices2.append(n)
					rh = (random.random() - 0.5) * delta * 1e-5
					hreg.append(rh)
					n += 1
			for i1 in range(0,nI1nodes):
				for i2 in range(1,nI2nodes-1):
					extras1 = 1 #/ ((i2- 0.5) / (nI2nodes))
					extras2 = 1 #/ ((i2 + 0.) / (nI2nodes))
					extras3 = 1 #/ ((i2 + 0.5) / (nI2nodes))
					o1 = self.dfgrid.index_to_orbitnr(i1,i2-1, 0) 
					o2 = self.dfgrid.index_to_orbitnr(i1,i2, 0)
					o3 = self.dfgrid.index_to_orbitnr(i1,i2+1, 0)
					# minimize -x_{i,j-1} + 2 x_{i,j} -x_{i,j+1} 
					r1 = 1 + (random.random() - 0.5) * rscale
					r2 = 1 + (random.random() - 0.5) * rscale
					r3 = 1 + (random.random() - 0.5) * rscale
					values.append(-1.*r1/weights[i1,i2-1]/scale_i2*extras1)
					indices1.append(o1)
					indices2.append(n)
					
					values.append(2.*r2/weights[i1,i2]/scale_i2*extras2)
					indices1.append(o2)
					indices2.append(n)
					
					values.append(-1.*r3/weights[i1,i2+1]/scale_i2*extras3)
					indices1.append(o3)
					indices2.append(n)
					rh = (random.random() - 0.5) * delta * 1e-5
					hreg.append(rh)
					n += 1
		Greg = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		hreg = array(hreg)
		return Greg, hreg
		
			
class RegularizationQP(object):
	def __init__(self, modelpath, regularization_delta, dfgrid, light_model, scale=8.):
		self.modelpath = modelpath
		self.regularization_delta = regularization_delta
		self.dfgrid = dfgrid
		self.light_model = light_model
		self.scale = scale
		
	def assemble_system(self):
		Erborders = self.dfgrid.r_borders
		density_Etarget = array([self.light_model.light_profile.cumdensityr((Erborders[i]), (Erborders[i+1]), M=1.) for i in range(self.dfgrid.n_I1)])
		
		n_nodes = self.dfgrid.get_dof()
		nI1nodes = self.dfgrid.n_I1
		nI2nodes = self.dfgrid.n_I2
		
		indices1 = []
		indices2 = []
		values = []
		hreg = []
		#delta = 10000.1
		logger.info("regularization delta: %.5f" % self.regularization_delta)
		delta =  self.regularization_delta #.regdelta #sqrt(delta / (2*orbits))
		#delta = 1.0
		#delta = sqrt(delta / (2*orbits))
		logger.info("internal regularization delta: %.5f" % delta)
		n = 0
		s = 1.
		scale_i1 = 1 * s #/ sqrt(self.dfgrid.n_I1/8.)
		scale_i2 = self.scale * s # / sqrt(self.dfgrid.n_I2/20.)
		rscale = 1e-4
		
		#s = 1
		#scale_i1 = 1 * s
		#scale_i2 = 8 * s
		#rscale = 1e-4
		
		#s = 0.01
		#scale_i1 = 1 * s
		#scale_i2 = 1 * s
		#rscale = 1e-4
		if 0:
			for i in range(4):
				print "=" * 10
				for vi in range(7):
					for ui in range(7):
						u = ui/6.
						v = vi/6.
						print "% 8.6e" % self.dfgrid.basis(i, u, v),
					print
				print "=" * 10
		#sys.exit(0)
		if self.dfgrid.order == 0:
			for i1 in range(1,nI1nodes-1):
				for i2 in range(nI2nodes):
					extras = 1 #/ ((i2 + 0.5) / (nI2nodes))
					o1 = self.dfgrid.index_to_orbitnr(i1-1,i2, 0) 
					o2 = self.dfgrid.index_to_orbitnr(i1+0,i2, 0)
					o3 = self.dfgrid.index_to_orbitnr(i1+1,i2, 0)
					# minimize -x_{i-1,j} + 2 x_{i,j} -x_{i+1,j} 
					r1 = 1 + (random.random() - 0.5) * rscale
					r2 = 1 + (random.random() - 0.5) * rscale
					r3 = 1 + (random.random() - 0.5) * rscale
					values.append(-1.*r1/density_Etarget[i1-1]/scale_i1*extras)
					assert o1 < n_nodes
					assert o2 < n_nodes
					assert o3 < n_nodes, "o3 = %d, n_nodes = %d" % (o3, n_nodes)
					indices1.append(o1)
					indices2.append(n)
					
					values.append(2.*r2/density_Etarget[i1]/scale_i1*extras)
					indices1.append(o2)
					indices2.append(n)
					
					values.append(-1.*r3/density_Etarget[i1+1]/scale_i1*extras)
					indices1.append(o3)
					indices2.append(n)
					rh = (random.random() - 0.5) * delta * 1e-5
					hreg.append(rh)
					n += 1
			for i1 in range(0,nI1nodes):
				for i2 in range(1,nI2nodes-1):
					extras1 = 1 #/ ((i2- 0.5) / (nI2nodes))
					extras2 = 1 #/ ((i2 + 0.) / (nI2nodes))
					extras3 = 1 #/ ((i2 + 0.5) / (nI2nodes))
					o1 = self.dfgrid.index_to_orbitnr(i1,i2-1, 0) 
					o2 = self.dfgrid.index_to_orbitnr(i1,i2, 0)
					o3 = self.dfgrid.index_to_orbitnr(i1,i2+1, 0)
					# minimize -x_{i,j-1} + 2 x_{i,j} -x_{i,j+1} 
					r1 = 1 + (random.random() - 0.5) * rscale
					r2 = 1 + (random.random() - 0.5) * rscale
					r3 = 1 + (random.random() - 0.5) * rscale
					values.append(-1.*r1/density_Etarget[i1]/scale_i2*extras1)
					indices1.append(o1)
					indices2.append(n)
					
					values.append(2.*r2/density_Etarget[i1]/scale_i2*extras2)
					indices1.append(o2)
					indices2.append(n)
					
					values.append(-1.*r3/density_Etarget[i1]/scale_i2*extras3)
					indices1.append(o3)
					indices2.append(n)
					rh = (random.random() - 0.5) * delta * 1e-5
					hreg.append(rh)
					n += 1
		Greg = array(cvxopt.matrix(cvxopt.spmatrix(values, indices1, indices2, size=(n_nodes, n))))
		hreg = array(hreg)
		return Greg, hreg
		
		