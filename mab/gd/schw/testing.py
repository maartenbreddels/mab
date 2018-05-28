from numpy import *
from kaplot import *
import mab.parallelize
import mab.utils.numpy

class Test(object):
	def __init__(self, modelpath, dfgrid, light_model, profile_model, solution, storage_3d, orbitweights):
		self.modelpath = modelpath
		self.dfgrid = dfgrid
		self.profile_model = profile_model
		self.light_model = light_model
		self.solution = solution
		self.storage_3d = storage_3d
		self.orbitweights = orbitweights
		
		
	def run(self, args):
		self.orbitweights.load()
		self.storage_3d.load()
		
		ini, scope = mab.gd.schw.configure.loadini(self.solution.modelpath, objectname="galaxy")
		galaxy = scope["galaxy"]
		jeans = galaxy.jeans()
		logr = arange(-3, 3, 0.1)
		r = 10**logr
		#r = arange(0, 1.25, 0.01)
		sigmar = jeans.sigmar(r)
		document(size="30cm,25cm")
		mozaic(3,3,box)
		select(0,0)
		graph(logr, sigmar)
		graph(logr, sigmar*sqrt(1+0.5))
		select(1,0)
		orbitweights = self.orbitweights.orbitweights.reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		df = self.orbitweights.df.reshape(self.dfgrid.n_I2, self.dfgrid.n_I1)
		indexedimage(df)
		
		
		select(2,1)
		for i in range(self.dfgrid.n_I1):#[::5]:
			red = float(i)/self.dfgrid.n_I1
			for j in range(self.dfgrid.n_I2):#[::4]:
				index = self.dfgrid.dof_index(i, j, 0)
				#print self.storage_3d.moments3d.shape
				#if (i == 15) and (j == 1):
				if 1:
					y = self.storage_3d.moments3d[index,4]
					n = self.storage_3d.moments3d[index,0]
					#y[n>0] /= n[n>0]
					#graph(self.storage_3d.x, y)#/(10**self.storage_3d.x)**3) 
				#self.storage_3d.moments[
		select(1,1)
		Epot = self.profile_model.potentialr(r)
		varr_tot = sigmar * 0
		vart_tot = sigmar * 0
		dens_tot = sigmar * 0
		varr_tot_check = self.storage_3d.moments3d[index,0] * 0
		counts = varr_tot_check * 0 
		Lmaxs = [self.profile_model.Lmax_at_E(E) for E in self.dfgrid.subgrid.E_borders]
		fs = []
		#self.modelpath = galaxy.modelpath
		filename = os.path.join(self.modelpath, "df", "fE.npy")
		#logger.info("loading f1(E) as: %s" % filename)
		dff = load(filename)
		
		filename = os.path.join(self.modelpath, "df", "fE-E.npy")
		#logger.info("loading corresponding energy as: %s" % filename)
		dfEs = load(filename)
		
		for i in range(self.dfgrid.n_I1): #[::5]:
			red = float(i)/self.dfgrid.n_I1
			print red
			for j in range(self.dfgrid.n_I2): #[::4]:
				green = float(j)/self.dfgrid.n_I1
				dither = self.dfgrid.dither
				I1 = i
				I2 = j
				E1 = self.dfgrid.E_borders[I1]
				E2 = self.dfgrid.E_borders[I1+1]
				L11 = self.dfgrid.l_borders[I2] * self.profile_model.Lmax_at_E(E1)
				L12 = self.dfgrid.l_borders[I2+1]*  self.profile_model.Lmax_at_E(E1)
				L21 = self.dfgrid.l_borders[I2] * self.profile_model.Lmax_at_E(E2)
				L22 = self.dfgrid.l_borders[I2+1]*  self.profile_model.Lmax_at_E(E2)
				for d1 in range(dither):
					for d2 in range(dither):
						Lmax1 = Lmaxs[i*dither+d1]
						Lmax2 = Lmaxs[(i)*dither+d1+1]
						E1 = self.dfgrid.subgrid.E_borders[i*dither+d1]
						E2 = self.dfgrid.subgrid.E_borders[(i)*dither+d1+1]
						L11 = self.dfgrid.subgrid.l_borders[j*dither+d2] * Lmax1
						L12 = self.dfgrid.subgrid.l_borders[(j)*dither+d2+1]*  Lmax1
						L21 = self.dfgrid.subgrid.l_borders[j*dither+d2] * Lmax2
						L22 = self.dfgrid.subgrid.l_borders[(j)*dither+d2+1]*  Lmax2
				
						def int_dens(r, L, E):
							return -2./3. * pi * (2.*E - 2. * Epot-L**2/r**2)**(3/2.) #* r**3 #* r**2#/r**3
							#return 2 * pi * L/r**2 / sqrt(2.*E - 2. * Epot-L**2/r**2)
						def int_varr(r, L, E):
							return -2./15. * pi * (2.*E - 2. * Epot-L**2/r**2)**(5/2.) #* r**3 #* r**2#/r**3
							#return 2 * pi * L/r**2 * sqrt(2.*E - 2. * Epot-L**2/r**2)
						def int_vart(r, L, E):
							return 2./15. * pi * (-4.*E + 4. * Epot - 3*L**2/r**2) * (2.*E - 2. * Epot-L**2/r**2)**(3/2.)
							#return 2 * pi * L/r**2 * sqrt(2.*E - 2. * Epot-L**2/r**2)
						
						Lmax = (Lmax1 + Lmax2)/2
						L1 = (L11 + L12)/2
						L2 = (L21 + L22)/2
						
						E = (E1+E2)/2
						L = (L1+L2)/2
						dE = E2 - E1
						dL = L2 - L1
						
						#mask = 
						#continue
						varr = ((int_varr(r, L2, E2) - int_varr(r, L1, E2)) - (int_varr(r, L2, E1) - int_varr(r, L1, E1)))
						#varr = int_varr(r, L, E) * dE * dL
						varr[isnan(varr)] = 0
						vart = ((int_vart(r, L2, E2) - int_vart(r, L1, E2)) - (int_vart(r, L2, E1) - int_vart(r, L1, E1)))
						#vart = int_vart(r, L, E) * dE * dL
						vart[isnan(vart)] = 0
						dens = (int_dens(r, L2, E2) - int_dens(r, L1, E2)) - (int_dens(r, L2, E1) - int_dens(r, L1, E1))
						#dens = int_dens(r, L, E) * dE * dL
						dens[isnan(dens)] = 0
						#if all(varr == 0):
						#	import pdb
						#	pdb.set_trace()
						#print ">", varr
						#if (i == 15) and (j == 1):
						#g = self.dfgrid.profile_model.gEL(E, L)
						if 1:
							y = (varr) * r ** 3  / (dE * dL)  # *r**2#* self.light_model.densityr(r)*r**3 
							#graph(logr, y/y.max(), color=Color(green, 0, 0))
							#graph(logr, y, color=Color(green, 0, 0))
							index = self.dfgrid.dof_index(i, j, 0)
							#y = self.storage_3d.moments3d[index,4]/self.storage_3d.moments3d[index,0]
							y = self.storage_3d.moments3d[index,4]
							y = y#/(10**self.storage_3d.x)**3
							#n = self.storage_3d.moments3d[index,0]
							#y[n>0] /= n[n>0]
							#graph(self.storage_3d.x, y/y.max(), color="red") 
							
							#print varr
						Ei = argmin(abs(dfEs-E))
						#print Ei,
						ddE = (dfEs[Ei+1] - dfEs[Ei])
						w = dff[Ei] * galaxy.fL(L) / (ddE) #/ Lmax # /ddE * dE #* L ** (-2*self.anisotropy_beta)
						#w = df[j,i] /(dE * dL)
						varr_tot += varr * w #df[j,i] /(dE * dL)
						vart_tot += vart * w
						dens_tot += dens * w #df[j,i] /(dE * dL)
						#varr_tot += varr * orbitweights[j,i] / (dE * dL * g) #df[j,i] #/ (dE * dL * g)
						#dens_tot += dens * orbitweights[j,i] / (dE * dL * g)#df[j,i]# / (dE * dL * g)
						fs.append(df[j,i])#/(dE * dL))
						index = self.dfgrid.index_to_orbitnr(i, j, 0)
						index = self.dfgrid.n_I2 * i + j
						index = self.dfgrid.n_I1 * j + i
						#print index
						#varr_tot_check += self.storage_3d.moments3d[index,4] * self.orbitweights.orbitweights[index]# * orbitweights[j,i]
						#counts += self.storage_3d.moments3d[index,0] * self.orbitweights.orbitweights[index]#orbitweights[j,i]
		select(0,1)
		#dens_tot /= sum(fs)
		#varr_tot /= sum(fs)
		print "tot f", sum(fs), sqrt(sum(fs)), sum(df)
		nu = self.light_model.densityr(r)/self.light_model.light_profile.M
		y = varr_tot * 1.0
		y[dens_tot>0] /= (dens_tot * r ** 0)[dens_tot>0] 
		#y[dens_tot>0] /= (dens_tot * r ** 0)[dens_tot>0]
		yr = y
		graph(logr, (y)**0.5)
		
		y = vart_tot * 1.0
		y[dens_tot>0] /= (dens_tot * r ** 0)[dens_tot>0] 
		yt = y
		yt /= 2
		print y, vart_tot
		graph(logr, (y)**0.5, color="blue")
		#graph(logr, (varr_tot/(nu))**0.5)
		
		#y = varr_tot_check
		#y[counts>0] /= counts[counts>0]
		#print y**0.5, self.storage_3d.x
		
		select(0,0) 
		graph(logr, (yt)**0.5, color="green")
		graph(logr, (yr)**0.5, color="green")
		 
		#print self.orbitweights.orbitweights.shape, self.storage_3d.moments3d.shape
		m = tensordot(self.orbitweights.orbitweights, self.storage_3d.moments3d, axes=[(0,), (0,)])
		m0 = m[0,:]
		m2r = m[4,:]
		m2r[m0>0] /= m0[m0>0]
		m2t1 = m[5,:]
		m2t1[m0>0] /= m0[m0>0]
		m2t2 = m[6,:]
		m2t2[m0>0] /= m0[m0>0]
		m2t = (m2t1 + m2t2)/2.
		#print m2
		graph(self.storage_3d.x, m2r**0.5, color="blue")
		graph(self.storage_3d.x, m2t**0.5, color="blue")
		#graph(self.storage_3d.x, m2t2**0.5, color="blue", linestyle="dash")
		select(2,0)
		graph(self.storage_3d.x, m0/m0.max(), color="blue")
		graph(self.storage_3d.x, counts/counts.max(), color="red")
		y = dens_tot*r**3
		graph(logr, y/y.max(), color="green")
		beta = 1-yt/yr
		mask = ~isnan(beta)
		graph(logr[mask], beta[mask])
		ylim(-2,1)
		hline(-0.5)
		print dens_tot
		
		r2 = 10**self.storage_3d.x
		y = self.light_model.densityr(r2)*r2**3
		graph(self.storage_3d.x, y/y.max())
		
		select(0,2)
		graph(logr, dens_tot)
		select(1,2)
		graph(logr, nu, color="red")
		select(2,2)
		graph(logr, dens_tot/dens_tot.max())
		graph(logr, nu/nu.max(), color="red")
		
				
				
				
		
		#solution.solve()
	
		
		draw()
		
		
	def run_(self, args):
		def pot(r):
			return self.profile_model.potentialr(r)
		def f(r, E, L):
			a = sqrt(2*E - L**2/r**2 - 2*pot(r))
			return 1/r**2 * (E * r * arctan2(L * r * a , (-L**2 + 2*E*r**2 - 2*r**2*pot(r)) ) \
			+ r * arctan2(r * a,  L) * pot(r) \
			 + 0.5 * L * a)
		def f1(r, E, L):
			a = sqrt(2*E - L**2/r**2 - 2*pot(r))
			return 1/r**2 * 1 / a
		def dens1(r, E, L):
			return f1(r, E, L)
		def dens(r, E1, E2, L11, L12, L21, L22):
			L1 = 0.5 * (L11 + L21)
			L2 = 0.5 * (L12 + L22)
			L1 = (L11)
			L2 = (L12)
			return (f(r, E2, L2) - f(r, E2, L1)) - (f(r, E1, L2) - f(r, E1, L1))
			#return (f(r, E2, L12) - f(r, E2, L11)) - (f(r, E1, L12) - f(r, E1, L11))
			#return (f(r, E2, L22) - f(r, E2, L21)) - (f(r, E1, L12) - f(r, E1, L11))
			#return (f(r, E2, L22) - f(r, E1, L21)) - (f(r, E2, L12) - f(r, E1, L11))
		def dens(r, E1, E2, L11, L12, L21, L22):
			dL1 = L12 - L11
			dL2 = L22 - L21
			L1 = L11
			def t(e):
				return ((2*e-L1**2/r**2-2*pot(r))**(3./2) - (-(dL2**2+2*dL2*L1+L1**2-2*e*r**2+2*r**2*pot(r))/r**2)**(3./2)) / (3 * dL2)
			return (t(E2)-t(E1)) * 2 * pi * (dL1 + (dL2-dL1)) 
		# / r**2
			#dL1 = 0.5 * (L11 + L21)
			#L2 = 0.5 * (L12 + L22)
			#L1 = (L11)
			#L2 = (L12)
			#return (f(r, E2, L2) - f(r, E2, L1)) - (f(r, E1, L2) - f(r, E1, L1))
		print "args", args
		box()
		#vsplit(box)
		self.storage_3d.load()
		if 0:
			I1 = int(args[0]) 
			I2 = int(args[1])
			
			E1 = self.dfgrid.E_borders[I1]
			E2 = self.dfgrid.E_borders[I1+1]
			L11 = self.dfgrid.l_borders[I2] * self.profile_model.Lmax_at_E(E1)
			L12 = self.dfgrid.l_borders[I2+1]*  self.profile_model.Lmax_at_E(E1)
			L21 = self.dfgrid.l_borders[I2] * self.profile_model.Lmax_at_E(E2)
			L22 = self.dfgrid.l_borders[I2+1]*  self.profile_model.Lmax_at_E(E2)
			#print E1, E2, L1, L2
			#logrs = arange(-3, 3, 0.1)
			#logrs = arange(-3, 3, 0.1)
			#rs = 10**logrs
			rs = self.storage_3d.x
			rho_total = rs * 0
			dither = self.dfgrid.dither
			for d1 in range(dither):
				for d2 in range(dither):
					E = self.dfgrid.subgrid.Es[I1*dither+d1]
					L = self.dfgrid.subgrid.ls[I2*dither+d2] * self.profile_model.Lmax_at_E(E)
					#rho = dens(rs, E1, E2, L11, L12, L21, L22)
					rho = dens1(rs, E, L)
					rho[isnan(rho)] = 0
					rho_total += rho * rs**2
			rho = rho_total
			#print 2*E1, L1**2/rs**2, 2*pot(rs)
			#print "->", sqrt(2*E1 - L1**2/rs**2 - 2*pot(rs)) 
			print rho
			#print rs[mask]
			#y = rho[mask]*rs[mask]**2
			y = rho
			y /= sum(y)# * 1.2
			graph(rs, y)
			xlim(min(rs), max(rs))
			ylim(0, max(y))
			orbitnr = self.dfgrid.index_to_orbitnr(I1, I2, 0)
			print orbitnr
			print self.storage_3d.moments3d.shape
			rho = self.storage_3d.moments3d[orbitnr, 0]
			#rho /= sum(rho)
			graph(self.storage_3d.x, rho, color="red", linestyle="dash")
		else:
			logrs = arange(-3, 3, 0.01/4)
			rs = 10**logrs
			#rs = rs = self.storage_3d.x
			#rho_total = rs* 0
			global rho_total
			rho_totals = mab.utils.numpy.mmapzeros((self.dfgrid.n_I1, self.dfgrid.n_I2, len(rs)))
			weights = self.solution.findsolution()
			weights = weights.reshape((self.dfgrid.n_I2, self.dfgrid.n_I1))
			#indexedimage(weights)
			#draw()
			cores = 48
			@mab.parallelize.parallelize(cores=cores, info=info)
			def do(i, j):
					global rho_total
			#for i in range(self.dfgrid.n_I1):
			#	for j in range(self.dfgrid.n_I2):
					if 1:
						E1 = self.dfgrid.E_borders[i]
						E2 = self.dfgrid.E_borders[i+1]
						L11 = self.dfgrid.l_borders[j] * self.profile_model.Lmax_at_E(E1)
						L12 = self.dfgrid.l_borders[j+1]*  self.profile_model.Lmax_at_E(E1)
						L21 = self.dfgrid.l_borders[j] * self.profile_model.Lmax_at_E(E2)
						L22 = self.dfgrid.l_borders[j+1]*  self.profile_model.Lmax_at_E(E2)
						dE = E2-E1
						dL = L22-L21
						
						dither = self.dfgrid.dither
						for d1 in range(dither):
							for d2 in range(dither):
								E1 = self.dfgrid.subgrid.E_borders[i*dither+d1]
								E2 = self.dfgrid.subgrid.E_borders[(i)*dither+d1+1]
								L11 = self.dfgrid.subgrid.l_borders[j*dither+d2] * self.profile_model.Lmax_at_E(E1)
								L12 = self.dfgrid.subgrid.l_borders[(j)*dither+d2+1]*  self.profile_model.Lmax_at_E(E1)
								L21 = self.dfgrid.subgrid.l_borders[j*dither+d2] * self.profile_model.Lmax_at_E(E2)
								L22 = self.dfgrid.subgrid.l_borders[(j)*dither+d2+1]*  self.profile_model.Lmax_at_E(E2)
								dEd = E2-E1
								dLd = L22-L21
								rho = dens(rs, E1, E2, L11, L12, L21, L22)  #* dLd # *( dE * dL)
								#E = self.dfgrid.subgrid.Es[i*dither+d1]
								#L = self.dfgrid.subgrid.ls[j*dither+d2] * self.profile_model.Lmax_at_E(E)
								#rho = dens1(rs, E, L) # * 1 / (dE * dL / (dEd * dLd))
								rho[isnan(rho)] = 0
								rho = rho
								if sum(rho) > 0:
									rho /= sum(rho)
								rho_totals[i,j] += rho  * rs**1 * weights[j,i]
					else:
						E1 = self.dfgrid.E_borders[i]
						E2 = self.dfgrid.E_borders[i+1]
						L11 = self.dfgrid.l_borders[j] * self.profile_model.Lmax_at_E(E1)
						L12 = self.dfgrid.l_borders[j+1]*  self.profile_model.Lmax_at_E(E1)
						L21 = self.dfgrid.l_borders[j] * self.profile_model.Lmax_at_E(E2)
						L22 = self.dfgrid.l_borders[j+1]*  self.profile_model.Lmax_at_E(E2)
						#dE = E2-E1
						#dL = L22-L21
						#E = self.dfgrid.Es[i]
						#L = self.dfgrid.ls[j] * self.profile_model.Lmax_at_E(E)
						#rho = dens1(rs, E, L) * weights[j,i]# * dE * dL
						#d = quad(lambda r: dens(rs, E1, E2, L11, L12, L21, L22), 
						d = dens(rs, E1, E2, L11, L12, L21, L22)
						rho = d * weights[j,i]# *( dE * dL)
						rho[isnan(rho)] = 0
						rho_totals[i,j] += rho * rs**1
						
			nE = self.dfgrid.n_I1
			nL = self.dfgrid.n_I2
			I1s = [i for i in range(nE) for j in range(nL)] #[:1]
			I2s = [j for i in range(nE) for j in range(nL)] #[:1]
			#print len(I1s), len(dfgrid.subgrid.E_borders), len(dfgrid.subgrid.Es)
			do(I1s, I2s)
			rho = sum(sum(rho_totals, axis=0), axis=0)
			print rho.shape
			rho = rho.reshape((len(rho)/8, 8))
			print rho.shape
			rho = sum(rho, axis=1)
			print rho.shape
					
			y = rho
			y /= sum(y)*8
			logdr = logrs[1] - logrs[0] 
			graph(logrs[::8]+logdr*4, y)
			y = self.light_model.densityr(rs)*rs**3/1e4
			y /= sum(y)
			graph(logrs, y, color="red")
		#xlim(min(logrs), max(logrs))
		#select(1)
		
		draw()
		