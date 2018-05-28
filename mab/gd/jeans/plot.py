from kaplot import *
import scipy.optimize
import numpy.linalg

class Test(object):
	def __init__(self, jeans):
		pass
	
	def run(self, args, opts, scope):
		mozaic(3,2,box)
		modelpath = scope["modelpath"]
		#names = ["jeans.beta", "dm_density_twoslope.alpha", "dm_density_twoslope.beta", "dm_density_twoslope.M1kpc", "dm_density_twoslope.rs"]
		names = ["dm_density_twoslope.alpha", "jeans.beta", "dm_density_twoslope.beta", "dm_density_twoslope.M1kpc", "dm_density_twoslope.rs"]
		
		Rs = load(os.path.join(modelpath, "R.npy"))
		deviations = zeros((len(Rs), len(names))) 
		for i, name in enumerate(names):
			filename = os.path.join(modelpath, name+".npy")
			x = load(filename)
			deviations[:,i] = x
			graph(Rs, x)
		select(1,0)
		if 0:
			n = len(Rs)
			print "n=", n
			chisqmap = zeros((n,n))
			for i in range(n):
				for j in range(n):
					#print i, j
					x2 = zeros(len(names))
					i, j = 1,-1
					def f(x):
						x2[0] = 1
						x2[1:] = x
						#x2[1] = -0.02
						return sum(x2*deviations[i] + x2*deviations[j])**2
						#print i
						#return sum((x2*deviations[i]))**2# + x2*deviations[j])**2)
					x0 = ones(len(names)-1)
					#print deviations[i], deviations[i]**2
					ret = scipy.optimize.fmin(f, x0, full_output=1, disp=1,xtol=0.000000000001, ftol=0.00000000001,)
					chisq = ret[1]
					chisqmap[i,j] = chisq
					print chisq,
					print ret
					print sum(x2*deviations[i]), sum(x2*deviations[j])
					sys.exit(0)
					#dot(ret, 
					#print ret
					#dsa
			#sekect
		indexedimage(deviations, colormap="whiteblack")
		u,s,v = numpy.linalg.svd(deviations, full_matrices = 1)
		print u.shape, s.shape, v.shape
		print s
		select(1,1)
		for i, color in zip(range(5), "black red green blue orange".split()):
			s = [0, 0, 0, 0, 0]
			s[i] = 1
			#S = numpy.diag(s)
			#import pdb
			#pdb.set_trace()
			#print S
			#print dot(u, dot(S, v))
			nv = dot(s, v)
			#print nv.shape
			dev = dot(nv, deviations.T)
			print dev
			print nv
			#print dot(u, s)
			graph(dev, color=color)
		for i, color in zip(range(5), "black red green blue orange".split()):
			dev = u[i]
			#graph(dev, color=color, linestyle="dot")
			dev = u[:,i]
			#graph(dev, color=color, linestyle="dash")
			
		select(0, 1)
		indexedimage(abs(u), colormap="whiteblack")
		select(2,1)
		indexedimage(abs(v), colormap="whiteblack")
		
				 
		draw()


class PlotJeansDerivatives(object):
	def __init__(self, jeans_derivatives):
		self.jeans_derivatives = jeans_derivatives
		
	def run(self, args, opts, scope):
		self.jeans_derivatives.load(scope)
		R = self.jeans_derivatives.R
		sigmaR = self.jeans_derivatives.jeans.sigma_los(R)
		mozaic(2,2,box)
		graph(R, sigmaR)
		ylim(0, sigmaR.max()*1.1)
		n = len(self.jeans_derivatives.parameters)
		select(0, 1)
		for i in range(n):
			graph(R, self.jeans_derivatives.deriv[:,i], color=nicecolors[i])
			scatter(R, self.jeans_derivatives.deriv[:,i], color=nicecolors[i], symbolsize="15pt")
			print self.jeans_derivatives.deriv[:,i]
		hline(0)
		select(1, 1)
		linestyles = "normal dash dot dashdot".split()
		for i in range(n):
			for j in range(n):
				#if i == j:
				if 0:
					graph(R, self.jeans_derivatives.hessian[:,i,j], color=nicecolors[i], linestyle=linestyles[j])
					print i, j, self.jeans_derivatives.hessian[:,i,j]
		hline(0)
		
		#H = self.jeans_derivatives.hessian
		D = self.jeans_derivatives.deriv
		#D[:,0] /= 5
		y = sigmaR
		N = len(R)
		def mat(ri):
			m = zeros((n,n))
			for i in range(n):
				for j in range(n):
					m[i,j] = D[ri,i] * D[ri,j]
			return m
					# + (y[ri]-
					
		matrices = [mat(k) for k in range(N)]
		print D
		#print D[4]
		print "-" * 10
		print matrices[0]
		print matrices[5]
		m = matrices[5]
		#print numpy.matrix(m).I
		U, S, Vh = numpy.linalg.svd(m)
		S[S<1e-5] = 0
		Si = S
		Si[S>1e-5] = 1/Si[S>1e-5]
		print Si
		Sd = numpy.diag(S)
		Si = numpy.diag(Si)
		print "svd",
		print U
		print Sd
		print Vh
		print "inv"
		print U
		print Si
		print Vh
		print "-" * 10
		M = Vh.T * Si * U.T
		print M
		#das
		info = 1.
		def fopt(w):
			w = abs(w)
			#w *= 1e-6
			#w[4] = 1000.
			weights = w# / sum(w)# * 10
			mat = sum([m*w*info for m,w in zip(matrices, weights)], axis=0)
			det = numpy.linalg.det(mat)
			#print mat
			mat = numpy.matrix(mat).I
			mat = array(mat)
			sigmas = [sqrt(mat[i,i]) for i in range(n)]
			for i in range(n):
				for j in range(n):
					mat[i,j] /= sigmas[i]*sigmas[j]
			
			#print mat
			#d = numpy.linalg.det(mat)
			rho1 = mat[0,1]
			#rho2 = mat[0,2]
			#rho3 = mat[1,2]
			#print 1/det, rho1, rho2, rho3, sigmas
			print 1/det, rho1, sigmas
			return (0.66-sigmas[0])**2 + (1-sigmas[1])**2 + abs(rho1*10)
			return abs(rho1*100) + 1/det
			return 1/det
			#return abs(rho1) + abs(rho2)
			#return sigmas[0] #abs(1/det)
			#return abs(rho3)
			#return mat[2,2] + mat[1,1]
			#return abs(rho1) + abs(rho2)-abs(det)# + abs(rho3)
			#return abs(det) + abs(rho2)
			#return abs(rho2)*sigmas[1]*sigmas[2]# + sigmas[0] 
			#return abs(det/100) + abs(rho1*10) + abs(rho2*20)# + sigmas[0] + sigmas[1]
			#return -sqrt(d)
			#u1, u2 = mat[0]
			#length = sqrt(u1**2+u2**2)
			lengths = [sqrt(sum([k**2 for k in mat[i]])) for i in range(n)]
			us = [mat[i,i]/lengths[i] for i in range(n)]
			#return -mat[
			return -us[1]
			u1 /= length
			u2 /= length
			m1 = u1 
			u1, u2 = mat[1]
			length = sqrt(u1**2+u2**2)
			u1 /= length
			u2 /= length
			m2 = u2
			#print m1, m2
			#print mat
			#return -abs(mat[0,0])
			#return -m2
			return -(m2-m1*0.5)
			#print "d[", d, "]"
			#return -d
			#return mat[0,1] - mat[0,0]
			return -(mat[0,1]*10 + mat[0,0])
		
		w = ones(len(y)) * 1.0
		w = numpy.random.random(len(w))
		if 0:
			w *= 0
			w[0] = 1.
			w[3] = 1.
			w[4] = 1.
		#w = 1-abs(D[:,1])**2
		
		import scipy.optimize
		w = scipy.optimize.fmin_powell(fopt, w, maxfun=100000, maxiter=10000)
		w = abs(w)
		w += 1e-10
		if 0:
			w *= 0
			#w[0] = 0.5
			#w[3] = 0.5
			w[4] = 0.5
			#w[0] = 0.5
			w += 1e-10
		weights = w# / sum(w)
		select(0,1)
		graph(R, weights*4, color="purple")
		scatter(R, weights*4, color="purple", symbolsize="25pt")
		#ylim(0, 1)
		#print weights
		for w in weights:
			print "%.3f" % w,
		weights = weights / sum(weights)
		#sigmas = 
		print
		mat = sum([m*w for m,w in zip(matrices, weights)], axis=0)
		print "Matrix:\n", mat
		mat = numpy.matrix(mat)
		if 0:
			print "SVD",
			def pseudo_inverse(M, eps=1e-9):
				U, S, Vh = numpy.linalg.svd(M)
				Smax = abs(S[0])
				mask = abs(S/Smax) >= eps
				Sinv = zeros_like(S)
				print S,
				Sinv[mask] = 1/S[mask]
				print Sinv
				Sinv = numpy.diag(Sinv)
				#import pdb; pdb.set_trace()
				return Vh.T * Sinv * U.T
			U, S, Vh = numpy.linalg.svd(mat)
			print "U:\n",U
			print "S:",S
			print "V:\n",Vh
			mati = mat.I
			print "inverse"
			print mati
			print "pseudo inverse"
			print pseudo_inverse(mat)
			print "dsa"
			print matrices[4]
			print pseudo_inverse(matrices[4])
			print mat
			print pseudo_inverse(mat)
			print "pk"
			print numpy.matrix(matrices[4]).I
			print pseudo_inverse(matrices[4])
			print pseudo_inverse(matrices[5])
			print "check"
			print mat.I.I
			mati = pseudo_inverse(mat)
		mati = mat.I
		print "Inverse:\n",mati
		sigmas = [sqrt(mati[i,i]) for i in range(n)]
		for i in range(n):
			for j in range(n):
				mati[i,j] /= sigmas[i]*sigmas[j]
		print mati
		
					
					
		#dsa
		
		draw()
	 
class JeansDerivatives(object):
	def __init__(self, jeans, parameters, rmin=0.01, rmax=1.5, nr=11):
		self.jeans = jeans
		self.parameters = parameters
		
		u = (arange(nr) +0.5) / (nr+1) * (rmax - rmin) + rmin
		#self.logr = u
		self.R = u
		
	def load(self, scope):
		modelpath = scope["modelpath"]
		n = len(self.parameters)
		self.hessian = zeros((len(self.R), n, n))
		self.deriv = zeros((len(self.R), n))
		print "loading"
		
		for i, param1 in enumerate(self.parameters):
			print param1
			filename = os.path.join(modelpath, "d_"+param1 +".npy")
			print "load", filename
			self.deriv[:,i] = load(filename)
			if 0:
				for j, param2 in enumerate(self.parameters):
					if i <= j:
						print param2
						filename = os.path.join(modelpath, "dd_"+param1+"_"+param2+".npy")
						print "load", filename
						self.hessian[:,i,j] = load(filename)
					else: # i > j, but H[i,j] == H[j,i]
						self.hessian[:,i,j] = self.hessian[:,j,i]
	
	def run(self, args, opts, scope):
		modelpath = scope["modelpath"]
		filename = os.path.join(modelpath, "R.npy")
		print filename
		save(filename, self.R)
		
		scope.reset()
		scope.re_readfiles()
		scope.init()
		cleanscope = scope.clone()
		
		#self.sigmaR0 = self.jeans.sigma_los(self.R)
		for i, param1 in enumerate(self.parameters):
			print param1
			filename = os.path.join(modelpath, "d_"+param1 +".npy")
			print "save to", filename
			x = self.dparam(cleanscope.clone(), param1, param1, 1e-3)
			#print x
			save(filename, x)
			if 0:
				for j, param2 in enumerate(self.parameters):
					if i <= j:
						print param2
						x = self.ddparam(cleanscope.clone(), param1, param2, 1e-3)
						filename = os.path.join(modelpath, "dd_"+param1+"_"+param2+".npy")
						print "save to", filename
						save(filename, x)
			
	def ddparam(self, scope, paramname1, paramname2, delta):
		#orig = scope.clone()
		#x0 = scope[paramname1]
		d1 = self.dparam(scope.clone(), paramname2, paramname1, delta)
		#scope.reset()
		#scope.re_readfiles()
		#scope.init()
		logarithmic = "M1kpc" in paramname1 or "rs" in paramname1
		logarithmic = True
		if logarithmic:
			scope[paramname1] = scope[paramname1]*(1 + delta)
		else:
			scope[paramname1] = scope[paramname1] + delta
		check = scope[paramname1]
		d2 = self.dparam(scope.clone(), paramname2, paramname1, delta)
		assert scope[paramname1] == check
		if logarithmic:
			print "log"
			ddparam = (d2-d1)/(delta)
		else:
			ddparam = (d2-d1)/delta
		print ddparam
		#scope.reset()
		#scope.re_readfiles()
		#scope.init()
		#scope[paramname1] = x0
			
		return ddparam
		
	def dparam(self, scope, paramname, paramname2, delta, logarithmic=False, **kwargs):
		orig = scope.clone()
		#print orig.dict.keys()
		jeans = scope["jeans"]
		#print orig.dict.keys()
		sigmaR0 = jeans.sigma_los(self.R)
		scope = orig
		#x0 = scope[paramname]
		#y0 = scope[paramname2]
		
		#scope.reset()
		#scope.re_readfiles()
		#scope.init()
		logarithmic = "M1kpc" in paramname or "rs" in paramname
		#logarithmic = True
		if logarithmic:
			scope[paramname] = scope[paramname]*(1+delta)
		else:
			scope[paramname] = scope[paramname]+delta
		jeans = scope["jeans"]
		
		sigmaR = jeans.sigma_los(self.R)
		dparam = (sigmaR-sigmaR0)/delta
		print "d", dparam
		
		#scope.reset()
		#scope.re_readfiles()
		#scope.init()
		#scope[paramname] = x0
		#scope[paramname2] = y0
		return dparam
		#graph(self.rlinear, dparam, addlegend=False, **kwargs)
		#filename = os.path.join(scope["modelpath"], paramname+".npy")
		#print filename
		#save(filename, dparam)
		
		#scope.reset()
		#scope.re_readfiles()
		

	
class Profile(object):
	def __init__(self, jeans, parameterset_iter, logrmin=-2, logrmax=2, nr=25, rmin=0.01, rmax=1.5):
		self.jeans = jeans
		self.parameterset_iter = parameterset_iter
		self.logrmin = logrmin
		self.logrmax = logrmax
		u = (arange(nr) +0.5) / (nr+1) * (logrmax - logrmin) + logrmin
		self.logr = u
		self.r = 10**u
		u = (arange(nr+1)) / (nr+0.) * (logrmax - logrmin) + logrmin
		self.logr_borders = 10**u
		self.r_borders = 10**u
		
		u = (arange(nr) +0.5) / (nr+1) * (rmax - rmin) + rmin
		#self.logr = u
		self.rlinear = u
		
		
	def run(self, args, opts, scope):
		#box()
		filename = os.path.join(scope["modelpath"], "R.npy")
		print filename
		save(filename, self.rlinear)
		 
		mozaic(3,2,box)
		
		select(0,0)
		sigmar = self.jeans.sigmar(self.r)
		graph(self.logr, sigmar)
		labels("log r/kpc", "&sigma;<sub>r</sub>")
		ylim(0, sigmar.max() * 1.1)
		#draw()
		
		select(1,0)
		sigmar = self.jeans.sigmar(self.rlinear)
		graph(self.rlinear, sigmar)
		labels("r/kpc", "&sigma;<sub>r</sub>")
		ylim(0, sigmar.max() * 1.1)
		
		select(2,0)
		sigmaR = self.jeans.sigma_los(self.rlinear)
		graph(self.rlinear, sigmaR)
		labels("r/kpc", "&sigma;<sub>R</sub>")
		ylim(0, sigmaR.max() * 1.1)
		draw()
		
		
		light_scale = scope["plummer"].b
		select(0, 1)
		labels("log r/kpc", "&Delta;")
		hline(0)
		select(1,1)
		labels("r/kpc", "&Delta;")
		hline(0)
		select(2,1)
		labels("r/kpc", "&Delta;")
		hline(0)
		select(2,1)
		
		select(2,0)
		vline(light_scale)
		select(2,1)
		vline(light_scale)
		
		clearautolegend()
		self.dparam(scope, "jeans.beta", 1e-4)
		self.dparam(scope, "dm_density_twoslope.alpha", 1e-4, color="red")
		self.dparam(scope, "dm_density_twoslope.beta", 1e-4, color="blue")
		self.dparam(scope, "dm_density_twoslope.M1kpc", 1, logarithmic=True, color="orange", linestyle="dash")
		self.dparam(scope, "dm_density_twoslope.rs", 1, color="green", linestyle="dot")
		autolegend("ani", "alpha", "beta", "log M1kpc", "rs")
		select(1,0)
		xlim(0, self.rlinear.max())
		select(1,1)
		xlim(0, self.rlinear.max())
		select(2,0)
		xlim(0, self.rlinear.max())
		select(2,1)
		xlim(0, self.rlinear.max())
		#draw()
		if 0:
			r0 = 0.1
			r1 = 1.0
			sigmar0 = self.jeans.sigmar(r0)
			sigmar1 = self.jeans.sigmar(r1)
			
			dbeta = 1e-4
			self.jeans.beta0 += dbeta
			sigmar = self.jeans.sigmar(r0)
			self.jeans.beta0 -= dbeta
			dsigma = sigmar-sigmar0
			a = dsigma/dbeta
			print a
			
			dalpha = 1e-3
			dm_density = scope["dm_density_twoslope"]
			dm_profile = scope["dm_profile"]
			dm_density.alpha += dalpha
			dm_profile.update()
			sigmar = self.jeans.sigmar(r0)
			dm_density.alpha -= dalpha
			dm_profile.update()
			dsigma = sigmar - sigmar0
			b = dsigma/dalpha
			print b
			#a*dbeta + b*dalpha = 0
			dbetatest = -b*dalpha/a
			print "dbeta", dbetatest
			print "0 = ", a*dbetatest + b*dalpha
			
			dm_density.alpha += dalpha
			self.jeans.beta0 += dbetatest
			dm_profile.update()
			sigmar = self.jeans.sigmar(r0)
			print sigmar-sigmar0
			sigmar = self.jeans.sigmar(r1)
			print sigmar-sigmar1
			
			sys.exit(0)
			
			
			print self.parameterset_iter
			for scope in self.parameterset_iter.iter(scope):
				print scope
				jeans = scope["jeans"]
				print jeans.beta(1), scope["dm_density_twoslope.alpha"]
	
				sigmar = jeans.sigmar(self.r)
				graph(self.logr, sigmar, color="red")
				#labels("log r/kpc", "&sigma;<sub>r</sub>")
				ylim(0, sigmar.max() * 1.1)
		draw()
		
	def dparam(self, scope, paramname, delta, logarithmic=False, **kwargs):
		sigmar0 = self.jeans.sigmar(self.r)
		sigmar0linear = self.jeans.sigmar(self.rlinear)
		sigmaR0 = self.jeans.sigma_los(self.rlinear)
		x0 = scope[paramname]
		
		scope.reset()
		scope.re_readfiles()
		scope.init()
		#for name, value in parametervalue.items:
		if logarithmic:
			scope[paramname] = x0+delta
		else:
			scope[paramname] = x0+delta
		jeans = scope["jeans"]
		
		sigmar1 = jeans.sigmar(self.r)
		sigmar1linear = jeans.sigmar(self.rlinear)
		
		select(0, 1)
		dparam = (sigmar1-sigmar0)/delta
		if logarithmic:
			dparam = dparam * x0
		#print dparam
		graph(self.logr, dparam, **kwargs)
		select(1, 1)
		
		dparam = (sigmar1linear-sigmar0linear)/delta
		if logarithmic:
			dparam = dparam * x0
		graph(self.rlinear, dparam, addlegend=False, **kwargs)
		
		select(2, 1)
		
		sigmaR = jeans.sigma_los(self.rlinear)
		dparam = (sigmaR-sigmaR0)/delta
		print dparam
		if logarithmic:
			dparam = dparam * x0
		graph(self.rlinear, dparam, addlegend=False, **kwargs)
		filename = os.path.join(scope["modelpath"], paramname+".npy")
		print filename
		save(filename, dparam)
		
		#scope.reset()
		#scope.re_readfiles()
		
		