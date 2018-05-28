from kaplot import *
import numpy

class Metropolis(object):
	def __init__(self, dim, f):
		self.dim = dim
		self.f = f
		
		
	def run(self, argv, opts, scope):
		#self.f.load()
		
		parameterset_plot = scope["command_plot_parameter_set"]
		parameterset = self.f
		dim = parameterset.dimension
		w = 7 * dim
		document("%fcm,%fcm" % (w,w))
		
		parameterset_plot.load()
		mozaic(dim, dim, box)
		
		for i in range(dim):
			for j in range(dim):
				if i >= j:
					select(i, dim-1-j)
					parameterset_plot.plot_parameters(i, j)
					
		x = array((0.5,) * dim)
		L = parameterset.logp_uniform(*x)
		#print L
		print "alala"
		Nmax = 1000
		n = 0
		sigma = 0.01
		chain = []
		Naccept = 0
		Nreject = 0
		sigmas = [0.05 for k in range(dim)]
		#sigmas  = [0.01, 0.01, 0.1]
		while n < Nmax:
			#print x, L
			n += 1
			jump = numpy.random.normal(0, sigmas, dim)
			xnew = x + jump
			xnew = maximum(xnew, 0)
			xnew = minimum(xnew, 1)
			Lnew = parameterset.logp_uniform(*xnew)
			#print "\t", xnew, Lnew
			u = exp(Lnew-L)
			#print u
			if u > numpy.random.uniform():
				#print "accept"
				x = xnew
				L = Lnew
				Naccept += 1
				chain.append(x)
			else:
				#print "reject"
				Nreject += 1
		chain = chain[100:]
		#print chain
		print "accepted", Naccept
		print "rejected", Nreject
		print "accepted percentage", Naccept * 100. / (Naccept + Nreject)  
		chain = array(chain)
		select(0, 0)
		scatter(chain[:,1],chain[:,0])
		xlim(0, 1)
		ylim(0, 1)
		
		select(0, 1)
		scatter(chain[:,2],chain[:,0])
		xlim(0, 1)
		ylim(0, 1)
		
		select(1, 0)
		scatter(chain[:,2],chain[:,1])
		xlim(0, 1)
		ylim(0, 1)
		
		draw()
		
		