# -*- coding: utf-8 -*-
from kaplot import *
import mab
import scipy

class TestForeground(object):
	def __init__(self, membership_test, foreground_distribution, object_distribution, N, fraction_foreground=0.5):
		self.membership_test = membership_test
		self.foreground_distribution = foreground_distribution
		self.object_distribution = object_distribution
		self.N = N
		self.fraction_foreground = fraction_foreground
		self.fraction_members = (1-self.fraction_foreground)
		
	def run(self, args, opts, scope):
		M = 100
		statistics = mab.utils.numpy.mmapzeros((M, 7))#, dtype=dtype_initial_conditions)
		
		@mab.parallelize.parallelize(cores=opts.cores, info=scope["progressbar"])
		#for i in range(10):
		def do(i):
			numpy.random.seed(i)
			vs = []
			v_members = []
			for j in range(self.N):
				if random.random() < self.fraction_foreground:
					vs.append(self.foreground_distribution.rvs())
				else:
					v = self.object_distribution.rvs()
					vs.append(v)
					v_members.append(v)
			vs = array(vs)
			v_members = array(v_members)
			mean_sigma, p_members = self.membership_test.calc_p_member(vs, vs * 0, self.fraction_members)
			statistics[i,0] = v_members.var()
			statistics[i,1] = mean_sigma**2
			statistics[i,2] = vs[p_members>0.80].var()
			statistics[i,3] = vs.var()
			vs = vs[p_members>0.90]
			statistics[i,4] = scipy.stats.moment(vs, 4)/scipy.stats.moment(vs, 2)**2
			def moment(x, k):
				m = 0 #x.mean()
				return sum((x-m)**k)/len(x)
			statistics[i,5] = moment(v_members, 4)/moment(v_members, 2)**2
			statistics[i,6] = scipy.stats.moment(v_members, 4)/scipy.stats.moment(v_members, 2)**2
			#print v_members.std()
			#print vs[p_members>0.90].std(),vs.std()
			
		do(range(M))
		print statistics.mean(axis=0)
		print statistics.mean(axis=0)[:3]**0.5
		print statistics.mean(axis=0)[4:]
		print statistics.std(axis=0)
		return
		dsa
		#box()
		mozaic(2,2,box)
		histogram(vs, datamin=-200, datamax=200, binwidth=1., normalize=True)
		x = arange(-200, 200, 0.5)
		y = self.fraction_foreground * self.foreground_distribution.pdf(x) + self.fraction_members * self.object_distribution.pdf(x)
		
		graph(x, y, color="red")
		p_members = self.membership_test.calc_p_member(vs, vs * 0, self.fraction_members)
		#print p_members
		#print p_members>0.90
		
		select(0, 1)
		scatter(vs, p_members, color="red")
		grow(1.1)
		xlim(-200, 200)
		
		draw()