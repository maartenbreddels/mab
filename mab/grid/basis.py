import scipy.integrate
from numpy import *
import numpy
import scipy.linalg


class FunctionApproximation(object):
	def __init__(self, f, basis):
		self.f = f
		self.basis = basis
		self.cp = [fe*f for fe in basis]
		self.M = array([[fi*fj for fi in basis] for fj in basis])
		for fi in basis:
			for fj in basis:
				assert fi*fj == fj*fi
		d = (self.M - transpose(self.M)).flatten()
		print max(abs(d))
		#self.c = numpy.linalg.solve(self.M, self.cp)
		self.c = scipy.linalg.solve(self.M, self.cp)
		#self.c = dot(Minv, self.cp)
		
	def __call__(self, x):
		return sum([c*fe(x) for c, fe in zip(self.c, self.basis)])
		

class Grid1d(object):
	def __init__(self, basis):
		self.basis = basis
		self.cp = zeros(len(self.basis))
		self.w = zeros(len(self.basis))
		self.M = array([[fj*fi for fi in self.basis] for fj in self.basis])
		self.N = 0
		
	def finish(self):
		self.c = numpy.linalg.solve(self.M, self.cp/(self.w))
		
	def add(self, x, y):
		self.N += len(x)
		for xi, yi in zip(x, y):
			for i, fi in enumerate(self.basis):
				self.cp[i] += fi(xi) * yi
				self.w[i] += fi(xi) * 0 + 0.1
		
	def __call__(self, x):
		return sum([c*fe(x) for c, fe in zip(self.c, self.basis)])
		

class BasisFunction(object):
	def __init__(self, f, x1, x2, x=None, n=1.0, level=0, w=1.):
		self.f = f
		self.x1 = x1
		self.x2 = x2
		self.dx = x2-x1
		self.w = w
		if 0:
			if level == 0:
				self.w = self.dx/ 2
			else:
				self.w = self.dx/2
		if x is None:
			self.x = 0.5 * (x1+x2)
		else:
			self.x = x
		self.n = n
		self.level = level
		
	def overlap(self, x):
		return (self.level <= 1) or ((x >= self.x1) and (x <= self.x2))
		#return (x <= self.x2)
		
	def __call__(self, *args):
		return self.n * self.f(*args)
	
	def __mul__(self, other):
		if isinstance(other,  BasisFunction):
			x1 = max(self.x1, other.x1)
			x2 = min(self.x2, other.x2)
			i, err = scipy.integrate.quad(lambda x: self(x) * other(x), x1, x2)
		else:
			i, err = scipy.integrate.quad(lambda x: self(x) * other(x), self.x1, self.x2)
		#i, err = scipy.integrate.quad(lambda x: self.n * self.f(x) * other(x), 0., 10.)
		return i
	
	
class BasisBlock(object):
	def __init__(self, x1, x2, N):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.N = N
		self.fe = []
		dx = (x2 - x1)/N
		for i in range(N):
			def bf(x, i=i):
				return (x >= dx*i+x1) * (x < dx*(i+1)+x1) #/ dx # * sqrt(2)
			f = BasisFunction(bf, x1+dx*i, x1+dx*(i+1)) 
			self.fe.append(f)
			
	def __len__(self):
		return self.N
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in self.fe:
			p = fi * fi
			fi.n = (1/p)
						
class BasisTri(object):
	def __init__(self, x1, x2, N):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.N = N
		self.fe = []
		dx = (x2 - x1)/(N)
	
		def trileft(x):
			return (x+1) * (x >= -1) * (x < 0)
		def triright(x):
			return (1-x) * (x >= 0) * (x < 1)
		
		for i in range(N+1):
			x1 = self.x1 + dx*i
			x2 = self.x1 + dx*(i+1)
			def bf(x, i=i, x1=x1, x2=x2):
				u = (x-x1)/(x2-x1)
				if i == 0:
					return triright(u)
				elif i == N:
					return trileft(u)
				else:
					return trileft(u) + triright(u)
			if i == 0:
				f = BasisFunction(bf, x1, x2)
			elif i == N: 
				f = BasisFunction(bf, x1-dx, x1)
			else: 
				f = BasisFunction(bf, x1-dx, x2) 
			self.fe.append(f)
			
	def __len__(self):
		return self.N+1
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)						

class BasisTriH(object):
	def __init__(self, x1, x2, levels, scale=True):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.fe = []
	
		def trileft(x):
			return (x+1) * (x >= -1) * (x < 0)
		def triright(x):
			return (1-x) * (x >= 0) * (x < 1)
		
		def bf(x):
			u = (x-self.x1)/(self.x2-self.x1)-1
			return trileft(u)
		f = BasisFunction(bf, self.x1, self.x2, level=-1)
		self.fe.append(f)
		
		def bf(x):
			u = (x-self.x1)/(self.x2-self.x1)
			return triright(u)
		f = BasisFunction(bf, self.x1, self.x2, level=-1) 
		self.fe.append(f)
		
		for level in range(levels):
			count = 2**level
			#print "level", level, count
			dx = (self.x2 - self.x1)/(count*2)
			for i in range(count):
				x1 = self.x1 + dx*(i*2+1)#*2
				x2 = self.x1 + dx*(i*2+2)#*2
				#print x1, x2, (i*2+1), (i*2+2)
				scaling = 1 if not scale else 1 - level*0.15
				def bf(x, i=i, x1=x1, x2=x2, level=level, scaling=scaling):
					u = (x-x1)/(x2-x1)
					return (trileft(u) + triright(u))/(level+1)
				f = BasisFunction(bf, x1-dx, x2, level=level) 
				self.fe.append(f)
			
	def __len__(self):
		return len(self.fe)
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)
			
class BasisTriLevels(object):
	def __init__(self, x1, x2, levels):
		self.basis_level = [BasisTriHLevel(x1, x2, k) for k in range(levels+1)]
		
	def __len__(self):
		return len(self.basis_level)
	
	def __getitem__(self, i):
		return self.basis_level[i]
	

class BasisTriHLevel(object):
	def __init__(self, x1, x2, level, scale=True):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.fe = []
	
		def trileft(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x > -1) * (x <= 0)
		def triright(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x > 0) * (x <= 1)
		def trilefti(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x >= -1) * (x <= 0)
		def tririghti(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x >= 0) * (x <= 1)
		#def trilefti(x): # x (-1, 0) y -> (0, 1) /0
		#	return (x+1) * (x > -1) * (x <= 0)
		#def tririghti(x):# x (0, 1) y -> (1, 0)   0\
		
		#	return (1-x) * (x >= 0) * (x < 1)
		
		if level == 0:
			self.overlap = True
			def bf(x):
				u = (x-self.x1)/(self.x2-self.x1)
				return tririghti(u)
			f = BasisFunction(bf, self.x1, self.x2, x1, level=0, w=0.5) 
			self.fe.append(f)
			
			def bf(x):
				u = (x-self.x1)/(self.x2-self.x1)-1
				return trilefti(u)
			f = BasisFunction(bf, self.x1, self.x2, x2, level=0, w=0.5)
			self.fe.append(f)
			
			
		else:
			self.overlap = False
			count = 2**(level-1)
			#print "level", level, count
			#dx = (self.x2 - self.x1)/(count*2)
			dx = (self.x2 - self.x1)/count
			for i in range(count):
				x1 = self.x1 + dx*i
				x2 = self.x1 + dx*(i+1)
				#print x1, x2, (i*2+1), (i*2+2)
				scaling = 1 if not scale else 1 - level*0.15
				def bf(x, i=i, x1=x1, x2=x2, level=level, scaling=scaling):
					u = (x-x1-dx/2)/(dx/2)
					return (trileft(u) + triright(u))#s/(level+1)
				#print ">>", level, count, x1, x2, dx
				f = BasisFunction(bf, x1, x2, x1+dx/2, level=level, w=dx/2) 
				self.fe.append(f)
			
	def __len__(self):
		return len(self.fe)
	
	def indices(self, x):
		count = len(self)
		if self.overlap:
			return range(count)
		else:
			normalized = (x - self.x1)/(self.x2-self.x1)
			if (normalized < 0) or (normalized >= 1):
				return []
			else:
				return int(normalized * count)
			#dx = (self.x2-self.x1)/count
			#index = int( 
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)
			
class BasisTriModLevels(object):
	def __init__(self, x1, x2, levels):
		self.basis_level = [BasisTriModHLevel(x1, x2, k) for k in range(levels+1)]
		
	def __len__(self):
		return len(self.basis_level)
	
	def __getitem__(self, i):
		return self.basis_level[i]
	

class BasisTriModHLevel(object):
	def __init__(self, x1, x2, level, scale=True):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.fe = []
	
		def trileft(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x > -1) * (x <= 0)
		def triright(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x > 0) * (x <= 1)
		def trilefti(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x >= -1) * (x <= 0)
		def tririghti(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x >= 0) * (x <= 1)
		#def trilefti(x): # x (-1, 0) y -> (0, 1) /0
		#	return (x+1) * (x > -1) * (x <= 0)
		#def tririghti(x):# x (0, 1) y -> (1, 0)   0\
		
		#	return (1-x) * (x >= 0) * (x < 1)
		
		if level == 0:
			def cf(x):
				return x * 0 + 1
			f = BasisFunction(cf, self.x1, self.x2, (x1+x2)/2, level=0, w=1.0)
			self.fe.append(f)
		else:
			count = 2**(level)
			dx = (self.x2 - self.x1)/count
			for i in range(count):
				x1 = self.x1 + dx*i
				x2 = self.x1 + dx*(i+1)
				#print x1, x2, (i*2+1), (i*2+2)
				if i == 0:
					def left(x, x1=x1, x2=x2, i=i, level=level):
						return (2 - 2**level * x - 1) * (x >= 0) * (x < dx) 
					f = BasisFunction(left, x1, x2, x1+dx/2, level=level, w=dx/2) 
					self.fe.append(f)
				elif i == count-1:
					def right(x, x1=x1, x2=x2, i=i, level=level):
						#return (2**level * x - i - 1) * (x >= x1) * (x <= x2) 
						return (x/dx - i ) * (x >= x1) * (x <= x2)
					f = BasisFunction(right, x1, x2, x1+dx/2, level=level, w=dx/2) 
					self.fe.append(f)
				else:
					def fb(x, x1=x1, x2=x2, i=i, level=level):
						u = x/(dx) - i - 0.5
						#return (0 * x + 0.1)  * (x >= x1) * (x <= x2) 
						return (0.5 - abs(u)) * (x >= x1) * (x <= x2) 
					f = BasisFunction(fb, x1, x2, x1+dx/2, level=level, w=dx/4) 
					self.fe.append(f)
			
	def __len__(self):
		return len(self.fe)
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)
			

class BasisTriCCLevels(object):
	def __init__(self, x1, x2, levels):
		self.basis_level = [BasisTriCCHLevel(x1, x2, k) for k in range(levels+1)]
		
	def __len__(self):
		return len(self.basis_level)
	
	def __getitem__(self, i):
		return self.basis_level[i]
	

class BasisTriCCHLevel(object):
	def __init__(self, x1, x2, level, scale=True):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.fe = []
	
		def trileft(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x > -1) * (x <= 0)
		def triright(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x > 0) * (x < 1)
		def trilefti(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x >= -1) * (x <= 0)
		def tririghti(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x >= 0) * (x <= 1)
		#def trilefti(x): # x (-1, 0) y -> (0, 1) /0
		#	return (x+1) * (x > -1) * (x <= 0)
		#def tririghti(x):# x (0, 1) y -> (1, 0)   0\
		
		#	return (1-x) * (x >= 0) * (x < 1)
		
		def hat_center(x):
			return (1-abs((x-0.5)*2)) * (x >= 0) * (x <= 1)
		if level == 0:
			def bf(x):
				u = (x-self.x1)/(self.x2-self.x1)
				return hat_center(u)
			f = BasisFunction(bf, self.x1, self.x2, (x1+x2)/2, level=0, w=0.5) 
			self.fe.append(f)
		else:
			count = 2**(level+1)
			dx = (self.x2 - self.x1)/count
			n = 2**level+1
			#theta = pi * arange(n)/(n-1.);
			print "level", level, count, dx
			xs = arange(n)/(n-1.) #(cos(theta)[1::2] + 1) / 2
			print xs
			if level == 1:
				xs = xs[::2]
				w = 0.25
			else: #elif level == 2:
				xs = xs[1::2]
				w = 0.5**(len(xs)*2)
			print xs
			if level == 2:
				#dsa
				pass
			#dx = 1./n
			dx = (xs[1] - xs[0])#/2
			count = len(xs)
			if level == 1:
				w = 0.25 #dx/2
			else:
				w = 0.5**level
				w = dx/2
			print "w/dx", w, dx
			print xs, count 
			for i in range(len(xs)):
				#x1 = self.x1 - dx/2+ dx*i
				#x2 = self.x1 + dx/2 + dx*i
				print i
				x1 = (xs[i]-dx/2) #self.x1 - dx/2+ dx*i
				x2 = (xs[i]+dx/2) #self.x1 + dx/2 + dx*i
				#print x1, x2, (i*2+1), (i*2+2)
				scaling = 1 if not scale else 1 - level*0.15
				def bf(x, i=i, x1=x1, x2=x2, level=level, scaling=scaling):
					u = (x-x1-dx/2)/(dx/2)
					return (trileft(u) + triright(u))#s/(level+1)
				#print ">>", level, count, x1, x2, dx
				f = BasisFunction(bf, x1, x2, xs[i], level=level, w=w) 
				self.fe.append(f)
			
	def __len__(self):
		return len(self.fe)
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)
			
			
class BasisTriSwapLevels(object):
	def __init__(self, x1, x2, levels):
		self.basis_level = [BasisTriSwapHLevel(x1, x2, k) for k in range(levels+1)]
		
	def __len__(self):
		return len(self.basis_level)
	
	def __getitem__(self, i):
		return self.basis_level[i]
	

class BasisTriSwapHLevel(object):
	def __init__(self, x1, x2, level, scale=True):
		self.x1 = float(x1)
		self.x2 = float(x2)
		self.fe = []
	
		def trileft(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x > -1) * (x <= 0)
		def triright(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x > 0) * (x <= 1)
		def trilefti(x): # x (-1, 0) y -> (0, 1) /0
			return (x+1) * (x >= -1) * (x <= 0)
		def tririghti(x):# x (0, 1) y -> (1, 0)   0\
			return (1-x) * (x >= 0) * (x <= 1)
		#def trilefti(x): # x (-1, 0) y -> (0, 1) /0
		#	return (x+1) * (x > -1) * (x <= 0)
		#def tririghti(x):# x (0, 1) y -> (1, 0)   0\
		
		#	return (1-x) * (x >= 0) * (x < 1)
		
		if level == 1:
			def bf(x):
				u = (x-self.x1)/(self.x2-self.x1)
				return tririghti(u)
			f = BasisFunction(bf, self.x1, self.x2, x1, level=1) 
			self.fe.append(f)
			
			def bf(x):
				u = (x-self.x1)/(self.x2-self.x1)-1
				return trilefti(u)
			f = BasisFunction(bf, self.x1, self.x2, x2, level=1)
			self.fe.append(f)
			
			
		else:
			if level == 0:
				count = 1
			else:
				count = 2**(level-1)
			#print "level", level, count
			#dx = (self.x2 - self.x1)/(count*2)
			dx = (self.x2 - self.x1)/count
			for i in range(count):
				x1 = self.x1 + dx*i
				x2 = self.x1 + dx*(i+1)
				#print x1, x2, (i*2+1), (i*2+2)
				scaling = 1 if not scale else 1 - level*0.15
				def bf(x, i=i, x1=x1, x2=x2, level=level, scaling=scaling):
					u = (x-x1-dx/2)/(dx/2)
					return (trileft(u) + triright(u))#s/(level+1)
				#print ">>", level, count, x1, x2, dx
				f = BasisFunction(bf, x1, x2, x1+dx/2, level=level) 
				self.fe.append(f)
			
	def __len__(self):
		return len(self.fe)
	
	def __getitem__(self, i):
		return self.fe[i]
	
	def normalize(self):
		for fi in basis:
			p = fi * fi
			fi.n = (1/p)
			
			