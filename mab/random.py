from numpy import random, arccos, pi, matrix, cos, sin, sqrt, array, dot, transpose, ravel, linalg, identity
import numpy, math

def randomdir(N=1):
	"""generate N unrandom phi and theta angles"""
	#r = p1.sampler(rmax)
	phi = random.random(size=N) * 2 * pi
	#theta = arccos(random.random()#
	costheta = random.random(size=N) * 2 - 1
	#sintheta = sqrt(1-costheta**2)
	return phi, arccos(costheta)

def random_spherical(R, N=10000, R0=0):
	"""generate N uniform random points between R0 and R0+R"""
	u1 = numpy.random.random(size=N)
	r = u1 ** (1./3.) * R + R0
	u2 = numpy.random.random(size=N) * 2 -1
	phi = numpy.random.random(size=N) * 2 * math.pi
	x = numpy.sqrt(1-u2**2) * numpy.cos(phi) * r
	y = numpy.sqrt(1-u2**2) * numpy.sin(phi) * r
	z = u2 * r
	return x, y, z

def randomSO3():
	"""return random rotatation matrix, algo by James Arvo"""
	u1 = random.random()
	u2 = random.random()
	u3 = random.random()
	R = array([[cos(2*pi*u1), sin(2*pi*u1), 0], [-sin(2*pi*u1), cos(2*pi*u1), 0], [0, 0, 1]])
	v = array([cos(2*pi*u2)*sqrt(u3), sin(2*pi*u2)*sqrt(u3), sqrt(1-u3)])
	H = identity(3)-2*v*transpose([v])
	#print "v", v
	#print "vvT", v*transpose([v])
	#print "H", H
	#print linalg.det(R), linalg.det(H)
	#print H, v * transpose([v])
	return - dot(H,  R)

if __name__ == "__main__":
	from numpy import *
	print "x,y,z"
	for i in range(1):
		M = randomSO3()
		x = ravel(dot(M, [1, 0, 0]))
		y = ravel(dot(M, [0, 1, 0]))
		z = ravel(dot(M, [0, 0, 1]))
		#print v
		#print v[0]
		print "-" * 70
		print x[0], x[1], x[2], sqrt(sum(x**2))
		print y[0], y[1], y[2]
		print z[0], z[1], z[2]
		print cross(x, y)
		print dot(x, y)
		print dot(x, z)
		print dot(y, z)
		#print 
		
		