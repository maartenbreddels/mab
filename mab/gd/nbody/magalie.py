from numpy import *
import os
from mab.gd.nbody.snapshot import Snapshot, Component

class Run(object):
	def __init__(self, inputfile, cwd=None, binary="magalie.exe"):
		self.inputfile = inputfile
		self.binary = binary
		self.cwd = cwd
		
	def run(self, args, opts, scope):
		oldcwd = os.getcwd()
		try:
			cmd = "%s < %s" % (self.binary, self.inputfile)
			if self.cwd:
				print self.cwd
				os.chdir(self.cwd)
			print cmd
			os.system(cmd)
		finally:
			os.chdir(oldcwd)
			
class IC(Snapshot):
	def __init__(self, filename, configuration, mass, scale):
		super(IC, self).__init__()
		self.filename = filename
		self.configuration = configuration
		self.massscale = mass
		self.scale = scale
		
	def load(self):
		#self.components = {}
		self.masses = {}
		lines = open(self.filename).readlines()
		grid = []
		for line in lines[1:]:
			#print line
			values = [float(k) for k in line.strip().split()]
			grid.append(values)
		grid = array(grid)
		print grid.shape
		
		# for scaling, see hernquist 1993
		Mdisk = 5.6e10
		rdisk = 3.5
		vscale = 262.333498 * (self.massscale/Mdisk * rdisk/self.scale)**0.5 ##* 200#**0.5
		if 0:
			self.massscale = 1.
			vscale = 1.
			self.scale = 1.
			
		#vscale = 220.0 * (self.mass/Mdisk * rdisk/self.scale)**0.5
		# format is : mass, x, y, z, vx, vy, vz
		self.p = grid[:,4:].T * vscale
		self.mass = grid[:,0] * self.massscale
		self.q = grid[:,1:4].T * self.scale
		x, y, z = self.q
		r = sqrt(x**2+y**2+z**2)
		mask = r < 8.
		
		print self.mass.shape, self.q.shape, self.p.shape
		print "total mass", sum(self.mass), sum(self.mass)/Mdisk
		i = 0
		self.add("disk", i, self.configuration.disk_N)
		i += self.configuration.disk_N
		self.add("bulge", i, i+self.configuration.bulge_N)
		i += self.configuration.bulge_N
		self.add("halo", i, i+self.configuration.halo_N)
		i += self.configuration.halo_N
		print "enc mass", "%g" % sum(self.mass[mask])
		
	def add(self, name, i1, i2):
		q = self.q[:,i1:i2]
		mass = self.mass[i1:i2]
		p = self.p[:,i1:i2]# / (mass)**0.5
		x, y, z = q
		r = sqrt(x**2+y**2+z**2)
		mask = r < 8.
		print "enc mass for", name, "%g" % sum(mass[mask])
		#self.components[name] = (q, p)
		if i1 == i2:
			self.masses[name] = 0
			mass = 0
		else:
			self.masses[name] = mass[0] # ASSUME MASSES ARE EQUAL
			mass = mass[0]
		if i2 > i1:
			self.add_component(Component(name, q, p, mass))
		
		
		#dsa