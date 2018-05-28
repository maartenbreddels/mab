import snapshot
import os
from numpy import *
velocity_unit = 100# km/s

mass_unit = 2.325e9# M_sun

class IC(snapshot.Snapshot):
	def __init__(self, dirname, read_disk=True, read_halo=True, read_bulge=True):
		super(IC, self).__init__()
		self.dirname = dirname
		self.read_disk = read_disk
		self.read_halo = read_halo
		self.read_bulge = read_bulge
		
	def load(self):
		if self.read_disk:
			self.add_ics("disk")
		if self.read_halo:
			self.add_ics("halo")
		if self.read_bulge:
			self.add_ics("bulge")

	def add_ics(self, name):
		filename = os.path.join(self.dirname, name)
		lines = file(filename).readlines()
		n = int(lines[0].split()[0])
		grid = zeros((7, n))
		i = 0
		for line in lines[1:]:
			values = [float(k) for k in line.split()]
			grid[:,i] = values
			i += 1
		mass = grid[0]
		q = grid[1:1+3]
		x, y, z = q
		r = sqrt(x**2 + y**2 +z**2)
		mask = r < 1
		p = grid[4:4+3]
		# assume all masses are equal
		print "total mass", name, "%g" % sum(mass*mass_unit)
		print "total mass within 1kpc", name, "% .3g" % sum(mass[mask]*mass_unit)
		c = snapshot.Component(name, q, p*velocity_unit, mass=mass[0]*mass_unit)
		self.add_component(c)
		