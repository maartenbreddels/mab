from kaplot import *
from mab.astrounits import *


class Orbit(object):
	def __init__(self, profile_model, orbit_integrator):
		self.profile_model = profile_model
		self.orbit_integrator = orbit_integrator
		
	def run(self, args, opts, scope):
		x = float(args[1])
		print "dsadsadas"
		print "x=",x
		E0 = self.profile_model.light_profile.potentialxz(x, 0.01)
		E = E0
		print "E0", E0
		#Einf = self.profile_model.light_profile.potentialxz(1e4, 0.)
		#E0 -= Einf
		xLzmax, Lzmax = self.profile_model.x_and_Lz_at_Lzmax(E0)
		vc = Lzmax/xLzmax
		Lz = Lzmax * float(args[2])
		print "Ekin", self.profile_model.light_profile.potentialxz_eff(x, 0., Lz)
		Eleft = self.profile_model.light_profile.potentialxz_eff(x, 0., Lz)
		vz = sqrt(Eleft)
		#saSS
		#E
		print "E0", E0
		#vzmax = sqrt(-2*E0)
		#print "vz_max", vzmax
		#vz = float(args[2]) * vzmax
		print "vz=", vz
		#E = E0+Einf
		
		E = self.profile_model.light_profile.potentialxz_eff(x, 0., Lz) + vz**2/2
		timescale = xLzmax/vc
		kpc_to_km = (1*KPC).asNumber(KM)
		s_to_gyr = (S/GYR).asNumber()
		
		timescale *= kpc_to_km
		print "t", timescale
		print "s_gyr", s_to_gyr
		print timescale*s_to_gyr
		
		
		x0, z0 = x, 0. #find_xz_zvc(theta, Lz, xLzmax)
		#dt = timescale * 0.01
		Lzs = array([Lz])
		s = 150
		dt = array([timescale/self.orbit_integrator.orbital_points*10.1*s])
		q0 = array([[x0, z0]])
		eps = 1e-16
		v = vz
		angle = math.radians(float(args[3]))
		p0 = array([[v*sin(angle), v*cos(angle)]])
		print dt, Lzs, q0, p0
		q1, p1 = self.orbit_integrator.integrate(dt, Lzs, q0, p0)
		x1 = ravel(q1[0,:])
		z1 = ravel(q1[1,:])
		vx1 = ravel(p1[0,:])
		vz1 = ravel(p1[1,:])
		box()
		scatter(x1, z1)
		xlim(0, z1.max())
		
		
		
		xmin = self.profile_model.x_min(E, Lz) * 1.001
		xmax = self.profile_model.x_max(E, Lz) * 0.999
		vline(xmin, color="red")
		vline(xmax, color="red")
		print xmin, xmax
		N = 100
		xs = arange(N) / (N-1.) * (xmax-xmin) + xmin
		print xs
		zmaxs = [self.profile_model.z_max_at_x(E, Lz, x) for x in xs]
		graph(xs, zmaxs, color="blue")
		
		print x1.max()
		print z1.max()
		xlim(0, xmax*1.1)
		
		draw()
		