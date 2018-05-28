# -*- coding: utf-8 -*-
from numpy import *
import numpy
from kapteyn import celestial
from mab.gd import gdfast
# in a right handed coordinates system! component 0 points to GC, 1 in direction of rotation, 2 to NGP
v_sun_lsr = array([10.0, 5.2, 7.2])
v_lsr_gsr = array([0., 220, 0])
v_sun_gsr = v_lsr_gsr + v_sun_lsr
v_sun_gsr_ = v_sun_gsr # alias

pc_to_meter = 1/3.24077649e-17
kpc_to_km = pc_to_meter
year_to_second = 31556926.
mas_to_degree = 1./(60*60*1000)

theta0 = 122.932
#d_NGP = 27.128336111111111
#al_NGP = 192.10950833333334
al_NGP = 192.85948
d_NGP = 27.12825

#theta0 = 123.
#al_NGP = 192.25
#d_NGP = 27.4


k = 4.74057

U_sol = 10.0
V_sol = 5.25
W_sol = 7.17
v_lsr = 220.

import mab

zero_ref_frame = gdfast.ZeroReferenceFrame()
cartesian = gdfast.Cartesian()
spherical = gdfast.SphericalGalactic()
cylindrical = gdfast.Cylindrical()

eq_reference_frame = gdfast.EqReferenceFrame(radians(theta0), radians(al_NGP), radians(d_NGP), zero_ref_frame)


def vlos_helio_to_sys(ra, dec, vlos_helio, correction):
	vlos_helio = stars.vlos_helio
	#e_vlos = stars.e_vlos
	l, b = mab.astrometry.eq_to_gal(ra, dec)
	e_x = cos(radians(l)) * cos(radians(b))
	e_y = sin(radians(l)) * cos(radians(b))
	e_z = sin(radians(b))
	e_los = transpose(array([e_x, e_y, e_z]))
	
	v_sun_gsr = self.correction.v_sun_gsr
	v_sys_gsr = self.correction.get_vsys_gsr()
	logger.info("systemic velocity of object: %r" % v_sys_gsr)
	
	# TODO: double check with paper!
	vlos_gsr = vlos_helio - dot(e_los, v_sys_gsr-v_sun_gsr)
	print vlos_gsr
	import pdb
	pdb.set_trace()
	#stars = stars.clone()
	#stars.vlos = vlos_gsr
	#stars[0].attributes.append("vlos")

class ProperMotion(object):
	def __init__(self, ra, dec, pm_alpha, pm_delta, distance, sigma_pm_alpha, sigma_pm_delta, sigma_vr, v_sun_gsr=None, vr_helio=None):
		self.ra = ra
		self.dec = dec
		self.pm_alpha = pm_alpha
		self.pm_delta = pm_delta
		self.distance = distance * 1000 # kpc to pc
		self.vr_helio = vr_helio
		self.v_sun_gsr = v_sun_gsr
		self.sigma_vr = sigma_vr
		self.sigma_pm_alpha = sigma_pm_alpha
		self.sigma_pm_delta = sigma_pm_delta
		if v_sun_gsr is None:
			self.v_sun_gsr = v_sun_gsr_
		self.v_sun_gsr = array(self.v_sun_gsr)
			
		self.alpha = radians(mab.astrometry.read_ra_hour(self.ra))
		self.delta = radians(mab.astrometry.read_dec_deg(self.dec))
		
	def load(self):
		pass
		
	def get_vsys_gsr(self):
		return self.vrpm_eq_to_eq(self.vr_helio, self.pm_alpha, self.pm_delta)
	
	def vrpm_eq_to_eq(self, vr, pm_alpha, pm_delta):
		# convert vr pm_a, pm_d to velocity wrt 
		coordinate = gdfast.Coordinate(1, self.alpha, -self.delta+pi/2, spherical)
		position = gdfast.Position(coordinate, eq_reference_frame)
		
		velocity_coordinate_scl_radial = gdfast.VelocityCoordinate(vr, pm_alpha*1e-5 * k *self.distance, -pm_delta*1e-5 * k *self.distance, spherical)
		velocity = gdfast.Velocity(velocity_coordinate_scl_radial, eq_reference_frame, position, eq_reference_frame)
		#print velocity.to_coordinate_system(zero_ref_frame, spherical)
		#v = velocity.to_coordinate_system(zero_ref_frame, spherical)
		v = velocity.to_coordinate_system(zero_ref_frame, cartesian) + self.v_sun_gsr
		return v
		
	def convert(self):
		#pm_a, pm_d = 9., 2.
		#sigma_pma, sigma_pmd = 13., 13.
		#vr_scl = 109.9
		#vr = 109.9
		
		if 0:
			T = mab.astrometry.getT()
			print "T=", T
			#alpha = radians(mab.astrometry.read_ra_hour("1 00 09"))
			#delta = radians(mab.astrometry.read_dec_deg("-33 42 30"))
			print degrees(self.alpha), degrees(self.delta)
			x = array([cos(self.delta) * cos(self.alpha), cos(self.delta) * sin(self.alpha), sin(self.delta)])
			y = ravel(dot(T, x))
			print y
			sinb = y[2]
			b = arcsin(sinb)
			cosl = y[0]/cos(b)
			l = arccos(cosl)
			sinl = y[1]/cos(b)
			l = arcsin(sinl)
			l = (l + 2 * pi) % (2*pi)
			print l, b
			print "l,b",degrees(l), degrees(b)
			
		v = self.vrpm_eq_to_eq(self.vr_helio, self.pm_alpha, self.pm_delta)
		if 0:
			coordinate = gdfast.Coordinate(1, alpha, -delta+pi/2, spherical)
			position = gdfast.Position(coordinate, eq_reference_frame)
			
			velocity_coordinate_scl_radial = gdfast.VelocityCoordinate(vr, pm_a*1e-5 * k *d, -pm_d*1e-5 * k *d, spherical)
			velocity = gdfast.Velocity(velocity_coordinate_scl_radial, eq_reference_frame, position, eq_reference_frame)
			#print velocity.to_coordinate_system(zero_ref_frame, spherical)
			#v = velocity.to_coordinate_system(zero_ref_frame, spherical)
			v = velocity.to_coordinate_system(zero_ref_frame, cartesian) + v_sun_gsr
		print v
		
	
		# finite differenceing function
		def fd(f, i, *args):
			x = array(args)
			org = x* 1.0
			dy = 10.1#1e-2
			x1 = x
			x2 = x * 1.0
			x2[i] += dy
			l = list((f(*x2) - f(*x1)) / dy)
			#print ["%e" % k for k in l]
			return l		
		
		v1 = fd(self.vrpm_eq_to_eq, 0, self.vr_helio, self.pm_alpha, self.pm_delta)
		v2 = fd(self.vrpm_eq_to_eq, 1, self.vr_helio, self.pm_alpha, self.pm_delta)
		v3 = fd(self.vrpm_eq_to_eq, 2, self.vr_helio, self.pm_alpha, self.pm_delta)
		M = matrix(array([v1, v2, v3])).T
		print "M", M
		print "M", matrix(M) * matrix(M).T
		
		
		cov = array([[self.sigma_vr**2, 0, 0],  [0, self.sigma_pm_alpha**2, 0], [0, 0, self.sigma_pm_delta**2]])
		covN = M * cov * M.T
		print "cov ", cov
		print "covN", covN
		self.covN = covN
		
		
	def run(self, opts, args, scope):
		self.convert()
		
		
 
def gal_to_eq(l, b):
	M = celestial.sky2sky(celestial.gal, (celestial.eq, celestial.fk5), l, b)
	a, d = ravel(M[:,0]), ravel(M[:,1])
	return a, d

def eq_to_gal(a, d):
	M = celestial.sky2sky((celestial.eq, celestial.fk5), celestial.gal, a, d)
	try:
		len(a)
		l, b = ravel(M[:,0]), ravel(M[:,1])
	except:
		l = M[0,0]
		b = M[0,1]
	return l, b

def _lbr_to_xyz(l, b, r):
	return  r * cos(radians(l)) * cos(radians(b)),\
			r * sin(radians(l)) * cos(radians(b)),\
			r* sin(radians(b))
			
def xyz_to_lbr(x, y, z):
	r = sqrt(x*x+y*y+z*z)
	b = degrees(arcsin(z/r))
	l = ((degrees(arctan2(y, x))+360)%360)
	return l, b, r

def lbr_to_xyz(l, b, r):
	return  r * cos(radians(l)) * cos(radians(b)),\
			r * sin(radians(l)) * cos(radians(b)),\
			r * sin(radians(b))
			
def _xyz_to_lbr(x, y, z):
	r = sqrt(x*x+y*y+z*z)
	b = degrees(arcsin(z/r))
	l = ((degrees(arctan2(y, x))+360)%360)
	return l, b, r


N = 5000

#degra = 4.0 * numpy.arctan(1.0)/180
cosd = lambda x: math.cos(math.radians(x))
sind = lambda x: math.sin(math.radians(x))

def getT(theta0=theta0, d_NGP=d_NGP, al_NGP=al_NGP):
	
	# numpy matrices: row,column format
	c = numpy.matrix([
			[cosd(al_NGP),  sind(al_NGP), 0],
			[sind(al_NGP), -cosd(al_NGP), 0],
			[0, 0, 1]])
	b = numpy.matrix([
			[-sind(d_NGP), 0, cosd(d_NGP)],
			[0, -1, 0],
			[cosd(d_NGP), 0, sind(d_NGP)]])
	a = numpy.matrix([
			[cosd(theta0),  sind(theta0), 0],
			[sind(theta0), -cosd(theta0), 0],
			[0, 0, 1]])
	T = a*b*c
	return T

def lbrv_to_proper_motion(l, b, distance, vrel, dt = 1e-3):
	l0, b0 = l, b
	x0, y0, z0 = lbr_to_xyz(l, b, distance)
	p0 = array([x0, y0, z0])
	
	vrel_kpc_per_year = vrel/kpc_to_km*year_to_second
	vrel_kpc_per_century = vrel_kpc_per_year * 100
	p1 = p0 + vrel_kpc_per_century*dt
	l1, b1, r = xyz_to_lbr(*p1)
	
	a0, d0 = gal_to_eq(l0, b0)
	a1, d1 = gal_to_eq(l1, b1)
	cosd = cos(radians(d0))
	cosb = cos(radians(b0))
	
	return (l1-l0)/dt/mas_to_degree*cosb, (b1-b0)/dt/mas_to_degree, (a1-a0)/dt/mas_to_degree*cosd, (d1-d0)/dt/mas_to_degree

def lbrv_to_proper_motion_xi_eta(l, b, distance, vrel, dt = 1e-3):
	x0, y0, z0 = lbr_to_xyz(l, b, distance)
	p0 = array([x0, y0, z0])
	
	vrel_kpc_per_year = vrel/kpc_to_km*year_to_second
	vrel_kpc_per_century = vrel_kpc_per_year * 100
	p1 = p0 + vrel_kpc_per_century*dt
	l0, b0 = l, b
	l1, b1, r = xyz_to_lbr(*p1)
	
	a0, d0 = gal_to_eq(l0, b0)
	a1, d1 = gal_to_eq(l1, b1)
	
	xi0, eta0 = ra_dec_to_xi_eta(a0, d0, a0, d0)
	xi1, eta1 = ra_dec_to_xi_eta(a1, d1, a0, d0)
	
	#cosd = cos(radians(d0))
	#cosb = cos(radians(b0))
	
	return (xi1-xi0)/dt/mas_to_degree, (eta1-eta0)/dt/mas_to_degree

#def vrpm_eq_to_gal(vr, pm_a, pm_d, a, d):


def read_ra_hour(ra):
	hour, minute, sec = ra.split()
	sign = 1
	if ra.strip()[0] == "-":
		sign = -1
	#sign = 1
	deg = float(hour)/24*360 + sign*float(minute)/60*360/24 + sign*float(sec)/(60*60)*360/24
	return deg
	
def read_dec_deg(dec):
	deg, minute, sec = dec.split()
	sign = 1
	if dec.strip()[0] == "-":
		sign = -1
	#sign = 1
	deg = float(deg) + sign*float(minute)/60 + sign*float(sec)/(60*60)
	return deg
	
def r__a_dec_to_re(ra, dec, ra_center, dec_center, position_angle, ellipticity):
	#ra_c = 15.0375
	#dec_c = -33.7083333

	#pa = 90. - 99.
	#ell = 0.32

	#ra = ra_scl
	#dec = dec_scl
	radeg = 180/pi

	den = sin(dec_center) * tan(dec) + cos(dec_center) * cos(ra - ra_center)
	
	xi = sin((ra - ra_center))/den
	eta = (cos(dec_center) * tan(dec) - sin(dec_center) * cos(ra - ra_center))/den
	
	xin = xi * cos(position_angle) + eta * sin(position_angle)
	etan = (-xi * sin(position_angle) + eta * cos(position_angle))*(1 - ellipticity)
	
	eps = sqrt(xin**2 + etan**2)*radeg
	return eps
"""
den = sin(dec_c/!radeg) * tan(dec/!radeg) + $
  cos(dec_c/!radeg) * cos((ra - ra_c)/!radeg)

xi = sin((ra - ra_c)/!radeg)/den
eta = (cos(dec_c/!radeg) * tan(dec/!radeg) - $
      sin(dec_c/!radeg) * cos((ra - ra_c)/!radeg))/den

xin = xi * cos(pa/!radeg) + eta * sin(pa/!radeg)
etan = (-xi * sin(pa/!radeg) + eta * cos(pa/!radeg))/(1 - ell)




pa = 90. - 99.
ell = 0.32

ra = ra_scl
dec = dec_scl

den = sin(dec_c/!radeg) * tan(dec/!radeg) + $
  cos(dec_c/!radeg) * cos((ra - ra_c)/!radeg)

xi = sin((ra - ra_c)/!radeg)/den
eta = (cos(dec_c/!radeg) * tan(dec/!radeg) - $
      sin(dec_c/!radeg) * cos((ra - ra_c)/!radeg))/den

xin = xi * cos(pa/!radeg) + eta * sin(pa/!radeg)
etan = (-xi * sin(pa/!radeg) + eta * cos(pa/!radeg))/(1 - ell)

eps = sqrt(xin^2 + etan^2)*!radeg
"""
	
def xi_eta_to_re(xi, eta, position_angle, ellipticity):
	xin = (xi * cos(radians(position_angle)) + eta * sin(radians(position_angle)))
	etan = (-xi * sin(radians(position_angle)) + eta * cos(radians(position_angle)))/(1 - ellipticity)
	eps = sqrt(xin**2 + etan**2)
	return eps
	
def ra_dec_to_xi_eta(ra, dec, ra_center, dec_center):
	den = sin(radians(dec_center)) * tan(radians(dec)) + cos(radians(dec_center)) * cos(radians(ra - ra_center))
	
	xi = sin(radians(ra - ra_center))/den* 180/pi
	eta = (cos(radians(dec_center)) * tan(radians(dec)) - sin(radians(dec_center)) * cos(radians(ra - ra_center)))/den * 180/pi
	return xi, eta

def xi_eta_to_xin_etan(xi, eta, position_angle, ellipticity):
	xin = xi * cos(radians(position_angle)) + eta * sin(radians(position_angle))
	etan = (-xi * sin(radians(position_angle)) + eta * cos(radians(position_angle))) / (1 - ellipticity)
	return xin, etan

def __lalala(xin, etan, ra_center, dec_center, position_angle, ellipticity):
	
	x = xin * cos(radians(-position_angle)) + etan * sin(-radians(-position_angle))
	y = (-xin * sin(radians(-position_angle)) + etan * cos(-radians(-position_angle)))
	#return x, y
	return xin, etan

