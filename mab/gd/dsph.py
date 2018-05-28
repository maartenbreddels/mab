from math import *
import mab.astrometry
map = {}

class Sculptor(object):
	l0 = 287.5
	b0 = -83.2 
	ra0, dec0 = mab.astrometry.gal_to_eq(l0, b0)
	ra0 = ra0[0]
	dec0 = dec0[0]
	distance = 79. #kpc
	stellar_mass = 6.4e6
	core_radius_arcmin = 5.8
	core_radius_arcsec = core_radius_arcmin * 60
	core_radius_rad = core_radius_arcsec/3600 * pi/180
	core_radius_kpc = distance * tan(core_radius_rad)

class Fornax(object):
	l0 = 237.1
	b0 = -65.7
	ra0, dec0 = mab.astrometry.gal_to_eq(l0, b0)
	ra0 = ra0[0]
	dec0 = dec0[0]
	distance = 138. #kpc
	stellar_mass = 68.e6
	core_radius_arcmin = 13.8
	core_radius_arcsec = core_radius_arcmin * 60
	core_radius_rad = core_radius_arcsec/3600 * pi/180
	core_radius_kpc = distance * tan(core_radius_rad)

class Sextans(object):
	l0 = 243.5
	b0 = 42.3
	ra0, dec0 = mab.astrometry.gal_to_eq(l0, b0)
	ra0 = ra0[0]
	dec0 = dec0[0]
	distance = 86. #kpc
	stellar_mass = 19.e6
	core_radius_arcmin = 16.6
	core_radius_arcsec = core_radius_arcmin * 60
	core_radius_rad = core_radius_arcsec/3600 * pi/180
	core_radius_kpc = distance * tan(core_radius_rad)
	
class Carina(object):
	l0 = 260.1
	b0 = -22.2
	ra0, dec0 = mab.astrometry.gal_to_eq(l0, b0)
	ra0 = ra0[0]
	dec0 = dec0[0]
	distance = 101. # kpc
	stellar_mass = 13.e6
	core_radius_arcmin = 8.8 
	core_radius_arcsec = core_radius_arcmin * 60
	core_radius_rad = core_radius_arcsec/3600 * pi/180
	core_radius_kpc = distance * tan(core_radius_rad)
	
map["car"] = Carina
map["scl"] = Sculptor
map["fnx"] = Fornax
map["sex"] = Sextans


if __name__ == "__main__":
	for name, dsph in map.items():
		print name
		print "\t", dsph.core_radius_kpc