from math import *
if 0:
	a = 7.5657e-15
	k = 1.3807e-16
	eVtoerg = 1.6022e-12
	ec = 4.803204e-10 # element. charge
	h = 6.626068e-27
	hbar = h / (2*pi)
	amu = 1.66053886e-24
	barn2cm2 = 1e-24
	c = 3e10
	G = 6.67428e-8
	mh = 1.00794 * amu
	Rg = k/mh
	Msun = 1.9891e33
	Lsun = 3.827e33
	Rsun = 6.960e10
# SI units
from unum.units import *
#from unum import Unum


class solar(object):
        mass = 1.98892e30


amu = 1.66053886e-27 * KG
mh = 1.007825 * amu
#secondsperyear = 31556926.

G = 6.67428e-11 * M ** 3 / (KG * S**2)

epsilon0 = 8.854187817e-12 * A**2 * S**4 * KG**-1 * M**-3
e_charge = 1.602176487e-19 * C
me = 9.1093821545e-31 * KG
plank = 6.626068e-34 * M**2 * KG / S
kb = 1.3806503e-23* M**2 * KG * S**-2 * K**-1
c = 299792458 * M / S

sigma_SB = 5.6704e-8 * J/S/M**2/K**4