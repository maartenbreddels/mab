from unum.units import *
from unum import Unum

MSOL = Unum.unit("MSOL", 1.98892e30 * KG)
RSOL = Unum.unit("RSOL", 6.955e8 * M)
PC = Unum.unit("pc", 3.08568025e16 * M)
KPC = Unum.unit("Kpc", PC * 1e3)
MPC = Unum.unit("Mpc", PC * 1e6)

CM = Unum.unit("cm", 1e-2 * M)
KM = Unum.unit("km", 1e3 * M)

YR = Unum.unit("yr", 31556926. * S)
MYR = Unum.unit("Myr", YR * 1e6)
GYR = Unum.unit("Gyr", YR * 1e9)


AMU = Unum.unit("amu", 1.66053886e-27 * KG)
MH = Unum.unit("mh", 1.007825 * AMU)
MHe = Unum.unit("mhe", 4.002602 * AMU)

hubbleH = Unum.unit("hubbleH")

LSOL = 3.839e26 *J/S