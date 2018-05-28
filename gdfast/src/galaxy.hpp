#pragma once
#include "math.h"
#include <pyublas/numpy.hpp>
#include "math.h"
#include <tuple>
#include "rootfinder.hpp"
#include "optimization.hpp"
#include "profile.hpp"
#include "integration.hpp"

namespace gd {

void py_export_galaxy();

//! Represents a galaxy
/** 
Defines the potential
*/
class Galaxy {
public:
	virtual double dphidr(double r) = 0;
	virtual double kpc_to_arcsec(double x) = 0;
	virtual Profile* getStellarProfile() = 0;
};

class Galaxy1C_constant_anisotropy : public Galaxy {
public:
	Galaxy1C_constant_anisotropy(Profile* stellar_profile, double distance, double beta) : stellar_profile(stellar_profile), distance(distance), beta(beta) {
	}
	
	virtual double dphidr(double r) {
		return stellar_profile->dphidr(r);
	}
	double potentialr(double r) {
		//printf("%p\n", stellar_profile);
		return stellar_profile->potentialr(r);
	}
	virtual double kpc_to_arcsec(double x) {
		return x * 1.0/(distance) / (M_PI/180/60/60);
	}
	virtual Profile* getStellarProfile() {
		return stellar_profile;
	};
private:
	double energy_z(double z, double R, double energy) {
		double r = sqrt(R*R+z*z);
		double Epot = potentialr(r);
		//printf("-- %f %f\n", energy, Epot);
		return energy-Epot; // leftover energy
	}
	/*double _ekinrad(double r, double L, double E) {
		return (-L*L/(2*r*r) - (potentialr(r) - E));
	}*/
	double neg_vrad_mag(double r, double E) {
		r = fabs(r);
		return -r*r*(2*E-2*potentialr(r));
	}
public:		
	double findzmax(double R, double energy, double rmax) {
		std::function<double(double)> f(bind(&Galaxy1C_constant_anisotropy::energy_z, this, _1, R, energy));
		RootFinderGSL<> rootFinder(f); // the root finder
		//printf("%f %f\n", energy_z(0, R, energy), 0);
		//printf("%f %f\n", energy_z(rmax*1.1, R, energy), rmax*1.1);
		double zmax = rootFinder.findRoot(0, rmax*1.01);
		return zmax;
	}

	double Tr(double E, double L, double ra, double rp) ;

	double rmax_at_E(double E) {
		std::function<double(double)> f(bind(&Galaxy1C_constant_anisotropy::neg_vrad_mag, this, _1, E));
		MinimizerGSL<> minimizer(f);
		double r = minimizer.findMinimum(1e-5, 1);
		return r; }
/*
		def f(r, E=E):
			r = abs(r)
			return -r**2*(2*(E-self.potentialr(r)))
		rLmax, Lmaxsq = fmin(f, 1. ,full_output=True, disp=False)[:2]
		return sqrt(abs(Lmaxsq)), abs(rLmax[0])*/
		
		
	/*boost::python::tuple get_apo_peri(double E, double L, double rtry=1e-5) {
	}
	
	double findR_at_EL(double E, double L, double rtry=1) {
		def f(r, E=E, L=L):
			r = abs(r)
			return L**2/(2*r**2) + (self.potentialr(r) - E)
		r = fsolve(f, rtry)
		#Lmax = extra[0]
		#print r
		return abs(r)
	*/

	Profile* stellar_profile;
	double distance;
	double beta; 

};

class Galaxy2C_constant_anisotropy : public Galaxy1C_constant_anisotropy {
public:
	Galaxy2C_constant_anisotropy(Profile* stellar_profile, Profile* dm_profile, double distance, double beta) : Galaxy1C_constant_anisotropy(stellar_profile, distance, beta), dm_profile(dm_profile) {
	}
	
	double potentialr(double r) {
		//printf("%p\n", stellar_profile);
		return stellar_profile->potentialr(r) + dm_profile->potentialr(r);
	}
	virtual double dphidr(double r) {
		return stellar_profile->dphidr(r) + dm_profile->dphidr(r);
	}
	Profile* dm_profile;
};

}