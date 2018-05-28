#pragma once
#include "math.h"
#include "math.h"
#include <tuple>
#include "integration.hpp"
//#include "datatypes.hpp"
#include <boost/python.hpp>
#include <boost/math/special_functions.hpp>
#include "gsl/gsl_sf_bessel.h"

namespace gd {
using namespace std;

void py_export_profile();

class Density3d {
public:
	virtual double density(double x, double y, double z) = 0;  
};

class Profile3d : public Density3d {
public:
	virtual double density(double x, double y, double z) = 0;  
	virtual double potential(double x, double y, double z) = 0;
	virtual double dphidx(double x, double y, double z) = 0;
	virtual double dphidy(double x, double y, double z) = 0;
	virtual double dphidz(double x, double y, double z) = 0;

	
	/*double enclosed_mass(double r) {
		auto dmass = [this](double rp){ return this->densityr(rp) * rp*rp * 4 * M_PI; };
		IntegratorGSL<> integratorGSL(dmass); // the integrator
		double integral = integratorGSL.integrate(0,r);
		return integral;
	}
	double total_mass() {
		auto dmass = [this](double rp){ return this->densityr(rp) * rp*rp * 4 * M_PI; };
		IntegratorGSL<> integratorGSL(dmass); // the integrator
		double integral = integratorGSL.integrate_to_inf(0);
		return integral;
	}*/
};

class ProfileModel3d {
public:
	virtual double density(double x, double y, double z) = 0;  
	virtual double potential(double x, double y, double z) = 0;
	virtual double dphidx(double x, double y, double z) = 0;
	virtual double dphidy(double x, double y, double z) = 0;
	virtual double dphidz(double x, double y, double z) = 0;
};

class ProfileModel3d1C : public ProfileModel3d {
public:
	ProfileModel3d1C(Profile3d* p) : p(p) {} 
	virtual double density(double x, double y, double z) { return p->density(x, y, z); }
	virtual double potential(double x, double y, double z) { return p->potential(x, y, z); }
	virtual double dphidx(double x, double y, double z) { return p->dphidx(x, y, z); }
	virtual double dphidy(double x, double y, double z) { return p->dphidy(x, y, z); }
	virtual double dphidz(double x, double y, double z) { return p->dphidz(x, y, z); }
	Profile3d* p;
};

class ProfileModel3d2C : public ProfileModel3d {
public:
	ProfileModel3d2C(Profile3d* p1, Profile3d* p2) : p1(p1), p2(p2) {} 
	virtual double density(double x, double y, double z) { return p1->density(x, y, z) +  p2->density(x, y, z); }
	virtual double potential(double x, double y, double z) { return p1->potential(x, y, z) + p2->potential(x, y, z); }
	virtual double dphidx(double x, double y, double z) { return p1->dphidx(x, y, z) + p2->dphidx(x, y, z); }
	virtual double dphidy(double x, double y, double z) { return p1->dphidy(x, y, z) + p2->dphidy(x, y, z); }
	virtual double dphidz(double x, double y, double z) { return p1->dphidz(x, y, z) + p2->dphidz(x, y, z); }
	Profile3d* p1;
	Profile3d* p2;
};


class NFW3d : public Profile3d {
public:
	NFW3d(double mass200, double rs, double G, double rho_crit) : mass200(mass200), rs(rs), G(G) {
		r200 = pow(mass200/200*3/(4*M_PI)/rho_crit, 1./3);
		c = r200/rs;
		double x = r200 / rs;
		rho0 = mass200 / (4*M_PI*pow(rs,3)*(log(1+x)-x/(1+x)));
	}
	virtual double density(double x, double y, double z) {
		return this->densityr(sqrt(x*x + y*y + z*z));
	}
	virtual double potential(double x, double y, double z) {
		return this->potentialr(sqrt(x*x + y*y + z*z));
	}
	virtual double dphidx(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double dphidr = this->dphidr(r);
		return -dphidr*x/r;
	}
	virtual double dphidy(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double dphidr = this->dphidr(r);
		return -dphidr*y/r;
	}
	virtual double dphidz(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double dphidr = this->dphidr(r);
		return -dphidr*z/r;
	}
	double potentialr(double r) {
		double x = r/rs;
		return - 4 * M_PI * G * rho0 * rs*rs * log(1+x)/x;
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (x * pow(1+x, 2));
	}
	double dphidr(double r) {
		double x = r/rs;
		double f = -4.0*M_PI*G*rho0*pow(rs, 3)/pow(r,2);
		double t1 = log(1.0 + x);
		double t2 = x/(1.0+x);
		//cout << x << ", " << f << ", " << t1 << ", " << t2 << endl;
		return -f*(t1 - t2);
	}
	double mass200, rho0, rs, r200, c, G;
};


class NFW3dTriax : public Profile3d {
public:
	NFW3dTriax(double mass200, double rs, double G, double rho_crit, double a, double b, double c, double ra) : mass200(mass200), rs(rs), G(G), a(a), b(b), c(c), ra(ra)  {
		r200 = pow(mass200/200*3/(4*M_PI)/rho_crit, 1./3);
		concentration = r200/rs;
		double x = r200 / rs;
		rho0 = mass200 / (4*M_PI*pow(rs,3)*(log(1+x)-x/(1+x)));
	}
	virtual double density(double x, double y, double z) {
		return 0 ;//this->densityr(sqrt(x*x + y*y + z*z));
	}
	virtual double potential(double x, double y, double z) {
		return this->potentialr( this->rtilde(x, y, z));
	}
	double rtilde(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double re = sqrt(x*x / (a*a) + y*y/ (b*b) + z*z/(c*c) );
		return (ra+r)*re / (ra+re);
	}
	virtual double dphidx(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double re = sqrt(x*x / (a*a) + y*y/ (b*b) + z*z/(c*c) );
		double rtilde = this->rtilde(x, y, z);
		double dphidrtilde = this->dphidr(rtilde);
		double drtildedx = x * (ra * ra * r + ra * r*r + a*a*pow(re,3)  + ra*a*a*re*re) / (a*a*r*re*pow(ra+re, 2));
		return -dphidrtilde*drtildedx;
	}
	virtual double dphidy(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double re = sqrt(x*x / (a*a) + y*y/ (b*b) + z*z/(c*c) );
		double rtilde = this->rtilde(x, y, z);
		double dphidrtilde = this->dphidr(rtilde);
		double drtildedy = y * (ra * ra * r + ra * r*r + b*b*pow(re,3)  + ra*b*b*re*re) / (b*b*r*re*pow(ra+re, 2));
		return -dphidrtilde*drtildedy;
	}
	virtual double dphidz(double x, double y, double z) {
		double r = sqrt(x*x + y*y + z*z);
		double re = sqrt(x*x / (a*a) + y*y/ (b*b) + z*z/(c*c) );
		double rtilde = this->rtilde(x, y, z);
		double dphidrtilde = this->dphidr(rtilde);
		double drtildedz = z * (ra * ra * r + ra * r*r + c*c*pow(re,3)  + ra*c*c*re*re) / (c*c*r*re*pow(ra+re, 2));
		return -dphidrtilde*drtildedz;
	}
	double potentialr(double r) {
		double x = r/rs;
		return - 4 * M_PI * G * rho0 * rs*rs * log(1+x)/x;
	}
	double densityr(double r) {
		double x = r/rs;
		return rho0 / (x * pow(1+x, 2));
	}
	double dphidr(double r) {
		double x = r/rs;
		double f = -4.0*M_PI*G*rho0*pow(rs, 3)/pow(r,2);
		double t1 = log(1.0 + x);
		double t2 = x/(1.0+x);
		//cout << x << ", " << f << ", " << t1 << ", " << t2 << endl;
		return -f*(t1 - t2);
	}
	double mass200, rho0, rs, r200, concentration, G, a, b, c, ra;
};


class  MiyamotoNagai3d : public Profile3d {
public:
	MiyamotoNagai3d(double mass, double G, double a, double b) : mass(mass), G(G), a(a), b(b)  {
	}
	virtual double density(double x, double y, double z) {
		return 0 ;//this->densityr(sqrt(x*x + y*y + z*z));
	}
	virtual double potential(double x, double y, double z) {
		return -G*mass*pow(x*x + y*y +pow(a + pow(z*z+b*b, .5), 2.), -.5);
		//return this->potentialr( this->rtilde(x, y, z));
	}
	virtual double dphidx(double x, double y, double z) {
		double part = G*mass*pow( x*x + y*y + pow(a+ pow(z*z+b*b,.5) , 2.), -1.5);
		return -x*part;
	}
	virtual double dphidy(double x, double y, double z) {
		double part = G*mass*pow( x*x + y*y + pow(a+ pow(z*z+b*b,.5) , 2.), -1.5);
		return -y*part;
	}
	virtual double dphidz(double x, double y, double z) {
		double part = G*mass*pow( x*x + y*y + pow(a+ pow(z*z+b*b,.5) , 2.), -1.5);
		return -z*part*(a/sqrt(b*b+z*z) + 1);
	}
	double mass, G, a, b;
};


}
