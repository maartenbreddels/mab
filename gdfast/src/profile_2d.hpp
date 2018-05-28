#pragma once
#include "math.h"
//#include <tuple>
//#include "integration.hpp"
//#include "datatypes.hpp"

namespace gd {
using namespace std;

void py_export_profile_axi();

class Density2d {
public:
	virtual double densityxy(double x, double y) = 0;
};

class Potential2d {
public:
	virtual double dphidx(double x, double y) = 0;
	virtual double dphidy(double x, double y) = 0;
	virtual double potentialxy(double x, double y) = 0;
	virtual double potentialxy_eff(double x, double z, double Lz) {
		return this->potentialxy(x, z);// + Lz*Lz/(2*x*x);
	}
	virtual double dphidx_eff(double x, double y, double Lz) {
		return dphidx(x, y);// -2* Lz * Lz / ( 2* x*x*x);
	}
	virtual double dphidy_eff(double x, double y, double) {
		return dphidy(x, y); 
	}
};

class Profile2d : public Density2d, public Potential2d {
};

class ProfileModel2d  {
public:
	virtual double densityxy(double x, double y) = 0;
	virtual double dphidx(double x, double y) = 0;
	virtual double dphidy(double x, double y) = 0;
	virtual double potentialxy(double x, double y) = 0;
	virtual double potentialxy_eff(double x, double y, double Lz) = 0;
	virtual double dphidx_eff(double x, double y, double Lz) = 0;
	virtual double dphidy_eff(double x, double y, double Lz) = 0;
};

class ProfileModel2d1C : public ProfileModel2d {
public:
	ProfileModel2d1C(Profile2d* p) : p(p) {} 
	virtual double densityxy(double x, double z) { return p->densityxy(x, z); }
	virtual double dphidx(double x, double z) { return p->dphidx(x, z); }
	virtual double dphidy(double x, double y) { return p->dphidy(x, y); }
	virtual double potentialxy(double x, double y) { return p->potentialxy(x, y); }
	virtual double potentialxy_eff(double x, double y, double Lz) { return p->potentialxy_eff(x, y, Lz); }
	virtual double dphidx_eff(double x, double y, double Lz) { return p->dphidx_eff(x, y, Lz); }
	virtual double dphidy_eff(double x, double y, double Lz) { return p->dphidy_eff(x, y, Lz); }
	Profile2d* p;
};

/*class ProfileModel2d2C : public ProfileModel2d {
public:
ProfileModel2C(Profile* p1, Profile* p2) : p1(p1), p2(p2) {} 
	virtual double densityr(double r) { return p1->densityr(r) + p2->densityr(r); } 
	virtual double dphidr(double r) { return p1->dphidr(r) + p2->dphidr(r); }
	virtual double potentialr(double r) { return p1->potentialr(r) + p2->potentialr(r); }
	Profile *p1, *p2;
};*/


class Logarithmic2d : public Profile2d {
public:
	/* comment: q >= 1/sqrt(2), otherwise the density can become negative (BT) */
	Logarithmic2d(double v0, double q, double Rc, double G) : v0(v0), q(q), Rc(Rc), G(G) {}
	virtual double densityxy(double x, double z) {
		return v0*v0 / (4*M_PI*G) * (1.); // implement this
	}
	virtual double potentialxy(double x, double y) {
		return 0.5*v0*v0 * log(Rc*Rc + x*x + y*y/(q*q));
	}
	virtual double dphidx(double x, double y) {
		return v0*v0 * x /(Rc*Rc + x*x + y*y/(q*q));
	}
	virtual double dphidy(double x, double y) {
		double qsq = q * q;
		return v0*v0 * y/qsq / (Rc*Rc + x*x + y*y/qsq);
	}
protected:
	double v0, q, Rc, G;
};


}