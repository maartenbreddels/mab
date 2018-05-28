#pragma once
#include "math.h"
//#include <tuple>
#include "integration.hpp"
//#include "datatypes.hpp"

namespace gd {
using namespace std;

void py_export_profile_axi();

class DensityAxi {
public:
	virtual double densityxz(double x, double z) = 0;
};

class PotentialAxi {
public:
	virtual double dphidx(double x, double z) = 0;
	virtual double dphidz(double x, double z) = 0;
	virtual double potentialxz(double x, double z) = 0;
	virtual double potentialxz_eff(double x, double z, double Lz) {
		return this->potentialxz(x, z) + Lz*Lz/(2*x*x);
	}
	virtual double dphidx_eff(double x, double z, double Lz) {
		return dphidx(x, z) -2* Lz * Lz / ( 2* x*x*x);
	}
	virtual double dphidz_eff(double x, double z, double) {
		return dphidz(x, z); 
	}
};

class ProfileAxi : public DensityAxi, public PotentialAxi {
};

class ProfileModelAxi  {
public:
	virtual double densityxz(double x, double z) = 0;
	virtual double dphidx(double x, double z) = 0;
	virtual double dphidz(double x, double z) = 0;
	virtual double potentialxz(double x, double z) = 0;
	virtual double potentialxz_eff(double x, double z, double Lz) = 0;
	virtual double dphidx_eff(double x, double z, double Lz) = 0;
	virtual double dphidz_eff(double x, double z, double Lz) = 0;
};

class ProfileModelAxi1C : public ProfileModelAxi {
public:
	ProfileModelAxi1C(ProfileAxi* p) : p(p) {} 
	virtual double densityxz(double x, double z) { return p->densityxz(x, z); }
	virtual double dphidx(double x, double z) { return p->dphidx(x, z); }
	virtual double dphidz(double x, double z) { return p->dphidz(x, z); }
	virtual double potentialxz(double x, double z) { return p->potentialxz(x, z); }
	virtual double potentialxz_eff(double x, double z, double Lz) { return p->potentialxz_eff(x, z, Lz); }
	virtual double dphidx_eff(double x, double z, double Lz) { return p->dphidx_eff(x, z, Lz); }
	virtual double dphidz_eff(double x, double z, double Lz) { return p->dphidz_eff(x, z, Lz); }
	ProfileAxi* p;
};

/*class ProfileModelAxi2C : public ProfileModelAxi {
public:
ProfileModel2C(Profile* p1, Profile* p2) : p1(p1), p2(p2) {} 
	virtual double densityr(double r) { return p1->densityr(r) + p2->densityr(r); } 
	virtual double dphidr(double r) { return p1->dphidr(r) + p2->dphidr(r); }
	virtual double potentialr(double r) { return p1->potentialr(r) + p2->potentialr(r); }
	Profile *p1, *p2;
};*/


class LogarithmicAxi : public ProfileAxi {
public:
	/* comment: q >= 1/sqrt(2), otherwise the density can become negative (BT) */
	LogarithmicAxi(double v0, double q, double xc, double G) : v0(v0), q(q), xc(xc), G(G) {}
	virtual double densityxz(double x, double z) {
		return v0*v0 / (4*M_PI*G) * (1.); // implement this
	}
	virtual double potentialxz(double x, double z) {
		return 0.5*v0*v0 * log(xc*xc + x*x + z*z/(q*q));
	}
	virtual double dphidx(double x, double z) {
		return v0*v0 * x /(xc*xc + x*x + z*z/(q*q));
	}
	virtual double dphidz(double x, double z) {
		double qsq = q * q;
		return v0*v0 * z/qsq / (xc*xc + x*x + z*z/qsq);
	}
protected:
	double v0, q, xc, G;
};

class ProfileStackelOblate : public ProfileAxi {
public:
	ProfileStackelOblate(double gamma, double alpha, double G) : gamma(gamma), alpha(alpha), G(G) {\
	}
	virtual double psi_tau(double tau)  = 0;

	virtual double densityxz(double x, double z) {
		double nu, la;
		xy_to_conf(x, z, nu, la);
		//printf("x=%f y=%f nu=%f la=%f\n", x, z, nu, la);
		return density_conf(nu, la);
	}
	virtual double potentialxz(double x, double z) {
		return 0;
	}
	double dphi_dla(double nu, double la) {
		double a = (la + gamma) * G_tau(la) - (nu + gamma) * G_tau(nu);
		return a/ pow(la - nu, 2) - (G_tau(la)+(gamma+la) * G_tau_prime(la))/(la-nu);
	} 
	double dphi_dnu(double nu, double la) {
		double a = (la + gamma) * G_tau(la) - (nu + gamma) * G_tau(nu);
		return -a/ pow(la - nu, 2) - (-G_tau(nu)-(gamma+nu) * G_tau_prime(nu))/(la-nu);
	} 

	virtual double dphidx(double x, double z) {
		double nu, la;
		xy_to_conf(x, z, nu, la);
		double Psq = (la - nu) / (4*(la+alpha) * (la + gamma));
		double Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma));
		double a = x / (2*(la+alpha)*Psq);
		double b = x / (2*(nu+alpha)*Qsq);
		return a * dphi_dla(nu, la) + b * dphi_dnu(nu, la);
	}
	virtual double dphidz(double x, double z) {
		double nu, la;
		xy_to_conf(x, z, nu, la);
		double Psq = (la - nu) / (4*(la+alpha) * (la + gamma));
		double Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma));
		double c = z / (2*(la+gamma)*Psq);
		double d = z / (2*(nu+gamma)*Qsq);
		return c * dphi_dla(nu, la) + d * dphi_dnu(nu, la);
	}


	void xy_to_conf(double x, double y, double& nu, double& la) {
		double a1 = pow(x, 2) + pow(y, 2) -alpha - gamma;
		double a2 = pow(x, 4) + pow(y*y + alpha - gamma, 2) + 2 * x*x*(y*y-alpha + gamma);
		//printf("%f %f | %f %f %f\n", gamma, alpha, a1, a2, sqrt(a2));
		la = 0.5 * (a1 + sqrt(a2));
		nu = 0.5 * (a1 - sqrt(a2));
	}

/*
	def rho_conf(self, nu, la):
		g_la = (la + self.alpha) / (la - nu)
		g_nu = (nu + self.alpha) / (nu - la) 
		return g_la**2 * self.psi_tau(la) + g_nu**2 * self.psi_tau(nu) + 2 * g_la * g_nu * (self.Psi_tau(la) - self.Psi_tau(nu))/(la - nu)
	
	def rho_xy(self, x, y):
		nus, las = self.xy_to_conf(x, y)
		return self.rho_conf(nus, las)
*/

	virtual double density_conf(double nu, double la) {
		//printf("nu = %f la = %f\n", nu, la);
		double g_la = (la + alpha) / (la - nu);
		double g_nu = (nu + alpha) / (nu - la);
		return pow(g_la,2) * psi_tau(la) + pow(g_nu, 2) * psi_tau(nu) + 2 * g_la * g_nu * (Psi_tau(la) - Psi_tau(nu))/(la - nu);
	} 

	double Psi_tau(double tau) {
		//printf("tau = %f\n", tau);
		auto dpsi = [this](double sigma){ return this->psi_tau(sigma); };
		IntegratorGSL<> integratorGSL(dpsi); // the integrator
		double integral = integratorGSL.integrate(-gamma, tau);
		return integral;
	}
	virtual double G_tau_prime(double tau) {
		return 0;
	}
	virtual double G_tau(double tau) {
		/*C = 2 * pi * G / sqrt(abs(tau+self.alpha))
		#C = -2 * pi * G / sqrt(abs(tau+self.gamma))
		def f(sigma):
			return sqrt(abs(sigma+self.alpha))/(2*(sigma+self.gamma)) * self.Psi_tau(sigma)
			#return (sigma+self.alpha)/(2*(sigma+self.gamma)**(3./2)) * self.Psi_tau(sigma)
		I, err = scipy.integrate.quad(f, -self.alpha, tau)
		#I, err = scipy.integrate.quad(f, -self.gamma, tau)
		U = - 2 * pi * G * self.psiinf + C * I 
		return -U*/
		return 0;
	}
		
	double gamma, alpha, G;
};

class PerfectSpheroidOblate : public ProfileStackelOblate {
public:
	/* comment: q >= 1/sqrt(2), otherwise the density can become negative (BT) */
	PerfectSpheroidOblate(double rho0, double gamma, double alpha, double G) : ProfileStackelOblate(gamma, alpha, G), s(4), rho0(rho0) {
		c = sqrt(-gamma); 
	}

	double psi_tau(double tau) {
		return rho0 * pow(c, s) / pow(tau, s/2.);
	}

	double G_tau(double tau) {
		return 2 * M_PI * G * rho0 * -alpha * sqrt(-gamma/(tau+gamma)) * atan(sqrt((tau+gamma)/-gamma));
	}
	double G_tau_prime(double tau) {
		double x = (gamma+tau)/gamma;
		double D = pow(-1./x, 3./2) * (gamma*sqrt(-x) + tau * atan(sqrt(-x)))/(2*gamma * tau);
		return 2 * M_PI * G * rho0 * -alpha * D;
	}

	void dphi_dxdy(double x, double y, double& dphidx, double& dphidy) {
		double nu, la;
		xy_to_conf(x, y, nu, la);
		double Psq = (la - nu) / (4*(la+alpha) * (la + gamma));
		double Qsq = (nu - la) / (4*(nu+alpha) * (nu + gamma));
		double a = x / (2*(la+alpha)*Psq);
		double b = x / (2*(nu+alpha)*Qsq);
		double c = y / (2*(la+gamma)*Psq);
		double d = y / (2*(nu+gamma)*Qsq);
		double dphi_dla_ = dphi_dla(nu, la);
		double dphi_dnu_ = dphi_dnu(nu, la);
		dphidx = a * dphi_dla_ + b * dphi_dnu_;
		dphidx = c * dphi_dla_ + d * dphi_dnu_;
	}
	/*double density_conf(double nu, double la) {
		return rho0 * pow((alpha * gamma / (nu * la) ),  2);
	}*/
protected:
	double c, s, rho0;
};


}