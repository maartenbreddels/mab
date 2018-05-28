#pragma once
#include <math.h>


namespace gd {


class Transformation1d_in_3d {
public:
	virtual double transform(double r) = 0;
	virtual double inverse_transform(double r) = 0;
	virtual double d3xdu(double u) = 0;
	virtual double drdu(double u) = 0;
	virtual double laplace_u1_1(double u) = 0;
	virtual double laplace_u1_2(double u) = 0;
	virtual double laplace_u1_1_times_d3xdu(double u) {
		return this->laplace_u1_1(u) * this->d3xdu(u); 
	}
};

class TransformationSpherical_in_3d : public Transformation1d_in_3d {
public:
	virtual double transform(double u) {
		return u;
	}
	virtual double inverse_transform(double r) {
		return r;
	}
	virtual double d3xdu(double u) {
		return 4 * M_PI * u * u;
	}
	virtual double drdu(double) {
		return 1;
	}
	virtual double laplace_u1_1(double u) {
		return 1/(u*u);
	}
	virtual double laplace_u1_2(double u) {
		return u*u;
	}
	virtual double laplace_u1_1_times_d3xdu(double) {
		return 4 * M_PI;
	}
};

class TransformationSphericalLog_in_3d : public Transformation1d_in_3d {
public:
	double base, logbase;
	TransformationSphericalLog_in_3d(double base) : base(base) {
		logbase = log(base);
	}
	virtual double transform(double u) {
		return pow(base, u);
	}
	virtual double inverse_transform(double r) {
		return log(r)/logbase;
	}
	virtual double d3xdu(double u) {
		double r = transform(u);
		return 4 * M_PI * pow(r, 3) * logbase;
	}
	virtual double drdu(double u) {
		double r = transform(u);
		return r * logbase;
	}
	virtual double laplace_u1_1(double u) {
		double r = transform(u);
		return 1/(pow(r,3)*logbase*logbase);
	}
	virtual double laplace_u1_2(double u) {
		double r = transform(u);
		return r;
	}
	virtual double laplace_u1_1_times_d3xdu(double) {
		return 4 * M_PI / logbase;
	}
};

class TransformationSphericalTan_in_3d : public Transformation1d_in_3d {
public:
	virtual double transform(double u) {
		return tan(u * M_PI/2);
	}
	virtual double inverse_transform(double r) {
		return atan(r) * 2 / M_PI;
	}
	virtual double drdu(double u) {
		double c = cos(u*M_PI/2);
		return M_PI/2 / (c*c);
	}
	virtual double d3xdu(double u) {
		//double t = tan(u*M_PI/2);
		double s = sin(u*M_PI/2);
		double c = cos(u*M_PI/2);
		return 2 * M_PI * M_PI * pow(s, 2) / pow(c, 4);
	}
	virtual double laplace_u1_1(double u) {
		double s = tan(u*M_PI/2);
		double c = sin(u*M_PI/2);
		return 4 / (M_PI*M_PI) * pow(c, 4) / pow(s, 4);
	}
	virtual double laplace_u1_2(double u) {
		double t = tan(u*M_PI/2);
		double c = cos(u*M_PI/2);
		return t * t * c * c;
	}
	virtual double laplace_u1_1_times_d3xdu(double) {
		return 8; 
	}
};


// return 16./2 * basis1.dfdx((x-xleft)/dx) / dx * basis2.dfdx((x-xleft)/dx) / dx * t * t * c * c;
}