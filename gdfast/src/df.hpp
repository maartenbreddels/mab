#pragma once
#include "galaxy.hpp"
#include "datatypes.hpp"

namespace gd {
	
void py_export_df();

class DF {
public:
	DF(ProfileModel* profile_model, double_vector rs, double_vector rborders) : profile_model(profile_model), rs(rs), rborders(rborders)
	{}
	double velocities(double_vector density, double_vector vr, double_vector vt, double E, double L, double rp, double ra);
	double moments(double_vector density, double_vector varvr, double_vector varvt, double E, double L, double rp, double ra);
	double momentsL(double_vector density, double_vector varvr, double_vector varvt, double E, double L1, double L2, double rp1, double ra1, double rp2, double ra2, double L, double rp, double ra);
	double momentsL2(double_vector density, double_vector varvr, double_vector varvt, double_vector v22, double_vector v40, double_vector v04, double E, double L1, double L2, double rp, double ra, bool debug);
	double momentsE(double_vector density, double_vector varvr, double_vector varvt, double E1, double E2, double L, double rp1, double ra1, double rp2, double ra2, double E, double rp, double ra);
	//double moments_projected(double_vector mass3d, double_vector density, double_vector v, double E, double L, double rp, double ra, double N);
protected:
	ProfileModel* profile_model;
	double_vector rs;
	double_vector rborders;
	//double_vector Rs;
	//double_vector Rborders;
};

};
