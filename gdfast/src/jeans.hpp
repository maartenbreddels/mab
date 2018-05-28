#pragma once
#include "profile.hpp"

namespace gd {
	
void py_export_jeans();

class JeansAnisotropicConstant {
public:
	JeansAnisotropicConstant(ProfileModel& stellar_profile, ProfileModel& halo_profile, double beta, double rmax);
	double sigmar(double r);

	double moment4r(double r);
	double dmoment4r(double r);

	double sigma_los(double R);
	double dsigma_los(double v, double R);
	double sigma_los_approx(double R);
	//double dsigmar(double r);
	//double dsigma_los_approx(double r, double R, double surface_density);
	double prob_losv(double R, double v, double sigma_obs);
	double dprob_losv(double r, double R, double v, double sigma_obs, double surface_density);
	ProfileModel &stellar_profile, &halo_profile;
	
	double beta, rmax;
	double integral_rmax_inf; 

};

};
