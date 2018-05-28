#include "profile.hpp"
#include "profile_python.hpp"
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>
#include "optimization.hpp"
#include "rootfinder.hpp"

namespace gd {
using namespace boost::python;

// in km/s^2 Msun/kpc
const double GravConst = 4.30200406461e-06;

boost::python::tuple ProfileModel::Lmax_and_rcirc_at_E_(double E) {
	double Lmax, rcirc;
	Lmax_and_rcirc_at_E(E, Lmax, rcirc);
	return boost::python::make_tuple(Lmax, rcirc);
}
void ProfileModel::Lmax_and_rcirc_at_E(double E, double& Lmax, double& rcirc) {
	/*	def f(r, E=E):
			r = abs(r)
			return -r**2*(2*(E-self.potentialr(r)))
		rLmax, Lmaxsq = scipy.optimize.fmin(f, rtry, full_output=True, disp=False)[:2]
		return sqrt(abs(Lmaxsq)), abs(rLmax[0])
	*/
	auto f = [&](double logr) -> double {
		//r = fabs(r);
		double r = pow(10, logr);
		//printf("r = %f logr\n ", r);
		return -r*r*(2*(E-this->potentialr(r)));
	};
	MinimizerNLopt2<> opt(f);
	//double lower_bounds[1] = { 0};
	//nlopt_set_lower_bounds(opt.opt, lower_bounds);
	//opt.x = 1.;
	rcirc = pow(10, opt.optimize());
	//printf("r = %f Ekin = %f\n", r, 2*(E-this->potentialr(r)));
	Lmax = sqrt(2*(E-this->potentialr(rcirc)))*rcirc;
}
double ProfileModel::Lmax_at_E(double E) {
	double Lmax, rcirc;
	Lmax_and_rcirc_at_E(E, Lmax, rcirc);
	return Lmax;
}
double ProfileModel::rcirc_at_E(double E) {
	double Lmax, rcirc;
	Lmax_and_rcirc_at_E(E, Lmax, rcirc);
	return rcirc;
}
double ProfileModel::rmax_at_E(double E, double rcirc) {
	auto f = [&](double r) -> double {
		return E - this->potentialr(r);
	};
	double scale = 1.5;
	double r_far_way = rcirc * scale;
	assert(f(rcirc) >= 0);
	while(f(r_far_way) > 0) {
		r_far_way *= scale;
	}
	RootFinderGSL<> rootfinder(f);
	double rmax = rootfinder.findRoot(rcirc, r_far_way);
	return rmax;
	
	/*auto f = [&](double r) -> double {
		return E - this->potentialr(r);
	};
	auto root = [&](double x, double * y, double *dy) -> void {
		double r = pow(10, x) + rcirc;
		// x = log10(r-rcirc)
		// dx = 1/(r-rcirc) / log(10) dr
		*y = f(r);
		double gradient = -this->dphidr(r) * (r-rcirc) * log(10);
		*dy = gradient;
		//*
		double r1 = r;
		double dx = 1e-4;
		double r2 = pow(10, x+dx)+rcirc;
		double y1 = f(r1);
		double y2 = f(r2);
		printf("rmax grad: logr=%f r=%10.5f: erad=%10.5f [%10.5f %10.5f] (%f) \n", x, r, *y, (y2-y1)/dx, gradient, this->potentialr(r));
		//return g
		/**/
	/*};
	
	// find apo or peri center
	RootFinderGSLDerivative<> rootfinderFirst(root);
	double rstart = pow(10, rootfinderFirst.findRoot(1, 0, 1e-5))+rcirc; //pow(10, rootfinderFirst.findRoot(0));
	double scale = 1.+1e-4;
	assert(f(rstart*scale) < 0);
	assert(f(rstart/scale) > 0);
	return rstart;*/
}


boost::python::tuple ProfileModel::get_apo_peri(double E, double L, double rmin, double rcirc, double rmax) {

	auto ekinrad = [&](double r) -> double { // kinetic energy in radial direction
		return (-L*L/(2*r*r) - (this->potentialr(r) - E));
	};
	auto root = [&](double x, double * y, double *dy) -> void {
		if(x < 0) {
			//*y =FP_INFINITE;// inf;
			//*dy = FP_INFINITE;//inf;
			//return;
			//printf("grad: x=%f\n", x);
			//x = -x;
		}
		double r = pow(10, x);
		*y = ekinrad(r);
		double gradient = (-this->dphidr(r) + L*L/pow(r, 3)) * r * log(10);
		*dy = gradient;
		//*
		double r1 = r;
		double dx = 1e-4;
		double r2 = pow(10, x+dx);
		double y1 = ekinrad(r1);
		double y2 = ekinrad(r2);
		//printf("grad: logr=%f r=%10.5f: erad=%10.5f [%10.5f %10.5f] (%f) \n", x, r, *y, (y2-y1)/dx, gradient, this->potentialr(r));
		//return g
		/**/
	};

	// find apo or peri center
	//RootFinderGSLDerivative<> rootfinderFirst(root);
	//double rstart = pow(10, rootfinderFirst.findRoot(1)); //pow(10, rootfinderFirst.findRoot(0));

	double scale = (1+1e-4);
	double apo, peri;
	RootFinderGSL<> rootfinder(ekinrad);
	//printf("found a valid r %f (E=%f L=%f Eleft=%f (%f %f))\n", rstart, E, L, ekinrad(rstart), ekinrad(rstart*scale), ekinrad(rstart/scale));
	//printf("rmin : %f ekinrad: %f %f %f\n", rmin, ekinrad(rmin), ekinrad(rmin*scale), ekinrad(rmin/scale));
	//printf("rcirc: %f ekinrad: %f %f %f\n", rcirc, ekinrad(rcirc), ekinrad(rcirc*scale), ekinrad(rcirc/scale));
	//printf("rmax : %f ekinrad: %f %f %f\n", rmax, ekinrad(rmax), ekinrad(rmax*scale), ekinrad(rmax/scale));
	int n;
	n = 0;
	assert(ekinrad(rcirc) >= 0);
	while((n < 10) & (ekinrad(rmin) > 0)) {
		rmin /= scale;
		n++;
	}
	if(ekinrad(rmin) > 0) {
		throw std::range_error("rmin cannot be properly found");
	}

	n = 0;
	while((n < 10) & (ekinrad(rmax) > 0)) {
		rmax *= scale;
		n++;
	}
	if(ekinrad(rmax) > 0) {
		throw std::range_error("rmax cannot be properly found");
	}


	peri = rootfinder.findRoot(rmin, rcirc);
	apo = rootfinder.findRoot(rcirc, rmax);
	
		//printf("finding apocenter [%f]\n", rstart);
		//printf("found apocenter [%f]\n", apo);
		//printf("finding pericenter [%f] (%f %f %f)\n", rstart, ekinrad(1e-6), ekinrad(rstart/scale/scale), ekinrad(rstart*scale*scale));
		//printf("found pericenter [%f]\n", peri);
	scale = (1+1e-5);
	while(ekinrad(peri) < 0) {
		peri *= scale;
	}
	while(ekinrad(apo) < 0) {
		apo /= scale;
	}
	return boost::python::make_tuple(apo/scale, peri*scale);
	
	/*
	# just find apo or pericenter
	rstart = self.findR_at_EL(E, L, rtry)
	try:
		rstart = rstart[0]
	except:
		pass
	s = (1+1e-5) # scale factor for testing apo/peri
	r = rstart
	#print rstart, ekinrad(rstart/s), ekinrad(rstart*s), -L**2/(2*r**2) - (self.potentialr(r) - E), -L**2/(2*r**2) , (self.potentialr(r) - E)
	if (ekinrad(rstart/s) < 0) and (ekinrad(rstart*s) > 0): # we found pericenter
		rp = rstart
		ra = brentq(ekinrad, rstart*s, 1e9)
	else: # we found apocenter
		ra = rstart
		rp = brentq(ekinrad, 1e-9, rstart/s)
	
	# sanity checks
	assert ekinrad(ra*s) < 0, "available energy ar r > r_apo should be negative"
	assert ekinrad(rp/s) < 0, "available energy ar r < r_peri should be negative"
	assert ekinrad(ra/s) > 0, "available energy ar r < r_apo should be positive"
	assert ekinrad(rp*s) > 0, "available energy ar r > r_peri should be positive"
	assert ra > rp, "apocenter should be larger than pericenter" 
	return ra/s, rp*s*/
}



void py_export_profile() {
	class_< ProfileModel, boost::noncopyable >("ProfileModel", no_init)
		.def("get_apo_peri", (&ProfileModel::get_apo_peri))
		.def("Lmax_at_E", (&ProfileModel::Lmax_at_E))
		.def("Lmax_and_rcirc_at_E", (&ProfileModel::Lmax_and_rcirc_at_E_))
		.def("rcirc_at_E", (&ProfileModel::rcirc_at_E))
		.def("rmax_at_E", (&ProfileModel::rmax_at_E))
		
		;
	class_<ProfileModel1C, bases<ProfileModel>, boost::noncopyable  >("ProfileModel1C", init<Profile*>())
		.def("densityr", (&ProfileModel1C::densityr))
		.def("dphidr", (&ProfileModel1C::dphidr))
		.def("potentialr", (&ProfileModel1C::potentialr))
		;
	class_<ProfileModel2C, bases<ProfileModel>, boost::noncopyable  >("ProfileModel2C", init<Profile*,Profile*>())
		.def("densityr", (&ProfileModel2C::densityr))
		.def("dphidr", (&ProfileModel2C::dphidr))
		.def("potentialr", (&ProfileModel2C::potentialr))
		;
	
	class_< Profile, boost::noncopyable >("Profile", no_init)
		.def("densityr", pure_virtual(&Profile::densityr))
		;
	class_< Density, boost::noncopyable >("Density", no_init)
		.def("densityr", pure_virtual(&Density::densityr))
		;
	//py_export_profile_profile<Plummer, init<double, double, double>>("Plummer");
	//init<double, double, double>((boost::python::arg("b"), boost::python::arg("M"), boost::python::arg("G")))
	py_export_profile_profile_kw<Plummer>("Plummer", init<double, double, double>((boost::python::arg("M"), boost::python::arg("b"), boost::python::arg("G")=GravConst)));
	/*class_<Plummer, bases<Profile, Density> >("Plummer", init<double, double, double>())
		.def("potentialr", (&Plummer::potentialr))
		.def("densityr", (&Plummer::densityr))
		.def("densityR", (&Plummer::densityR))
		.def("dphidr", (&Plummer::dphidr))
		//.def("dphidr2", (&Plummer::dphidr2))
		;
	*/
	//py_export_profile_profile_kw<ProjectedExponential>("ProjectedExponential", init<double, double, double>((boost::python::arg("M"), boost::python::arg("scale"), boost::python::arg("G")=GravConst)));
	py_export_profile_profile_kw<ProjectedExponential>("ProjectedExponential", init<double, double, double>());
	
	
	class_<TestCase, bases<Profile, Density> >("TestCase", init<double, double, double>())
		.def("potentialr", (&TestCase::potentialr))
		.def("densityr", (&TestCase::densityr))
		.def("densityR", (&TestCase::densityR))
		.def("dphidr", (&TestCase::dphidr))
		//.def("dphidr2", (&TestCase::dphidr2))
		;
	class_<LogarithmicProfile, bases<Profile, Density> >("LogarithmicProfile", init<double, double>())
		.def("densityr", (&LogarithmicProfile::densityr))
		.def("densityR", (&LogarithmicProfile::densityR))
		.def("dphidr", (&LogarithmicProfile::dphidr))
		//.def("dphidr2", (&LogarithmicProfile::dphidr2))
	;
	class_<Isochrone, bases<Profile, Density> >("Isochrone", init<double, double, double>())
		.def("densityr", (&Isochrone::densityr))
		.def("densityR", (&Isochrone::densityR))
		.def("dphidr", (&Isochrone::dphidr))
		//.def("dphidr2", (&Isochrone::dphidr2))
		.def("potentialr", (&Isochrone::potentialr))
		;
	class_<NullProfile, bases<Profile, Density> >("NullProfile", init<>())
		.def("densityr", (&NullProfile::densityr))
		.def("densityR", (&NullProfile::densityR))
		.def("dphidr", (&NullProfile::dphidr))
		//.def("dphidr2", (&NullProfile::dphidr2))
	;
	class_<Hernquist, bases<Profile, Density> >("Hernquist", init<double, double, double>())
		.def("densityr", (&Hernquist::densityr))
		.def("densityR", (&Hernquist::densityR))
		.def("dphidr", (&Hernquist::dphidr))
		.def("potentialr", (&Hernquist::potentialr))
		;
	class_<Jaffe, bases<Profile, Density> >("Jaffe", init<double, double, double>())
		.def("densityr", (&Jaffe::densityr))
		.def("densityR", (&Jaffe::densityR))
		.def("dphidr", (&Jaffe::dphidr))
		.def("potentialr", (&Jaffe::potentialr))
		;
	class_<Einasto, bases<Profile, Density> >("Einasto", init<double, double, double, double>())
		.def("densityr", (&Einasto::densityr))
		.def("densityR", (&Einasto::densityR))
		.def("dphidr", (&Einasto::dphidr))
		.def("potentialr", (&Einasto::potentialr))
		;
	class_<Burkert, bases<Profile, Density> >("Burkert", init<double, double, double>())
		.def("densityr", (&Burkert::densityr))
		.def("densityR", (&Burkert::densityR))
		.def("dphidr", (&Burkert::dphidr))
		.def("potentialr", (&Burkert::potentialr))
		;
	class_<NFW, bases<Profile, Density> >("NFW", init<double, double, double, double>())
		.def("densityr", (&NFW::densityr))
		.def("potentialr", (&NFW::potentialr))
		.def("dphidr", (&NFW::dphidr))
		.def_readonly("c", &NFW::c)
		.def_readonly("mass200", &NFW::mass200)
		.def_readonly("r200", &NFW::r200)
		.def_readonly("rs", &NFW::rs)
		.def_readonly("rho0", &NFW::rho0);
	;
	class_<NFWCut, bases<Density> >("NFWCut", init<double, double, double>())
		.def("densityr", (&NFWCut::densityr))
		.def_readwrite("rs", &NFWCut::rs)
		.def_readwrite("rte", &NFWCut::rte)
		.def_readwrite("rho0", &NFWCut::rho0)
		;
	class_<TwoSlopeDensity, bases<Density> >("TwoSlopeDensity", init<double, double, double, double, double>())
		.def("densityr", (&TwoSlopeDensity::densityr))
		.def_readwrite("alpha", &TwoSlopeDensity::alpha)
		.def_readwrite("beta", &TwoSlopeDensity::beta)
		.def_readwrite("rho0", &TwoSlopeDensity::rho0)
		.def_readwrite("rs", &TwoSlopeDensity::rs)
		.def_readwrite("gamma", &TwoSlopeDensity::gamma)
		;
	class_<TwoSlopeDensityCut, bases<Density> >("TwoSlopeDensityCut", init<double, double, double, double, double, double>())
		.def("densityr", (&TwoSlopeDensityCut::densityr))
		.def_readwrite("alpha", &TwoSlopeDensityCut::alpha)
		.def_readwrite("beta", &TwoSlopeDensityCut::beta)
		.def_readwrite("rho0", &TwoSlopeDensityCut::rho0)
		.def_readwrite("rs", &TwoSlopeDensityCut::rs)
		.def_readwrite("gamma", &TwoSlopeDensityCut::gamma)
		.def_readwrite("rte", &TwoSlopeDensityCut::rte)
		;
	class_<BrokenPowerLawDensitySoft3, bases<Density> >("BrokenPowerLawDensitySoft3", init<double, double, double, double, double, double, double, double>())
		.def("densityr", (&BrokenPowerLawDensitySoft3::densityr))
		.def_readwrite("rho0", &BrokenPowerLawDensitySoft3::rho0)
		.def_readwrite("s1", &BrokenPowerLawDensitySoft3::s1)
		.def_readwrite("s2", &BrokenPowerLawDensitySoft3::s2)
		.def_readwrite("s3", &BrokenPowerLawDensitySoft3::s3)
		.def_readwrite("gamma1", &BrokenPowerLawDensitySoft3::gamma1)
		.def_readwrite("gamma2", &BrokenPowerLawDensitySoft3::gamma2)
		.def_readwrite("rs1", &BrokenPowerLawDensitySoft3::rs1)
		.def_readwrite("rs2", &BrokenPowerLawDensitySoft3::rs2)
		;
	
}


}
