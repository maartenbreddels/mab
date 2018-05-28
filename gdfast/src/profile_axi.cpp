#include "profile_axi.hpp"
//#include "profile_python.hpp"
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>

namespace gd {
using namespace boost::python;

template<class ProfileType, class Ctor>
	void py_export_profile_potential_axi(string name) {
		class_<ProfileType, bases<PotentialAxi> >(name.c_str(), Ctor())
			.def("potentialrhoz", (&ProfileType::potentialxz))
			.def("potentialxz_eff", (&ProfileType::potentialxz_eff))
			;
	}
template<class ProfileType, class Ctor>
	class_<ProfileType, bases<ProfileAxi> > py_export_profile_profile_axi(string name) {
		return class_<ProfileType, bases<ProfileAxi> >(name.c_str(), Ctor())
			.def("densityxz", (&ProfileType::densityxz))
			.def("potentialxz", (&ProfileType::potentialxz))
			.def("potentialxz_eff", (&ProfileType::potentialxz_eff))
			.def("dphidx", (&ProfileType::dphidx))
			.def("dphidz", (&ProfileType::dphidz))
			.def("dphidx_eff", (&ProfileType::dphidx_eff))
			.def("dphidz_eff", (&ProfileType::dphidz_eff))
			;
	}


void py_export_profile_axi() {
	class_< PotentialAxi, boost::noncopyable >("PotentialAxi", no_init);
	class_< ProfileAxi, boost::noncopyable >("ProfileAxi", no_init);
	class_< ProfileModelAxi, boost::noncopyable >("ProfileModel", no_init)
		;
	class_<ProfileModelAxi1C, bases<ProfileModelAxi>, boost::noncopyable  >("ProfileModelAxi1C", init<ProfileAxi*>())
		.def("densityxz", (&ProfileModelAxi1C::densityxz))
		.def("potentialxz", (&ProfileModelAxi1C::potentialxz))
		.def("potentialxz_eff", (&ProfileModelAxi1C::potentialxz_eff))
		.def("dphidx", (&ProfileModelAxi1C::dphidx))
		.def("dphidz", (&ProfileModelAxi1C::dphidz))
		.def("dphidx_eff", (&ProfileModelAxi1C::dphidx_eff))
		.def("dphidz_eff", (&ProfileModelAxi1C::dphidz_eff))
		;
	/*class_<ProfileModel2C, bases<ProfileModel>, boost::noncopyable  >("ProfileModel2C", init<Profile*,Profile*>())
		.def("densityr", (&ProfileModel2C::densityr))
		.def("dphidr", (&ProfileModel2C::dphidr))
		.def("potentialr", (&ProfileModel2C::potentialr))
		;
	
	*/
		//.def("densityr", pure_virtual(&Profile::densityr))
	/*class_< DensityAxi, boost::noncopyable >("DensityAxi", no_init)
		//.def("densityr", pure_virtual(&Density::densityr))
		;*/
	py_export_profile_profile_axi<LogarithmicAxi, init<double, double, double, double>>("LogarithmicAxi");
	py_export_profile_profile_axi<PerfectSpheroidOblate, init<double, double, double, double>>("PerfectSpheroidOblate")
		.def("density_conf", (&PerfectSpheroidOblate::density_conf))
		.def("G_tau", (&PerfectSpheroidOblate::G_tau))
		.def("G_tau_prime", (&PerfectSpheroidOblate::G_tau_prime))
		.def("dphi_dla", (&PerfectSpheroidOblate::dphi_dla))
		.def("dphi_dnu", (&PerfectSpheroidOblate::dphi_dnu))
		.def("dphi_dx", (&PerfectSpheroidOblate::dphidx))
		.def("dphi_dx", (&PerfectSpheroidOblate::dphidz))
		;

}

}
	