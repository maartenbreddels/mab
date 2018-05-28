#include "profile_2d.hpp"
//#include "profile_python.hpp"
#include <boost/python.hpp>
#include <pyublas/numpy.hpp>

namespace gd {
using namespace boost::python;

template<class ProfileType, class Ctor>
	void py_export_profile_potential_axi(string name) {
		class_<ProfileType, bases<Potential2d> >(name.c_str(), Ctor())
			.def("potentialrhoz", (&ProfileType::potentialxy))
			.def("potentialxy_eff", (&ProfileType::potentialxy_eff))
			;
	}
template<class ProfileType, class Ctor>
	void py_export_profile_profile_axi(string name) {
		class_<ProfileType, bases<Profile2d> >(name.c_str(), Ctor())
			.def("densityxy", (&ProfileType::densityxy))
			.def("potentialxy", (&ProfileType::potentialxy))
			.def("potentialxy_eff", (&ProfileType::potentialxy_eff))
			.def("dphidx", (&ProfileType::dphidx))
			.def("dphidy", (&ProfileType::dphidy))
			.def("dphidx_eff", (&ProfileType::dphidx_eff))
			.def("dphidy_eff", (&ProfileType::dphidy_eff))
			;
	}


void py_export_profile_2d() {
	class_< Potential2d, boost::noncopyable >("Potential2d", no_init);
	class_< Profile2d, boost::noncopyable >("Profile2d", no_init);
	class_< ProfileModel2d, boost::noncopyable >("ProfileModel", no_init)
		;
	class_<ProfileModel2d1C, bases<ProfileModel2d>, boost::noncopyable  >("ProfileModel2d1C", init<Profile2d*>())
		.def("densityxy", (&ProfileModel2d1C::densityxy))
		.def("potentialxy", (&ProfileModel2d1C::potentialxy))
		.def("potentialxy_eff", (&ProfileModel2d1C::potentialxy_eff))
		.def("dphidx", (&ProfileModel2d1C::dphidx))
		.def("dphidy", (&ProfileModel2d1C::dphidy))
		.def("dphidx_eff", (&ProfileModel2d1C::dphidx_eff))
		.def("dphidy_eff", (&ProfileModel2d1C::dphidy_eff))
		;
	/*class_<ProfileModel2C, bases<ProfileModel>, boost::noncopyable  >("ProfileModel2C", init<Profile*,Profile*>())
		.def("densityr", (&ProfileModel2C::densityr))
		.def("dphidr", (&ProfileModel2C::dphidr))
		.def("potentialr", (&ProfileModel2C::potentialr))
		;
	
	*/
		//.def("densityr", pure_virtual(&Profile::densityr))
	/*class_< Density2d, boost::noncopyable >("Density2d", no_init)
		//.def("densityr", pure_virtual(&Density::densityr))
		;*/
	py_export_profile_profile_axi<Logarithmic2d, init<double, double, double, double>>("Logarithmic2d");

}

}
	